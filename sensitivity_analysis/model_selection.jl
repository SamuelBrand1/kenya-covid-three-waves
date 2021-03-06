## This script gives a sensitivity analysis of models that fit to Nairobi data

using StatsPlots,Dates,JLD2,Statistics,Optim,Parameters,Distributions,DataFrames,CSV
using DynamicHMC,OrdinaryDiffEq,DiffEqCallbacks
using NamedArrays,TransformVariables,LogDensityProblems,LinearAlgebra,BlackBoxOptim,Random
using DynamicHMC.Diagnostics,MCMCDiagnostics,MCMCChains,ForwardDiff,LogExpFunctions
using Plots.PlotMeasures,LogExpFunctions
using CountTimeSeries
import KenyaCoVSD
##Load data for all Kenya
@load("data/linelist_data_with_pos_neg_20feb_to_27apr_c_age.jld2")
@load("data/cleaned_linelist20210521_c.jld2")
@load("data/deaths_data20210521.jld2")
@load("data/serological_data_with_20feb_to_31dec_age.jld2")
@load("data/serological_data_with_20feb_to_10Mar2021_age.jld2")
@load("data/N_kenya.jld2")

##Create a `CovAreaModel` object for Kenya which we will used with alternative models explaining the data
#`basic_nai_model` contains the PCR, serological and death data for Nairobi.

basic_nai_model = KenyaCoVSD.CoVAreaModel("Nairobi",θ -> 0.,θ -> 0.0;
                                            case_data = linelist_data_with_pos_neg,
                                            sero_data = serological_data,
                                            death_data = deaths_data,
                                            pop_data = N_kenya)


##Create a table of model score estimates

modelscores = DataFrame(model = String[],DIC = Union{Float64,String}[],ΔDIC = Union{Float64,String}[])



## Model 1: paper model
@load("modelfits/Nairobi_model.jld2")
nai_two_group = deepcopy(model)
DIC_two_group = KenyaCoVSD.modeldic(nai_two_group)
push!(modelscores,["Two-group SES",DIC_two_group,DIC_two_group-DIC_two_group])


## Model 2: One group SEIR with fitted R(t)
##Load code dependencies -- These define the transmission model for the one-group (daily varying R(t))
#They also define a log-likelihood function, priors and code for the EM algorithm using a combination 
#of HMC inference (E-Step) and ADAM optimisation for fitting contact rate (ct) (M-step).
# R(t) = R₀ * ct (t)

include("one_group_model.jl");
include("one_group_inference.jl");
include("one_group_fitting_ct.jl");

## 
# Set the ADAM optimisation parameters
baseline_adam_params = (β₁ = 0.9,
                        β₂ = 0.9,
                        η = 0.05,
                        ϵ = 1e-8,
                        total_num_steps = 500,
                        averaging_num_steps=100,
                        λ = 1e6)


                        ## Get first draw from MCMC with baseline ct estimate
## NB: The EM algorithm approach is slow, there are saved results accessed in the next code
inferparameters!(nai_one_group,2000,trans_one_group,0.05,D,cts₀)
# and get the initial E[log L]
final_LL = 0
for k = 1:size(nai_one_group.MCMC_results.chain,1)
        θ = NamedTuple{Tuple(keys(nai_one_group.MCMC_results.chain))}([nai_one_group.MCMC_results.chain[k,n,1] for n = 1:size(nai_one_group.MCMC_results.chain,2)])
        final_LL += ll_onegroup_newvariant_infboost_ct(θ,nai_one_group,0.0,cts₀)
end
final_LL = (final_LL/size(nai_one_group.MCMC_results.chain,1)) - ct_penalty(baseline_adam_params.λ,cts₀)

EM_steps = [(cts₀,final_LL)]

# Run alternate MCMC (for parameters not contact rate ct) and ADAM (for ct)
for k = 1:5
        log_ct_fit,mean_log_LL = ADAM_optim(log.(EM_steps[end][1][25:end]),
                                                random_grad_log_ct!,
                                                baseline_adam_params,
                                                nai_one_group,
                                                baseline_ct,
                                                trans_one_group,
                                                0.0)

        # Add each iteration to EM_steps                                        
        push!(EM_steps,(vcat(baseline_ct,exp.(log_ct_fit)),mean_log_LL))
        inferparameters!(nai_one_group,500,trans_one_group,0.05,D,EM_steps[end][1])
end

#Uncomment to save the sequence of fitted contact 
# @save("sensitivity_analysis/nai_EM_fit.jld2",nai_one_group,EM_steps)
## Plot the one-group model fit

@load("sensitivity_analysis/nairobi_fits/nai_EM_fit.jld2")
@load("forecasts/condensed_county_forecasts.jld2")
include("../analysis_scripts/fitting_methods.jl");
include("../analysis_scripts/plotting_methods.jl")


ct_fitted = EM_steps[end][1]
N = sum(N_kenya[:,"Nairobi"])
ct_plt = plot_ct(ct_fitted)


##Calculate predictions

fit_mats = gather_uncertainty_one_group(nai_one_group,ct_fitted)

#Join the % PCR pos and number PCR pos predictions based on whether each day has neg PCR test counts or not
nai_tests = vcat(nai_one_group.PCR_cases[1:(end-14),2],fill(-17,14+size(fit_mats.prop_PCR_pos_mat,1) - size(nai_one_group.PCR_cases,1)))
pred_num_PCR_pos_mat = fit_mats.prop_PCR_pos_mat.*nai_tests.*(nai_tests .>= 0) .+ fit_mats.no_neg_PCR_pos_mat.*(nai_tests .< 0)
#Get the posterior mean and credible intervals
pred_prop_sero_pos = get_credible_intervals(fit_mats.prop_sero_pos_mat)
pred_num_PCR_pos = get_credible_intervals(pred_num_PCR_pos_mat)
pred_incidence = get_credible_intervals(fit_mats.infection_mat)
#Compare to two-group fits
nai_fit = condensed_county_forecasts[[fit.name == "Nairobi" for fit in condensed_county_forecasts]][1]
sero_plt_compare = plot_pop_exposure(nai_fit,serological_data,serology_data,N_kenya);
PCR_plt_compare = plot_PCR(nai_fit,linelist_data_with_pos_neg,linelist_data);

plot!(PCR_plt_compare, 4:(length(pred_num_PCR_pos.pred)-3),
                        weekly_mv_av(pred_num_PCR_pos.pred),
                        title = "Nairobi PCR test, model comparison",
                        color = :green,lw = 3,
                        ribbon = (weekly_mv_av(pred_num_PCR_pos.lb),weekly_mv_av(pred_num_PCR_pos.ub)),
                        lab = "One group model: fit and forecast")

plot!(sero_plt_compare,pred_prop_sero_pos.pred./nai_one_group.N,
        title = "Nairobi population exposure, model comparison",
        lab = "One group model fit: seroposivity",lw = 3,color = :green,ls = :dash,
        ribbon = (pred_prop_sero_pos.lb./nai_one_group.N,pred_prop_sero_pos.ub./nai_one_group.N))

plot!(sero_plt_compare,cumsum(pred_incidence.pred,dims = 1)./nai_one_group.N,
        lab = "One group model fit: Overall population exposure",lw = 3,color = :red,ls = :dash,
        ribbon = (cumsum(pred_incidence.lb,dims = 1)./nai_one_group.N,cumsum(pred_incidence.ub,dims = 1)./nai_one_group.N )
        )

savefig(sero_plt_compare,"plots/sensitivity_analysis_plots/comparison_onegrp_vs_twogrp_pop_exposure.png")
savefig(PCR_plt_compare,"plots/sensitivity_analysis_plots/comparison_onegrp_vs_twogrp_PCR_pos.png")

## Fill table for model 2
@load("sensitivity_analysis/nairobi_fits/nai_EM_fit.jld2") #<---- Saved EM algorithm fits
DIC_one_group_model = KenyaCoVSD.modeldic(nai_one_group)

push!(modelscores,["One-group SES",DIC_one_group_model,DIC_one_group_model-DIC_two_group])


## Model 3: lower 25%, middle 50% and upper 25% wealth percentiles

#Include code for three group model, log-likelihood, prior
include("three_group_model.jl");

#Run MCMC for three group model
KenyaCoVSD.inferparameters!(nai_three_group,2000,trans_three_groups,0.05,D,q₀;serowaningrate = 0.0)
#uncomment to save
# @save("nai_three_group.jld2",nai_three_group)

##
@load("sensitivity_analysis/nairobi_fits/nai_three_group.jld2")

DIC_three_group_model = KenyaCoVSD.modeldic(nai_three_group)
push!(modelscores,["Three-group SES",DIC_three_group_model,DIC_three_group_model-DIC_two_group])

## Read-out the model-score table
display(modelscores)
