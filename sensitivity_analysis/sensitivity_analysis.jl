#Load dependencies
using StatsPlots,Dates,JLD2,Statistics,Optim,Parameters,Distributions,DataFrames,CSV
using DynamicHMC,OrdinaryDiffEq,DiffEqCallbacks
using NamedArrays,TransformVariables,LogDensityProblems,LinearAlgebra,BlackBoxOptim,Random
using DynamicHMC.Diagnostics,MCMCDiagnostics,MCMCChains,ForwardDiff,LogExpFunctions
using Plots.PlotMeasures,LogExpFunctions
using CountTimeSeries
import KenyaCoVSD
include("inference_variable_transmissibility.jl");
## Gather the fits for each scenatio
# waning_immunity_scenarios = [(σ = 0.0,ω = 1/180,scenario_name = "_no_waning_immunity"),
#               (σ = 0.5,ω = 1/180,scenario_name = "_med_susceptibility"),
#               (σ = 1.0,ω = 1/180,scenario_name = "_full_susceptibility"),
#               (σ = 0.16,ω = 1/90,scenario_name = "_fast_reversion"),
#               (σ = 0.16,ω = 1/365,scenario_name = "_slow_reversion"),
#               (σ = 1.0,ω = 1/(5*365),scenario_name = "_slow_reversion_to_full_susceptibility")]
# infectiousness_scenarios = [(σ = 0.16,ω = 1/180,ι = 0.5, scenario_name = "reduced_inf"),
#             (σ = 0.5,ω = 1/180,ι = 0.5,scenario_name = "_med_susceptibility_reduced_inf"),
#             (σ = 1.0,ω = 1/180,ι = 0.5,scenario_name = "_full_susceptibility_reduced_inf"),
#             (σ = 1.0,ω = 1/90,ι = 0.5,scenario_name = "_full_susceptibility_fast_wane_reduced_inf")]

## Gather the fits for the different sensitivity analyses
all_scenarios = []

@load("modelfits/Nairobi_model.jld2")
push!(all_scenarios,(σ = 0.16,ω = 1/180,ι = 1.0,scenario_name = "Baseline (mean 180 waning -> 16% sus. + 100% inf.)",chain = model.MCMC_results.chain,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_model_no_waning_immunity.jld2")
push!(all_scenarios,(σ = 0.0,ω = 1/180,ι = 1.0,scenario_name = "No waning immunity",chain = model.MCMC_results.chain,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_model_med_susceptibility.jld2")
push!(all_scenarios,(σ = 0.5,ω = 1/180,ι = 1.0,scenario_name = "Wane -> 50% sus.",chain = model.MCMC_results.chain,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_model_full_susceptibility.jld2")
push!(all_scenarios,(σ = 1.0,ω = 1/180,ι=1.0,scenario_name = "Wane -> 100% sus.",chain = model.MCMC_results.chain,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_model_fast_reversion.jld2")
push!(all_scenarios,(σ = 0.16,ω = 1/90,ι = 1.0, scenario_name = "mean 90 day waning",chain = model.MCMC_results.chain,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_model_slow_reversion.jld2")
push!(all_scenarios,(σ = 0.16,ω = 1/365,ι = 1.0, scenario_name = "mean 1 year waning",chain = model.MCMC_results.chain,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_model_slow_reversion_to_full_susceptibility.jld2")
push!(all_scenarios,(σ = 1.0,ω = 1/(5*365),ι = 1.0, scenario_name = "mean 5 year waning -> 100% sus.",chain = model.MCMC_results.chain,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_modelreduced_inf.jld2")
push!(all_scenarios,(σ = 0.16,ω = 1/180,ι = 0.5, scenario_name = "Wane -> 50% inf.",chain = model.MCMC_results.chain,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_model_med_susceptibility_reduced_inf.jld2")
push!(all_scenarios,(σ = 0.5,ω = 1/180,ι = 0.5,scenario_name = "Wane -> 50% inf. + 50% sus.",chain = model.MCMC_results.chain,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_model_full_susceptibility_reduced_inf.jld2")
push!(all_scenarios,(σ = 1.0,ω = 1/180,ι = 0.5,scenario_name = "Wane -> 50% inf. + 100% sus.",chain = model.MCMC_results.chain,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_model_full_susceptibility_fast_wane_reduced_inf.jld2")
push!(all_scenarios,(σ = 1.0,ω = 1/90,ι = 0.5,scenario_name = "mean 90 day Wane -> 50% inf. + 100% sus.",chain = model.MCMC_results.chain,DIC = KenyaCoVSD.modeldic(model)))

## Gather the posterior means and quantiles for each scenario and for each parameter
meanarray = zeros(size(all_scenarios[1].chain,2),11)
for k = 1:11
    meanarray[:,k] .= mean(all_scenarios[k].chain)[:,2]
end
lb_array = zeros(size(all_scenarios[1].chain,2),11)
for k = 1:11
    Q = quantile(all_scenarios[k].chain)[:,:]
    lb_array[:,k] .= Q[:,1]
end
ub_array = zeros(size(all_scenarios[1].chain,2),11)
for k = 1:11
    Q = quantile(all_scenarios[k].chain)[:,:]
    ub_array[:,k] .= Q[:,end]
end
variable_names = string.(keys(all_scenarios[1].chain))
scenario_name = [scenario.scenario_name for scenario in all_scenarios]
## For each of 15 parameters, make a scatter plot of posterior mean estimate for each sentivity scenario (95% CIs are
# xerror bars). These 15 plots are saved.
#Nicer title names for parameters
variable_title_name = Dict("ct_min1" => "% mobility of lower SES mid-April 2020",
                            "R₀" => "Baseline reproductive number",
                            "ϵ" => "In-group mixing preference",
                            "χ₁" => "Bias in favour of pos. PCR (lower SES)",
                            "χ₂" => "Bias in favour of pos. PCR (higher SES)",
                            "p_test₁" => "Per day rate of getting tested if positive (lower SES)",
                            "p_test₂" => "Per day rate of getting tested if positive (higher SES)",
                            "p_choose1" => "Prop. tests among lower SES",
                            "P_eff" => "Prop. population among lower SES",
                            "schooleffect" => "Reproductive number multiplier when schools shut",
                            "extra_transmissibility" => "Rel. trans new variants compared to baseline",
                            "influx_exposed_new_variant" => "Inc. number exposed Feb 2021 (new variants)",
                            "p_test_boost" => "Boost to testing rate (new variants)",
                            "χ_boost" => "Change in bias in favour of pos. PCR (new variants)",
                            "E₀" => "Number of initial exposed 20th Feb 2020")
for k = 1:15 #Loop over parameters
    plt = scatter(meanarray[k,:],1:11,
        xerror = (meanarray[k,:] .- lb_array[k,:],ub_array[k,:] .- meanarray[k,:] ),
        yticks = (1:11,scenario_name),
        title = "Parameter: "*variable_title_name[variable_names[k]],
        lab = "",
        left_margin = 5mm,right_margin = 5mm,
        size = (800,800),dpi = 250,
        guidefont = 12,ytickfont = 9,titlefont = 10,xtickfont = 9)
        savefig(plt,"plots/sensitivity_analysis_plots/$(variable_names[k]).png")
end

