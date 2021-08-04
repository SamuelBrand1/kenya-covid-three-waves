## This script gives a sensitivity analysis of models that fit to Nairobi data


## Active and precompile the KenyaCoVSD environment
using Pkg
Pkg.activate(".");Pkg.precompile()
using StatsPlots,Dates,JLD2,Statistics,Optim,Parameters,Distributions,DataFrames,CSV
using DynamicHMC,OrdinaryDiffEq,DiffEqCallbacks
using NamedArrays,TransformVariables,LogDensityProblems,LinearAlgebra,BlackBoxOptim,Random
using DynamicHMC.Diagnostics,MCMCDiagnostics,MCMCChains,ForwardDiff
using Plots.PlotMeasures
using CountTimeSeries
import KenyaCoVSD
#Load data for all Kenya
@load("data/linelist_data_with_pos_neg_20feb_to_27apr_c_age.jld2")
@load("data/cleaned_linelist20210521_c.jld2")
@load("data/deaths_data20210521.jld2")
@load("data/serological_data_with_20feb_to_31dec_age.jld2")
@load("data/serological_data_with_20feb_to_10Mar2021_age.jld2")
@load("data/N_kenya.jld2")

#Create a `CovAreaModel` object for Kenya which we will mutate with alternative models explaining the data
#`basic_nai_model` contains the PCR, serological and death data for Nairobi.

basic_nai_model = KenyaCoVSD.CoVAreaModel("Nairobi",θ -> 0.,θ -> 0.0;
                                            case_data = linelist_data_with_pos_neg,
                                            sero_data = serological_data,
                                            death_data = deaths_data,
                                            pop_data = N_kenya)


##Create a table of model score estimates

modelscores = DataFrame(model = String[],
                        BIC = Union{Float64,String}[],
                        AIC = Union{Float64,String}[],
                        DIC = Union{Float64,String}[],
                        lpd_cases = Float64[],
                        lpd_serology = Float64[])

## Model 1: Simple INGARCH model
#Can get find the best AIC model here, but what about lookahead predictions?
#Find the best (p,1) INGARCH model using AIC criterion
y = basic_nai_model.PCR_cases[1:(end-13),1]

function BIC_and_AIC_by_p(p)
        res = fit(y, Model(pastObs = 1:p, pastMean = 1,
                        distr = "Poisson",link = "Log"))
        BIC(res,55),AIC(res,55)
end

BICs_AICs = [BIC_and_AIC_by_p(p) for p = 1:30]
BICs = [x[1] for x in BICs_AICs]
AICs = [x[2] for x in BICs_AICs]

scatter(AICs)
findmin(AICs)
findmin(BICs)
findmin(AICs)

findmin(BICs)
findmin(AICs)

model = Model(pastObs = 1:24, pastMean = 1:1,
                distr = "Poisson",link = "Log")
# y = simulate(1000, model, [10, 0.5, 0.2])[1]
res = fit(y, model)
p = predict(res, 30,1000)
plot(y,lab = "All data",legend = :topleft)
# plot!(1:length(y),y,lab = "Fitting data")
plot!((length(y)+1):(length(y)+length(p[1])),p[1],lab = "Prediction",lw = 3,
        ribbon = (p[1] .- p[2][1,:],p[2][2,:] .- p[1]) )
# AIC(res,55)

plot(basic_nai_model.PCR_cases[:,1],lab = "All data",legend = :topleft)
plot!(res.λ)
plot(cumsum(basic_nai_model.PCR_cases[:,1])./basic_nai_model.N,lab = "All data",legend = :topleft)
plot!(cumsum(res.λ)./basic_nai_model.N)
p = predict(res, 30,1000)
plot(p[1])

## Model 2: One group SEIR with fitted R(t)

include("one_group_model.jl");
include("one_group_inference.jl");

inferparameters!(nai_one_group,2000,trans_one_group,0.05,D,cts₀)
# @save("sensitivity_analysis/nai_one_group_rd_1.jld2",nai_one_group,cts₀)

baseline_adam_params = (β₁ = 0.9,
                        β₂ = 0.9,
                        η = 0.05,
                        ϵ = 1e-8,
                        total_num_steps = 500,
                        averaging_num_steps=100,
                        λ = 1e6)

log_ct_fit,mean_log_LL = ADAM_optim(log.(cts₀[25:end]),
                                        random_grad_log_ct!,
                                        baseline_adam_params,
                                        nai_one_group,
                                        baseline_ct,
                                        trans_one_group,
                                        0.0)
