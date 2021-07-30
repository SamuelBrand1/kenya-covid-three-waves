## This script gives a sensitivity analysis of models that fit to Nairobi data


## Active and precompile the KenyaCoVSD environment
using Revise,Pkg
Pkg.activate(".");Pkg.precompile()
using StatsPlots,Dates,JLD2,Statistics,Optim,Parameters,Distributions,DataFrames,CSV
using Plots.PlotMeasures,OrdinaryDiffEq,DiffEqCallbacks,NamedArrays
using CountTimeSeries
import KenyaCoVSD

##

#Load data for all Kenya

@load("data/linelist_data_with_pos_neg_20feb_to_27apr_c_age.jld2")
@load("data/deaths_data20210521.jld2")
@load("data/serological_data_with_20feb_to_31dec_age.jld2")
@load("data/serological_data_with_20feb_to_10Mar2021_age.jld2")
@load("data/N_kenya.jld2")

#Create a `CovAreaModel` object for Kenya which we will mutate with alternative models explaining the data

basic_nai_model = KenyaCoVSD.CoVAreaModel("Nairobi",θ -> 0.,θ -> 0.0;
                                            case_data = linelist_data_with_pos_neg,
                                            sero_data = serological_data,
                                            death_data = deaths_data,
                                            pop_data = N_kenya)

#`basic_nai_model` contains the PCR, serological and death data for Nairobi.

scatter(basic_nai_model.PCR_cases[:,1])

##Create

## Model 1: Simple INGARCH model
#Can get find the best AIC model here, but what about lookahead predictions?
#Find the best (p,1) INGARCH model using AIC criterion
y = basic_nai_model.PCR_cases[1:410,1]

function AIC_by_p(p)
        res = fit(y, Model(pastObs = 1:p, pastMean = 1,
                        distr = "NegativeBinomial",link = "Log"))
        BIC(res,55)
end

AICs = [AIC_by_p(p) for p = 1:10]
scatter(AICs)
findmin(AICs)
model = Model(pastObs = 1, pastMean = 1,
                distr = "NegativeBinomial",link = "Log")
# y = simulate(1000, model, [10, 0.5, 0.2])[1]
res = fit(y, model)
p = predict(res, 100,1000)
plot(basic_nai_model.PCR_cases[1:410,1],lab = "All data",legend = :topleft)
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

## Model 2: One group model
