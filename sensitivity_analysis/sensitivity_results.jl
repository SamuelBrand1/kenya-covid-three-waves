using StatsPlots,Dates,JLD2,Statistics,Optim,Parameters,Distributions,DataFrames,CSV
using DynamicHMC,OrdinaryDiffEq,DiffEqCallbacks
using NamedArrays,TransformVariables,LogDensityProblems,LinearAlgebra,BlackBoxOptim,Random
using DynamicHMC.Diagnostics,MCMCDiagnostics,MCMCChains,ForwardDiff,LogExpFunctions
using Plots.PlotMeasures,LogExpFunctions
using CountTimeSeries
import KenyaCoVSD
include("inference_variable_transmissibility.jl");
## Gather the fits
waning_immunity_scenarios = [(σ = 0.0,ω = 1/180,scenario_name = "_no_waning_immunity"),
              (σ = 0.5,ω = 1/180,scenario_name = "_med_susceptibility"),
              (σ = 1.0,ω = 1/180,scenario_name = "_full_susceptibility"),
              (σ = 0.16,ω = 1/90,scenario_name = "_fast_reversion"),
              (σ = 0.16,ω = 1/365,scenario_name = "_slow_reversion"),
              (σ = 1.0,ω = 1/(5*365),scenario_name = "_slow_reversion_to_full_susceptibility")]
infectiousness_scenarios = [(σ = 0.16,ω = 1/180,ι = 0.5, scenario_name = "reduced_inf"),
            (σ = 0.5,ω = 1/180,ι = 0.5,scenario_name = "_med_susceptibility_reduced_inf"),
            (σ = 1.0,ω = 1/180,ι = 0.5,scenario_name = "_full_susceptibility_reduced_inf"),
            (σ = 1.0,ω = 1/90,ι = 0.5,scenario_name = "_full_susceptibility_fast_wane_reduced_inf")]

##
all_scenarios = []

@load("sensitivity_analysis/nairobi_fits/Nairobi_model_no_waning_immunity.jld2")
push!(all_scenarios,(σ = 0.0,ω = 1/180,ι = 1.0,scenario_name = "_no_waning_immunity",fit= model,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_model_med_susceptibility.jld2")
push!(all_scenarios,(σ = 0.5,ω = 1/180,ι = 1.0,scenario_name = "_med_susceptibility",fit= model,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_model_full_susceptibility.jld2")
push!(all_scenarios,(σ = 1.0,ω = 1/180,ι=1.0,scenario_name = "_full_susceptibility",fit= model,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_model_fast_reversion.jld2")
push!(all_scenarios,(σ = 0.16,ω = 1/90,ι = 1.0, scenario_name = "_fast_reversion",fit= model,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_model_slow_reversion.jld2")
push!(all_scenarios,(σ = 0.16,ω = 1/365,ι = 1.0, scenario_name = "_slow_reversion",fit= model,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_model_slow_reversion_to_full_susceptibility.jld2")
push!(all_scenarios,(σ = 1.0,ω = 1/(5*365),ι = 1.0, scenario_name = "_slow_reversion_to_full_susceptibility",fit= model,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_modelreduced_inf.jld2")
push!(all_scenarios,(σ = 0.16,ω = 1/180,ι = 0.5, scenario_name = "reduced_inf",fit= model,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_model_med_susceptibility_reduced_inf.jld2")
push!(all_scenarios,(σ = 0.5,ω = 1/180,ι = 0.5,scenario_name = "_med_susceptibility_reduced_inf",fit= model,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_model_full_susceptibility_reduced_inf.jld2")
push!(all_scenarios,(σ = 1.0,ω = 1/180,ι = 0.5,scenario_name = "_full_susceptibility_reduced_inf",fit= model,DIC = KenyaCoVSD.modeldic(model)))
@load("sensitivity_analysis/nairobi_fits/Nairobi_model_full_susceptibility_fast_wane_reduced_inf.jld2")
push!(all_scenarios,(σ = 1.0,ω = 1/90,ι = 0.5,scenario_name = "_full_susceptibility_fast_wane_reduced_inf",fit= model,DIC = KenyaCoVSD.modeldic(model)))
