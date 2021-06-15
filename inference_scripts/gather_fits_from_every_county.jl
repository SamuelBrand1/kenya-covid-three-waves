using Dates: include
using OrdinaryDiffEq,StatsPlots,Dates,Parameters,DiffEqCallbacks,Distributions,NamedArrays
using JLD2,DataFrames,CSV,LinearAlgebra,TransformVariables,LogDensityProblems,MCMCChains,Optim
import KenyaCoVSD

##Load data, delay distribution data and methods
@load("data/p_ID.jld2")
include("fitting_methods.jl")
include("plotting_methods.jl")
@load("modelfits/Nairobi_model.jld2")
relative_testing_rate = model.relative_testing_rate
@load("data/linelist_data_with_pos_neg_20feb_to_27apr_c_age.jld2")
@load("data/cleaned_linelist20210521_deaths_c__date_of_lab_confirmation.jld2")

##Loop over county fits and condense into one data file

fitfiles = readdir("modelfits",join=true)

condensed_county_forecasts = []

for filename in fitfiles
    fit_dict = load(filename)
    model = fit_dict[first(keys(fit_dict))]
    name = model.areaname
    println("Analysing and forecasting for $(model.areaname)")
    pred = get_predictions(model,Date(2021,8,1),p_ID)
    println("Condensing fit for county $(name)")
    uprname = uppercase(name)
    tests = sum(linelist_data_with_pos_neg.cases[:,linelist_data_with_pos_neg.areas .== uprname,:,2],dims = 3)[:,1,1,1]
    deaths = deaths_data.deaths[:,deaths_data.areas .== uprname][1:405]
    condensedpred = condense_prediction(pred,
                            tests,deaths,p_ID,relative_testing_rate)
    push!(condensed_county_forecasts,(name = name,pred=condensedpred))                            
end

#Save
@save("forecasts/condensed_county_forecasts.jld2",condensed_county_forecasts)

