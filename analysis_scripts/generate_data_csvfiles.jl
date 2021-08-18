using MCMCChains: promote_eltype_tuple_type
using Dates,JLD2,Statistics,Optim,Parameters,Distributions,DataFrames,CSV,MCMCChains
import KenyaCoVSD
include("fitting_methods.jl");
@load("forecasts/condensed_county_forecasts.jld2")
## Data S1
## Get all parameter names
modelfitfilenames = readdir("modelfits",join = true)
modeldict = load(modelfitfilenames[1])
model = modeldict[first(keys(modeldict))]
chn = model.MCMC_results.chain
parameternames = names(chn)

##Fill in a matrix of strings with the posterior mean and 95% CIs for each parameter
q = quantile(chn)
m = mean(chn)

parameterfits = fill("",47,length(parameternames))
mean_parameterfits = fill(0.0,47,length(parameternames))

for (i,filename) in enumerate(modelfitfilenames)
    modeldict = load(filename)
    model = modeldict[first(keys(modeldict))]
    chn = model.MCMC_results.chain
    parameternames = names(chn)
    q = quantile(chn)
    m = mean(chn)
    for j = 1:length(parameternames)
        str = "$(round(m.nt.mean[j],digits = 3)) ($(round(q.nt[2][j],digits = 3)),$(round(q.nt[end][j],digits = 3)))"
        parameterfits[i,j] = str
        mean_parameterfits[i,j] = m.nt.mean[j]
    end
end

##Fill in a matrix of IFR AND detection %s from the condensed fitting information
fatality_detection_percent = zeros(47,2)
for i = 1:47,j = 1:2
    fit = condensed_county_forecasts[i]
    fatality_detection_percent[i,j] = round(fit.pred.μ_fit[j]*100,sigdigits = 3)
end

##Get the population exposure estimates
@load("forecasts/county_pop_exposure_estimates_firstjune2021.jld2")
population_exposure_fits = fill("",47)
for i = 1:47
    fit = population_exposure_estimates[i]
    mean_est = round(fit.pred_pop_exposure*100,digits = 1)
    lest = round((fit.pred_pop_exposure - fit.lerr_pop_exposure)*100,digits = 1)
    uest = round((fit.pred_pop_exposure + fit.uerr_pop_exposure)*100,digits = 1)
    str = "$(mean_est) ($(lest),$(uest))"
    population_exposure_fits[i] = str
end

##Combine fits into one document
countynames = [fit.name for fit in condensed_county_forecasts]
data = [countynames parameterfits string.(fatality_detection_percent) population_exposure_fits]
_data = [countynames mean_parameterfits]
colnames = [:county_name;parameternames;[:fatality_detection_lower_SES,:fatality_detection_upper_SES,:population_exposure]]
_colnames = [:county_name;parameternames]
dataS1 = DataFrame()
parameter_post_means = DataFrame()
for (j,colname) in enumerate(colnames)
    dataS1[!,string(colname)] = data[:,j]
end

for (j,colname) in enumerate(_colnames)
    parameter_post_means[!,string(colname)] = _data[:,j]
end

CSV.write("opendatacsvs/dataS1.csv",dataS1)
CSV.write("forecasts/parameter_posterior_means.csv",parameter_post_means)

## Data S2, S3 and S4

dataS2 = DataFrame()
dataS3 = DataFrame()
dataS4 = DataFrame()

for (i,filename) in enumerate(modelfitfilenames)
    modeldict = load(filename)
    model = modeldict[first(keys(modeldict))]
    if i == 1
        dataS2[!,"Date"] = [Date(2020,2,20) + Day(k) for k = 1:size(model.PCR_cases,1)]
        dataS3[!,"Date"] = [Date(2020,2,20) + Day(k) for k = 1:size(model.sero_cases,1)]
        dataS4[!,"Date"] = [Date(2020,2,20) + Day(k) for k = 1:size(model.deaths,1)]
    end
    dataS2[!,model.areaname*"_PCR_pos"] = model.PCR_cases[:,1]
    dataS2[!,model.areaname*"_PCR_neg"] = model.PCR_cases[:,2]
    dataS3[!,model.areaname*"_sero_pos"] = model.sero_cases[:,1,1]
    dataS3[!,model.areaname*"_sero_neg"] = model.sero_cases[:,2,1]
    dataS4[!,model.areaname*"_deaths_pos"] = model.deaths
end
CSV.write("opendatacsvs/dataS2.csv",dataS2)
CSV.write("opendatacsvs/dataS3.csv",dataS3)
CSV.write("opendatacsvs/dataS4.csv",dataS4)
