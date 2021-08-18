using StatsPlots,Dates,JLD2,Statistics,Optim,Parameters,Distributions,DataFrames,CSV
using Plots.PlotMeasures,OrdinaryDiffEq,DiffEqCallbacks
import KenyaCoVSD
include("fitting_methods.jl");
include("plotting_methods.jl");
## Load the collated fits
@load("data/N_kenya.jld2")
@load("data/linelist_data_with_pos_neg_20feb_to_27apr_c_age.jld2")
@load("data/serological_data_with_20feb_to_31dec_age.jld2")
@load("data/serological_data_with_20feb_to_10Mar2021_age.jld2")
@load("data/cleaned_linelist20210521_deaths_c__date_of_lab_confirmation.jld2")
@load("data/p_ID.jld2")
@load("forecasts/condensed_county_forecasts.jld2")
@load("data/rel_sero_detection_after_infection.jld2")
@load("data/cleaned_linelist20210521_c.jld2")#<--- Positive tests only looking ahead of linelist we fitted to

## Loop over fitted counties to get county specific plots of incidence rates, population exposure, deaths and positive case rates

for fit in condensed_county_forecasts
    inc_plt = plot_group_incidence(fit);
    sero_plt = plot_pop_exposure(fit,serological_data,serology_data,N_kenya);
    deaths_plt = plot_deaths(fit,linelist_data_with_pos_neg);        
    PCR_plt = plot_PCR(fit,linelist_data_with_pos_neg,linelist_data);
    savefig(inc_plt,"plots/county_plots/infection_rates/$(fit.name)_infection_rate.png")
    savefig(sero_plt,"plots/county_plots/pop_exposure/$(fit.name)_pop_exposure.png")
    savefig(deaths_plt,"plots/county_plots/deaths_rates/$(fit.name)_deaths_rate.png")
    savefig(PCR_plt,"plots/county_plots/PCR_positives/$(fit.name)_pos_cases_rate.png")
end

##Loop over county fits to get R(t) plots for every county

fitfiles = readdir("modelfits",join=true)


for filename in fitfiles
    fit_dict = load(filename)
    model = fit_dict[first(keys(fit_dict))]
    plt_Rt = plot_Rt_both_SES_groups(model,Date(2021,6,1))
    savefig(plt_Rt,"plots/county_plots/Rt/$(model.areaname)_Rt.png")
end