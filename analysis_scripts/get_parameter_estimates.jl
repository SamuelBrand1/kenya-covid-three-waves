using Dates,JLD2,Statistics,Optim,Parameters,Distributions,DataFrames,CSV
import KenyaCoVSD
include("fitting_methods.jl");

## Load the collated fits
@load("data/N_kenya.jld2")
@load("data/linelist_data_with_pos_neg_20feb_to_27apr_c_age.jld2")
@load("data/serological_data_with_20feb_to_31dec_age.jld2")
@load("data/serological_data_with_20feb_to_10Mar2021_age.jld2")
@load("data/cleaned_linelist20210426_deaths_c__date_of_lab_confirmation.jld2")
@load("data/p_ID.jld2")
@load("data/smoothedPCR_detection_after_infection.jld2")


#schools effect

fits_schools_effect = []

for filename in readdir("modelfits",join = true)
    modeldict = load(filename)
    model = modeldict[first(keys(modeldict))]
    println("Gathering data from $(model.areaname)")
    params = get(model.MCMC_results.chain,:schooleffect)
    push!(fits_schools_effect,(name = model.areaname,
                                mean = mean(params.schooleffect[:]),
                                lwr = quantile(params.schooleffect[:],0.025),
                                upr = quantile(params.schooleffect[:],0.975)))
end

mean_school_effects = [1 .- fit.mean for fit in fits_schools_effect]
idx_nai = findfirst([fit.name .== "Nairobi" for fit in fits_schools_effect])
idx_mom = findfirst([fit.name .== "Mombasa" for fit in fits_schools_effect])
1 - fits_schools_effect[idx_nai].mean
1 - fits_schools_effect[idx_nai].upr
1 - fits_schools_effect[idx_nai].lwr
1 - fits_schools_effect[idx_mom].mean
1 - fits_schools_effect[idx_mom].upr
1 - fits_schools_effect[idx_mom].lwr

quantile(mean_school_effects,[0.25,0.5,0.75])

## Ct min
fits_ct_min1 = []

for filename in readdir("modelfits",join = true)
    modeldict = load(filename)
    model = modeldict[first(keys(modeldict))]
    println("Gathering data from $(model.areaname)")
    params = get(model.MCMC_results.chain,:ct_min1)
    push!(fits_ct_min1,(name = model.areaname,
                                mean = mean(params.ct_min1[:]),
                                lwr = quantile(params.ct_min1[:],0.025),
                                upr = quantile(params.ct_min1[:],0.975)))
end

mean_ct_min1= [1 - fit.mean for fit in fits_ct_min1]
quantile(mean_ct_min1,[0.25,0.5,0.75])

idx_nai = findfirst([fit.name .== "Nairobi" for fit in fits_ct_min1])
idx_mom = findfirst([fit.name .== "Mombasa" for fit in fits_ct_min1])
1 - fits_ct_min1[idx_nai].mean
1 - fits_ct_min1[idx_nai].upr
1 - fits_ct_min1[idx_nai].lwr
1 - fits_ct_min1[idx_mom].mean
1 - fits_ct_min1[idx_mom].upr
1 - fits_ct_min1[idx_mom].lwr

## P_eff
fits_P_eff = []

for filename in readdir("modelfits",join = true)
    modeldict = load(filename)
    model = modeldict[first(keys(modeldict))]
    println("Gathering data from $(model.areaname)")
    params = get(model.MCMC_results.chain,:P_eff)
    push!(fits_P_eff,(name = model.areaname,
                                mean = mean(params.P_eff[:]),
                                lwr = quantile(params.P_eff[:],0.025),
                                upr = quantile(params.P_eff[:],0.975)))
end

mean_P_eff = [1 - fit.mean for fit in fits_P_eff]
quantile(mean_P_eff,[0.25,0.5,0.75])

pop_weighted_peff = [(1 - fit.mean)*sum(N_kenya[:,fit.name])/sum(N_kenya) for fit in fits_P_eff]
sum(pop_weighted_peff)

idx_nai = findfirst([fit.name .== "Nairobi" for fit in fits_P_eff])
idx_mom = findfirst([fit.name .== "Mombasa" for fit in fits_P_eff])
1 - fits_P_eff[idx_nai].mean
1 - fits_P_eff[idx_nai].upr
1 - fits_P_eff[idx_nai].lwr
1 - fits_P_eff[idx_mom].mean
1 - fits_P_eff[idx_mom].upr
1 - fits_P_eff[idx_mom].lwr


## Rel detection rate
fits_rel_detect = []

for filename in readdir("modelfits",join = true)
    modeldict = load(filename)
    model = modeldict[first(keys(modeldict))]
    println("Gathering data from $(model.areaname)")
    params = get(model.MCMC_results.chain,[:p_test₁,:p_test₂])
    prob_detect_1 = [dailytestrate_to_probofdetections(p_test,PCR_array) for p_test in params.p_test₁[:]]
    prob_detect_2 = [dailytestrate_to_probofdetections(p_test,PCR_array) for p_test in params.p_test₂[:]]

    push!(fits_rel_detect,(name = model.areaname,
                                mean = mean(prob_detect_2./prob_detect_1),
                                lwr = quantile(prob_detect_2./prob_detect_1,0.025),
                                upr = quantile(prob_detect_2./prob_detect_1,0.975)))
end

idxs_P_eff = mean_P_eff .> 0.3
mean_rel_detect = [fit.mean for fit in fits_rel_detect]
fits_rel_detect[idxs_P_eff]

quantile(mean_rel_detect,[0.25,0.5,0.75])
idx_nai = findfirst([fit.name .== "Nairobi" for fit in fits_rel_detect])
idx_mom = findfirst([fit.name .== "Mombasa" for fit in fits_rel_detect])
fits_rel_detect[idx_nai].mean
fits_rel_detect[idx_nai].upr
fits_rel_detect[idx_nai].lwr
fits_rel_detect[idx_mom].mean
fits_rel_detect[idx_mom].upr
fits_rel_detect[idx_mom].lwr

## p_choose1
fits_p_choose1 = []

for filename in readdir("modelfits",join = true)
    modeldict = load(filename)
    model = modeldict[first(keys(modeldict))]
    println("Gathering data from $(model.areaname)")
    params = get(model.MCMC_results.chain,:p_choose1)
    push!(fits_p_choose1,(name = model.areaname,
                                mean = mean(params.p_choose1[:]),
                                lwr = quantile(params.p_choose1[:],0.025),
                                upr = quantile(params.p_choose1[:],0.975)))
end

mean_p_choose1 = [fit.mean for fit in fits_p_choose1]
quantile(mean_p_choose1,[0.25,0.5,0.75])
## p_choose1 bias
fits_p_choose2_bias = []

for filename in readdir("modelfits",join = true)
    modeldict = load(filename)
    model = modeldict[first(keys(modeldict))]
    println("Gathering data from $(model.areaname)")
    params = get(model.MCMC_results.chain,[:p_choose1,:P_eff])
    push!(fits_p_choose2_bias,(name = model.areaname,
                                mean = mean((1 .- params.p_choose1[:])./(1 .- params.P_eff[:])),
                                lwr = quantile((1 .- params.p_choose1[:])./(1 .- params.P_eff[:]),0.025),
                                upr = quantile((1 .- params.p_choose1[:])./(1 .- params.P_eff[:]),0.975)))
end

mean_p_choose2_bias = [fit.mean for fit in fits_p_choose2_bias]
quantile(mean_p_choose2_bias,[0.25,0.5,0.75])

idx_nai = findfirst([fit.name .== "Nairobi" for fit in fits_p_choose2_bias])
idx_mom = findfirst([fit.name .== "Mombasa" for fit in fits_p_choose2_bias])
fits_P_eff[idx_nai].mean
fits_p_choose2_bias[idx_nai].upr
fits_p_choose2_bias[idx_nai].lwr
fits_p_choose2_bias[idx_mom].mean
fits_p_choose2_bias[idx_mom].upr
fits_p_choose2_bias[idx_mom].lwr

## extra_transmissibility bias
fits_extra_transmissibility = []

for filename in readdir("modelfits",join = true)
    modeldict = load(filename)
    model = modeldict[first(keys(modeldict))]
    println("Gathering data from $(model.areaname)")
    params = get(model.MCMC_results.chain,:extra_transmissibility)
    push!(fits_extra_transmissibility,(name = model.areaname,
                                mean = mean(params.extra_transmissibility[:]),
                                lwr = quantile(params.extra_transmissibility[:],0.025),
                                upr = quantile(params.extra_transmissibility[:],0.975)))
end

mean_extra_transmissibility = [fit.mean - 1 for fit in fits_extra_transmissibility]
quantile(mean_extra_transmissibility,[0.25,0.5,0.75])
idx_nai = findfirst([fit.name .== "Nairobi" for fit in fits_extra_transmissibility])
idx_mom = findfirst([fit.name .== "Mombasa" for fit in fits_extra_transmissibility])
fits_extra_transmissibility[idx_nai].mean
fits_extra_transmissibility[idx_nai].upr
fits_extra_transmissibility[idx_nai].lwr
fits_extra_transmissibility[idx_mom].mean
fits_extra_transmissibility[idx_mom].upr
fits_extra_transmissibility[idx_mom].lwr