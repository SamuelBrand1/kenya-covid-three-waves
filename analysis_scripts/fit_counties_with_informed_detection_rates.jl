## Do inference in parallel

using Distributed

## add number of processors desired

addprocs(5) 

## Load relevant code into scope for each worker

@everywhere begin
    using Pkg;Pkg.activate(".")
end
@everywhere using JLD2,LinearAlgebra,TransformVariables,LogDensityProblems,MCMCChains,Dates
@everywhere import KenyaCoVSD

## Load data

@load("data/linelist_data_with_pos_neg_20feb_to_27apr_c_age.jld2")
@load("data/serological_data_with_20feb_to_31dec_age.jld2")
@load("data/deaths_data20210521.jld2")
@load("data/N_kenya.jld2")

## Parameter transformations and initial KE guess for HMC

trans = as((ct_min1 = as(Real, 0.4, 0.999),
        R₀ = as(Real, 0.0, 7.0),
        ϵ = as(Real, 0.7, 0.999),
        χ₁=as(Real, 0.1, 30.0),χ₂=as(Real, 0.1, 30.0),
        p_test₁=as(Real, 5e-6, 5e-3),p_test₂ = as(Real, 5e-6, 5e-3),
        p_choose1=as(Real, 0.1, 0.999),
        P_eff=as(Real, 0.5, 0.999),
        schooleffect = as(Real, 0.5, 0.999),
        extra_transmissibility = as(Real, 0.5, 2.5),
        influx_exposed_new_variant = as(Real, 0.0, 10e3),
        p_test_boost = as(Real, 1.0, 3.0),
        χ_boost= as(Real, 1.0, 3.0),
        E₀=as(Real, 0.0, 10e3)))

D = Diagonal(0.1*ones(15))

## Name for each county with informed priors for detection rates

semiurbancounties_lowserodata = ["Kajiado", "Kiambu", "Machakos","Isiolo" , "Kirinyaga", "Murang'a",
                                "Trans Nzoia","Nyandarua","Meru","Bungoma","Nyamira","Kakamega",
                                "Laikipia", "Taita Taveta", "Makueni","Elgeyo-Marakwet",
                                "Tharaka-Nithi", "Kericho","Bomet", "Homa Bay", "Busia","Migori",
                                "Nandi","Vihiga"]

ruralcounties_lowserodata = [ "Garissa", "Kitui",
                            "Lamu", "Marsabit","Narok", "Tana River", "Turkana", "Wajir",
                            "Baringo", "Mandera","Samburu", "West Pokot"]

counties_with_informed_priors = [semiurbancounties_lowserodata;ruralcounties_lowserodata]

priors =[fill(KenyaCoVSD.priors_twogroup_newvariant_semi_urbanrural_fitted_detection,length(semiurbancounties_lowserodata));
        fill(KenyaCoVSD.priors_twogroup_newvariant_rural_fitted_detection,length(ruralcounties_lowserodata))]

## parallelised loop for fitting each county                             
@distributed for i = 1:36
    name = counties_with_informed_priors[i]
    prior = priors[i]    
    model = KenyaCoVSD.CoVAreaModel(name,KenyaCoVSD.ll_twogroup_newvariant_infboost,prior;
                case_data = linelist_data_with_pos_neg,
                sero_data = serological_data,
                death_data = deaths_data,
                pop_data = N_kenya)
    
    KenyaCoVSD.inferparameters!(model,2000,trans,0.05,D;serowaningrate = 0.0)

    @save("modelfits/$(name)_model.jld2",model)
end


