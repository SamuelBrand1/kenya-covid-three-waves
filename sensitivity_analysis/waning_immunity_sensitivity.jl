## Do inference in parallel

using Distributed

## add number of processors desired

addprocs(6)

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

## Transformation of variables
D = Diagonal(0.1*ones(TransformVariables.dimension(trans)))

## Define scenarios 
scenarios = [(σ = 0.0,ω = 1/180,scenario_name = "_no_waning_immunity"),
              (σ = 0.5,ω = 1/180,scenario_name = "_med_susceptibility"),
              (σ = 1.0,ω = 1/180,scenario_name = "_full_susceptibility"),
              (σ = 0.16,ω = 1/90,scenario_name = "_fast_reversion"),
              (σ = 0.16,ω = 1/365,scenario_name = "_slow_reversion"),
              (σ = 1.0,ω = 1/(5*365),scenario_name = "_slow_reversion_to_full_susceptibility")]


##Fit each scenario
@distributed for i = 1:length(scenarios)
    scenario = scenarios[i]
    model = KenyaCoVSD.CoVAreaModel("Nairobi",KenyaCoVSD.ll_twogroup_newvariant_infboost,KenyaCoVSD.priors_twogroup_newvariant_cities;
                case_data = linelist_data_with_pos_neg,
                sero_data = serological_data,
                death_data = deaths_data,
                pop_data = N_kenya)
    
    KenyaCoVSD.inferparameters!(model,2000,trans,0.05,D;serowaningrate = 0.0,σ = scenario.σ,ω = scenario.ω)

    @save("sensitivity_analysis/Nairobi_model$(scenario.scenario_name).jld2",model)
end



