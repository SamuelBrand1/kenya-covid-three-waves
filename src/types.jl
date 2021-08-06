"""
struct MCMCResults

Struct for holding the MCMC chain (a vector of NamedTuples which represent the accepted parameter sets), and the MCMC tree statistics
"""
struct MCMCResults
    chain
    logL::Vector{Float64}
    treestatistics::Vector{DynamicHMC.TreeStatisticsNUTS}
end


"""
    mutable struct CoVAreaModel

Struct for holding the fixed data of the area and MCMC results.
"""
Base.@kwdef mutable struct CoVAreaModel
    areaname::String #Name of area
    PCR_cases::Matrix{Int64} = zeros(Int64,365,2) #Daily PCR cases: dim 1: day, dim 2: num pos and total
    sero_cases::Array{Int64,3} = zeros(Int64,365,2,3) #Daily serology tests: dim 1: day, dim 2: num pos and neg, dim 3: 1 = Overall, 2 = lower SES sero tests, 3 = higher SES sero tests
    deaths::Vector{Int64} = zeros(Int64,365) #Daily deaths
    dates::Vector{Date} = [Date(2020,2,20) + Day(k) for k = 1:365]#Dates of each day
    N::Float64 = 4.3e6 #Area's population size
    α::Float64 = 1/3.1 # 1/mean latent period
    γ::Float64 = 1/2.4#Recovery rate
    relative_testing_rate::Vector{Float64} = relative_testing_rate
    prob::ODEProblem = build_prob()
    PCR_array::Vector{Float64} = PCR_array #Relative Sensitivity of PCR test on days post-infection
    PCR_sensitivity::Float64 = 1. #Base sensitivity for RT-PCR test
    PCR_specificity::Float64 = 0.995 #Base specificity for RT_PCR test
    baseline_sero_array::Vector{Float64} = rel_sero_array_26days #Relative sensitivity of serology test on days post-infection
    sero_sensitivity::Float64 = 0.825 #Base sensitivity for serology test
    sero_specificity::Float64 = 0.992 #Base specificity for serology test
    M_BB::Float64 = 40.74 #Shrinkage factor for BetaBinomial distribution that incorporates the uncertainty in the sensitivity (fitted to serology group)
    log_likelihood::Function = θ -> 0.
    log_priors::Function = θ -> 0.
    MCMC_results::Union{MCMCResults,Nothing} = nothing #This field gets added after the MCMC run
end

"""
    struct CoVCountryModel

Struct for holding a collection of `CovAreaModel` objects, fitting to them simultaneously, and a sero-waning rate that covers all of them.
"""
Base.@kwdef struct CoVCountryModel
    name::String = "Kenya"
    areas::Vector{CoVAreaModel}
    p_serowaning::Float64 = 1/180.
    log_priors::Function = θ -> 0.
end

"""
    CoVAreaModel(name::String,loglk_func,log_prior_func;
                        case_data,sero_data,death_data,pop_data::NamedArray)

Constructor for a `CoVAreaModel` struct using country data objects.
"""
function CoVAreaModel(name::String,loglk_func,log_prior_func;
                        case_data,sero_data,death_data,pop_data::NamedArray)
        uprname = uppercase(name)
        PCR_cases = sum(case_data.cases[:,case_data.areas .== uprname,:,:],dims = [2,3])[:,1,1,:]
        sero_cases = sum(sero_data.serodata[:,sero_data.areas .== uprname,:,:],dims = [2,3])[:,1,1,:]
        sero_cases = cat(sero_cases,zeros(Int64,size(sero_cases)),zeros(Int64,size(sero_cases)),dims = 3)
        deaths = death_data.deaths[:,death_data.areas .== uprname][:]
        dates = [Date(2020,2,20) + Day(k) for k = 1:size(PCR_cases,1)]
        N = sum(pop_data[:,name])

        return KenyaCoVSD.CoVAreaModel(areaname = name,
                                        PCR_cases=PCR_cases,
                                        sero_cases=sero_cases,
                                        deaths=deaths,
                                        dates = dates,
                                        N=N,
                                        log_likelihood=loglk_func,
                                        log_priors = log_prior_func)

end

function (model::KenyaCoVSD.CoVAreaModel)(θ)
    model.log_likelihood(θ,model,1/180) + model.log_priors(θ)
end

"""
    function modeldic(areamodel::CoVAreaModel)

Calculate the deviance information criterion (DIC) for the posterior parameter draws in `areamodel`.
"""
function modeldic(areamodel::CoVAreaModel)
        D = -2*areamodel.MCMC_results.logL
        D̄ = mean(D)
        p_D = 0.5*var(D)
        return D̄ + p_D
end
