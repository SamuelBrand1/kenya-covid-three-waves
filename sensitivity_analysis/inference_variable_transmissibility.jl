"""
    function ll_twogroup_newvariant_infboost_var_trans(θ,model::KenyaCoVSD.CoVAreaModel,seroreversionrate;σ = 0.16,ω = 1/180,ct_min2 = 0.445)

Log-likelihood of the case and serology data in `model` assuming that new variant is more transmissible rather than immune evading.
"""
function ll_twogroup_newvariant_infboost_var_trans(θ,model::KenyaCoVSD.CoVAreaModel,seroreversionrate;σ = 0.16,ω = 1/180,ct_min2 = 0.445,ι=1.0)
    @unpack PCR_cases,sero_cases,baseline_sero_array,PCR_array,sero_sensitivity,sero_specificity,N,M_BB,prob,α,γ,relative_testing_rate = model
    @unpack ct_min1,R₀,ϵ,χ₁,χ₂,p_test₁,p_test₂,p_choose1,P_eff,schooleffect,extra_transmissibility,influx_exposed_new_variant,p_test_boost,χ_boost,E₀ = θ

    #Set variance scalers
    clustering_factor_PCR = 0.5
    M_PCR = 30.
    T = eltype(R₀)

    #Set transmission parameters
	N₁ = P_eff*N
	N₂ = N - N₁
    p = convert.(T,[R₀,α,γ,ϵ,σ,N₁,N₂,ω,schooleffect,ct_min1,ct_min2,ι])
    u0 = convert.(T,[N₁,E₀,0.,0.,0.,0.,
                     N₂,E₀,0.,0.,0.,0.,
                     0.0,0.0,
		     0.0,0.0,0.0,0.0])

    #Create new variant introduction
    function new_variant_effect!(integrator)
        integrator.p[1] *= extra_transmissibility
        integrator.u[2] += influx_exposed_new_variant #Lower SES
        integrator.u[8] += influx_exposed_new_variant #Higher SES
    end

    janendpoint = (Date(2021,1,30) - Date(2020,2,20)).value
	aprilendpoint = (Date(2021,4,30) - Date(2020,2,20)).value
    variant_cb = PresetTimeCallback([janendpoint],new_variant_effect!)

    #Sero-waning
    sero_array = vcat(baseline_sero_array[1:30],[(1-seroreversionrate)^k for k in 1:length(baseline_sero_array[31:end])])


    #Solve for daily incidence by age (not currently broken down by severity of disease)
    LL = T(0.)
    try
        sol = solve(prob, BS3();tspan = (0,aprilendpoint),reltol = 1e-3,
				callback = variant_cb,
				u0=u0,
				p=p,
				saveat = 1,
				verbose = true,
			        isoutofdomain=(u,p,t) -> any(x -> x < 0, u))

        ι₁ = diff(sol[:C₁])
        ι₂ = diff(sol[:C₂])
        ι_sero₁ = diff(sol[:C_sero₁])
        ι_sero₂ = diff(sol[:C_sero₂])
        PCR₁ = KenyaCoVSD.simple_conv(ι₁,PCR_array)
        PCR₂ = KenyaCoVSD.simple_conv(ι₂,PCR_array)
    	sero₁ = sero_sensitivity.*KenyaCoVSD.simple_conv(ι_sero₁,sero_array)
	    sero₂ = sero_sensitivity.*KenyaCoVSD.simple_conv(ι_sero₂,sero_array)

        #Calculate log-likelihood for PCR testing
        for t in 55:janendpoint
            #Convert from μ,α parameterisation to p,r parameterisation for using Neg. binomial model
            μ = relative_testing_rate[t]*(p_test₁*PCR₁[t] + p_test₂*PCR₂[t]) + 0.001
            σ² = μ + clustering_factor_PCR*μ^2
            p_negbin = 1 - (clustering_factor_PCR*μ^2/σ²)
            r_negbin  = 1/clustering_factor_PCR
            LL += logpdf(NegativeBinomial(r_negbin,p_negbin),PCR_cases[t,1])#likelihood contribution from PCR testing --- positive case detection

            if PCR_cases[t,2] >= 0 #Negative tests available
                p_PCR_pred = p_choose1*(χ₁*PCR₁[t]/((χ₁-1)*PCR₁[t] + sum(N₁))) + (1-p_choose1)*(χ₂*PCR₂[t]/((χ₂-1)*PCR₂[t] + sum(N₂))) + 0.001
                #Convert from p_PCR,inv_M_PCR parameterisation of Beta-Binomial to standard α,β parameterisation
                LL += logpdf(BetaBinomial(PCR_cases[t,2],p_PCR_pred*M_PCR,(1-p_PCR_pred)*M_PCR),PCR_cases[t,1])#likelihood contribution from PCR testing --- proportion postive
            end
        end

        #Calculate log-likelihood for PCR testing after new variant introduction
        for t in (janendpoint+1):min(size(PCR_cases,1),size(PCR₁,1))
            #Convert from μ,α parameterisation to p,r parameterisation for using Neg. binomial model
            μ = p_test_boost*relative_testing_rate[t]*(p_test₁*PCR₁[t] + p_test₂*PCR₂[t]) + 0.001
            σ² = μ + clustering_factor_PCR*μ^2
            p_negbin = 1 - (clustering_factor_PCR*μ^2/σ²)
            r_negbin  = 1/clustering_factor_PCR
            LL += logpdf(NegativeBinomial(r_negbin,p_negbin),PCR_cases[t,1])#likelihood contribution from PCR testing --- positive case detection

            if PCR_cases[t,2] >= 0 #Negative tests available
                p_PCR_pred = p_choose1*(χ_boost*χ₁*PCR₁[t]/((χ_boost*χ₁-1)*PCR₁[t] + sum(N₁))) + (1-p_choose1)*(χ_boost*χ₂*PCR₂[t]/((χ_boost*χ₂-1)*PCR₂[t] + sum(N₂))) + 0.001
                #Convert from p_PCR,inv_M_PCR parameterisation of Beta-Binomial to standard α,β parameterisation
                LL += logpdf(BetaBinomial(PCR_cases[t,2],p_PCR_pred*M_PCR,(1-p_PCR_pred)*M_PCR),PCR_cases[t,1])#likelihood contribution from PCR testing that day in that age group
            end
        end
        #Calculate log-likelihood contribution from serological testing
        for t in 1:min(size(sero_cases,1),size(sero₁,1))
            #Convert from the p_hat,M_BB parameterisation of the Betabinomial distribution to the α_BB,β_BB parameterisation

            #Overall serological testing
            if sero_cases[t,1,1] + sero_cases[t,2,1] > 0
                p_sero_pred  = (sero₁[t] + sero₂[t])./(N₁ + N₂)
                p_hat  = p_sero_pred + (1-p_sero_pred)*(1-sero_specificity)
                LL += logpdf(BetaBinomial(sero_cases[t,1,1] + sero_cases[t,2,1],M_BB*p_hat,M_BB*(1-p_hat)),sero_cases[t,1,1])#Likelihood contribution from sero testing
            end
            #Lower SES serological testing
            if sero_cases[t,1,2] + sero_cases[t,2,2] > 0
                p_sero_pred  = sero₁[t]./N₁
                p_hat  = p_sero_pred + (1-p_sero_pred)*(1-sero_specificity)
                LL += logpdf(BetaBinomial(sero_cases[t,1,2] + sero_cases[t,2,2],M_BB*p_hat,M_BB*(1-p_hat)),sero_cases[t,1,2])#Likelihood contribution from sero testing
            end
            #Higher SES serological testing
            if sero_cases[t,1,3] + sero_cases[t,2,3] > 0
                p_sero_pred  = sero₂[t]./N₂
                p_hat  = p_sero_pred + (1-p_sero_pred)*(1-sero_specificity)
                LL += logpdf(BetaBinomial(sero_cases[t,1,3] + sero_cases[t,2,3],M_BB*p_hat,M_BB*(1-p_hat)),sero_cases[t,1,3])#Likelihood contribution from sero testing
            end

        end
    catch
        LL = T(-Inf)
    end
    return LL::T
end


"""
        function inferparameters!(model,samples,trans,stepsize::Float64,D::Diagonal,q₀;serowaningrate = 1/180,num_chains = 1)

Infer the unknown transmission and observation parameters for the `k`th county in by `country_model`. Inference is performed by drawing `samples` number of replicates from the posterior distribution with log-density
    given in the model with sero-waning rate fixed. The parameters and their transformation are encoded in `trans`. Initial parameter search starts at `q₀` in the transformed domain defined by `trans`.
    whilst the initial kinetic energy search begins at `D` and the HMC step size is fixed to be `stepsize`. Inference output is stored in-place.

    Serowaning rate can be added, as well as additional chains for cross-comparison and validation.

"""
function inferparameters!(model,samples,trans,stepsize::Float64,D::Diagonal,q₀;serowaningrate = 1/365,num_chains = 1,σ = 0.16,ω = 1/180,ι = 1.0)
        println("Starting MCMC parameter inference for $(model.areaname)")
        n = length(q₀)
        l_area(x) = model.log_likelihood(x,model,serowaningrate;σ = σ,ω = ω,ι=ι) + model.log_priors(x)
        ℓ = TransformedLogDensity(trans, l_area)#transformed log-likelihood
        ∇ℓ = LogDensityProblems.ADgradient(:ForwardDiff, ℓ)#transformed log-likelihood gradient wrt the parameters

        results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇ℓ, samples,
                                initialization = (q = q₀,κ=GaussianKineticEnergy(D),ϵ = stepsize),
                                warmup_stages = fixed_stepsize_warmup_stages(M=Symmetric),
                                reporter = NoProgressReport())
        # results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇ℓ, samples,
        #                         initialization = (q = q₀,κ=GaussianKineticEnergy(D)))

        transformed_results = TransformVariables.transform.(trans,results.chain)
        val = zeros(length(transformed_results),length(transformed_results[1]),1)

        for i = 1:size(val,1),j = 1:size(val,2)
                val[i,j,1] = transformed_results[i][j]
        end

        for chain_reps = 2:num_chains
                println("Running additional chain $(chain_reps) for $(model.areaname)")
                results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇ℓ, samples,
                                        initialization = (q = q₀,κ=GaussianKineticEnergy(D),ϵ = stepsize),
                                        warmup_stages = fixed_stepsize_warmup_stages(M=Symmetric))#,
                                        # reporter = NoProgressReport())
                transformed_results = TransformVariables.transform.(trans,results.chain)
                for i = 1:size(val,1),j = 1:size(val,2)
                        val[i,j,chain_reps] = transformed_results[i][j]
                end
        end

        chn = Chains(val,[String(k) for k in keys(transformed_results[1])])

        model.MCMC_results = MCMCResults(chn,
                                        [l_area(transformed_results[i]) - model.log_priors(transformed_results[i]) for i = 1:length(transformed_results)],
                                        results.tree_statistics)

        return nothing
end


"""
        function inferparameters!(model,samples,trans,stepsize::Float64,D::Diagonal,searchrange;serowaningrate = 1/180,num_chains = 1)

Infer the unknown transmission and observation parameters for the `k`th county in by `country_model`. Inference is performed by drawing `samples` number of replicates from the posterior distribution with log-density
    given in the model with sero-waning rate fixed. The parameters and their transformation are encoded in `trans`. Initial parameter search has a pre-optimisation using the default adaptive differential evolution algorithm implemented by the `BlackBoxOptim` package,
    whilst the initial kinetic energy search begins at `D` and the HMC step size is fixed to be `stepsize`. Inference output is stored in-place.

    serowaning rate can be added, as well as additional chains for cross-comparison and validation.
"""
function inferparameters!(model,samples,trans,stepsize::Float64,D::Diagonal;serowaningrate = 1/365,num_chains = 1,σ = 0.16,ω = 1/180,ι = 1.0)
        println("Searching for a good initial condition for $(model.areaname)")
        l_area(x) = model.log_likelihood(x,model,serowaningrate;σ = σ,ω = ω,ι=ι) + model.log_priors(x)
        f(x) = -transform_logdensity(trans, l_area,x)
        searchrange = fill((-3.,3.),TransformVariables.dimension(trans))
        res = bboptimize(f; SearchRange = searchrange,PopulationSize=1000,MaxSteps=30000)#,TraceMode = :silent)
        q₀ = best_candidate(res)
        inferparameters!(model,samples,trans,stepsize,D,q₀;serowaningrate=serowaningrate,num_chains=num_chains,σ = σ,ω = ω,ι=ι)
        return nothing
end
