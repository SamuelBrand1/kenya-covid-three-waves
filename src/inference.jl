"""
        function inferparameters!(model,samples,trans,stepsize::Float64,D::Diagonal,q₀;serowaningrate = 1/180,num_chains = 1)

Infer the unknown transmission and observation parameters for the `k`th county in by `country_model`. Inference is performed by drawing `samples` number of replicates from the posterior distribution with log-density
    given in the model with sero-waning rate fixed. The parameters and their transformation are encoded in `trans`. Initial parameter search starts at `q₀` in the transformed domain defined by `trans`.
    whilst the initial kinetic energy search begins at `D` and the HMC step size is fixed to be `stepsize`. Inference output is stored in-place.

    Serowaning rate can be added, as well as additional chains for cross-comparison and validation.

"""
function inferparameters!(model,samples,trans,stepsize::Float64,D::Diagonal,q₀;serowaningrate = 1/365,num_chains = 1)
        println("Starting MCMC parameter inference for $(model.areaname)")
        n = length(q₀)
        l_area(x) = model.log_likelihood(x,model,serowaningrate) + model.log_priors(x)
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
                                        warmup_stages = fixed_stepsize_warmup_stages(M=Symmetric),
                                        reporter = NoProgressReport())
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
function inferparameters!(model,samples,trans,stepsize::Float64,D::Diagonal;serowaningrate = 1/365,num_chains = 1)
        println("Searching for a good initial condition for $(model.areaname)")
        l_area(x) = model.log_likelihood(x,model,serowaningrate) + model.log_priors(x)
        f(x) = -transform_logdensity(trans, l_area,x)
        searchrange = fill((-3.,3.),TransformVariables.dimension(trans))
        res = bboptimize(f; SearchRange = searchrange,PopulationSize=1000,MaxSteps=30000,TraceMode = :silent)
        q₀ = best_candidate(res)
        inferparameters!(model,samples,trans,stepsize,D,q₀;serowaningrate=serowaningrate,num_chains=num_chains)
        return nothing
end
