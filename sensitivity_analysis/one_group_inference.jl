#This file defines methods for using the EM algorithm to converge on a joint inference of basic parameters and the daily contact rate

#Define the initial 'guess' for the contact rate based on Google mobility data
cts₀ = [KenyaCoVSD.ct_kenya(t,0.55) for (t,d) in enumerate(Date(2020,2,21):Day(1):Date(2021,4,15))]

#Define variable transformation

trans_one_group = as((R₀ = as(Real, 0.0, 7.0),
        χ=as(Real, 0.1, 30.0),
        p_test=as(Real, 0.0, 50.0),
        extra_transmissibility = as(Real, 0.5, 2.5),
        influx_exposed_new_variant = as(Real, 0.0, 10e3),
        p_test_boost = as(Real, 1.0, 3.0),
        χ_boost= as(Real, 1.0, 3.0),
        E₀=as(Real, 0.0, 10e3)))
D = Diagonal(0.1*ones(TransformVariables.dimension(trans_one_group)))

#Modified MCMC methods for working with daily fitted ct value

function inferparameters!(model,samples,trans,stepsize::Float64,D::Diagonal,cts,q₀;serowaningrate = 1/365,num_chains = 1)
        println("Starting MCMC parameter inference for $(model.areaname)")
        n = length(q₀)
        l_area(x) = model.log_likelihood(x,model,serowaningrate,cts) + model.log_priors(x)
        ℓ = TransformedLogDensity(trans, l_area)#transformed log-likelihood
        ∇ℓ = LogDensityProblems.ADgradient(:ForwardDiff, ℓ)#transformed log-likelihood gradient wrt the parameters

        results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇ℓ, samples,
                                initialization = (q = q₀,κ=GaussianKineticEnergy(D),ϵ = stepsize),
                                warmup_stages = fixed_stepsize_warmup_stages(M=Symmetric))
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
                                        warmup_stages = fixed_stepsize_warmup_stages(M=Symmetric,stepsize_search = nothing),
                                        reporter = NoProgressReport())
                transformed_results = TransformVariables.transform.(trans,results.chain)
                for i = 1:size(val,1),j = 1:size(val,2)
                        val[i,j,chain_reps] = transformed_results[i][j]
                end
        end

        chn = Chains(val,[String(k) for k in keys(transformed_results[1])])

        model.MCMC_results = KenyaCoVSD.MCMCResults(chn,
                                        [l_area(transformed_results[i]) - model.log_priors(transformed_results[i]) for i = 1:length(transformed_results)],
                                        results.tree_statistics)

        return nothing
end

function inferparameters!(model,samples,trans,stepsize::Float64,D::Diagonal,cts;serowaningrate = 1/365,num_chains = 1)
        println("Searching for a good initial condition for $(model.areaname)")
        l_area(x) = model.log_likelihood(x,model,serowaningrate,cts) + model.log_priors(x)
        f(x) = -transform_logdensity(trans, l_area,x)
        searchrange = fill((-3.,3.),TransformVariables.dimension(trans))
        res = bboptimize(f; SearchRange = searchrange,PopulationSize=1000,MaxSteps=30000)
        q₀ = best_candidate(res)
        inferparameters!(model,samples,trans,stepsize,D,cts,q₀;serowaningrate=serowaningrate,num_chains=num_chains)
        return nothing
end
