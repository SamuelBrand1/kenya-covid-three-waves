#Define the baseline ct at beginning
baseline_ct = ones(24)

#Define a penalty for having too many turning points in ct fit
function ct_penalty(λ,cts)
        λ*sum([(cts[t-1] + cts[t+1] - 2*cts[t])^2 for t = 2:(length(cts)-1)])
end

"""
        function ct_neg_log_density_grad!(grad,log_cts,θ,model::KenyaCoVSD.CoVAreaModel,baseline_contact,λ,p_waning)

In-place calculation of gradient of contact rate negative log-posterior density (first `length(baseline_contact)` days are fixed as Google estimate).
With sero-waning `p_waning`. This is based on the log-ct value (to maintain positivity)
"""
function ct_neg_log_density_grad!(grad,log_cts,θ,model::KenyaCoVSD.CoVAreaModel,baseline_contact,λ,p_waning)
        l = log_cts -> -ll_onegroup_newvariant_infboost_ct(θ,model,p_waning,vcat(baseline_contact,exp.(log_cts))) + ct_penalty(λ,vcat(baseline_contact,exp.(log_cts)))
        ForwardDiff.gradient!(grad,l,log_cts)
        return nothing
end


"""
	function random_grad_log_ct!(grad,x₀,model_ct,baseline_contact,λ,trans_name,p_serowaning)

In-place calculation for gradient of the log-contact rate with sero-waning for a randomly chosen draw from the MCMC
"""
function random_grad_log_ct!(grad,log_cts,model::KenyaCoVSD.CoVAreaModel,baseline_contact,λ,trans_name,p_serowaning)
	k = rand(1:size(model.MCMC_results.chain,1))
	θ = NamedTuple{Tuple(keys(model.MCMC_results.chain))}([model.MCMC_results.chain[k,n,1] for n = 1:size(model.MCMC_results.chain,2)])
	ct_neg_log_density_grad!(grad,log_cts,θ,model,baseline_contact,λ,p_serowaning)
	return nothing
end

"""
	ADAM_step!(X_n::Vector{Float64},
						∇V::Vector{Float64},
						∇V_n::Vector{Float64},
						∇V_sq_n::Vector{Float64},
						n::Int64,
						ADAM_params::NamedTuple)

In-place calculation of parameter, gradient and squared gradient estimates using ADAM SGD method.
"""
function ADAM_step!(X_n::Vector{Float64},
					∇V::Vector{Float64},
					∇V_n::Vector{Float64},
					∇V_sq_n::Vector{Float64},
					hat_∇V_n::Vector{Float64},
					hat_∇V_sq_n::Vector{Float64},
					n::Int64,
					ADAM_params::NamedTuple)
	@unpack β₁,β₂,η,ϵ = ADAM_params
	∇V_n .= β₁.*∇V_n .+ (1-β₁).*∇V
	∇V_sq_n .= β₂.*∇V_sq_n .+ (1-β₂).*(∇V.^2)
	hat_∇V_n .= ∇V_n./(1 - β₁^(n+1))
	hat_∇V_sq_n .= ∇V_sq_n./(1 - β₂^(n+1))
	X_n .= X_n .- η.*hat_∇V_n./(sqrt.(hat_∇V_sq_n) .+ ϵ)
	return nothing
end

"""
	function ADAM_optim(X_opt::Vector{Float64},
                        grad_func!::Function,
                        ADAM_params::NamedTuple,
                        model::KenyaSerology.CoVAreaModel,
                        baseline_contact::Vector{Float64},
                        trans_name,
                        p_waning)

Optimize the vector initialized at `X_opt` using ADAM method with `grad_func!` gradient with sero-waning.
"""
function ADAM_optim(X_opt::Vector{Float64},
			grad_func!::Function,
			ADAM_params::NamedTuple,
			model::KenyaCoVSD.CoVAreaModel,
			baseline_contact::Vector{Float64},
			trans_name,
			p_waning)
	println("Starting ADAM SGD optimisation for effective contact rates.")
	@unpack total_num_steps,averaging_num_steps,λ = ADAM_params

	∇V_n = zeros(length(X_opt))
	∇V_sq_n = zeros(length(X_opt))
	hat_∇V_n = zeros(length(X_opt))
	hat_∇V_sq_n = zeros(length(X_opt))
	grad = zeros(length(X_opt))
	X̂ = zeros(length(X_opt))

	for n = 0:(total_num_steps-1)
		if n%10 == 0
			println("On ADAM step $(n)")
		end
		grad_func!(grad,X_opt,model,baseline_contact,λ,trans_name,p_waning)
		ADAM_step!(X_opt,
					grad,
					∇V_n,∇V_sq_n,
					hat_∇V_n,hat_∇V_sq_n,
					n,ADAM_params)
		if n > total_num_steps - averaging_num_steps #
			m = n - (total_num_steps - averaging_num_steps)
			X̂ .+= (X_opt - X̂)/m
		end

	end
	final_LL = 0
	for k = 1:size(model.MCMC_results.chain,1)
		θ = NamedTuple{Tuple(keys(model.MCMC_results.chain))}([model.MCMC_results.chain[k,n,1] for n = 1:size(model.MCMC_results.chain,2)])
		final_LL += ll_onegroup_newvariant_infboost_ct(θ,model,p_waning,vcat(baseline_contact,exp.(X̂)))
	end
	final_LL = (final_LL/size(model.MCMC_results.chain,1)) - ct_penalty(λ,vcat(baseline_contact,exp.(X̂)))

	return X̂,final_LL
end
