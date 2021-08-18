"""
	loglikelihood_contactratemodel_ct(contactrate,θ,model::KenyaSerology.CoVAreaModel)

Calculate the log-posterior density for `contactrate` for parameters fixed as `θ` and data from `model`.
"""
function loglikelihood_contactratemodel_ct(contactrate,θ,model::KenyaSerology.CoVAreaModel,λ)
    @unpack R,E₀,I₀,χ,M_PCR,α,p_test,P_eff = θ #Parameters to fit 1) initial conditions and R₀ 2) parameters for Neg. bin.
    @unpack PCR_cases,sero_cases,N,σ,γ,prob,PCR_array,sero_array,log_priors,PCR_sensitivity,PCR_specificity,sero_specificity,sero_sensitivity,M_BB,relative_testing_rate = model
    T = eltype(contactrate)
    N_eff = N*P_eff
    u0 = convert.(T, [N_eff,E₀,I₀,0.,0.])
    p = convert.(T,vcat([R,σ,γ,N_eff],contactrate))
	smoothness_penalty = T(0.)
	LL = T(0.)

    try
        sol = solve(prob, BS3(); u0=u0, p=p, saveat = 1)
        incidence = KenyaSerology.get_incidence(sol)
        num_PCR = KenyaSerology.simple_conv(incidence,PCR_array)
        p_sero = sero_sensitivity*KenyaSerology.simple_conv(incidence,sero_array)/N
        for t in 59:min(length(incidence),size(PCR_cases,1))
            d_PCR_pos = PCR_cases[t,1]#Detected number of PCR pos that day
            d_PCR_neg = PCR_cases[t,2]#Detected number of PCR neg that day, -1 === no negative tests available
            d_sero_pos = 0
            n_sero = 0
            if t<= size(sero_cases,1)
                d_sero_pos = sero_cases[t,1] #Detected number of sero pos that day
                n_sero = sero_cases[t,1] + sero_cases[t,2] #Number sero tested on that day
            end
            if d_PCR_neg >= 0
                #Convert from p_PCR,inv_M_PCR parameterisation of Beta-Binomial to standard α,β parameterisation
                p_PCR = χ*num_PCR[t]/((χ-1)*num_PCR[t] + N)
                LL += logpdf(BetaBinomial(d_PCR_pos+d_PCR_neg,p_PCR*M_PCR,(1-p_PCR)*M_PCR),d_PCR_pos)#likelihood contribution from PCR testing --- given
            end
            if d_PCR_neg == -1
                #Convert from μ,α parameterisation to p,r parameterisation
                μ = relative_testing_rate[t]*p_test*num_PCR[t] + 0.001
                σ² = μ + α*μ^2
                p_negbin = 1 - (α*μ^2/σ²)
                r_negbin = 1/α
                LL += logpdf(NegativeBinomial(r_negbin,p_negbin),d_PCR_pos)#likelihood contribution from PCR testing
            end
            #Convert from the p_hat,M_BB parameterisation of the Betabinomial distribution to the α_BB,β_BB parameterisation
            if n_sero > 0
                p_hat = p_sero[t] + (1- p_sero[t])*(1-sero_specificity)
                α_BB = M_BB*p_hat
                β_BB = M_BB*(1-p_hat)
                LL += logpdf(BetaBinomial(n_sero,α_BB,β_BB),d_sero_pos)#Likelihood contribution from sero testing
            end
        end
        # return LL + log_priors(θ)
		for t in 2:(length(contactrate)-1)
			smoothness_penalty += (contactrate[t-1] + contactrate[t+1] - 2*contactrate[t])^2
		end
		return LL + log_priors(θ) - (λ*smoothness_penalty/length(contactrate))

    catch
        return T(-Inf)
    end
end

"""
        function loglikelihood_contactratemodel_ct(contactrate,θ,model::KenyaSerology.CoVAreaModel,λ,p_waning)

Calculate the log-posterior density for `contactrate` for parameters fixed as `θ` and data from `model` with serological waning rate `p_waning`.
"""
function loglikelihood_contactratemodel_ct(contactrate,θ,model::KenyaSerology.CoVAreaModel,λ,p_waning)
    @unpack R,E₀,I₀,χ,M_PCR,α,p_test,P_eff = θ #Parameters to fit 1) initial conditions and R₀ 2) parameters for Neg. bin.
    @unpack PCR_cases,sero_cases,N,σ,γ,prob,PCR_array,sero_array,log_priors,PCR_sensitivity,PCR_specificity,sero_specificity,sero_sensitivity,M_BB,relative_testing_rate = model
    T = eltype(contactrate)
    N_eff = N*P_eff
    u0 = convert.(T, [N_eff,E₀,I₀,0.,0.])
    p = convert.(T,vcat([R,σ,γ,N_eff],contactrate))
    smoothness_penalty = T(0.)

    for t = 31:length(sero_array)
        sero_array[t] = (1-p_waning)^(t-30)
    end
    LL = T(0.)

    try
        sol = solve(prob, BS3(); u0=u0, p=p, saveat = 1)
        incidence = KenyaSerology.get_incidence(sol)
        num_PCR = KenyaSerology.simple_conv(incidence,PCR_array)
        p_sero = sero_sensitivity*KenyaSerology.simple_conv(incidence,sero_array)/N
        for t in 59:min(length(incidence),size(PCR_cases,1))
            d_PCR_pos = PCR_cases[t,1]#Detected number of PCR pos that day
            d_PCR_neg = PCR_cases[t,2]#Detected number of PCR neg that day, -1 === no negative tests available
            d_sero_pos = 0
            n_sero = 0
            if t<= size(sero_cases,1)
                d_sero_pos = sero_cases[t,1] #Detected number of sero pos that day
                n_sero = sero_cases[t,1] + sero_cases[t,2] #Number sero tested on that day
            end
            if d_PCR_neg >= 0
                #Convert from p_PCR,inv_M_PCR parameterisation of Beta-Binomial to standard α,β parameterisation
                p_PCR = χ*num_PCR[t]/((χ-1)*num_PCR[t] + N)
                LL += logpdf(BetaBinomial(d_PCR_pos+d_PCR_neg,p_PCR*M_PCR,(1-p_PCR)*M_PCR),d_PCR_pos)#likelihood contribution from PCR testing --- given
            end
            if d_PCR_neg == -1
                #Convert from μ,α parameterisation to p,r parameterisation
                μ = relative_testing_rate[t]*p_test*num_PCR[t] + 0.001
                σ² = μ + α*μ^2
                p_negbin = 1 - (α*μ^2/σ²)
                r_negbin = 1/α
                LL += logpdf(NegativeBinomial(r_negbin,p_negbin),d_PCR_pos)#likelihood contribution from PCR testing
            end
            #Convert from the p_hat,M_BB parameterisation of the Betabinomial distribution to the α_BB,β_BB parameterisation
            if n_sero > 0
                p_hat = p_sero[t] + (1- p_sero[t])*(1-sero_specificity)
                α_BB = M_BB*p_hat
                β_BB = M_BB*(1-p_hat)
                LL += logpdf(BetaBinomial(n_sero,α_BB,β_BB),d_sero_pos)#Likelihood contribution from sero testing
            end
        end
        # return LL + log_priors(θ)
		for t in 2:(length(contactrate)-1)
			smoothness_penalty += (contactrate[t-1] + contactrate[t+1] - 2*contactrate[t])^2
		end
		
		return LL + log_priors(θ) - (λ*smoothness_penalty/length(contactrate))

    catch
        return T(-Inf)
    end
end

"""
        function loglikelihood_contactratemodel_ct(contactrate,θ,model::KenyaSerology.CoVAreaModel,λ,p_waning,σ²_aboveone)

Calculate the log-posterior density for `contactrate` for parameters fixed as `θ` and data from `model` with serological waning rate `p_waning`.
`σ²_aboveone` scales the penalty for the Ct value going above 1.
"""
function loglikelihood_contactratemodel_ct(contactrate,θ,model::KenyaSerology.CoVAreaModel,λ,p_waning,σ²_aboveone)
    @unpack R,E₀,I₀,χ,M_PCR,α,p_test,P_eff = θ #Parameters to fit 1) initial conditions and R₀ 2) parameters for Neg. bin.
    @unpack PCR_cases,sero_cases,N,σ,γ,prob,PCR_array,sero_array,log_priors,PCR_sensitivity,PCR_specificity,sero_specificity,sero_sensitivity,M_BB,relative_testing_rate = model
    T = eltype(contactrate)
    N_eff = N*P_eff
    u0 = convert.(T, [N_eff,E₀,I₀,0.,0.])
    p = convert.(T,vcat([R,σ,γ,N_eff],contactrate))
    smoothness_penalty = T(0.)
	above_one_penalty = T(0.)

    for t = 31:length(sero_array)
        sero_array[t] = (1-p_waning)^(t-30)
    end
    LL = T(0.)

    try
        sol = solve(prob, BS3(); u0=u0, p=p, saveat = 1)
        incidence = KenyaSerology.get_incidence(sol)
        num_PCR = KenyaSerology.simple_conv(incidence,PCR_array)
        p_sero = sero_sensitivity*KenyaSerology.simple_conv(incidence,sero_array)/N
        for t in 59:min(length(incidence),size(PCR_cases,1))
            d_PCR_pos = PCR_cases[t,1]#Detected number of PCR pos that day
            d_PCR_neg = PCR_cases[t,2]#Detected number of PCR neg that day, -1 === no negative tests available
            d_sero_pos = 0
            n_sero = 0
            if t<= size(sero_cases,1)
                d_sero_pos = sero_cases[t,1] #Detected number of sero pos that day
                n_sero = sero_cases[t,1] + sero_cases[t,2] #Number sero tested on that day
            end
            if d_PCR_neg >= 0
                #Convert from p_PCR,inv_M_PCR parameterisation of Beta-Binomial to standard α,β parameterisation
                p_PCR = χ*num_PCR[t]/((χ-1)*num_PCR[t] + N)
                LL += logpdf(BetaBinomial(d_PCR_pos+d_PCR_neg,p_PCR*M_PCR,(1-p_PCR)*M_PCR),d_PCR_pos)#likelihood contribution from PCR testing --- given
            end
            if d_PCR_neg == -1
                #Convert from μ,α parameterisation to p,r parameterisation
                μ = relative_testing_rate[t]*p_test*num_PCR[t] + 0.001
                σ² = μ + α*μ^2
                p_negbin = 1 - (α*μ^2/σ²)
                r_negbin = 1/α
                LL += logpdf(NegativeBinomial(r_negbin,p_negbin),d_PCR_pos)#likelihood contribution from PCR testing
            end
            #Convert from the p_hat,M_BB parameterisation of the Betabinomial distribution to the α_BB,β_BB parameterisation
            if n_sero > 0
                p_hat = p_sero[t] + (1- p_sero[t])*(1-sero_specificity)
                α_BB = M_BB*p_hat
                β_BB = M_BB*(1-p_hat)
                LL += logpdf(BetaBinomial(n_sero,α_BB,β_BB),d_sero_pos)#Likelihood contribution from sero testing
            end
        end
        # return LL + log_priors(θ)
		for t in 2:(length(contactrate)-1)
			smoothness_penalty += (contactrate[t-1] + contactrate[t+1] - 2*contactrate[t])^2
		end
		for t in 1:(length(contactrate))
			if contactrate[t] > 1
				above_one_penalty += (log(contactrate[t])^2)/(2*σ²_aboveone) #Log-normal assumption above Ct = 1,0.174 chosen so that < 50% greater contacts are 99% probable
			end
		end
		return LL + log_priors(θ) - (λ*smoothness_penalty/length(contactrate)) - above_one_penalty

    catch
        return T(-Inf)
    end
end

"""
	C_t_neg_log_density_grad!(grad,ct,θ,model::KenyaSerology.CoVAreaModel,log_prior_ct::Function,baseline_contact)

In-place calculation of gradient of contact rate negative log-posterior density (first `length(baseline_contact)` days are fixed as Google estimate).
"""
function C_t_neg_log_density_grad!(grad,ct,θ,model::KenyaSerology.CoVAreaModel,baseline_contact,λ)
        l = x -> -loglikelihood_contactratemodel_ct(vcat(baseline_contact,x),θ,model,λ)
        ForwardDiff.gradient!(grad,l,ct)
        return nothing
end

"""
        function C_t_neg_log_density_grad!(grad,ct,θ,model::KenyaSerology.CoVAreaModel,baseline_contact,λ,p_waning)

In-place calculation of gradient of contact rate negative log-posterior density (first `length(baseline_contact)` days are fixed as Google estimate).
With sero-waning `p_waning`.
"""
function C_t_neg_log_density_grad!(grad,ct,θ,model::KenyaSerology.CoVAreaModel,baseline_contact,λ,p_waning)
        l = x -> -loglikelihood_contactratemodel_ct(vcat(baseline_contact,x),θ,model,λ,p_waning)
        ForwardDiff.gradient!(grad,l,ct)
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
	function random_grad_log_ct!(grad,x₀,model_ct,baseline_contact,λ)

In-place calculation for gradient of the log-contact rate.
"""
function random_grad_log_ct!(grad,x₀,model_ct,baseline_contact,λ)
    k = rand(1:size(model_ct.MCMC_results.chain,1))
    θ = (R=model_ct.MCMC_results.chain[k,1,1],
        E₀ =model_ct.MCMC_results.chain[k,2,1],
        I₀ =model_ct.MCMC_results.chain[k,3,1],
        χ = model_ct.MCMC_results.chain[k,4,1],
        M_PCR = model_ct.MCMC_results.chain[k,5,1],
        α = model_ct.MCMC_results.chain[k,6,1],
        p_test = model_ct.MCMC_results.chain[k,7,1],
        P_eff = model_ct.MCMC_results.chain[k,8,1])
    KenyaSerology.C_t_neg_log_density_grad!(grad,exp.(x₀),θ,model_ct,baseline_contact,λ)
    grad .= grad./exp.(x₀)
    return nothing
end

"""
	function random_grad_log_ct!(grad,x₀,model_ct,baseline_contact,λ,trans_name,p_serowaning)

In-place calculation for gradient of the log-contact rate with sero-waning.
"""
function random_grad_log_ct!(grad,x₀,model_ct,baseline_contact,λ,trans_name,p_serowaning)
    k = rand(1:size(model_ct.MCMC_results.chain,1))
	θ = TransformVariables.transform(trans_name,[model_ct.MCMC_results.chain[1,n,1] for n = 1:size(model_ct.MCMC_results.chain,2)])
    KenyaSerology.C_t_neg_log_density_grad!(grad,exp.(x₀),θ,model_ct,baseline_contact,λ,p_serowaning)
    grad .= grad./exp.(x₀)
    return grad
end

"""
	ADAM_optim(X_opt::Vector{Float64},
					grad_func!::Function,
					total_num_steps::Int64,
					averaging_num_steps::Int64,
					ADAM_params::NamedTuple,
					model::KenyaSerology.CoVAreaModel,
					baseline_contact::Vector{Float64},λ)

Optimize the vector initialized at `X_opt` using ADAM method with `grad_func!` gradient.
"""
function ADAM_optim(X_opt::Vector{Float64},
					grad_func!::Function,
					ADAM_params::NamedTuple,
                    model::KenyaSerology.CoVAreaModel,
                    baseline_contact::Vector{Float64})
	@unpack total_num_steps,averaging_num_steps,λ = ADAM_params
	∇V_n = zeros(length(X_opt))
	∇V_sq_n = zeros(length(X_opt))
    hat_∇V_n = zeros(length(X_opt))
    hat_∇V_sq_n = zeros(length(X_opt))
    grad = zeros(length(X_opt))
	X̂ = zeros(length(X_opt))

	for n = 0:(total_num_steps-1)
        grad_func!(grad,X_opt,model,baseline_contact,λ)
        KenyaSerology.ADAM_step!(X_opt,
                                grad,
                                ∇V_n,∇V_sq_n,
                                hat_∇V_n,hat_∇V_sq_n,
                                n,ADAM_params)
		if n > total_num_steps - averaging_num_steps #
			m = n - (total_num_steps - averaging_num_steps)
			X̂ .+= (X_opt - X̂)/m
		end
	end

	return X̂
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
					model::KenyaSerology.CoVAreaModel,
					baseline_contact::Vector{Float64},
					trans_name,
					p_waning)
	println("Starting ADAM SGD optimisation for effective contact rates.")
	@unpack total_num_steps,averaging_num_steps,λ,σ²_aboveone = ADAM_params

	∇V_n = zeros(length(X_opt))
	∇V_sq_n = zeros(length(X_opt))
	hat_∇V_n = zeros(length(X_opt))
	hat_∇V_sq_n = zeros(length(X_opt))
	grad = zeros(length(X_opt))
	X̂ = zeros(length(X_opt))

	for n = 0:(total_num_steps-1)
		grad_func!(grad,X_opt,model,baseline_contact,λ,trans_name,p_waning)
		KenyaSerology.ADAM_step!(X_opt,
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
		θ = TransformVariables.transform(trans_name,[model.MCMC_results.chain[k,n,1] for n = 1:size(model.MCMC_results.chain,2)])
		final_LL += KenyaSerology.loglikelihood_contactratemodel_ct(vcat(baseline_contact,exp.(X̂)),θ,model,λ,p_waning,σ²_aboveone) #NB this is without smoothing penalty or priors just predictive ability
	end
	final_LL = final_LL/size(model.MCMC_results.chain,1)
        
	return X̂,final_LL
end

"""
	function ADAM_optim(X_opt::Vector{Float64},
                        grad_func!::Function,
                        ADAM_params::NamedTuple,
                        model::KenyaSerology.CoVAreaModel,
                        baseline_contact::Vector{Float64},
                        trans_name,
                        p_waning,
                        plt::Plots.Plot)

Optimize the vector initialized at `X_opt` using ADAM method with `grad_func!` gradient.
Every 10 steps plots onto plt so convergence can be visualized.
"""

function ADAM_optim(X_opt::Vector{Float64},
                        grad_func!::Function,
                        ADAM_params::NamedTuple,
                        model::KenyaSerology.CoVAreaModel,
                        baseline_contact::Vector{Float64},
                        trans_name,
                        p_waning,
                        plt::Plots.Plot)

	println("Starting ADAM SGD optimisation for effective contact rates, with plotting output.")

	@unpack total_num_steps,averaging_num_steps,λ = ADAM_params
        ∇V_n = zeros(length(X_opt))
        ∇V_sq_n = zeros(length(X_opt))
        hat_∇V_n = zeros(length(X_opt))
        hat_∇V_sq_n = zeros(length(X_opt))
        grad = zeros(length(X_opt))
        X̂ = zeros(length(X_opt))

	for n = 0:(total_num_steps-1)
                grad_func!(grad,X_opt,model,baseline_contact,λ,trans_name,p_waning)
                KenyaSerology.ADAM_step!(X_opt,
                        grad,
                        ∇V_n,∇V_sq_n,
                        hat_∇V_n,hat_∇V_sq_n,
                        n,ADAM_params)
                if n > total_num_steps - averaging_num_steps #
                        m = n - (total_num_steps - averaging_num_steps)
                        X̂ .+= (X_opt - X̂)/m
                end
                if n%100 == 0
					LL = 0
					for k = 1:size(model.MCMC_results.chain,1)
                        θ = TransformVariables.transform(trans_name,[model.MCMC_results.chain[1,n,1] for n = 1:size(model.MCMC_results.chain,2)])
                        LL -= KenyaSerology.loglikelihood_contactratemodel_ct(vcat(baseline_contact,exp.(X_opt)),θ,model,λ,p_waning)
					end
					LL = LL/size(model.MCMC_results.chain,1)
                    scatter!(plt,[n],[LL],lab = "",color = :black)
					# display(plot(plt))
                end
	end
        plot!(plt,xlabel = "Steps",ylabel = "Mean Neg. LL",yscale = :log10)
	final_LL = 0
	for k = 1:size(model.MCMC_results.chain,1)
		θ = TransformVariables.transform(trans_name,[model.MCMC_results.chain[k,n,1] for n = 1:size(model.MCMC_results.chain,2)])
		final_LL += KenyaSerology.loglikelihood_contactratemodel_ct(vcat(baseline_contact,exp.(X̂)),θ,model,λ,p_waning)
	end
	final_LL = final_LL/size(model.MCMC_results.chain,1)
	return X̂,final_LL,plt
end



"""
	get_model_and_optimise_ADAM(filename,projected_contact_rate;
									varname = "model",
									enddate::Date,
									fixed_days::Int64,
									ADAM_params::NamedTuple)

Load modelfit from `filename`, and use ADAM SGD optimization to improve C(t).
"""
function get_model_and_optimise_ADAM(filename,projected_contact_rate;
									varname = "model",
									enddate::Date,
									fixed_days::Int64,
									ADAM_params::NamedTuple)

	model = load(filename)[varname]
	# model_ct = deepcopy(model)
	model.prob = make_odeproblemforinference_parameter_ct(projected_contact_rate;
									startdate=Date(2020,2,20),
									enddate=enddate)

	n = (enddate - Date(2020,2,20)).value
	baseline_contact = copy(projected_contact_rate.contactrate[1:fixed_days])
	log_ct = log.(copy(projected_contact_rate.contactrate[(fixed_days+1):end]))
	len_ct = length(projected_contact_rate.contactrate)
	if n > len_ct
		log_ct = vcat(log_ct,fill(log_ct[end],n-len_ct))
	end
	println("Got fit using Google mobility and starting first ADAM grad. descent for $(model.areaname)")
	log_ct_fit = ADAM_optim(copy(log_ct),
				random_grad_log_ct!,
				ADAM_params,
	            model,
	            baseline_contact)

    smoothed_ct = vcat(baseline_contact,exp.(log_ct_fit))
    return smoothed_ct,model
end

"""
	optimise_ct_model_ADAM(model::CoVAreaModel,projected_contact_rate;
							enddate::Date,
							fixed_days::Int64,
							ADAM_params::NamedTuple)

Fit a contact rate using ADAM SGD.
"""
function optimise_ct_model_ADAM(model::CoVAreaModel,projected_contact_rate;
									enddate::Date,
									fixed_days::Int64,
									ADAM_params::NamedTuple)

	model.prob = make_odeproblemforinference_parameter_ct(projected_contact_rate;
									startdate=Date(2020,2,20),
									enddate=enddate)

	n = (enddate - Date(2020,2,20)).value
	baseline_contact = copy(projected_contact_rate.contactrate[1:fixed_days])
	log_ct = log.(copy(projected_contact_rate.contactrate[(fixed_days+1):end]))
	len_ct = length(projected_contact_rate.contactrate)
	if n > len_ct
		log_ct = vcat(log_ct,fill(log_ct[end],n-len_ct))
	end
	println("Refit eff. contact rate using ADAM grad. descent for $(model.areaname)")

	log_ct_fit = ADAM_optim(copy(log_ct),
				random_grad_log_ct!,
				ADAM_params,
	            model,
	            baseline_contact)

	smoothed_ct = vcat(baseline_contact,exp.(log_ct_fit))

    return smoothed_ct
end


"""
	function EM_optimise_ADAM(filename,projected_contact_rate,trans;
							num_steps=3,
							varname="model",
							enddate::Date,
							fixed_days::Int64,
							ADAM_params::NamedTuple)

Perform `num_steps` EM algorithm iterations.
"""
function EM_optimise_ADAM(filename,projected_contact_rate,trans;
						num_steps=3,
						varname="model",
						enddate::Date,
						fixed_days::Int64,
						ADAM_params::NamedTuple)
	#Initial fit of C(t)
    smoothed_ct,model = get_model_and_optimise_ADAM(filename,projected_contact_rate;
													varname = varname,
													enddate=enddate,
													fixed_days=fixed_days,
													ADAM_params=ADAM_params)

	model.prob = make_odeproblemforinference_simple(smoothed_ct,#Method for defining the ODE problem underlying the inference
													startdate = Date(2020,2,21),#Don't change from Feb 21st as start date!
													enddate = enddate)
    model.contactrate_data = copy(smoothed_ct)
	if num_steps > 1
	    KenyaSerology.inferparameters!(model,1000,trans)
	end
	if num_steps == 1
		KenyaSerology.inferparameters!(model,10000,trans)
	end

	#Middle fits of C(t)
    for i = 2:(num_steps-1)
        smoothed_ct = optimise_ct_model_ADAM(model,projected_contact_rate;
											enddate=enddate,
											fixed_days=fixed_days,
											ADAM_params=ADAM_params)

        model.prob = make_odeproblemforinference_simple(smoothed_ct,#Method for defining the ODE problem underlying the inference
                                                        startdate = Date(2020,2,21),#Don't change from Feb 21st as start date!
                                                        enddate = enddate)

        model.contactrate_data = copy(smoothed_ct)
        KenyaSerology.inferparameters!(model,1000,trans)
    end

	#Final fit of C(t)
	if num_steps > 1
		smoothed_ct = optimise_ct_model_ADAM(model,projected_contact_rate;
											enddate=enddate,
											fixed_days=fixed_days,
											ADAM_params=ADAM_params)

		model.prob = make_odeproblemforinference_simple(smoothed_ct,#Method for defining the ODE problem underlying the inference
														startdate = Date(2020,2,21),#Don't change from Feb 21st as start date!
														enddate = enddate)

		model.contactrate_data = copy(smoothed_ct)
		KenyaSerology.inferparameters!(model,10000,trans)
	end
    return model
end
