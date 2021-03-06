#Code for defining the one social group model for Nairobi
# In this file are defined the ODE model, the log-likelihood function, and the log-prior function

"""
	function get_ct(t,daily_ct_values)

Get the contact rate value for time `t` from a generated time series of daily cts.
"""
function get_ct(t,daily_ct_values)
	ct = t >= length(daily_ct_values) ? daily_ct_values[end] : daily_ct_values[max(Int(floor(t))+1,1)]
	return ct
end

"""
    function onegrouptransmission(du,u,p,t)

In-place creation of vector field for the one group transmission model.
"""
function onegrouptransmission(du,u,p,t)
    #Get state
	S,E,I,R,W,C,C_sero  = u
    #Get parameters
	R₀,α,γ,σ,N,ω = p[1:6]
	#Get daily ct value
	cts = p[7:end]
    #Force of infection
	ct = get_ct(t,cts)
    λ = γ*R₀*ct*I
	#Dynamics
	du[1] = -ct*(S/N)*λ
	du[2] = ct*((S+σ*W)/N)*λ - α*E
	du[3] = α*E - γ*I
	du[4] = γ*I - ω*R
	du[5] = ω*R - σ*ct*(W/N)*λ
	du[6] = ct*((S+σ*W)/N)*λ
    du[7] = ct*(S/N)*λ
	return nothing
end

#Create basic ODEProblem for one group model
u0_onegroup = [basic_nai_model.N - 100.,100.,0.,0.,0.,0.,0.0]
f_onegroup = ODEFunction(onegrouptransmission;syms = [:S,:E,:I,:R,:W,:C,:C_sero])
prob_one_group = ODEProblem(f_onegroup,u0_onegroup,(0,500))

#Create log-likelihood function and log-prior function for one group model

"""
    function ll_onegroup_newvariant_infboost(θ,model::KenyaCoVSD.CoVAreaModel,seroreversionrate,cts;σ = 0.16,ω = 1/180,ct_min2 = 0.445)

Log-likelihood of the case and serology data in `model` assuming: 1) one group dynamics, and 2) that new variant is more transmissible rather than immune evading.
"""
function ll_onegroup_newvariant_infboost(θ,model::KenyaCoVSD.CoVAreaModel,seroreversionrate,cts;σ = 0.16,ω = 1/180)
    @unpack PCR_cases,sero_cases,baseline_sero_array,PCR_array,sero_sensitivity,sero_specificity,N,M_BB,prob,α,γ,relative_testing_rate = model
    @unpack R₀,χ,p_test,extra_transmissibility,influx_exposed_new_variant,p_test_boost,χ_boost,E₀ = θ

    #Set variance scalers
    clustering_factor_PCR = 0.5
    M_PCR = 30.
    T = eltype(R₀)
	# if T <: Real
	# 	T = eltype(cts)
	# end

    #Set transmission parameters with daily cts appended on end
    p = convert.(T,vcat([R₀,α,γ,σ,N,ω],cts))
    u0 = convert.(T,[N,E₀,0.,0.,0.,0.,0.])

    #Create new variant introduction
    function new_variant_effect!(integrator)
        integrator.p[1] *= extra_transmissibility
        integrator.u[2] += influx_exposed_new_variant
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
					verbose = false,
                    isoutofdomain=(u,p,t) -> any(x -> x < 0, u))

        ι = diff(sol[:C])
        ι_sero = diff(sol[:C_sero])
        PCR = KenyaCoVSD.simple_conv(ι,PCR_array)
    	sero = sero_sensitivity.*KenyaCoVSD.simple_conv(ι_sero,sero_array)

        #Calculate log-likelihood for PCR testing
        for t in 55:janendpoint
            #Convert from μ,α parameterisation to p,r parameterisation for using Neg. binomial model
            μ = relative_testing_rate[t]*p_test*PCR[t]*1e-4 + 0.001 #Covert approximately from total probability of detection in % to daily probability of detection
            σ² = μ + clustering_factor_PCR*μ^2
            p_negbin = 1 - (clustering_factor_PCR*μ^2/σ²)
            r_negbin  = 1/clustering_factor_PCR
            LL += logpdf(NegativeBinomial(r_negbin,p_negbin),PCR_cases[t,1])#likelihood contribution from PCR testing --- positive case detection

            if PCR_cases[t,2] >= 0 #Negative tests available
                p_PCR_pred = (χ*PCR[t]/((χ-1)*PCR[t] + N)) + 0.001
                #Convert from p_PCR,inv_M_PCR parameterisation of Beta-Binomial to standard α,β parameterisation
                LL += logpdf(BetaBinomial(PCR_cases[t,2],p_PCR_pred*M_PCR,(1-p_PCR_pred)*M_PCR),PCR_cases[t,1])#likelihood contribution from PCR testing --- proportion postive
            end
        end

        #Calculate log-likelihood for PCR testing after new variant introduction
        for t in (janendpoint+1):min(size(PCR_cases,1),size(PCR,1))
            #Convert from μ,α parameterisation to p,r parameterisation for using Neg. binomial model
            μ = p_test_boost*relative_testing_rate[t]*p_test*PCR[t]*1e-4 + 0.001
            σ² = μ + clustering_factor_PCR*μ^2
            p_negbin = 1 - (clustering_factor_PCR*μ^2/σ²)
            r_negbin  = 1/clustering_factor_PCR
            LL += logpdf(NegativeBinomial(r_negbin,p_negbin),PCR_cases[t,1])#likelihood contribution from PCR testing --- positive case detection

            if PCR_cases[t,2] >= 0 #Negative tests available
                p_PCR_pred = (χ_boost*χ*PCR[t]/((χ_boost*χ-1)*PCR[t] + N))  + 0.001
                #Convert from p_PCR,inv_M_PCR parameterisation of Beta-Binomial to standard α,β parameterisation
                LL += logpdf(BetaBinomial(PCR_cases[t,2],p_PCR_pred*M_PCR,(1-p_PCR_pred)*M_PCR),PCR_cases[t,1])#likelihood contribution from PCR testing that day in that age group
            end
        end
        #Calculate log-likelihood contribution from serological testing
        for t in 1:min(size(sero_cases,1),size(sero,1))
            #Convert from the p_hat,M_BB parameterisation of the Betabinomial distribution to the α_BB,β_BB parameterisation
            #Overall serological testing
            if sero_cases[t,1,1] + sero_cases[t,2,1] > 0
                p_sero_pred  = sero[t]/N
                p_hat  = p_sero_pred + (1-p_sero_pred)*(1-sero_specificity)
                LL += logpdf(BetaBinomial(sero_cases[t,1,1] + sero_cases[t,2,1],M_BB*p_hat,M_BB*(1-p_hat)),sero_cases[t,1,1])#Likelihood contribution from sero testing
            end

        end
    catch
        LL = T(-Inf)
    end
    return LL::T
end

function ll_onegroup_newvariant_infboost_ct(θ,model::KenyaCoVSD.CoVAreaModel,seroreversionrate,cts;σ = 0.16,ω = 1/180)
    @unpack PCR_cases,sero_cases,baseline_sero_array,PCR_array,sero_sensitivity,sero_specificity,N,M_BB,prob,α,γ,relative_testing_rate = model
    @unpack R₀,χ,p_test,extra_transmissibility,influx_exposed_new_variant,p_test_boost,χ_boost,E₀ = θ

    #Set variance scalers
    clustering_factor_PCR = 0.5
    M_PCR = 30.
    # T = eltype(R₀)
	# if T <: Real
    T = eltype(cts)
	# end

    #Set transmission parameters with daily cts appended on end
    p = convert.(T,vcat([R₀,α,γ,σ,N,ω],cts))
    u0 = convert.(T,[N,E₀,0.,0.,0.,0.,0.])

    #Create new variant introduction
    function new_variant_effect!(integrator)
        integrator.p[1] *= extra_transmissibility
        integrator.u[2] += influx_exposed_new_variant
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
					verbose = false,
                    isoutofdomain=(u,p,t) -> any(x -> x < 0, u))

        ι = diff(sol[:C])
        ι_sero = diff(sol[:C_sero])
        PCR = KenyaCoVSD.simple_conv(ι,PCR_array)
    	sero = sero_sensitivity.*KenyaCoVSD.simple_conv(ι_sero,sero_array)

        #Calculate log-likelihood for PCR testing
        for t in 55:janendpoint
            #Convert from μ,α parameterisation to p,r parameterisation for using Neg. binomial model
            μ = relative_testing_rate[t]*p_test*PCR[t]*1e-4 + 0.001 #Covert approximately from total probability of detection in % to daily probability of detection
            σ² = μ + clustering_factor_PCR*μ^2
            p_negbin = 1 - (clustering_factor_PCR*μ^2/σ²)
            r_negbin  = 1/clustering_factor_PCR
            LL += logpdf(NegativeBinomial(r_negbin,p_negbin),PCR_cases[t,1])#likelihood contribution from PCR testing --- positive case detection

            if PCR_cases[t,2] >= 0 #Negative tests available
                p_PCR_pred = (χ*PCR[t]/((χ-1)*PCR[t] + N)) + 0.001
                #Convert from p_PCR,inv_M_PCR parameterisation of Beta-Binomial to standard α,β parameterisation
                LL += logpdf(BetaBinomial(PCR_cases[t,2],p_PCR_pred*M_PCR,(1-p_PCR_pred)*M_PCR),PCR_cases[t,1])#likelihood contribution from PCR testing --- proportion postive
            end
        end

        #Calculate log-likelihood for PCR testing after new variant introduction
        for t in (janendpoint+1):min(size(PCR_cases,1),size(PCR,1))
            #Convert from μ,α parameterisation to p,r parameterisation for using Neg. binomial model
            μ = p_test_boost*relative_testing_rate[t]*p_test*PCR[t]*1e-4 + 0.001
            σ² = μ + clustering_factor_PCR*μ^2
            p_negbin = 1 - (clustering_factor_PCR*μ^2/σ²)
            r_negbin  = 1/clustering_factor_PCR
            LL += logpdf(NegativeBinomial(r_negbin,p_negbin),PCR_cases[t,1])#likelihood contribution from PCR testing --- positive case detection

            if PCR_cases[t,2] >= 0 #Negative tests available
                p_PCR_pred = (χ_boost*χ*PCR[t]/((χ_boost*χ-1)*PCR[t] + N))  + 0.001
                #Convert from p_PCR,inv_M_PCR parameterisation of Beta-Binomial to standard α,β parameterisation
                LL += logpdf(BetaBinomial(PCR_cases[t,2],p_PCR_pred*M_PCR,(1-p_PCR_pred)*M_PCR),PCR_cases[t,1])#likelihood contribution from PCR testing that day in that age group
            end
        end
        #Calculate log-likelihood contribution from serological testing
        for t in 1:min(size(sero_cases,1),size(sero,1))
            #Convert from the p_hat,M_BB parameterisation of the Betabinomial distribution to the α_BB,β_BB parameterisation
            #Overall serological testing
            if sero_cases[t,1,1] + sero_cases[t,2,1] > 0
                p_sero_pred  = sero[t]/N
                p_hat  = p_sero_pred + (1-p_sero_pred)*(1-sero_specificity)
                LL += logpdf(BetaBinomial(sero_cases[t,1,1] + sero_cases[t,2,1],M_BB*p_hat,M_BB*(1-p_hat)),sero_cases[t,1,1])#Likelihood contribution from sero testing
            end

        end
    catch
        LL = T(-Inf)
    end
    return LL::T
end

function priors_onegroup_newvariant(θ)
	@unpack R₀,χ,p_test,extra_transmissibility,influx_exposed_new_variant,p_test_boost,χ_boost,E₀ = θ
    LP = 0.
    LP += logpdf(Gamma(3,2.5/3),R₀)
    LP += logpdf(Gamma(3,2/3),χ)
    LP += logpdf(Gamma(3,1/3),p_test)#In the region of 1% detection rate
    LP += logpdf(Gamma(10,1.5/10),extra_transmissibility)
    LP += logpdf(Gamma(5,100/5),influx_exposed_new_variant)
	LP += logpdf(Normal(0,0.1),log(χ_boost))
	LP += logpdf(Normal(0,0.1),log(p_test_boost))
	LP += logpdf(Gamma(3,100/3),E₀)

    return LP
end

## Create one group model
nai_one_group = deepcopy(basic_nai_model)
nai_one_group.prob = prob_one_group
nai_one_group.log_likelihood = ll_onegroup_newvariant_infboost
nai_one_group.log_priors = priors_onegroup_newvariant


function gather_uncertainty_one_group(nai_one_group,ct_fitted)
    #Placeholder solve run to get length of matrix
    p = [[2.5,nai_one_group.α,nai_one_group.γ,0.16,N,1/180];ct_fitted]
    u0 = [N-100,100,0.0,0.0,0.0,0.0,0.0]

    function new_variant_effect!(integrator)
            integrator.p[1] *= 1.5
            integrator.u[2] += 0.0
    end

    janendpoint = (Date(2021,1,30) - Date(2020,2,20)).value
    aprilendpoint = (Date(2021,4,30) - Date(2020,2,20)).value
    variant_cb = PresetTimeCallback([janendpoint],new_variant_effect!)


    sol = solve(nai_one_group.prob, BS3();tspan = (0,(Date(2021,6,1) - Date(2020,2,20)).value),
                       callback = variant_cb, u0=u0, p=p, saveat = 1)
    ι = diff(sol[:C])
    #Define arrays
    ## Gather the uncertainty over the PCR and seropos predictions
    prop_PCR_pos_mat = zeros(length(ι),size(nai_one_group.MCMC_results.chain,1))
    no_neg_PCR_pos_mat = zeros(length(ι),size(nai_one_group.MCMC_results.chain,1))
    prop_sero_pos_mat = zeros(length(ι),size(nai_one_group.MCMC_results.chain,1))
    infection_mat = zeros(length(ι),size(nai_one_group.MCMC_results.chain,1))

    sero_array = vcat(nai_one_group.baseline_sero_array[1:30],[(1-0)^k for k in 1:500])

    for j = 1:size(nai_one_group.MCMC_results.chain,1)
        R₀ = nai_one_group.MCMC_results.chain[:R₀][j]
        extra_transmissibility = nai_one_group.MCMC_results.chain[:extra_transmissibility][j]
        influx_exposed_new_variant = nai_one_group.MCMC_results.chain[:influx_exposed_new_variant][j]
        
        p_test = nai_one_group.MCMC_results.chain[:p_test][j]
        χ = nai_one_group.MCMC_results.chain[:χ][j]
        E₀ = nai_one_group.MCMC_results.chain[:E₀][j]
        p_test_boost = nai_one_group.MCMC_results.chain[:p_test_boost][j]
        χ_boost = nai_one_group.MCMC_results.chain[:χ_boost][j]

        p = [[R₀,nai_one_group.α,nai_one_group.γ,0.16,N,1/180];ct_fitted]
        u0 = [N-E₀,E₀,0.0,0.0,0.0,0.0,0.0]

        function affect!(integrator)
                integrator.p[1] *= extra_transmissibility
                integrator.u[2] += influx_exposed_new_variant
        end

        janendpoint = (Date(2021,1,30) - Date(2020,2,20)).value
        aprilendpoint = (Date(2021,4,30) - Date(2020,2,20)).value
        variant_cb = PresetTimeCallback([janendpoint],affect!)
        
        sol = solve(nai_one_group.prob, BS3();tspan = (0,(Date(2021,6,1) - Date(2020,2,20)).value),
                        callback = variant_cb, u0=u0, p=p, saveat = 1)
        ι = diff(sol[:C])
        ι_sero = diff(sol[:C_sero])
        PCR = KenyaCoVSD.simple_conv(ι,nai_one_group.PCR_array)
        sero = nai_one_group.sero_sensitivity.*KenyaCoVSD.simple_conv(ι_sero,sero_array)

        PCR_pred = [nai_one_group.relative_testing_rate[t]*p_test*PCR[t]*1e-4 + 0.001 for t = 1:janendpoint]
        PCR_pred = vcat(PCR_pred,[p_test_boost*nai_one_group.relative_testing_rate[t]*p_test*PCR[t]*1e-4 + 0.001 for t = (janendpoint+1):length(PCR) ])
        PCR_prop = [(χ*PCR[t]/((χ-1)*PCR[t] + N)) + 0.001 for t = 1:janendpoint]
        PCR_prop = vcat(PCR_prop,[(χ_boost*χ*PCR[t]/((χ_boost*χ-1)*PCR[t] + N))  + 0.001 for t = (janendpoint+1):length(PCR)])

        infection_mat[:,j] .= ι_sero
        prop_PCR_pos_mat[:,j] .= PCR_prop
        no_neg_PCR_pos_mat[:,j] .= PCR_pred
        prop_sero_pos_mat[:,j] .=sero
    end
    return (;infection_mat,prop_PCR_pos_mat,no_neg_PCR_pos_mat,prop_sero_pos_mat)
end