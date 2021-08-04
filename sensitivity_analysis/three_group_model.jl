"""
    function threegrouptransmission(du,u,p,t)

In-place creation of vector field for the three group transmission model.
"""
function threegrouptransmission(du,u,p,t)
    #Get state
	S₁,E₁,I₁,R₁,W₁,C₁,C_sero₁,S₂,E₂,I₂,R₂,W₂,C₂,C_sero₂,S₃,E₃,I₃,R₃,W₃,C₃,C_sero₃ = u
    #Get parameters
	R₀,α,γ,ϵ,σ,N₁,N₂,N₃,ω,schooleffect,ct_min1,ct_min2,ct_min3 = p
	#Adjust for schools being open or shut
    Reff = R₀*(t < 24 || t ≥ 319) + R₀*schooleffect*(t≥ 24 && t< 319)
	#Contact rates for each group
    ct1 = KenyaCoVSD.ct_kenya(t,ct_min1)
    ct2 = KenyaCoVSD.ct_kenya(t,ct_min2)
    ct3 = KenyaCoVSD.ct_kenya(t,ct_min3)
    #Force of infection on each group
    λ₁ = γ*Reff*(ct1*ϵ*I₁ + ct2*(1-ϵ)*I₂*(N₁/(N₁ + N₃)) + ct3*(1-ϵ)*I₃*(N₁/(N₁ + N₂)))
    λ₂ = γ*Reff*(ct1*(1-ϵ)*I₁*(N₂/(N₂ + N₃)) + ct2*ϵ*I₂ + ct3*(1-ϵ)*I₃*(N₂/(N₂ + N₁)))
    λ₃ = γ*Reff*(ct1*(1-ϵ)*I₁*(N₃/(N₂ + N₃)) + ct2*(1-ϵ)*I₂*(N₃/(N₃ + N₁)) + ct3*ϵ*I₃)
	#Dynamics
	du[1] = -ct1*(S₁/N₁)*λ₁
	du[2] = ct1*((S₁+σ*W₁)/N₁)*λ₁ - α*E₁
	du[3] = α*E₁ - γ*I₁
	du[4] = γ*I₁ - ω*R₁
	du[5] = ω*R₁ - σ*ct1*(W₁/N₁)*λ₁
	du[6] = ct1*((S₁+σ*W₁)/N₁)*λ₁
    du[7] = ct1*(S₁/N₁)*λ₁

	du[8] = -ct2*(S₂/N₂)*λ₂
	du[9] = ct2*((S₂+σ*W₂)/N₂)*λ₂ - α*E₂
	du[10] = α*E₂ - γ*I₂
	du[11] = γ*I₂ - ω*R₂
	du[12] = ω*R₂ - σ*ct2*(W₂/N₂)*λ₂
	du[13] = ct2*((S₂+σ*W₂)/N₂)*λ₂
    du[14] = ct2*(S₂/N₂)*λ₂

	du[15] = -ct3*(S₃/N₃)*λ₃
	du[16] = ct3*((S₃+σ*W₃)/N₃)*λ₃ - α*E₃
	du[17] = α*E₃ - γ*I₃
	du[18] = γ*I₃ - ω*R₃
	du[19] = ω*R₃ - σ*ct3*(W₃/N₃)*λ₃
	du[20] = ct3*((S₃+σ*W₃)/N₃)*λ₃
    du[21] = ct3*(S₃/N₃)*λ₃
	return nothing
end


u0_threegroup = [4.3e6*0.25 - 100.,100.,0.,0.,0.,0.,0.,
        4.3e6*0.5 - 100.,100.,0.,0.,0.,0.,0.,
        4.3e6*0.25 - 100.,100.,0.,0.,0.,0.,0.]
f_threegroup = ODEFunction(threegrouptransmission;syms = [:S₁,:E₁,:I₁,:R₁,:W₁,:C₁,:C_sero₁,
															:S₂,:E₂,:I₂,:R₂,:W₂,:C₂,:C_sero₂,
															:S₃,:E₃,:I₃,:R₃,:W₃,:C₃,:C_sero₃])
prob_three_group = ODEProblem(f_threegroup,u0_threegroup,(0,365))

"""
    function ll_onegroup_newvariant_infboost(θ,model::KenyaCoVSD.CoVAreaModel,seroreversionrate,cts;σ = 0.16,ω = 1/180,ct_min2 = 0.445)

Log-likelihood of the case and serology data in `model` assuming: 1) one group dynamics, and 2) that new variant is more transmissible rather than immune evading.
"""
function ll_threegroup_newvariant_infboost(θ,model::KenyaCoVSD.CoVAreaModel,seroreversionrate;σ = 0.16,ω = 1/180)
    @unpack PCR_cases,sero_cases,baseline_sero_array,PCR_array,sero_sensitivity,sero_specificity,N,M_BB,prob,α,γ,relative_testing_rate = model
    @unpack ct_min1,ct_min2,ct_min3,R₀,ϵ,χ₁,χ₂,p_test₁,p_test₂,p_choose1,P_eff,schooleffect,extra_transmissibility,influx_exposed_new_variant,p_test_boost,χ_boost,E₀ = θ

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
