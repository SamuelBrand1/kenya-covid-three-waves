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
    function ll_threegroup_newvariant_infboost(θ,model::KenyaCoVSD.CoVAreaModel,seroreversionrate,cts;σ = 0.16,ω = 1/180,ct_min2 = 0.445)

Log-likelihood of the case and serology data in `model` assuming: 1) Three group dynamics, and 2) that new variant is more transmissible rather than immune evading.
"""
function ll_threegroup_newvariant_infboost(θ,model::KenyaCoVSD.CoVAreaModel,seroreversionrate;σ = 0.16,ω = 1/180)
    @unpack PCR_cases,sero_cases,baseline_sero_array,PCR_array,sero_sensitivity,sero_specificity,N,M_BB,prob,α,γ,relative_testing_rate = model
    @unpack R₀,ϵ,schooleffect,ct_min1,ct_min2,ct_min3,χ₁,χ₂,χ₃,p_test1,p_test2,p_test3,z_pch1,z_pch2,z_pch3,extra_transmissibility,influx_exposed_new_variant,p_test_boost,χ_boost,E₀ = θ

    N₁ = 0.25*N
    N₂ = 0.5*N
    N₃ = 0.25*N
    pch1,pch2,pch3 = softmax([z_pch1,z_pch2,z_pch3])
    
    #Set variance scalers
    clustering_factor_PCR = 0.5
    M_PCR = 30.
    T = eltype(R₀)
	# if T <: Real
	# 	T = eltype(cts)
	# end
    #Set transmission parameters with daily cts appended on end
    p = convert.(T,[R₀,α,γ,ϵ,σ,N₁,N₂,N₃,ω,schooleffect,ct_min1,ct_min2,ct_min3])
    u0 = convert.(T,[N₁-E₀,E₀,0.,0.,0.,0.,0.,
                    N₂-E₀,E₀,0.,0.,0.,0.,0.,
                    N₃-E₀,E₀,0.,0.,0.,0.,0.])

    #Create new variant introduction
    function new_variant_effect!(integrator)
        integrator.p[1] *= extra_transmissibility
        integrator.u[1] -= influx_exposed_new_variant
        integrator.u[2] += influx_exposed_new_variant
        integrator.u[1 + 7] -= influx_exposed_new_variant
        integrator.u[2 + 7] += influx_exposed_new_variant
        integrator.u[1 + 14] -= influx_exposed_new_variant
        integrator.u[2 + 14] += influx_exposed_new_variant
    end

    janendpoint = (Date(2021,1,30) - Date(2020,2,20)).value
	aprilendpoint = (Date(2021,4,30) - Date(2020,2,20)).value
    variant_cb = PresetTimeCallback([janendpoint],new_variant_effect!)

    #Sero-waning
    sero_array = vcat(baseline_sero_array[1:30],[(1-seroreversionrate)^k for k in 1:length(baseline_sero_array[31:end])])


    #Solve for daily incidence by age (not currently broken down by severity of disease)
    LL = T(0.)
    # try
        sol = solve(prob, BS3();tspan = (0,aprilendpoint),reltol = 1e-3,
					callback = variant_cb,
					u0=u0,
					p=p,
					saveat = 1,
					verbose = false,
                    isoutofdomain=(u,p,t) -> any(x -> x < 0, u))

        ι₁ = diff(sol[:C₁])
        ι₂ = diff(sol[:C₂])
        ι₃ = diff(sol[:C₃])
        ι_sero₁ = diff(sol[:C_sero₁])
        ι_sero₂ = diff(sol[:C_sero₂])
        ι_sero₃ = diff(sol[:C_sero₃])

        PCR₁ = KenyaCoVSD.simple_conv(ι₁,PCR_array)
        PCR₂ = KenyaCoVSD.simple_conv(ι₂,PCR_array)
        PCR₃ = KenyaCoVSD.simple_conv(ι₃,PCR_array)
        weighted_PCR_ptest = (p_test1.*PCR₁ .+ p_test2.*PCR₂ .+ p_test3.*PCR₃).*1e-4 
        weighted_prop_PCR = pch1.*(χ₁.*PCR₁./((χ₁-1).*PCR₁ .+ N₁)) .+ pch2.*(χ₂.*PCR₂./((χ₂-1).*PCR₂ .+ N₂)) .+ pch3.*(χ₃.*PCR₃./((χ₃-1).*PCR₃ .+ N₃))
        boosted_weighted_prop_PCR = pch1.*(χ_boost.*χ₁.*PCR₁./((χ_boost.*χ₁-1).*PCR₁ .+ N₁)) .+ pch2.*(χ_boost.*χ₂.*PCR₂./((χ_boost.*χ₂-1).*PCR₂ .+ N₂)) .+ pch3.*(χ_boost.*χ₃.*PCR₃./((χ_boost.*χ₃-1).*PCR₃ .+ N₃))

    	sero = sero_sensitivity.*KenyaCoVSD.simple_conv(ι_sero₁.+ι_sero₂.+ι_sero₃,sero_array)

        #Calculate log-likelihood for PCR testing
        for t in 55:janendpoint
            #Convert from μ,α parameterisation to p,r parameterisation for using Neg. binomial model
            μ = relative_testing_rate[t]*weighted_PCR_ptest[t] + 0.001 #Covert approximately from total probability of detection in % to daily probability of detection
            σ² = μ + clustering_factor_PCR*μ^2
            p_negbin = 1 - (clustering_factor_PCR*μ^2/σ²)
            r_negbin  = 1/clustering_factor_PCR
            LL += logpdf(NegativeBinomial(r_negbin,p_negbin),PCR_cases[t,1])#likelihood contribution from PCR testing --- positive case detection

            if PCR_cases[t,2] >= 0 #Negative tests available
                p_PCR_pred = weighted_prop_PCR[t] + 0.001
                #Convert from p_PCR,inv_M_PCR parameterisation of Beta-Binomial to standard α,β parameterisation
                LL += logpdf(BetaBinomial(PCR_cases[t,2],p_PCR_pred*M_PCR,(1-p_PCR_pred)*M_PCR),PCR_cases[t,1])#likelihood contribution from PCR testing --- proportion postive
            end
        end

        #Calculate log-likelihood for PCR testing after new variant introduction
        for t in (janendpoint+1):min(size(PCR_cases,1),size(PCR₁,1))
            #Convert from μ,α parameterisation to p,r parameterisation for using Neg. binomial model
            μ = p_test_boost*relative_testing_rate[t]*weighted_PCR_ptest[t] + 0.001 #Covert approximately from total probability of detection in % to daily probability of detection
            σ² = μ + clustering_factor_PCR*μ^2
            p_negbin = 1 - (clustering_factor_PCR*μ^2/σ²)
            r_negbin  = 1/clustering_factor_PCR
            LL += logpdf(NegativeBinomial(r_negbin,p_negbin),PCR_cases[t,1])#likelihood contribution from PCR testing --- positive case detection

            if PCR_cases[t,2] >= 0 #Negative tests available
                p_PCR_pred = boosted_weighted_prop_PCR[t] + 0.001
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
    # catch
    #     LL = T(-Inf)
    # end
    return LL::T
end

function threegroup_priors(θ)
	@unpack R₀,ϵ,schooleffect,ct_min1,ct_min2,ct_min3,χ₁,χ₂,χ₃,p_test1,p_test2,p_test3,z_pch1,z_pch2,z_pch3,extra_transmissibility,influx_exposed_new_variant,p_test_boost,χ_boost,E₀ = θ
    LP = 0.
    LP += logpdf(Gamma(3,2.5/3),R₀)
    LP += logpdf(Beta(7,3),ϵ)
    LP += logpdf(Beta(8,2),schooleffect)
    LP += logpdf(Beta(7,3),ct_min1)
    LP += logpdf(Beta(7,3),ct_min2)
    LP += logpdf(Beta(55,45),ct_min3)
    LP += logpdf(Gamma(3,4/3),χ₁)
    LP += logpdf(Gamma(3,4/3),χ₂)
    LP += logpdf(Gamma(3,1.5/3),χ₃)
    LP += logpdf(Gamma(3,1/3),p_test1)#In the region of 1% detection rate
    LP += logpdf(Gamma(3,1/3),p_test2)#In the region of 1% detection rate
    LP += logpdf(Gamma(3,5/3),p_test3)#In the region of 5% detection rate
    LP += logpdf(Normal(0,1.0),z_pch1)
	LP += logpdf(Normal(0,1.0),z_pch2)
	LP += logpdf(Normal(0,1.0),z_pch3)
    LP += logpdf(Gamma(10,1.5/10),extra_transmissibility)
    LP += logpdf(Gamma(5,100/5),influx_exposed_new_variant)
	LP += logpdf(Normal(0,0.1),log(χ_boost))
	LP += logpdf(Normal(0,0.1),log(p_test_boost))
	LP += logpdf(Gamma(3,100/3),E₀)
    return LP
end


## Create three group model
nai_three_group = deepcopy(basic_nai_model)
nai_three_group.prob = prob_three_group
nai_three_group.log_likelihood = ll_threegroup_newvariant_infboost
nai_three_group.log_priors = threegroup_priors


#Define variable transformation for three group model
trans_three_groups = as((R₀ = as(Real, 0.0, 7.0),
        ϵ = as(Real,0.5,0.999),
        schooleffect = as(Real,0.25,1.0),
        ct_min1 = as(Real, 0.1, 0.999),
        ct_min2 = as(Real, 0.1, 0.999),
        ct_min3 = as(Real, 0.1, 0.999),
        χ₁=as(Real, 0.1, 30.0),χ₂=as(Real, 0.1, 30.0),χ₃=as(Real, 0.1, 30.0),
        p_test1=as(Real, 0.0, 50.0),p_test2=as(Real, 0.0, 50.0),p_test3=as(Real, 0.0, 50.0),
        z_pch1 = asℝ,z_pch2 = asℝ,z_pch3 = asℝ,
        extra_transmissibility = as(Real, 0.5, 2.5),
        influx_exposed_new_variant = as(Real, 0.0, 10e3),
        p_test_boost = as(Real, 0.5, 3.0),
        χ_boost= as(Real, 0.5, 3.0),
        E₀=as(Real, 0.0, 10e3)))

D = Diagonal(0.1*ones(TransformVariables.dimension(trans_three_groups)))
