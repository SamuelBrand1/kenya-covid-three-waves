function priors_twogroup_newvariant_cities(θ)
	@unpack ct_min1,R₀,ϵ,χ₁,χ₂,p_test₁,p_test₂,p_choose1,P_eff,schooleffect,extra_transmissibility,influx_exposed_new_variant,p_test_boost,χ_boost,E₀ = θ

    LP = 0.
    LP += logpdf(Beta(8,2),ct_min1)
    LP += logpdf(Gamma(3,2.5/3),R₀)
    LP += logpdf(Beta(45,5),ϵ)
    LP += logpdf(Gamma(10,4.5/10),χ₁)
    LP += logpdf(Gamma(10,1.5/10),χ₂)
    LP += logpdf(Gamma(3,1e-4/3),p_test₁)#In the region of 1% detection rate
    LP += logpdf(Gamma(3,5e-4/3),p_test₂)
    LP += logpdf(Beta(40,60),p_choose1)
	LP += logpdf(Beta(35,15),P_eff)
    LP += logpdf(Gamma(100,1.5/100),extra_transmissibility)
    LP += logpdf(Gamma(5,100/5),influx_exposed_new_variant)
	LP += logpdf(Normal(0,0.1),log(χ_boost))
	LP += logpdf(Normal(0,0.1),log(p_test_boost))
	LP += logpdf(Gamma(3,100/3),E₀)

    return LP
end

function priors_twogroup_newvariant_semi_urbanrural(θ)
	@unpack ct_min1,R₀,ϵ,χ₁,χ₂,p_test₁,p_test₂,p_choose1,P_eff,schooleffect,extra_transmissibility,influx_exposed_new_variant,p_test_boost,χ_boost,E₀ = θ

    LP = 0.
    LP += logpdf(Beta(8,2),ct_min1)
    LP += logpdf(Gamma(3,2/3),R₀)
    LP += logpdf(Beta(45,5),ϵ)
    LP += logpdf(Gamma(5,4.5/5),χ₁)
    LP += logpdf(Gamma(10,1.5/10),χ₂)
    LP += logpdf(Gamma(3,2e-5/3),p_test₁)
    LP += logpdf(Gamma(3,1e-4/3),p_test₂)
    LP += logpdf(Beta(80,20),p_choose1)
	LP += logpdf(Beta(95,5),P_eff)
    LP += logpdf(Gamma(100,1.5/100),extra_transmissibility)
    LP += logpdf(Gamma(5,100/5),influx_exposed_new_variant)
	LP += logpdf(Normal(0,0.1),log(χ_boost))
	LP += logpdf(Normal(0,0.1),log(p_test_boost))
	LP += logpdf(Gamma(3,1/3),E₀)

    return LP
end

function priors_twogroup_newvariant_rural(θ)
	@unpack ct_min1,R₀,ϵ,χ₁,χ₂,p_test₁,p_test₂,p_choose1,P_eff,schooleffect,extra_transmissibility,influx_exposed_new_variant,p_test_boost,χ_boost,E₀ = θ

    LP = 0.
    LP += logpdf(Beta(8,2),ct_min1)
    LP += logpdf(Gamma(3,1.5/3),R₀)
    LP += logpdf(Beta(45,5),ϵ)
    LP += logpdf(Gamma(3,4.5/3),χ₁)
    LP += logpdf(Gamma(10,1.5/10),χ₂)
    LP += logpdf(Gamma(3,1e-5/3),p_test₁)
    LP += logpdf(Gamma(3,5e-5/3),p_test₂)
    LP += logpdf(Beta(90,10),p_choose1)
	LP += logpdf(Beta(95,5),P_eff)
    LP += logpdf(Gamma(100,1.5/100),extra_transmissibility)
    LP += logpdf(Gamma(5,100/5),influx_exposed_new_variant)
	LP += logpdf(Normal(0,0.1),log(χ_boost))
	LP += logpdf(Normal(0,0.1),log(p_test_boost))
	LP += logpdf(Gamma(3,0.1/3),E₀)

    return LP
end

function priors_twogroup_newvariant_semi_urbanrural_fitted_detection(θ)
	@unpack ct_min1,R₀,ϵ,χ₁,χ₂,p_test₁,p_test₂,p_choose1,P_eff,schooleffect,extra_transmissibility,influx_exposed_new_variant,p_test_boost,χ_boost,E₀ = θ

    LP = 0.
    LP += logpdf(Beta(8,2),ct_min1)
    LP += logpdf(Gamma(3,2/3),R₀)
    LP += logpdf(Beta(45,5),ϵ)
    LP += logpdf(semi_urban_rural_fits.d_χ1,χ₁)
    LP += logpdf(semi_urban_rural_fits.d_χ2,χ₂)
    LP += logpdf(semi_urban_rural_fits.d_ptest1,p_test₁)
    LP += logpdf(semi_urban_rural_fits.d_ptest2,p_test₂)
    LP += logpdf(Beta(80,20),p_choose1)
	LP += logpdf(Beta(95,5),P_eff)
    LP += logpdf(Gamma(100,1.5/100),extra_transmissibility)
    LP += logpdf(Gamma(5,100/5),influx_exposed_new_variant)
	LP += logpdf(Normal(0,0.1),log(χ_boost))
	LP += logpdf(Normal(0,0.1),log(p_test_boost))
	LP += logpdf(Gamma(3,1/3),E₀)

    return LP
end


function priors_twogroup_newvariant_rural_fitted_detection(θ)
	@unpack ct_min1,R₀,ϵ,χ₁,χ₂,p_test₁,p_test₂,p_choose1,P_eff,schooleffect,extra_transmissibility,influx_exposed_new_variant,p_test_boost,χ_boost,E₀ = θ

    LP = 0.
    LP += logpdf(Beta(8,2),ct_min1)
    LP += logpdf(Gamma(3,1.5/3),R₀)
    LP += logpdf(Beta(45,5),ϵ)
    LP += logpdf(rural_fits.d_χ1,χ₁)
    LP += logpdf(rural_fits.d_χ2,χ₂)
    LP += logpdf(rural_fits.d_ptest1,p_test₁)
    LP += logpdf(rural_fits.d_ptest2,p_test₂)
    LP += logpdf(Beta(90,10),p_choose1)
	LP += logpdf(Beta(95,5),P_eff)
    LP += logpdf(Gamma(100,1.5/100),extra_transmissibility)
    LP += logpdf(Gamma(5,100/5),influx_exposed_new_variant)
	LP += logpdf(Normal(0,0.1),log(χ_boost))
	LP += logpdf(Normal(0,0.1),log(p_test_boost))
	LP += logpdf(Gamma(3,0.1/3),E₀)

    return LP
end
