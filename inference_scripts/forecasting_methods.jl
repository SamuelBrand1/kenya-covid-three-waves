function get_predictions_deltavariant(δ_var,model::KenyaCoVSD.CoVAreaModel,projection_date::Date,p_ID;seroreversionrate=0.0,σ = 0.16,ω = 1/180,ct_min2 = 0.445,fitlength = 400)
    @unpack PCR_cases,sero_cases,baseline_sero_array,PCR_array,sero_sensitivity,sero_specificity,N,M_BB,prob,α,γ,relative_testing_rate = model
    MCMCdraws = get(model.MCMC_results.chain,[:ct_min1,:R₀,:ϵ,:χ₁,:χ₂,:p_test₁,:p_test₂,:p_choose1,:P_eff,:schooleffect,:extra_transmissibility,:influx_exposed_new_variant,:p_test_boost,:χ_boost,:E₀])
    
    relative_testing_rate = [relative_testing_rate;ones(1000)]
    
    incidence₁ = zeros((projection_date - Date(2020,2,20)).value,length(first(MCMCdraws)))
    incidence₂ = zeros((projection_date - Date(2020,2,20)).value,length(first(MCMCdraws)))
    serocoverted₁ = zeros((projection_date - Date(2020,2,20)).value,length(first(MCMCdraws)))
    serocoverted₂ = zeros((projection_date - Date(2020,2,20)).value,length(first(MCMCdraws)))
    susceptible₁ = zeros((projection_date - Date(2020,2,20)).value,length(first(MCMCdraws)))
    susceptible₂ = zeros((projection_date - Date(2020,2,20)).value,length(first(MCMCdraws)))
    Rt₁ = zeros((projection_date - Date(2020,2,20)).value,length(first(MCMCdraws)))
    Rt₂ = zeros((projection_date - Date(2020,2,20)).value,length(first(MCMCdraws)))
    PCR_pred = zeros((projection_date - Date(2020,2,20)).value,length(first(MCMCdraws)))
    prop_PCR_pred = zeros((projection_date - Date(2020,2,20)).value,length(first(MCMCdraws)))
    forecast_testing_rate = zeros((projection_date - Date(2020,2,20)).value,length(first(MCMCdraws)))
    pred_deaths = zeros((projection_date - Date(2020,2,20)).value,length(first(MCMCdraws)))
    #Calculate prop_pos to test ratio
    test_fit = perc_pos_to_num_tests_scale(model)
    println("incidence fitting")
    for k = 1:length(first(MCMCdraws))
        N₁ = MCMCdraws.P_eff[k]*N
        N₂ = N - N₁
        p = [MCMCdraws.R₀[k],model.α,model.γ,MCMCdraws.ϵ[k],σ,N₁,N₂,ω,MCMCdraws.schooleffect[k],MCMCdraws.ct_min1[k],ct_min2]
        u0 = [N₁,MCMCdraws.E₀[k],0.,0.,0.,0.,N₂,MCMCdraws.E₀[k],0.,0.,0.,0.,0.0,0.0 ]

        function new_variant_effect!(integrator)
            integrator.p[1] *= MCMCdraws.extra_transmissibility[k]
            integrator.u[2] += MCMCdraws.influx_exposed_new_variant[k] #Lower SES
            integrator.u[8] += MCMCdraws.influx_exposed_new_variant[k] #Higher SES
        end

        function δ_variant_effect!(integrator)
            integrator.p[1] *= δ_var
        end

        janendpoint = (Date(2021,1,30) - Date(2020,2,20)).value
        maystartpoint = (Date(2021,5,1) - Date(2020,2,20)).value

        variant_cb = PresetTimeCallback([janendpoint],new_variant_effect!,save_positions=(false,false))
        δ_variant_cb = PresetTimeCallback([maystartpoint],δ_variant_effect!,save_positions=(false,false))
        cbs = CallbackSet(variant_cb,δ_variant_cb)
        #Sero-waning
        sero_array = vcat(baseline_sero_array[1:30],[(1-seroreversionrate)^k for k in 1:length(baseline_sero_array[31:end])])
        #solve model
        sol = solve(prob, BS3();tspan = (0,(projection_date - Date(2020,2,20)).value),
                        callback = cbs,
                        u0=u0, 
                        p=p, 
                        saveat = 1)
        ι₁ = diff(sol[:C₁])
        ι₂ = diff(sol[:C₂])
        ι_sero₁ = diff(sol[:C_sero₁])
        ι_sero₂ = diff(sol[:C_sero₂])

        PCR₁ = KenyaCoVSD.simple_conv(ι₁,PCR_array)
        PCR₂ = KenyaCoVSD.simple_conv(ι₂,PCR_array)
        incidence₁[:,k] = ι₁
        incidence₂[:,k] = ι₂
        susceptible₁[:,k] = sol[:S₁][2:end]
        susceptible₂[:,k] = sol[:S₂][2:end]
        serocoverted₁[:,k] .= sero_sensitivity.*KenyaCoVSD.simple_conv(ι_sero₁,sero_array)
        serocoverted₂[:,k] .= sero_sensitivity.*KenyaCoVSD.simple_conv(ι_sero₂,sero_array)

        #Calculate effective R
        Reff = [MCMCdraws.R₀[k]*(t < 24 || t ≥ 319) + MCMCdraws.R₀[k]*MCMCdraws.schooleffect[k]*(t≥ 24 && t< 319) for t = 1:length(ι₁)]
        new_variant_multiplier = [1.0*(t < janendpoint) + MCMCdraws.extra_transmissibility[k]*(t >= janendpoint) for t = 1:length(ι₁)]
        ct1 = [KenyaCoVSD.ct_kenya(t,MCMCdraws.ct_min1[k]) for t = 1:length(ι₁)]
        ct2 = [KenyaCoVSD.ct_kenya(t,ct_min2) for t = 1:length(ι₁)]
        Rt₁[:,k] .= new_variant_multiplier.*Reff.*(((ct1.*sol[:S₁][2:end]./N₁).*MCMCdraws.ϵ[k].*ct1) .+ ((ct2.*sol[:S₂][2:end]./N₂).*(1 .- MCMCdraws.ϵ[k]).*ct1))
        Rt₂[:,k] .= new_variant_multiplier.*Reff.*(((ct1.*sol[:S₁][2:end]./N₁).*(1 .- MCMCdraws.ϵ[k]).*ct2) .+ ((ct2.*sol[:S₂][2:end]./N₂).*MCMCdraws.ϵ[k].*ct2))

        #PCR test predictions
        new_variant_multiplier_detection_multiplier = [1.0*(t < janendpoint) + MCMCdraws.p_test_boost[k]*(t >= janendpoint) for t = 1:length(ι₁)]
        PCR_pred[:,k] = new_variant_multiplier_detection_multiplier.*relative_testing_rate[1:length(ι₁)].*(MCMCdraws.p_test₁[k].*PCR₁ .+ MCMCdraws.p_test₂[k].*PCR₂)


        prop_PCR_prenv = MCMCdraws.p_choose1[k].*(MCMCdraws.χ₁[k].*PCR₁[1:janendpoint]./((MCMCdraws.χ₁[k]-1).*PCR₁[1:janendpoint] .+ sum(N₁)))
        prop_PCR_prenv .+=  (1 .- MCMCdraws.p_choose1[k]).*(MCMCdraws.χ₂[k].*PCR₂[1:janendpoint]./((MCMCdraws.χ₂[k]-1).*PCR₂[1:janendpoint] .+ sum(N₂)))
        prop_PCR_postnv = MCMCdraws.p_choose1[k].*(MCMCdraws.χ_boost[k]*MCMCdraws.χ₁[k].*PCR₁[(janendpoint+1):end]./((MCMCdraws.χ_boost[k]*MCMCdraws.χ₁[k]-1).*PCR₁[(janendpoint+1):end] .+ sum(N₁)))
        prop_PCR_postnv .+=  (1 .- MCMCdraws.p_choose1[k]).*(MCMCdraws.χ_boost[k]*MCMCdraws.χ₂[k].*PCR₂[(janendpoint+1):end]./((MCMCdraws.χ_boost[k]*MCMCdraws.χ₂[k]-1).*PCR₂[(janendpoint+1):end].+ sum(N₂)))
        prop_PCR_pred[:,k] .= vcat(prop_PCR_prenv,prop_PCR_postnv)
        

        #PCR forecast accounting for test rates in last 60 days of period
        forecast_testing_rate[:,k] .= [ones(length(prop_PCR_pred[:,k])) prop_PCR_pred[:,k]]*test_fit.β
        
        

    end
                             

    return (incidence₁ = incidence₁,incidence₂ = incidence₂,
            serocoverted₁ = serocoverted₁,serocoverted₂ = serocoverted₂,
            susceptible₁ = susceptible₁,susceptible₂ = susceptible₂,
            Rt₁ = Rt₁,Rt₂ = Rt₂,
            PCR_pred = PCR_pred,prop_PCR_pred = prop_PCR_pred,
            forecast_testing_rate = forecast_testing_rate,
            pred_deaths = pred_deaths,
            test_fit=test_fit)
end