"""
    function perc_pos_to_num_tests_scale(model)

Fit a linear realationship between the % PCR positive (7 day running average) and the number of PCR tests performed on each day (7-day running average),
over the last 60 days of available data. Returns a `NamedTuple` with the regression coefficients and the std over error.
"""
function perc_pos_to_num_tests_scale(model)
    prop_pos = model.PCR_cases[:,1]./model.PCR_cases[:,2]
    prop_pos_mv_av = [mean(prop_pos[(t-3):t+3]) for t = 4:(length(prop_pos)-3)]
    number_pos_mv_av = [mean(model.PCR_cases[:,1][(t-3):t+3]) for t = 4:(length(model.PCR_cases[:,1])-3)]
    number_tests_mv_av = [mean(model.PCR_cases[:,2][(t-3):t+3]) for t = 4:(length(model.PCR_cases[:,2])-3)]

    idxs = (.~isnan.(prop_pos_mv_av)).&(prop_pos_mv_av .> 0).&(1:length(prop_pos_mv_av) .>= length(prop_pos_mv_av) - 60).&(number_tests_mv_av .> 0)
    xs = prop_pos_mv_av[idxs];ys = number_tests_mv_av[idxs];
    β = [ones(length(xs)) xs] \ ys
    pred_ys = [ones(length(xs)) xs]*β
    std_ys = std(ys .- pred_ys)
    return (β = β, σ = std_ys )
end

"""
    function get_unscaled_predictions(ι₁,ι₂,p_ID,rel_test_rate;fitlength=430)

Useful method for collating a matrix of of unnormalised (i.e. IFR not applied) estimates of the daily death rate from the model incidence estimates.
"""
function get_unscaled_predictions(ι₁,ι₂,p_ID,rel_test_rate;fitlength=430)
    unscaled_deaths1 = similar(ι₁)[1:fitlength,:]
    unscaled_deaths2 = similar(ι₂)[1:fitlength,:]
    for j = 1:size(unscaled_deaths1,2)
        unscaled_deaths1[:,j] .= KenyaCoVSD.simple_conv(ι₁[:,j],p_ID)[1:fitlength].*rel_test_rate[1:fitlength]
        unscaled_deaths2[:,j] .= KenyaCoVSD.simple_conv(ι₂[:,j],p_ID)[1:fitlength].*rel_test_rate[1:fitlength]
    end
    return unscaled_deaths1,unscaled_deaths2
end

"""
    function get_predictions(model::KenyaCoVSD.CoVAreaModel,projection_date::Date;seroreversionrate=0.0,σ = 0.16,ω = 1/180,ct_min2 = 0.445,fitlength = 400)

Main method for gathering predictions from a fitted `CoVAreaModel`.
"""
function get_predictions(model::KenyaCoVSD.CoVAreaModel,projection_date::Date,p_ID;seroreversionrate=0.0,σ = 0.16,ω = 1/180,ct_min2 = 0.445,fitlength = 400)
    @unpack PCR_cases,sero_cases,baseline_sero_array,PCR_array,sero_sensitivity,sero_specificity,N,M_BB,prob,α,γ,relative_testing_rate = model
    MCMCdraws = get(model.MCMC_results.chain,[:ct_min1,:R₀,:ϵ,:χ₁,:χ₂,:p_test₁,:p_test₂,:p_choose1,:P_eff,:schooleffect,:extra_transmissibility,:influx_exposed_new_variant,:p_test_boost,:χ_boost,:E₀])

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

        janendpoint = (Date(2021,1,30) - Date(2020,2,20)).value
        variant_cb = PresetTimeCallback([janendpoint],new_variant_effect!,save_positions=(false,false))
        #Sero-waning
        sero_array = vcat(baseline_sero_array[1:30],[(1-seroreversionrate)^k for k in 1:600])
        #solve model
        sol = solve(prob, BS3();tspan = (0,(projection_date - Date(2020,2,20)).value),callback = variant_cb,u0=u0, p=p, saveat = 1)
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


function condense_prediction(prediction,num_tests,deaths,p_ID,relative_testing_rate)
    mean_incidence₁ = mean(prediction.incidence₁,dims = 2)[:]
    std_incidence₁ = std(prediction.incidence₁,dims = 2)[:]
    mean_incidence₂ = mean(prediction.incidence₂,dims = 2)[:]
    std_incidence₂ = std(prediction.incidence₂,dims = 2)[:]
    mean_serocoverted₁ = mean(prediction.serocoverted₁,dims = 2)[:]
    std_serocoverted₁ = std(prediction.serocoverted₁,dims = 2)[:]
    mean_serocoverted₂ = mean(prediction.serocoverted₂,dims = 2)[:]
    std_serocoverted₂ = std(prediction.serocoverted₂,dims = 2)[:]
    mean_susceptible₁ = mean(prediction.susceptible₁,dims = 2)[:]
    std_susceptible₁ = std(prediction.susceptible₁,dims = 2)[:]
    mean_susceptible₂ = mean(prediction.susceptible₂,dims = 2)[:]
    std_susceptible₂ = std(prediction.susceptible₂,dims = 2)[:]
    mean_Rt₁ = mean(prediction.Rt₁,dims = 2)[:]
    std_Rt₁ = std(prediction.Rt₁,dims = 2)[:]
    mean_Rt₂ = mean(prediction.Rt₂,dims = 2)[:]
    std_Rt₂ = std(prediction.Rt₂,dims = 2)[:]
    mean_PCR_pred_no_test_rate = mean(prediction.PCR_pred,dims = 2)[:]
    std_PCR_pred_no_test_rate = std(prediction.PCR_pred,dims = 2)[:]
    mean_prop_PCR_pred = mean(prediction.prop_PCR_pred,dims = 2)[:]
    std_prop_PCR_pred = std(prediction.prop_PCR_pred,dims = 2)[:]
   
    # mean_PCR_forecast = mean(prediction.prop_PCR_pred.*prediction.forecast_testing_rate,dims = 2)[:]
    # std_PCR_forecast = std(prediction.prop_PCR_pred.*prediction.forecast_testing_rate,dims = 2)[:]

    n = length(num_tests)
    idxs_neg_tests = num_tests .>= 0

    PCR_pred_mat = similar(prediction.prop_PCR_pred)

    for j = 1:size(PCR_pred_mat,2)
        past = prediction.PCR_pred[1:n,j].*(.~idxs_neg_tests) .+ prediction.prop_PCR_pred[1:n,j].*num_tests.*(idxs_neg_tests)
        PCR_pred_mat[:,j] .= vcat(past[1:(n-8)],
                                    prediction.prop_PCR_pred[(n-7):end,j].*prediction.forecast_testing_rate[(n-7):end,j])
    end

    mean_PCR_forecast = mean(PCR_pred_mat,dims = 2)
    std_PCR_forecast = std(PCR_pred_mat,dims = 2)

    #Fit deaths

    unscaled_deaths1,unscaled_deaths2 = get_unscaled_predictions(prediction.incidence₁,
                                        prediction.incidence₂,
                                        p_ID,
                                        relative_testing_rate;
                                        fitlength=size(prediction.incidence₁,1))
    
    mean_unscaled_deaths1 = mean(unscaled_deaths1,dims = 2)[:]    
    mean_unscaled_deaths2 = mean(unscaled_deaths2,dims = 2)[:] 
    println("fitting deaths")
    function neg_ll(μ)
        LL = 0.
        for (t,D) in enumerate(deaths) 
            LL += logpdf(Poisson(exp(μ[1])*mean_unscaled_deaths1[t] + exp(μ[2])*mean_unscaled_deaths2[t]),D)
        end
        return -LL
    end

    res = optimize(neg_ll,log.([0.0001,0.0001]),LBFGS())    
    μ_fit = exp.(res.minimizer)

    
    pred_deaths = μ_fit[1].*unscaled_deaths1 .+ μ_fit[2].*unscaled_deaths2 
    mean_deaths = mean(pred_deaths,dims = 2)[:]
    std_deaths = std(pred_deaths,dims = 2)[:]

    return (mean_incidence₁ = mean_incidence₁,
            std_incidence₁ = std_incidence₁,
            mean_incidence₂ = mean_incidence₂,
            std_incidence₂ = std_incidence₂,
            mean_serocoverted₁ = mean_serocoverted₁,
            std_serocoverted₁ = std_serocoverted₁,
            mean_serocoverted₂=mean_serocoverted₂,
            std_serocoverted₂=std_serocoverted₂,
            mean_susceptible₁=mean_susceptible₁,
            std_susceptible₁=std_susceptible₁,
            mean_susceptible₂=mean_susceptible₂,
            std_susceptible₂=std_susceptible₂,
            mean_Rt₁=mean_Rt₁,
            std_Rt₁=std_Rt₁,
            mean_Rt₂=mean_Rt₂,
            std_Rt₂=std_Rt₂,
            mean_PCR_pred_no_test_rate=mean_PCR_pred_no_test_rate,
            std_PCR_pred_no_test_rate = std_PCR_pred_no_test_rate,
            mean_prop_PCR_pred=mean_prop_PCR_pred,
            std_prop_PCR_pred=std_prop_PCR_pred,
            mean_PCR_forecast=mean_PCR_forecast,
            std_PCR_forecast=std_PCR_forecast,
            mean_deaths = mean_deaths,
            std_deaths=std_deaths,
            test_fit = prediction.test_fit,
            μ_fit = μ_fit)

end

function weekly_mv_av(x)
     return [mean(x[(t-3):(t+3)]) for t = 4:(length(x)-3)]
end


function get_Rt(model::KenyaCoVSD.CoVAreaModel,projection_date::Date;seroreversionrate=1/365,σ = 0.16,ω = 1/180,ct_min2 = 0.445)
    @unpack PCR_cases,sero_cases,baseline_sero_array,PCR_array,sero_sensitivity,sero_specificity,N,M_BB,prob,α,γ,relative_testing_rate = model
    MCMCdraws = get(model.MCMC_results.chain,[:ct_min1,:R₀,:ϵ,:χ₁,:χ₂,:p_test₁,:p_test₂,:p_choose1,:P_eff,:schooleffect,:extra_transmissibility,:influx_exposed_new_variant,:p_test_boost,:χ_boost,:E₀])

    Rt₁ = zeros((projection_date - Date(2020,2,20)).value+1,length(first(MCMCdraws)))
    Rt₂ = zeros((projection_date - Date(2020,2,20)).value+1,length(first(MCMCdraws)))

    #Calculate prop_pos to test ratio

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

        janendpoint = (Date(2021,1,30) - Date(2020,2,20)).value
        variant_cb = PresetTimeCallback([janendpoint],new_variant_effect!)
        #solve model
        sol = solve(prob, BS3();tspan = (0,(projection_date - Date(2020,2,20)).value),callback = variant_cb,u0=u0, p=p, saveat = 1)
        #Calculate effective R
        Reff = [MCMCdraws.R₀[k]*(t < 24 || t ≥ 319) + MCMCdraws.R₀[k]*MCMCdraws.schooleffect[k]*(t≥ 24 && t< 319) for t = 1:length(sol[:S₁][2:end])]
        new_variant_multiplier = [1.0*(t < janendpoint) + MCMCdraws.extra_transmissibility[k]*(t >= janendpoint) for t = 1:length(sol[:S₁][2:end])]
        ct1 = [KenyaCoVSD.ct_kenya(t,MCMCdraws.ct_min1[k]) for t = 1:length(sol[:S₁][2:end])]
        ct2 = [KenyaCoVSD.ct_kenya(t,ct_min2) for t = 1:length(sol[:S₁][2:end])]
        Rt₁[:,k] .= new_variant_multiplier.*Reff.*(((ct1.*sol[:S₁][2:end]./N₁).*MCMCdraws.ϵ[k].*ct1) .+ ((ct2.*sol[:S₂][2:end]./N₂).*(1 .- MCMCdraws.ϵ[k]).*ct1))
        Rt₂[:,k] .= new_variant_multiplier.*Reff.*(((ct1.*sol[:S₁][2:end]./N₁).*(1 .- MCMCdraws.ϵ[k]).*ct2) .+ ((ct2.*sol[:S₂][2:end]./N₂).*MCMCdraws.ϵ[k].*ct2))
        
    end


    return (Rt₁_fit = get_credible_intervals(Rt₁),Rt₂_fit =get_credible_intervals(Rt₂))
end

function get_credible_intervals(X)
    pred = mean(X,dims = 2)
    lb =  pred .- [quantile(X[t,:],0.025) for t = 1:size(X,1)]
    ub = [quantile(X[t,:],0.975) for t = 1:size(X,1)] .- pred
    return (pred=pred,lb=lb,ub=ub)
end

function dailytestrate_to_probofdetections(p_test,PCR_array)
        1 - prod(1 .- p_test.*PCR_array)
end
