"""
    plot_Rt_both_SES_groups(model::KenyaCoVSD.CoVAreaModel,forecastdate::Date)

Plot R(t) for both SES groups from a fitted `CoVAreaModel` object.    
"""
function plot_Rt_both_SES_groups(model::KenyaCoVSD.CoVAreaModel,forecastdate::Date)
        june1day = (Date(2021,6,1) - Date(2020,2,24)).value

        xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,24)).value for k = 1:18 ]
        xticklabs = [monthname(k)[1:3]*"/20" for k = 3:12]
        xticklabs = vcat(xticklabs,[monthname(k)[1:3]*"/21" for k = 1:8])
        qs = quantile(model.MCMC_results.chain)
        ms = mean(model.MCMC_results.chain)
        idx = findfirst(ms[:,1].== :P_eff)
        mean_Peff_perc = round(ms[:,2][idx]*100,digits = 1)
        lower_CI = round(qs[:,2][idx]*100,digits = 1)
        higher_CI = round(qs[:,end][idx]*100,digits = 1)

        Rt_fits = get_Rt(model,forecastdate)
        plt = plot(Rt_fits.Rt₁_fit.pred,
                lab = "Lower SES group, $(mean_Peff_perc)% of county ($(lower_CI)% - $(higher_CI)%)",
                ribbon = (Rt_fits.Rt₁_fit.lb,Rt_fits.Rt₁_fit.ub),
                xticks = (xticktimes[1:2:(end-2)],xticklabs[1:2:(end-2)]),
                title = "Effective reproductive number by SES group: $(model.areaname)",
                size = (700,500),dpi = 250,
                xlims = (-5,june1day),
                ylabel = "R(t)",
                legendfont = 10,titlefont = 15,xtickfontsize=10,ytickfontsize=13,guidefont = 18,
                left_margin = 10mm,right_margin = 7.5mm)

        plot!(plt,Rt_fits.Rt₂_fit.pred,
                ribbon = (Rt_fits.Rt₂_fit.lb,Rt_fits.Rt₂_fit.ub),
                lab = "Higher SES group")
        return plt
end

function plot_PCR(fit,linelist_data_with_pos_neg,linelist_data)
        june1day = (Date(2021,6,1) - Date(2020,2,24)).value

        xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,24)).value for k = 1:18 ]
        xticklabs = [monthname(k)[1:3]*"/20" for k = 3:12]
        xticklabs = vcat(xticklabs,[monthname(k)[1:3]*"/21" for k = 1:8])

        n = size(fit.pred.mean_PCR_forecast,1)
        n_1 = size(linelist_data_with_pos_neg.cases,1) - 14
        county_pcr_forecast = fit.pred.mean_PCR_forecast[:]
        std_county_pcr_forecast =fit.pred.std_PCR_forecast[:]
        county_pcr_forecast_mv_av = weekly_mv_av(county_pcr_forecast)
        kenya_pcr_forecast_mv_av_std = sqrt.(weekly_mv_av(std_county_pcr_forecast.^2))
        uprname = uppercase(fit.name)
        county_pos = sum(linelist_data_with_pos_neg.cases[:,linelist_data_with_pos_neg.areas .== uprname,
                :,1],dims = [2,3])[:]
        county_pos_mv_av = weekly_mv_av(county_pos)

        PCR_plt = scatter(county_pos[1:(end)],
                ms = 4,markerstrokewidth = 0,color = :grey,alpha = 0.5,
                xticks = (xticktimes,xticklabs),
                lab = "Daily cases: Kenyan linelist (used for fitting)",legend = :topleft,
                title = "$(fit.name) PCR test positives",
                xlims = (-5,june1day),
                ylabel = "Daily PCR-confirmed cases",
                size = (1100,500),dpi = 250,
                legendfont = 13,titlefont = 24,tickfontsize=10,guidefont = 18,
                left_margin = 10mm,right_margin = 7.5mm)

        plot!(PCR_plt,(1+3):(length(county_pos_mv_av)+3),county_pos_mv_av,
                color = :black,lw = 3,
                lab = "Daily cases: 7 day mv-av (used in fitting)")

        smooth_cases_lookahead = weekly_mv_av(sum(linelist_data.cases[(n_1-2):end,linelist_data.areas .== uprname],dims = 2)[:])

        plot!(PCR_plt,(length(county_pos_mv_av)+4):(length(county_pos_mv_av)+3+length(smooth_cases_lookahead)),smooth_cases_lookahead,
                color = :black,lw = 3,ls = :dot,
                lab = "Daily cases: 7 day mv-av (not used in fitting)")

        plot!(PCR_plt,4:(n-3),county_pcr_forecast_mv_av,ribbon = min.(3*kenya_pcr_forecast_mv_av_std,county_pcr_forecast_mv_av),
                color = :red, lw = 2,lab = "Model fit and forecast (7 day mv-av)")

        return PCR_plt

end
        
function plot_deaths(fit,linelist_data_with_pos_neg)
        june1day = (Date(2021,6,1) - Date(2020,2,24)).value

        n = size(fit.pred.mean_PCR_forecast,1)
        n_1 = size(linelist_data_with_pos_neg.cases,1) - 14

        uprname = uppercase(fit.name)
        county_deaths_forecast = fit.pred.mean_deaths[:]
        std_county_deaths_forecast = fit.pred.std_deaths[:]
        county_deaths_forecast_mv_av = weekly_mv_av(county_deaths_forecast)
        county_deaths_forecast_mv_av_std = sqrt.(weekly_mv_av(std_county_deaths_forecast.^2))
        county_deaths = sum(deaths_data.deaths[:,deaths_data.areas .== uprname],dims = 2)[:]
        county_deaths_mv_av = weekly_mv_av(county_deaths)
        xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,24)).value for k = 1:18 ]
        xticklabs = [monthname(k)[1:3]*"/20" for k = 3:12]
        xticklabs = vcat(xticklabs,[monthname(k)[1:3]*"/21" for k = 1:8])

        deaths_plt = scatter(county_deaths,
                ms = 4,markerstrokewidth = 0,color = :grey,alpha = 0.5,
                xticks = (xticktimes,xticklabs),
                lab = "Daily deaths: Kenyan linelist",legend = :topleft,
                title = "$(fit.name) PCR-confirmed deaths",
                xlims = (-5,june1day),
                ylabel = "Daily PCR-confirmed deaths",
                size = (1100,500),dpi = 250,
                legendfont = 13,titlefont = 24,xtickfontsize=10,ytickfontsize=13,guidefont = 18,
                left_margin = 10mm,right_margin = 7.5mm)
        plot!(deaths_plt,(1+3):(length(county_deaths_mv_av)+3),county_deaths_mv_av,
                color = :black,lw = 3,
                lab = "Daily deaths: 7 day mv-av")
        plot!(deaths_plt,4:(n-3),county_deaths_forecast_mv_av,ribbon = 3.0.*county_deaths_forecast_mv_av_std,
                color = :green, lw = 5, ls = :dot,lab = "Model fit and forecast (7 day mv-av)")

        plot!(deaths_plt,(1+3):(length(county_deaths_mv_av)+3),cumsum(county_deaths_mv_av),
                xticks = (xticktimes[4:4:end],xticklabs[4:4:end]),
                xlims = (-5,june1day),
                color = :black,lw = 3,
                grid = nothing,
                inset = (1,bbox(0.35, -0.25, 0.25, 0.25, :center)),
                lab="",
                subplot = 2,
                title = "Cumulative confirmed deaths",
                bg_inside = nothing)
        plot!(deaths_plt,4:(n-3),cumsum(county_deaths_forecast_mv_av),
                xlims = (-5,june1day),color = :green, lw = 5, ls = :dot,lab = "",subplot=2)
        return deaths_plt
end

function plot_pop_exposure(fit,serological_data,serology_data,N_kenya)
        june1day = (Date(2021,6,1) - Date(2020,2,24)).value

        n = size(fit.pred.mean_PCR_forecast,1)
        n_1 = size(linelist_data_with_pos_neg.cases,1) - 14

        uprname = uppercase(fit.name)

        xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,24)).value for k = 1:18 ]
        xticklabs = [monthname(k)[1:3]*"/20" for k = 3:12]
        xticklabs = vcat(xticklabs,[monthname(k)[1:3]*"/21" for k = 1:8])

        zeropadtomonday = dayofweek(Date(2020,2,20)) - 1
        county_sero_pos_rnd1_2 = vcat(zeros(zeropadtomonday),
                                sum(serological_data.serodata[:,serological_data.areas .== uprname,:,1],dims = [2,3])[:])
        county_sero_neg_rnd1_2 = vcat(zeros(zeropadtomonday),
                                sum(serological_data.serodata[:,serological_data.areas .== uprname,:,2],dims = [2,3])[:])
        county_weekly_sero_pos_rnd1_2 = [sum(grp) for grp in Iterators.partition(county_sero_pos_rnd1_2,7)]
        county_weekly_sero_total_rnd1_2 = [sum(grp) for grp in Iterators.partition(county_sero_pos_rnd1_2.+county_sero_neg_rnd1_2,7)]

        county_sero_pos_rnd3 = vcat(zeros(zeropadtomonday),
                                sum(serology_data.sero[:,serology_data.areas .== uprname,:,1],dims = [2,3])[:])
        county_sero_neg_rnd3 = vcat(zeros(zeropadtomonday),
                                sum(serology_data.sero[:,serology_data.areas .== uprname,:,2],dims = [2,3])[:])
        county_weekly_sero_pos_rnd3 = [sum(grp) for grp in Iterators.partition(county_sero_pos_rnd3,7)]
        county_weekly_sero_total_rnd3 = [sum(grp) for grp in Iterators.partition(county_sero_pos_rnd3.+county_sero_neg_rnd3,7)]

        #Jeffery intervals
        seroidxs = county_weekly_sero_total_rnd3 .> 0
        # rnd_1_2idx =
        uerr = [invlogcdf(Beta(pos + 0.5,county_weekly_sero_total_rnd3[k] - pos + 0.5),log(0.975)) - pos/county_weekly_sero_total_rnd3[k] for (k,pos) in enumerate(county_weekly_sero_pos_rnd3) ]
        lerr = [pos/county_weekly_sero_total_rnd3[k] - invlogcdf(Beta(pos + 0.5,county_weekly_sero_total_rnd3[k] - pos + 0.5),log(0.025)) for (k,pos) in enumerate(county_weekly_sero_pos_rnd3) ]

        xs_mondays = [-3 + (k-1)*7 for k = 1:length(county_weekly_sero_total_rnd3)]
        rnd1_2_idxs = xs_mondays .< 300
        rnd3_idxs = .~rnd1_2_idxs

        county_sero_pos_nw = fit.pred.mean_serocoverted₁ .+ fit.pred.mean_serocoverted₂
        N = sum(N_kenya[:,fit.name])
        county_serology_forecast_nw = county_sero_pos_nw./N
        std_kenya_serology_forecast = sqrt.(((1.0/sum(N_kenya[:,fit.name]))^2 ).*(fit.pred.std_serocoverted₁.^2 .+ fit.pred.std_serocoverted₂.^2))
        county_infections_forecast = cumsum(fit.pred.mean_incidence₁[1:524] .+ fit.pred.mean_incidence₂[1:524])
        std_county_infections_forecast = sqrt.(cumsum((fit.pred.std_incidence₁[1:524] .+ fit.pred.std_incidence₂[1:524]).^2))

        plt_sero = scatter(xs_mondays[seroidxs.*rnd1_2_idxs],county_weekly_sero_pos_rnd3[seroidxs.*rnd1_2_idxs]./county_weekly_sero_total_rnd3[seroidxs.*rnd1_2_idxs],
                lab = "Weekly KNBTS: rounds 1 and 2 (used in fitting)",
                legend = :topleft,
                yerr = (lerr[seroidxs.*rnd1_2_idxs],uerr[seroidxs.*rnd1_2_idxs]),
                xticks = (xticktimes,xticklabs),
                title = "$(fit.name) population exposure",
                size = (1100,500),dpi = 250,
                xlims = xlims = (-5,june1day), ylims = (-0.025,1),
                ylabel = "Proportion of population",
                legendfont = 10,titlefont = 24,xtickfontsize=10,ytickfontsize=13,guidefont = 18,
                left_margin = 10mm,right_margin = 7.5mm)
        scatter!(plt_sero,xs_mondays[seroidxs.*rnd3_idxs],county_weekly_sero_pos_rnd3[seroidxs.*rnd3_idxs]./county_weekly_sero_total_rnd3[seroidxs.*rnd3_idxs],
                yerr = (lerr[seroidxs.*rnd3_idxs],uerr[seroidxs.*rnd3_idxs]),
                lab = "Weekly KNBTS: round 3 (not used in fitting)")
        plot!(plt_sero,county_serology_forecast_nw,lw = 2,color = :green,
                ribbon = min.(3*std_kenya_serology_forecast,county_serology_forecast_nw),
                lab = "Model fit: seroposivity")

        plot!(plt_sero,county_infections_forecast./N,
                ribbon = min.(9*std_county_infections_forecast./N,county_infections_forecast./N),
                lab = "Model fit: Overall population exposure",
                color = :red)

        return plt_sero
end

function plot_group_incidence(fit)
        june1day = (Date(2021,6,1) - Date(2020,2,24)).value

        xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,24)).value for k = 1:18 ]
        xticklabs = [monthname(k)[1:3]*"/20" for k = 3:12]
        xticklabs = vcat(xticklabs,[monthname(k)[1:3]*"/21" for k = 1:8])
        county_group1_incidence = fit.pred.mean_incidence₁
        std_county_group1_incidence = fit.pred.std_incidence₁
        county_group2_incidence = fit.pred.mean_incidence₂
        std_county_group2_incidence = fit.pred.std_incidence₂

        deleteat!(county_group1_incidence,county_group1_incidence.==0)
        deleteat!(std_county_group1_incidence,std_county_group1_incidence.==0)
        deleteat!(county_group2_incidence,county_group2_incidence.==0)
        deleteat!(std_county_group2_incidence,std_county_group2_incidence.==0)

        if length(std_county_group2_incidence) < length(county_group2_incidence)
                m = length(county_group2_incidence) - length(std_county_group2_incidence)
                std_county_group2_incidence = [std_county_group2_incidence;fill(0.0,m)]
        end

         if length(std_county_group2_incidence) > length(county_group2_incidence)
                m = length(county_group2_incidence) - length(std_county_group2_incidence)
                county_group2_incidence = [county_group2_incidence;fill(0.0,abs(m))]
        end

        plt_inc = plot(county_group1_incidence./1e5,
                lab = "Daily incidence: Lower SES",
                legend = :topleft,
                ribbon = min.(3*std_county_group1_incidence./1e5,county_group1_incidence./1e5),
                xticks = (xticktimes,xticklabs),
                title = "$(fit.name) transmission rates by SES group",
                size = (1100,500),dpi = 250,
                xlims = (-5,june1day),
                ylabel = "Daily infections (100,000s)",
                legendfont = 13,titlefont = 24,xtickfontsize=10,ytickfontsize=13,guidefont = 18,
                left_margin = 10mm,right_margin = 7.5mm)

        plot!(plt_inc,county_group2_incidence./1e5,
                lab = "Daily incidence: Higher SES",
                legend = :topleft,
                ribbon = min.(3*std_county_group2_incidence./1e5,county_group2_incidence./1e5) )
        return plt_inc
end