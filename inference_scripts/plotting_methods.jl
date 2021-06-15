"""
    plot_Rt_both_SES_groups(model::KenyaCoVSD.CoVAreaModel,forecastdate::Date)

Plot R(t) for both SES groups from a fitted `CoVAreaModel` object.    
"""
function plot_Rt_both_SES_groups(model::KenyaCoVSD.CoVAreaModel,forecastdate::Date)
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
