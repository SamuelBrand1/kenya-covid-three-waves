using JLD2,MCMCChains
import KenyaCoVSD
using StatsPlots,Distributions,Optim

## Group fitted models by semiurban/rural and rural

Semi_urban_filenames = ["modelfits/Nakuru_model.jld2",
                        "modelfits/Uasin Gishu_model.jld2",
                        "modelfits/Embu_model.jld2",
                        "modelfits/Kisumu_model.jld2",
                        "modelfits/Nyeri_model.jld2",
                        "modelfits/Siaya_model.jld2",
                         "modelfits/Kisii_model.jld2"]

rural_filenames = ["modelfits/Kilifi_model.jld2",
                        "modelfits/Kwale_model.jld2"]

## Combine posterior distribution draws across groups
fit_dict = load(Semi_urban_filenames[1])
model = fit_dict[first(keys(fit_dict))]
semi_urban_chain = model.MCMC_results.chain

for filename in Semi_urban_filenames[2:end]
    fit_dict = load(filename)
    model = fit_dict[first(keys(fit_dict))]
    semi_urban_chain = vcat(semi_urban_chain,model.MCMC_results.chain)
end


fit_dict = load(rural_filenames[1])
model = fit_dict[first(keys(fit_dict))]
rural_chain = model.MCMC_results.chain
for filename in rural_filenames[2:end]
    fit_dict = load(filename)
    model = fit_dict[first(keys(fit_dict))]
    rural_chain = vcat(rural_chain,model.MCMC_results.chain)
end


## Define a MLE method for creating a univariate fit for detection rate parameters across groups
"""
    fit_univariate_distributions_for_detection(post)

Fit a collection of univariate Gamma distributions to posterior samples of detection rate parameters.
"""
function fit_univariate_distributions_for_detection(post)
    v = get(post,[:p_test₁,:p_test₂,:χ₁,:χ₂])
    d_ptest1 = fit_mle(Gamma,v.p_test₁)
    d_ptest2 = fit_mle(Gamma,v.p_test₂)
    d_χ1 = fit_mle(Gamma,v.χ₁)
    d_χ2 = fit_mle(Gamma,v.χ₂)

    return (d_ptest1 = d_ptest1,
            d_ptest2= d_ptest2,
            d_χ1 = d_χ1,
            d_χ2 = d_χ2)
end



semi_urban_rural_fits = fit_univariate_distributions_for_detection(semi_urban_chain)
rural_fits = fit_univariate_distributions_for_detection(rural_chain)

@save("data/semi_urban_rural_detection_fits.jld2",semi_urban_rural_fits)
@save("data/rural_fits.jld2",rural_fits)

