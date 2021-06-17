using StatsPlots,Dates,JLD2,Statistics,Optim,Parameters,Distributions,DataFrames,CSV
using Plots.PlotMeasures,OrdinaryDiffEq,DiffEqCallbacks
import KenyaCoVSD
include("fitting_methods.jl");
include("plotting_methods.jl");

## Load the collated fits

@load("data/N_kenya.jld2")
@load("data/linelist_data_with_pos_neg_20feb_to_27apr_c_age.jld2")
@load("data/serological_data_with_20feb_to_31dec_age.jld2")
@load("data/serological_data_with_20feb_to_10Mar2021_age.jld2")
@load("data/cleaned_linelist20210521_deaths_c__date_of_lab_confirmation.jld2")
@load("data/p_ID.jld2")
@load("forecasts/condensed_county_forecasts.jld2")
@load("data/rel_sero_detection_after_infection.jld2")
@load("modelfits/Nairobi_model.jld2")
nai_model = model
@load("data/cleaned_linelist20210521_c.jld2")#<--- Positive tests only looking ahead of linelist we fitted to

##Group the serology by week
gr()
zeropadtomonday = dayofweek(Date(2020,2,24)) - 1
june1day = (Date(2021,6,1) - Date(2020,2,24)).value

nai_sero_pos_rnd1_2 = vcat(zeros(zeropadtomonday),
                        sum(serological_data.serodata[:,serological_data.areas .== "NAIROBI",:,1],dims = [2,3])[:])
nai_sero_neg_rnd1_2 = vcat(zeros(zeropadtomonday),
                        sum(serological_data.serodata[:,serological_data.areas .== "NAIROBI",:,2],dims = [2,3])[:])
nai_weekly_sero_pos_rnd1_2 = [sum(grp) for grp in Iterators.partition(nai_sero_pos_rnd1_2,7)]
nai_weekly_sero_total_rnd1_2 = [sum(grp) for grp in Iterators.partition(nai_sero_pos_rnd1_2.+nai_sero_neg_rnd1_2,7)]

nai_sero_pos_rnd3 = vcat(zeros(zeropadtomonday),
                        sum(serology_data.sero[:,serology_data.areas .== "NAIROBI",:,1],dims = [2,3])[:])
nai_sero_neg_rnd3 = vcat(zeros(zeropadtomonday),
                        sum(serology_data.sero[:,serology_data.areas .== "NAIROBI",:,2],dims = [2,3])[:])
nai_weekly_sero_pos_rnd3 = [sum(grp) for grp in Iterators.partition(nai_sero_pos_rnd3,7)]
nai_weekly_sero_total_rnd3 = [sum(grp) for grp in Iterators.partition(nai_sero_pos_rnd3.+nai_sero_neg_rnd3,7)]

#Jeffery intervals
seroidxs_nai = nai_weekly_sero_total_rnd3 .> 0
# rnd_1_2idx =
uerr_nai = [invlogcdf(Beta(pos + 0.5,nai_weekly_sero_total_rnd3[seroidxs_nai][k] - pos + 0.5),log(0.975)) - pos/nai_weekly_sero_total_rnd3[seroidxs_nai][k] for (k,pos) in enumerate(nai_weekly_sero_pos_rnd3[seroidxs_nai]) ]
lerr_nai = [pos/nai_weekly_sero_total_rnd3[seroidxs_nai][k] - invlogcdf(Beta(pos + 0.5,nai_weekly_sero_total_rnd3[seroidxs_nai][k] - pos + 0.5),log(0.025)) for (k,pos) in enumerate(nai_weekly_sero_pos_rnd3[seroidxs_nai]) ]

xs_mondays = [-3 + (k-1)*7 for k = 1:length(nai_weekly_sero_total_rnd3)]
rnd1_2_idxs = xs_mondays .< 300

## Nairobi sero predictions

xticktimes = [((Date(2020,2,1) + Month(k))- Date(2020,2,24)).value for k = 1:18]
xticklabs = [monthname(k)[1:3]*"/20" for k = 3:12]
xticklabs = vcat(xticklabs,[monthname(k)[1:3]*"/21" for k = 1:8])

nai_pred = get_predictions(nai_model,Date(2021,8,1),p_ID)
N1draws = nai_model.N.*get(nai_model.MCMC_results.chain,[:P_eff]).P_eff
N2draws = nai_model.N .- N1draws
sero1 = get_credible_intervals(nai_pred.serocoverted₁./N1draws')
sero2 = get_credible_intervals(nai_pred.serocoverted₂./N2draws')
sero_tot =  get_credible_intervals((nai_pred.serocoverted₁.+nai_pred.serocoverted₂)./nai_model.N)
popexposure1 = get_credible_intervals(cumsum(nai_pred.incidence₁,dims = 1)./N1draws')
popexposure2 = get_credible_intervals(cumsum(nai_pred.incidence₂ ,dims = 1)./N2draws')
popexposure_tot = get_credible_intervals(cumsum(nai_pred.incidence₁ .+nai_pred.incidence₂ ,dims = 1)./nai_model.N)


## Other survey data
num_pos_random_trial = 384
tot_pos_random_trial = 1164
num_neg_random_trial = tot_pos_random_trial-num_pos_random_trial
uerr_all_random_trial = invlogcdf(Beta(num_pos_random_trial + 0.5,num_neg_random_trial + 0.5),log(0.975)) - num_pos_random_trial/tot_pos_random_trial
lerr_all_random_trial = num_pos_random_trial/tot_pos_random_trial - invlogcdf(Beta(num_pos_random_trial + 0.5,num_neg_random_trial + 0.5),log(0.025))
random_trial_middle_age_grp_num = 244+265 + 241 + 134 + 61
random_trial_middle_age_grp_pos = 74+100+95+51+21
random_trial_middle_age_grp_neg = random_trial_middle_age_grp_num-random_trial_middle_age_grp_pos
uerr_middle_age_random_trial = invlogcdf(Beta(random_trial_middle_age_grp_pos + 0.5,random_trial_middle_age_grp_neg + 0.5),log(0.975)) - random_trial_middle_age_grp_pos/random_trial_middle_age_grp_num
lerr_middle_age_random_trial =random_trial_middle_age_grp_pos/random_trial_middle_age_grp_num - invlogcdf(Beta(random_trial_middle_age_grp_pos + 0.5,random_trial_middle_age_grp_neg + 0.5),log(0.025))

num_pos_kibera = 37
tot_pos_kibera = 87
num_neg_kibera = tot_pos_kibera-num_pos_kibera
uerr_all_kibera = invlogcdf(Beta(num_pos_kibera + 0.5,num_neg_kibera + 0.5),log(0.975)) - num_pos_kibera/tot_pos_kibera
lerr_all_kibera = num_pos_kibera/tot_pos_kibera - invlogcdf(Beta(num_pos_kibera + 0.5,num_neg_kibera + 0.5),log(0.025))

num_pos_roysambu = 8 + 19 + 30
tot_pos_roysambu = 59 + 61 + 98
num_neg_roysambu = tot_pos_roysambu-num_pos_roysambu
uerr_all_roysambu = invlogcdf(Beta(num_pos_roysambu + 0.5,num_neg_roysambu + 0.5),log(0.975)) - num_pos_roysambu/tot_pos_roysambu
lerr_all_roysambu = num_pos_roysambu/tot_pos_roysambu - invlogcdf(Beta(num_pos_roysambu + 0.5,num_neg_roysambu + 0.5),log(0.025))


# kibera_middle_age_grp_pos = 191
# kibera_middle_age_grp_num = 384
# kibera_middle_age_grp_neg = kibera_middle_age_grp_num-kibera_middle_age_grp_pos
# uerr_middle_age_kibera = invlogcdf(Beta(kibera_middle_age_grp_pos + 0.5,kibera_middle_age_grp_neg + 0.5),log(0.975)) - kibera_middle_age_grp_pos/kibera_middle_age_grp_num
# lerr_middle_age_kibera =kibera_middle_age_grp_pos/kibera_middle_age_grp_num - invlogcdf(Beta(kibera_middle_age_grp_pos + 0.5,kibera_middle_age_grp_neg + 0.5),log(0.025))


##
plt_nai_sero = plot(4:(length(sero1.pred)+3),sero1.pred,color = 1,
        lab = "Model prediction: lower SES group",legend = :topleft,
        title = "Nairobi seropositivity: Model predictions and surveillance data",
        xticks = (xticktimes,xticklabs),
        ribbon = (sero1.lb,sero1.ub),
        xlims = (0,june1day),
        size = (1100,700),dpi = 250,
        ylabel = "Seroprevalence",
        legendfont = 10,titlefont = 20,xtickfontsize=10,
        ytickfontsize=13,guidefont = 18,
        left_margin = 10mm,right_margin = 7.5mm)
# plot!(4:(length(popexposure1.pred)+3),popexposure1.pred,
#         color = 1,lw =2,ls = :dot,lab = "")
plot!(sero2.pred,lab = "Model prediction: upper SES group",ribbon = (sero2.lb,sero2.ub),color = 2)
# plot!(4:(length(popexposure2.pred)+3),popexposure2.pred,
#         color = 2,lw =2,ls = :dot,lab = "")
plot!(sero_tot.pred,lab = "Model prediction: overall",color = 3,
        ribbon = (sero_tot.lb,sero_tot.ub))
# plot!(4:(length(popexposure_tot.pred)+3),popexposure_tot.pred,color = 3,
#         lw =2,ls = :dot,lab = "")

scatter!(xs_mondays[seroidxs_nai.*rnd1_2_idxs],nai_weekly_sero_pos_rnd3[seroidxs_nai.*rnd1_2_idxs]./nai_weekly_sero_total_rnd3[seroidxs_nai.*rnd1_2_idxs],
        lab = "KNBTS: rounds 1 and 2 (used in fitting)",
        yerr = (lerr_nai[1:11],uerr_nai[1:11]))
scatter!(xs_mondays[seroidxs_nai.*(.~rnd1_2_idxs)],nai_weekly_sero_pos_rnd3[seroidxs_nai.*(.~rnd1_2_idxs)]./nai_weekly_sero_total_rnd3[seroidxs_nai.*(.~rnd1_2_idxs)],
        lab = "KNBTS: round 3 (not used in fitting)",
        yerr = (lerr_nai[12:end],uerr_nai[12:end]))

scatter!([mean([(Date(2020,11,2) - Date(2020,2,24)).value,(Date(2020,11,23) - Date(2020,2,24)).value])],
        [num_pos_random_trial/tot_pos_random_trial],lab = "Nairobi randomised survey: all ages (not used in fitting)",
        ms = 8,yerr = ([lerr_all_random_trial],[uerr_all_random_trial]),xerr = ([10.5],[10.5]),
        markershape = :square,mc = 3)

scatter!([mean([(Date(2020,11,2) - Date(2020,2,24)).value,(Date(2020,11,23) - Date(2020,2,24)).value])],
        [random_trial_middle_age_grp_pos/random_trial_middle_age_grp_num],lab = "Nairobi randomised survey: 10-60 year olds",
        ms = 8,yerr = ([lerr_middle_age_random_trial],[uerr_middle_age_random_trial]),xerr = ([10.5],[10.5]),
        markershape = :diamond,mc = 3)

scatter!([mean([(Date(2020,11,2) - Date(2020,2,24)).value,(Date(2020,11,23) - Date(2020,2,24)).value])],
        [num_pos_kibera/tot_pos_kibera],
        ms = 8,yerr = ([lerr_all_kibera],[uerr_all_kibera]),xerr = ([10.5],[10.5]),lab = "Nairobi randomised survey: Kibera",
        markershape = :square,mc = 1)

scatter!([mean([(Date(2020,11,2) - Date(2020,2,24)).value,(Date(2020,11,23) - Date(2020,2,24)).value])],
        [num_pos_roysambu/tot_pos_roysambu],
        ms = 8,yerr = ([lerr_all_roysambu],[uerr_all_roysambu]),xerr = ([10.5],[10.5]),lab = "Nairobi randomised survey: low poverty/informal settlement sub-counties",
        markershape = :square,mc = 2)        

# scatter!([mean([(Date(2020,11,27) - Date(2020,2,24)).value,(Date(2020,12,5) - Date(2020,2,24)).value])],[kibera_middle_age_grp_pos/kibera_middle_age_grp_num],
#         ms = 8,yerr = ([lerr_middle_age_kibera],[uerr_middle_age_kibera]),xerr = ([4],[4]),lab = " Kibera random household survey : 10-60 year olds (not used in fitting)",
#         markershape = :diamond,mc = 1)

# savefig("plots/nairobi_seroplot.png")


