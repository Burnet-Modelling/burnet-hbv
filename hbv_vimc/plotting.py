import atomica as at
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import math
from matplotlib.backends.backend_pdf import PdfPages
import hbv_vimc as hbv
import seaborn as sns
import itertools

def calib_diag_plots(country:str, cal, res_cal):
    
    """ Produces a combined pdf of calibration diagnostic plots
    Author: Chris Seaman (Updated 05/10/23)
    Page 1: Population sizes [done]
    Page 2: Demographic validations (all cause mortality, births, birth-sex ratio, age of births, migration) [done]
    Page 3: Transmission dynamics [done]
    Page 4: HBsAg positive (n) by population [done]
    Page 5: HBsAg prevalence by population [done]
    Page 6: HBeAg positive by population [done]
    Page 7: Disease progression rates by age/sex
    Page 8: Disease state distribution over time (IT, ICL, ICT, IE, CC, DC, HCC) in each age group [done]
    Page 9: Compensated Cirrhosis prevalence by population [done]
    Page 10: Decompensated Cirrhosis prevalence by population [done]
    Page 11: HCC incidence by population [done]
    Page 12: HCC prevalence by population [done]
    Page 13: Mortality by population [done]
    Page 14: DALYs (YLL, YLD) by population
    Page 15: Birth flow dynamics
    """

    # Output name and PDF initiation
    filename = hbv.root/"outputs"/"calibrations"/f"{country}_calplot_{res_cal.name}.pdf"
    pp = PdfPages(filename)

    # Required data and settings
    sns.set_theme(style="whitegrid")
    pop_bins, birth_bins, cons_pars, var_pars=hbv.age_bins("vimc_test")
    pops = list(cons_pars.columns)
    s_yr = hbv.s_yr
    e_yr = hbv.e_yr

    # Figure 1: Population size
    fig1=plt.figure(figsize=(20,20))

    for idx,pop in enumerate(pops):
        f1=fig1.add_subplot(int(np.ceil(len(pops)/4)), 4, idx+1)
        f1.scatter(np.arange(s_yr, e_yr, 1), cal.get_par("alive").ts[pop].vals, alpha=0.2, label="Data")
        f1.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="alive", pops=pop, t_bins=1).series[0].vals,
                  color="red", label="Model Projection")  # Model
        f1.set_title(pop)
        f1.set_ylabel("Population")
        f1.set_ylim(bottom=0)

        if idx==0:
            f1.legend(loc="best")

    fig1.suptitle(f"{country} Population: Model Projections", fontsize=25, style="normal", color="blue")
    fig1.tight_layout(pad=2)
    fig1.subplots_adjust(top=0.95)
    pp.savefig(fig1)  # Attaches plot to PDF output

    # Figure 2: Population parameters summary

    # Births (stacked by maternal age)
    births = np.zeros(((e_yr-s_yr, len(birth_bins)+1)))
    births[:, 0] = at.PlotData(res_cal, outputs="b_rate", pops="total", t_bins=1).series[0].vals  # Total births

    for idx, b_bin in enumerate(birth_bins):
        births [:, idx+1] = at.PlotData(res_cal, outputs="pregs", pops=f'{b_bin}F', t_bins=1).series[0].vals  # Moth Age

    # Birth Sex Ratio (male to female)
    b_sex = np.zeros((e_yr-s_yr, 1))
    b_sex[:, 0] = at.PlotData(res_cal, outputs="b_rate", pops="0-0M", t_bins=1).series[0].vals\
        / at.PlotData(res_cal, outputs="b_rate", pops="0-0F", t_bins=1).series[0].vals

    # All Cause Mortality
    acm = np.zeros((e_yr-s_yr, len(pops)))
    for idx, pop in enumerate(pops):
        acm[:, idx] = at.PlotData(res_cal, outputs="acm", pops=pop, t_bins=1).series[0].vals

    # TODO: Be mindful that changes in approach may need to be factored in here
    # Net Migration Rate (total)

    # Inputs
    net_m_in = np.zeros((len(cal.get_par("mig_rate").ts["0-0M"].t),2))
    net_m_in[:, 0] = cal.get_par("mig_rate").ts["0-0M"].t
    net_m_in[:, 1] = cal.get_par("mig_rate").ts["0-0M"].vals

    # Modelled
    net_m_mod = np.zeros((e_yr-s_yr, 3))
    net_m_mod[:, 0] = np.nan_to_num(at.PlotData(res_cal, outputs="im_rate", pops="total", pop_aggregation="weighted", t_bins=1).series[0].vals)
    net_m_mod[:, 1] = np.nan_to_num(at.PlotData(res_cal, outputs="em_rate", pops="total", pop_aggregation="weighted", t_bins=1).series[0].vals)
    net_m_mod[:, 1] = net_m_mod[:, 1] * -1
    net_m_mod[:, 2] = net_m_mod[:, 0] + net_m_mod[:, 1]

    #net_m_in = np.zeros((e_yr-s_yr, 6))
    #net_m[:, 0] = cal.get_par("mig_rate").ts["0-0M"].vals
    #net_m[:, 1] = cal.get_par("em_rate").ts["0-0M"].vals
    #net_m[:, 1] = net_m[:, 1] * -1
    #net_m[:, 2] = net_m[:, 0] + net_m[:, 1]


    # Plot
    fig2 = plt.figure(figsize=(20, 20))

    # Births
    f2a = fig2.add_subplot(2,2,1)
    f2a.plot(np.arange(s_yr, e_yr, 1), births[:, 0], color="blue", linestyle="dashed", alpha=0.6)
    f2a.bar(np.arange(s_yr, e_yr, 1), births[:, 1], label="10-19")
    f2a.bar(np.arange(s_yr, e_yr, 1), births[:, 2], bottom=births[:, 1], label="20-29")
    f2a.bar(np.arange(s_yr, e_yr, 1), births[:, 3], bottom=births[:, 1] + births[:, 2], label="30-39")
    f2a.bar(np.arange(s_yr, e_yr, 1), births[:, 4], bottom=births[:, 1] + births[:, 2] + births[:, 3],label="40-49")
    f2a.bar(np.arange(s_yr, e_yr, 1), births[:, 5], bottom=births[:, 1] + births[:, 2] + births[:, 3]+  births[:, 4],label="50-54")
    f2a.set_ylim(bottom=0)
    f2a.set_xlim(left=1990, right=2100)
    f2a.set_title("Births")
    f2a.set_ylabel("Annual Births")
    f2a.legend(loc="best")

    # Birth-Sex Ratio
    f2b = fig2.add_subplot(2,2,2)
    f2b.plot(np.arange(s_yr, e_yr, 1), b_sex[:,0], color="blue")
    f2b.hlines(xmin=1990, xmax=2100, y=1, linestyles="dashed", alpha=0.6, color="red")
    f2b.set_ylim(bottom=0.5, top=1.5)
    f2b.set_xlim(left=1990, right=2100)
    f2b.set_title("Male to Female Birth Sex Ratio")
    f2b.set_ylabel("Male Births/Female Births")

    # All-cause mortality rate
    acm_col=['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9',
             '#bc80bd', '#ccebc5', '#ffed6f']
    f2c = fig2.add_subplot(2,2,3)
    for idx, pop in enumerate(pops):
        if (idx % 2) == 0:
            f2c.scatter(np.arange(s_yr, e_yr, 1), acm[:, idx], marker="x", color=acm_col[int(idx/2)], label=pop)
        else:
            f2c.scatter(np.arange(s_yr, e_yr, 1), acm[:, idx], marker="*", color=acm_col[int(np.ceil(idx/2))-1], label=pop)
    f2c.set_xlim(left=1990, right=2100)
    f2c.set_ylim(bottom=0)
    f2c.set_ylabel("Annual mortality probability")
    f2c.set_title("All cause mortality rate")
    f2c.legend(loc="best", ncol=4)

    # Net migration rate
    f2d = fig2.add_subplot(2,2,4)
    f2d.scatter(net_m_in[:, 0], net_m_in[:, 1], color="blue", alpha=0.2, label="Data")
    f2d.plot(np.arange(s_yr, e_yr, 1), net_m_mod[:, 2], color="red", label="Model (Calibrated)")
    f2d.set_xlim(left=1990, right=2100)
    f2d.hlines(xmin=1990, xmax=2100, y=0, linestyles="dashed", alpha=0.6, color="red")
    f2d.set_ylabel("Annual migration rate")
    f2d.set_title("Net migration rate")
    f2d.legend(loc="best")

    fig2.suptitle(f"{country} Demographic Inputs", fontsize=25, style="normal", color="blue")
    fig2.tight_layout(pad=2)
    fig2.subplots_adjust(top=0.95)

    pp.savefig(fig2)  # Attaches plot to PDF output

    # Figure 3: Transmission dynamics
    fig3 = plt.figure(figsize=(20, 20))

    for idx, pop in enumerate(pops):
        trans = np.zeros((e_yr-s_yr, 2))
        trans[:, 0] = at.PlotData(res_cal, outputs="tot_inc", pops=pop, t_bins=1).series[0].vals  # CHB incidence
        trans[:, 1] = at.PlotData(res_cal, outputs={"trans": "sus:j_acu+mtct_inf+mtct_rec"}, pops=pop, t_bins=1).series[0].vals  # total transmission
        trans[:, 1] = trans[:, 1] - trans[:, 0]  # transmission, not chronic

        f3 = fig3.add_subplot(int(np.ceil(len(pops) / 4)), 4, idx + 1)
        f3.bar(np.arange(s_yr, e_yr, 1), trans[:, 0], label="Chronic Hepatitis B")
        f3.bar(np.arange(s_yr, e_yr, 1), trans[:, 1], bottom=trans[:, 0], label="Non-Chronic (Transmission)")
        f3.set_title(pop)
        f3.set_ylabel("HBV transmissions")
        f3.set_ylim(bottom=0)
        f3.set_xlim(left=1990, right=2100)

        if idx == 0:
            f3.legend(loc="best")

    fig3.suptitle(f"{country} Transmission Dynamics", fontsize=25, style="normal", color="blue")
    fig3.tight_layout(pad=2)
    fig3.subplots_adjust(top=0.95)
    pp.savefig(fig3)  # Attaches plot to PDF

    # Figure 4: HBsAg positive population size
    fig4 = plt.figure(figsize=(20, 20))

    for idx, pop in enumerate(pops):
        f4 = fig4.add_subplot(int(np.ceil(len(pops) / 4)), 4, idx + 1)
        f4.scatter(cal.get_par("chb_pop").ts[pop].t, cal.get_par("chb_pop").ts[pop].vals, alpha=0.2, label="Data")
        f4.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="chb_pop", pops=pop, t_bins=1).series[0].vals,
                color="red", label="Model Projection")  # Model
        f4.set_title(pop)
        f4.set_ylabel("CHB Positive")
        f4.set_ylim(bottom=0)
        f4.set_xlim(left=1990, right=2030)

        if idx == 0:
            f4.legend(loc="best")

    fig4.suptitle(f"{country} CHB Positive: Model Projections", fontsize=25, style="normal", color="blue")
    fig4.tight_layout(pad=2)
    fig4.subplots_adjust(top=0.95)
    pp.savefig(fig4)  # Attaches plot to PDF output

    # Figure 5: HBsAg prevalence
    fig5 = plt.figure(figsize=(20, 20))

    who_prev = pd.read_csv(hbv.root/"who_prev.csv")
    who_prev = who_prev[who_prev["country_code"] == country]

    for idx, pop in enumerate(pops):
        f5 = fig5.add_subplot(int(np.ceil(len(pops) / 4)), 4, idx + 1)
        f5.scatter(cal.get_par("prev").ts[pop].t, cal.get_par("prev").ts[pop].vals, alpha=0.2, label="GBD Estimate")
        # TODO: Finish this
        f5.scatter(who_prev["year"], who_prev[pop], alpha=0.8, marker="*", color="purple", label="WHO Estimate")
        f5.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="prev", pops=pop, t_bins=1).series[0].vals,
                color="red", label="Model Projection")  # Model
        f5.set_title(pop)
        f5.set_ylabel("CHB Prevalence")
        f5.set_ylim(bottom=0)
        f5.set_xlim(left=1990, right=2030)

        if idx == 0:
            f5.legend(loc="best")

    fig5.suptitle(f"{country} CHB Prevalence: Model Projections", fontsize=25, style="normal", color="blue")
    fig5.tight_layout(pad=2)
    fig5.subplots_adjust(top=0.95)
    pp.savefig(fig5)  # Attaches plot to PDF output

    # Figure 6: HBeAg prevalence
    fig6 = plt.figure(figsize=(20, 20))

    for idx, pop in enumerate(pops):
        f6 = fig6.add_subplot(int(np.ceil(len(pops) / 4)), 4, idx + 1)
        f6.scatter(cal.get_par("eag_ott").ts[pop].t, cal.get_par("eag_ott").ts[pop].vals, alpha=0.2, label="Data")
        f6.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="eag_ott", pops=pop, t_bins=1).series[0].vals,
                color="red", label="Model Projection")  # Model
        f6.set_title(pop)
        f6.set_ylabel("HBeAg Prevalence")
        f6.set_ylim(bottom=0)
        f6.set_xlim(left=1990, right=2030)

        if idx == 0:
            f6.legend(loc="best")

    fig6.suptitle(f"{country} HBeAg prevalence: Model Projections", fontsize=25, style="normal", color="blue")
    fig6.tight_layout(pad=2)
    fig6.subplots_adjust(top=0.95)
    pp.savefig(fig6)  # Attaches plot to PDF output

    # Figure 7: Disease Progression Rates (by age) post calibration

    col_names = ["pops", "it_hcc", "icl_hcc", "ict_hcc", "ie_hcc", "cc_hcc", "dc_hcc", "icl_cc", "ie_cc", "cc_dc",
            "m_acu", "m_dc", "m_hcc", "it_icl", "icl_ict", "ict_ie"]

    pars = ["it_hcc", "icl_hcc", "ict_hcc", "ie_hcc", "cc_hcc", "dc_hcc", "icl_cc", "ie_cc", "cc_dc",
            "m_acu", "m_dc", "m_hcc", "it_icl", "icl_ict", "ict_ie"]

    par_df = pd.DataFrame(columns=col_names)
    par_df["pops"] = pops
    par_df = par_df.set_index("pops")

    for i, pop in enumerate(pops):
        for par in pars:
            par_df.at[pop, par] = at.PlotData(res_cal, outputs=par, pops=pop, t_bins=1).series[0].vals[0]

    par_df=par_df.reset_index()

    fig,ax = plt.subplots()
    fig.patch.set_visible(False)
    ax.axis('off')
    ax.axis('tight')
    the_table = ax.table(cellText=par_df.values, colLabels=par_df.columns, loc="center", fontsize=20)
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(2)
    fig.tight_layout()
    pp.savefig(fig)# Attaches plot to PD

    # Figure 8: Disease State distribution over time (each population)
    fig8 = plt.figure(figsize=(20, 20))

    for idx, pop in enumerate(pops):
        dsd = np.zeros((e_yr-s_yr, 8))
        dsd[:, 0] =  at.PlotData(res_cal, outputs={"it": "it+it_dx+it_tx"}, pops=pop, t_bins=1).series[0].vals
        dsd[:, 1] =  at.PlotData(res_cal, outputs={"icl": "icl+icl_dx+icl_tx"}, pops=pop, t_bins=1).series[0].vals
        dsd[:, 2] =  at.PlotData(res_cal, outputs={"ict": "ict+ict_dx+ict_tx+ict_tx_i"}, pops=pop, t_bins=1).series[0].vals
        dsd[:, 3] =  at.PlotData(res_cal, outputs={"ie": "ie+ie_dx+ie_tx"}, pops=pop, t_bins=1).series[0].vals
        dsd[:, 4] =  at.PlotData(res_cal, outputs="comp", pops=pop, t_bins=1).series[0].vals
        dsd[:, 5] =  at.PlotData(res_cal, outputs="decomp", pops=pop, t_bins=1).series[0].vals
        dsd[:, 6] =  at.PlotData(res_cal, outputs="hepcc", pops=pop, t_bins=1).series[0].vals
        dsd[:, 7] =  at.PlotData(res_cal, outputs="chb_pop", pops=pop, t_bins=1).series[0].vals

        dsd[:, 0] = dsd[:, 0] / dsd[:, 7]
        dsd[:, 1] = dsd[:, 1] / dsd[:, 7]
        dsd[:, 2] = dsd[:, 2] / dsd[:, 7]
        dsd[:, 3] = dsd[:, 3] / dsd[:, 7]
        dsd[:, 4] = dsd[:, 4] / dsd[:, 7]
        dsd[:, 5] = dsd[:, 5] / dsd[:, 7]
        dsd[:, 6] = dsd[:, 6] / dsd[:, 7]

        f8 = fig8.add_subplot(int(np.ceil(len(pops) / 4)), 4, idx + 1)
        f8.bar(np.arange(s_yr, e_yr, 1), dsd[:, 0], label="Immune Tolerant")
        f8.bar(np.arange(s_yr, e_yr, 1), dsd[:, 1], bottom=dsd[:, 0], label="HBeAg positive active")
        f8.bar(np.arange(s_yr, e_yr, 1), dsd[:, 2], bottom=dsd[:, 0]+dsd[:, 1], label="HBeAg negative inactive")
        f8.bar(np.arange(s_yr, e_yr, 1), dsd[:, 3], bottom=dsd[:, 0]+dsd[:, 1]+dsd[:, 2], label="HBeAg negative active")
        f8.bar(np.arange(s_yr, e_yr, 1), dsd[:, 4], bottom=dsd[:, 0]+dsd[:, 1]+dsd[:, 2]+dsd[:, 3], label="Compensated Cirrhosis")
        f8.bar(np.arange(s_yr, e_yr, 1), dsd[:, 5], bottom=dsd[:, 0]+dsd[:, 1]+dsd[:, 2]+dsd[:, 3]+dsd[:, 4], label="Decompensated Cirrhosis")
        f8.bar(np.arange(s_yr, e_yr, 1), dsd[:, 6], bottom=dsd[:, 0]+dsd[:, 1]+dsd[:, 2]+dsd[:, 3]+dsd[:, 4]+dsd[:, 5], label="Hepatocellular Carcinoma")
        f8.set_title(pop)
        f8.set_ylabel("Proportion of CHB population")
        f8.set_ylim(bottom=0, top=1)
        f8.set_xlim(left=1990, right=2100)

        if idx == 0:
            f8.legend(loc="best")

    fig8.suptitle(f"{country} Disease State Distribution", fontsize=25, style="normal", color="blue")
    fig8.tight_layout(pad=2)
    fig8.subplots_adjust(top=0.95)
    pp.savefig(fig8)  # Attaches plot to PD

    # Figure 9: Compensated Cirrhosis prevalence by population
    fig9 = plt.figure(figsize=(20, 20))

    for idx, pop in enumerate(pops):
        f9 = fig9.add_subplot(int(np.ceil(len(pops) / 4)), 4, idx + 1)
        f9.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="comp", pops=pop, t_bins=1).series[0].vals,
                color="red")  # Model
        f9.set_title(pop)
        f9.set_ylabel("Prevalent Cases")
        f9.set_ylim(bottom=0)
        f9.set_xlim(left=1990, right=2100)

    fig9.suptitle(f"{country} Compensated Cirrhosis prevalence: Model Projections", fontsize=25, style="normal", color="blue")
    fig9.tight_layout(pad=2)
    fig9.subplots_adjust(top=0.95)
    pp.savefig(fig9)

    # Figure 10: Decompensated Cirrhosis prevalence by population
    fig10 = plt.figure(figsize=(20, 20))

    for idx, pop in enumerate(pops):
        f10 = fig10.add_subplot(int(np.ceil(len(pops) / 4)), 4, idx + 1)
        f10.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="decomp", pops=pop, t_bins=1).series[0].vals,
                color="red")  # Model
        f10.set_title(pop)
        f10.set_ylabel("Prevalent Cases")
        f10.set_ylim(bottom=0)
        f10.set_xlim(left=1990, right=2100)

    fig10.suptitle(f"{country} Decompensated Cirrhosis prevalence: Model Projections", fontsize=25, style="normal", color="blue")
    fig10.tight_layout(pad=2)
    fig10.subplots_adjust(top=0.95)
    pp.savefig(fig10)

    # Figure 11: HCC incidence by population
    fig11 = plt.figure(figsize=(20, 20))

    for idx, pop in enumerate(pops):
        f11 = fig11.add_subplot(int(np.ceil(len(pops) / 4)), 4, idx + 1)
        f11.scatter(cal.get_par("flw_hcc").ts[pop].t, cal.get_par("flw_hcc").ts[pop].vals, alpha=0.2, label="Data")
        f11.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="flw_hcc", pops=pop, t_bins=1).series[0].vals,
                color="red")  # Model
        f11.set_title(pop)
        f11.set_ylabel("Incident Cases")
        f11.set_ylim(bottom=0)
        f11.set_xlim(left=1990, right=2100)

    fig11.suptitle(f"{country} HCC Incidence: Model Projections", fontsize=25, style="normal", color="blue")
    fig11.tight_layout(pad=2)
    fig11.subplots_adjust(top=0.95)
    pp.savefig(fig11)

    # Figure 12: HCC prevalence by population
    fig12 = plt.figure(figsize=(20, 20))

    for idx, pop in enumerate(pops):
        f12 = fig12.add_subplot(int(np.ceil(len(pops) / 4)), 4, idx + 1)
        f12.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="hepcc", pops=pop, t_bins=1).series[0].vals,
                color="red")  # Model
        f12.set_title(pop)
        f12.set_ylabel("Prevalent Cases")
        f12.set_ylim(bottom=0)
        f12.set_xlim(left=1990, right=2100)

    fig12.suptitle(f"{country} HCC prevalence: Model Projections", fontsize=25, style="normal", color="blue")
    fig12.tight_layout(pad=2)
    fig12.subplots_adjust(top=0.95)
    pp.savefig(fig12)

    # Figure 13: Acute, Cirrhosis and HCC mortality by population
    fig13 = plt.figure(figsize=(20, 20))

    for idx, pop in enumerate(pops):
        f13 = fig13.add_subplot(int(np.ceil(len(pops) / 4)), 4, idx + 1)
        # Acute
        # f13.scatter(cal.get_par("cl_acu").ts[pop].t, cal.get_par("cl_acu").ts[pop].vals, alpha=0.2, color="green")  # Data
        # f13.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="cl_acu", pops=pop, t_bins=1).series[0].vals,
        #          color="green", label="Acute")  # Model
        # Cirrhosis
        f13.scatter(cal.get_par("cl_cir").ts[pop].t, cal.get_par("cl_cir").ts[pop].vals, alpha=0.2,
                    color="blue")  # Data
        f13.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="cl_cir", pops=pop, t_bins=1).series[0].vals,
                 color="blue", label="Cirrhosis")  # Model
        # HCC
        f13.scatter(cal.get_par("cl_hcc").ts[pop].t, cal.get_par("cl_hcc").ts[pop].vals, alpha=0.2,
                    color="red")  # Data
        f13.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="cl_hcc", pops=pop, t_bins=1).series[0].vals,
                 color="red", label="HCC")  # Model

        f13.set_title(pop)
        f13.set_ylabel("Deaths")
        f13.set_ylim(bottom=0)
        f13.set_xlim(left=1990, right=2100)

        if idx == 0:
            f13.legend(loc="best")

    fig13.suptitle(f"{country} HBV Deaths: Model Projections", fontsize=25, style="normal", color="blue")
    fig13.tight_layout(pad=2)
    fig13.subplots_adjust(top=0.95)
    pp.savefig(fig13)

    # Figure 14: DALYs (yld, yll) by population
    fig14 = plt.figure(figsize=(20, 20))

    for idx, pop in enumerate(pops):
        f14 = fig14.add_subplot(int(np.ceil(len(pops) / 4)), 4, idx + 1)
        # YLL
        f14.bar(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="yll", pops=pop, t_bins=1).series[0].vals,
                label="YLL")

        # YLD
        f14.bar(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="yld", pops=pop, t_bins=1).series[0].vals,
                bottom=at.PlotData(res_cal, outputs="yll", pops=pop, t_bins=1).series[0].vals, label="YLD")

        # DALYs (total)
        # f14.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="dalys", pops=pop, t_bins=1).series[0].vals,
        #          linestyle="dashed", color="blue")

        f14.set_title(pop)
        f14.set_ylabel("DALYs")
        f14.set_ylim(bottom=0)
        f14.set_xlim(left=1990, right=2100)

        if idx == 0:
            f14.legend(loc="best")

    fig14.suptitle(f"{country} HBV DALYs: Model Projections", fontsize=25, style="normal", color="blue")
    fig14.tight_layout(pad=2)
    fig14.subplots_adjust(top=0.95)
    pp.savefig(fig14)

    # Figure 15: Birth Junction Flow Dynamics
    fig15 = plt.figure(figsize=(20, 20))

    bflow = fig15.add_subplot(4,1,1)

    bflow.bar(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="b_sus",
              pops={"0-0": ["0-0M", "0-0F"]}, t_bins=1).series[0].vals, label="Susceptible")

    bflow.bar(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="b_vax",
              pops={"0-0": ["0-0M", "0-0F"]}, t_bins=1).series[0].vals,
              bottom=at.PlotData(res_cal, outputs="b_sus",
              pops={"0-0": ["0-0M", "0-0F"]}, t_bins=1).series[0].vals,
              label="Vaccinated")

    bflow.bar(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="mtct_inf",
              pops={"0-0": ["0-0M", "0-0F"]}, t_bins=1).series[0].vals,
              bottom=at.PlotData(res_cal, outputs="b_sus",
              pops={"0-0": ["0-0M", "0-0F"]}, t_bins=1).series[0].vals + at.PlotData(res_cal, outputs="b_vax",
              pops={"0-0": ["0-0M", "0-0F"]}, t_bins=1).series[0].vals,
              label="Infected")

    bflow.bar(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="mtct_rec",
              pops={"0-0": ["0-0M", "0-0F"]}, t_bins=1).series[0].vals,
              bottom=at.PlotData(res_cal, outputs="b_sus", pops={"0-0": ["0-0M", "0-0F"]}, t_bins=1).series[0].vals +
              at.PlotData(res_cal, outputs="b_vax", pops={"0-0": ["0-0M", "0-0F"]}, t_bins=1).series[0].vals +
              at.PlotData(res_cal, outputs="mtct_inf", pops={"0-0": ["0-0M", "0-0F"]}, t_bins=1).series[0].vals,
              label="Recovered")

    bflow.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="b_rate",
               pops={"0-0": ["0-0M", "0-0F"]}, t_bins=1).series[0].vals, linestyle="dashed", alpha=0.6, color="red",
               label="Total Births")

    bflow.set_title("Birth Flow")
    bflow.set_ylim(bottom=0)
    bflow.set_xlim(left=1990, right=2023)
    bflow.legend(loc="best")

    inc = fig15.add_subplot(4,1,2)

    inc.bar(np.arange(s_yr, e_yr, 1),at.PlotData(res_cal, outputs="mtct_inf",
               pops={"0-0": ["0-0M", "0-0F"]}, t_bins=1).series[0].vals, label="MTCT")

    inc.bar(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="j_acu:it", pops={"0-0": ["0-0M", "0-0F"]}, t_bins=1).series[0].vals,
            bottom=at.PlotData(res_cal, outputs="mtct_inf", pops={"0-0": ["0-0M", "0-0F"]}, t_bins=1).series[0].vals,
            label="ECHT")

    inc.set_title("0y transmission dynamics")
    inc.set_ylim(bottom=0)
    inc.set_xlim(left=1990, right=2023)
    inc.legend(loc="best")

    tot_inc = fig15.add_subplot(4,1,3)

    tot_inc.bar(np.arange(s_yr, e_yr, 1),at.PlotData(res_cal, outputs="mtct_inf",
               pops={"0-0": ["0-0M", "0-0F"]}, t_bins=1).series[0].vals, label="MTCT")

    tot_inc.bar(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="j_acu:it", pops={"0-4": ["0-0M", "0-0F", "1-4M", "1-4F"]}, t_bins=1).series[0].vals,
                bottom=at.PlotData(res_cal, outputs="mtct_inf", pops={"0-0": ["0-0M", "0-0F"]}, t_bins=1).series[0].vals,
                label="ECHT")

    tot_inc.bar(np.arange(s_yr, e_yr, 1),at.PlotData(res_cal, outputs="j_acu:icl", pops="total", t_bins=1).series[0].vals ,
                bottom=at.PlotData(res_cal, outputs="mtct_inf", pops={"0-0": ["0-0M", "0-0F"]}, t_bins=1).series[0].vals +
                at.PlotData(res_cal, outputs="j_acu:it", pops={"0-4": ["0-0M", "0-0F", "1-4M", "1-4F"]}, t_bins=1).series[0].vals,
                label=">5 years")

    tot_inc.set_title("Total Incidence")
    tot_inc.set_ylim(bottom=0)
    tot_inc.set_xlim(left=1990, right=2023)
    tot_inc.legend(loc="best")

    vax = fig15.add_subplot(4,1,4)
    vax.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="bd",
            pops={"0-0": ["0-0M", "0-0F"]}, pop_aggregation="weighted", t_bins=1).series[0].vals,
            alpha=0.6, color="red", label = "HepB-BD")
    vax.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="hb3",
            pops={"0-0": ["0-0M", "0-0F"]}, pop_aggregation="weighted", t_bins=1).series[0].vals,
            alpha=0.6, color="blue", label="HepB3")
    vax.set_title("Vaccination Coverage")
    vax.set_ylim(bottom=0)
    vax.set_xlim(left=1990, right=2023)
    vax.legend(loc="best")

    pp.savefig(fig15)

    # Figure 16: Weighted (interaction) prevalence and transmission values
    fig16 = plt.figure(figsize=(20, 20))

    # Mother to child transmission (weighted should sum to interaction prev)
    prv_nint = np.zeros((e_yr-s_yr, 13))

    #  Extract data
    prv_nint[:, 0] = at.PlotData(res_cal, outputs="prev", pops="10-19F", t_bins=1).series[0].vals
    prv_nint[:, 1] = at.PlotData(res_cal, outputs="pregs", pops="10-19F", t_bins=1).series[0].vals
    prv_nint[:, 2] = at.PlotData(res_cal, outputs="prev", pops="20-29F", t_bins=1).series[0].vals
    prv_nint[:, 3] = at.PlotData(res_cal, outputs="pregs", pops="20-29F", t_bins=1).series[0].vals
    prv_nint[:, 4] = at.PlotData(res_cal, outputs="prev", pops="30-39F", t_bins=1).series[0].vals
    prv_nint[:, 5] = at.PlotData(res_cal, outputs="pregs", pops="30-39F", t_bins=1).series[0].vals
    prv_nint[:, 6] = at.PlotData(res_cal, outputs="prev", pops="40-49F", t_bins=1).series[0].vals
    prv_nint[:, 7] = at.PlotData(res_cal, outputs="pregs", pops="40-49F", t_bins=1).series[0].vals
    prv_nint[:, 8] = at.PlotData(res_cal, outputs="prev", pops="50-59F", t_bins=1).series[0].vals
    prv_nint[:, 9] = at.PlotData(res_cal, outputs="pregs", pops="50-59F", t_bins=1).series[0].vals

    # Calculate weighted outcomes
    prv_nint[:, 10] = prv_nint[:, 0] * prv_nint[:, 1] + prv_nint[:, 2] * prv_nint[:, 3] + prv_nint[:, 4] *\
                      prv_nint[:, 5] + prv_nint[:, 6] * prv_nint[:, 7] + prv_nint[:, 8] * prv_nint[:, 9]

    prv_nint[:, 11] = prv_nint[:, 1] + prv_nint[:, 3] + prv_nint[:, 5] + prv_nint[:, 7] + prv_nint[:, 9]

    prv_nint[:, 12] = prv_nint[:, 10] / prv_nint[:, 11]


    mtct = fig16.add_subplot(2,2,1)
    mtct.plot(np.arange(s_yr, e_yr, 1), prv_nint[:, 12], linestyle="dashed", alpha=0.6, color="red",
               label="Weighted Prevalence")
    mtct.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs={"int_prv":"e_preg+s_preg"}, pops="0-0M", t_bins=1).series[0].vals,
              color="black", label="Interaction Prevalence")
    mtct.legend(loc="best")

    # <5 year prevalence (from interaction matrix) should equal weighted population prevalence
    echt = fig16.add_subplot(2,2,2)
    echt.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="prev", pops={"Under 5":["0-0M", "0-0F", "1-4M", "1-4F"]},
                                                    pop_aggregation="weighted", t_bins=1).series[0].vals, linestyle="dashed",
                                                    alpha=0.6, color="red", label="Weighted Prevalence")
    echt.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="foi_in_ec", pops="0-0M", t_bins=1).series[0].vals,
              color="black", label="Interaction Prevalence")
    echt.legend(loc="best")

    # >5 year prevalence (from interaction matrix) should equal weighted population prevalence
    adt = fig16.add_subplot(2,2,3)
    adt.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="prev", pops={"Over 5":["5-9M", "5-9F", "10-19M", "10-19F",
                                                    "20-29M", "20-29F", "30-39M", "30-39F", "40-49M", "40-49F", "50-59M",
                                                    "50-59F", "60-69M", "60-69F", "70-79M", "70-79F", "80-89M", "80-89F",
                                                                                             "90+M", "90+F"]},
                                                    pop_aggregation="weighted", t_bins=1).series[0].vals, linestyle="dashed",
                                                    alpha=0.6, color="red", label="Weighted Prevalence")
    adt.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs="foi_in_horiz", pops="10-19M", t_bins=1).series[0].vals,
              color="black", label="Interaction Prevalence")
    adt.legend(loc="best")

    pp.savefig(fig16)

    pp.close()


def stochastic_outcome_plots(country:str, central, st_plot):
    """Basic visualization of outcomes from modelled scenarios with
    uncertainity intervals (95% CrI)"""
    
    ## Pre-requisites for plots
    s_yr=hbv.s_yr
    e_yr=hbv.e_yr
    scens=hbv.get_scenarios(country)
    sns.set_theme(style="whitegrid")

    ## Outcomes of interest (note -order needs to match for this to work)
    outcomes=[["prev", "total", "weighted", "total HBV prevalence"],
              ["prev", {"0-4":["0-0F", "0-0M","1-4F", "1-4M"]}, "weighted", "<5y HBV prevalence"],
              ["flw_hcc", "total", "sum", "HCC incidence"],  # Incident HCC
              [":dd_hbv", "total", "sum", "HBV mortality"]]  # DALYs
    
    plot_names = ["t_prev", "c_prev", "i_hcc", "d_hbv"]

    
    ## Extract data for plotting
    cent_plot={}
    for scen in scens.values():
            # Initialize the dictionary by scenario name
            cent_plot[scen.name]=np.repeat(np.zeros((e_yr-s_yr,1)), len(outcomes), axis=0)
            # Extract the central results
            for idx,val in enumerate(outcomes):
                cent_plot[scen.name][int(idx*(e_yr-s_yr)):int((idx*(e_yr-s_yr))+(e_yr-s_yr)), 0]=np.nan_to_num(at.PlotData(central[scen.name], val[0], pops=val[1], pop_aggregation=val[2], t_bins=1).series[0].vals)
                
    ## Produce plots (central estimate with uncertainty from stochastic runs)
    fig = plt.figure(figsize=(15, 20))
    cols = 2
    rows = int(np.ceil(len(outcomes)/cols))
    
    for idx,val in enumerate(outcomes):
        ax=fig.add_subplot(cols,rows,idx+1)
        for scen in scens.values():
            ax.plot(np.arange(s_yr, e_yr,1),cent_plot[scen.name][int(idx*(e_yr-s_yr)):int((idx*(e_yr-s_yr))+(e_yr-s_yr)),0], label=scen.name)
            ax.fill_between(np.arange(s_yr, e_yr,1),np.percentile(st_plot[scen.name][plot_names[idx]], 2.5, axis=1),np.percentile(st_plot[scen.name][plot_names[idx]], 97.5, axis=1), alpha=0.2)
            ax.set_title(val[3])
            ax.set_xlim(s_yr+10, e_yr-1)
            
            if idx==0:
                ax.legend(loc="best")
    
    fig.tight_layout()
    fig.savefig(hbv.root/"outputs"/"uncert_plots"/f"{country}.png", dpi=300)
    
    ## Central vs median vs mean for each scenario
    colors = ["#8dd3c7", "#ffffb3", "#bebada"]
    
    filename = hbv.root/"outputs"/"sampling"/f"{country}_central comp.pdf"
    pp = PdfPages(filename)
    
    for idx, scen in enumerate(scens.values()):
        fig_b = plt.figure(figsize=(15, 20))
        
        for sp, out in enumerate(outcomes):
            ax = fig_b.add_subplot(cols,rows,sp+1)
            ax.plot(np.arange(s_yr, e_yr,1),cent_plot[scen.name][int(sp*(e_yr-s_yr)):int((sp*(e_yr-s_yr))+(e_yr-s_yr)),0],  label='central', color = colors[0])
            ax.plot(np.arange(s_yr, e_yr,1),np.percentile(st_plot[scen.name][plot_names[sp]], 50, axis=1), label='median', color = colors[1])
            ax.plot(np.arange(s_yr, e_yr,1),np.mean(st_plot[scen.name][plot_names[sp]], axis=1), label='mean', color = colors[2])
            ax.set_title(out[3])
            ax.set_xlim(s_yr+10, e_yr-1)
            if sp==0:
                ax.legend(loc="best")
        
        fig_b.suptitle(scen.name, color="blue", size=12)
        fig_b.tight_layout()
        pp.savefig(fig_b)
    
    pp.close()


def pop_calibs(country: str, cal, res_cal):
    """Fit of model outputs compared to population data"""

    # Output name and PDF initiation
    filename = hbv.root/"outputs"/"calibrations"/f"{country}_agg cal_{res_cal.name}.pdf"
    pp = PdfPages(filename)

    # Required data and settings
    sns.set_theme(style="whitegrid")
    pop_bins, birth_bins, cons_pars, var_pars = hbv.age_bins("vimc_test")
    pops = list(cons_pars.columns)
    s_yr = hbv.s_yr
    e_yr = hbv.e_yr

    # Total Population Size
    fig1 = plt.figure(figsize=(20,20))

    tpop = np.zeros((e_yr-s_yr, len(pops)+1))

    for idx, pop in enumerate(pops):
        tpop[:, idx] = cal.get_par("alive").ts[pop].vals
    tpop[:, len(pops)] = np.sum(tpop[:, 0:len(pops)], axis=1)

    pop = fig1.add_subplot(1,1,1)
    pop.scatter(np.arange(s_yr, e_yr, 1), tpop[:, len(pops)], alpha=0.2, label="Data")
    pop.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, pops="total", outputs="alive", pop_aggregation="sum",
               t_bins=1).series[0].vals, color="red", label="Model Output")

    pop.legend(loc="best")
    pop.set_ylim(bottom=0)
    pop.set_ylabel("Population")

    fig1.suptitle(f"{country} Total Population", fontsize=25, style="normal", color="blue")
    pp.savefig(fig1)


    # Prevalence (population and under 5)
    fig2 = plt.figure(figsize=(20,20))
    u5_pops = ["year", "0-0M", "0-0F", "1-4M", "1-4F"]

    alive = pd.DataFrame(columns=u5_pops)
    prev = pd.DataFrame(columns=u5_pops)

    # Extract the data
    for pop in u5_pops:
        if pop == "year":
            alive[pop] = cal.get_par("alive").ts["0-0M"].t
            prev[pop] = cal.get_par("prev").ts["0-0M"].t
        else:
            alive[pop] = cal.get_par("alive").ts[pop].vals
            prev[pop] = cal.get_par("prev").ts[pop].vals

    # Keep only matching years and reset indexes
    years = list(prev["year"])
    alive = alive[alive["year"].isin(years)]
    alive = alive.reset_index()

    # Calculate >5y prevalence
    for pop in u5_pops:
        if pop!="year":
            prev[pop] = prev[pop] * alive[pop]

    prev["total"] = np.sum(prev.iloc[:, 1:], axis=1)
    alive["total"] = np.sum(alive.iloc[:, 1:], axis=1)
    prev["total"] = prev["total"]/alive["total"]

    child = fig2.add_subplot(2, 1, 1)
    child.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, "prev", pops={"U5": ["0-0M", "0-0F", "1-4M", "1-4F"]},
               pop_aggregation="weighted", t_bins=1).series[0].vals, color="red", label="Model Output")
    child.scatter(prev["year"], prev["total"], alpha=0.2, label="GBD Estimates")
    child.legend(loc="best")
    child.set_title("Under 5y")
    child.set_ylabel("HBsAg seroprevalence")

    pops.insert(0, "year")
    alive = pd.DataFrame(columns=pops)
    prev = pd.DataFrame(columns=pops)

    # Extract the data
    for pop in pops:
        if pop == "year":
            alive[pop] = cal.get_par("alive").ts["0-0M"].t
            prev[pop] = cal.get_par("prev").ts["0-0M"].t
        else:
            alive[pop] = cal.get_par("alive").ts[pop].vals
            prev[pop] = cal.get_par("prev").ts[pop].vals

    # Keep only matching years and reset indexes
    years = list(prev["year"])
    alive = alive[alive["year"].isin(years)]
    alive = alive.reset_index()

    # Calculate >5y prevalence
    for pop in pops:
        if pop != "year":
            prev[pop] = prev[pop] * alive[pop]

    prev["total"] = np.sum(prev.iloc[:, 1:], axis=1)
    alive["total"] = np.sum(alive.iloc[:, 1:], axis=1)
    prev["total"] = prev["total"] / alive["total"]

    total = fig2.add_subplot(2, 1, 2)
    total.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, "prev", pops="total",
               pop_aggregation="weighted", t_bins=1).series[0].vals, color="red", label="Model Output")
    total.scatter(prev["year"], prev["total"], alpha=0.2, label="GBD Estimates")
    total.legend(loc="best")
    total.set_title("Population")
    total.set_ylabel("HBsAg seroprevalence")

    pp.savefig(fig2)

    # HCC incidence
    fig3 = plt.figure(figsize=(20,20))
    pops = list(cons_pars.columns)
    thcc = np.zeros((30, len(pops) + 1))

    for idx, pop in enumerate(pops):
        thcc[:, idx] = cal.get_par("flw_hcc").ts[pop].vals
    thcc[:, len(pops)] = np.sum(thcc[:, 0:len(pops)], axis=1)

    hcc = fig3.add_subplot(1, 1, 1)
    hcc.scatter(np.arange(s_yr, s_yr+30, 1), thcc[:, len(pops)], alpha=0.2, label="Data")
    hcc.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, pops="total", outputs="flw_hcc", pop_aggregation="sum",
                                                   t_bins=1).series[0].vals, color="red", label="Model Output")

    hcc.legend(loc="best")
    hcc.set_ylim(bottom=0)
    hcc.set_ylabel("Incident Cases")
    fig3.suptitle(f"{country} Total HCC incidence", fontsize=25, style="normal", color="blue")

    pp.savefig(fig3)

    # Mortality (total, acute, cirrhosis, hcc)
    fig4 = plt.figure(figsize=(20,20))

    pops.insert(0, "year")
    acute = pd.DataFrame(columns=pops)
    cirr = pd.DataFrame(columns=pops)
    canc = pd.DataFrame(columns=pops)

    for pop in pops:
        if pop == "year":
            acute[pop] = cal.get_par("cl_acu").ts["0-0M"].t
            cirr[pop] = cal.get_par("cl_cir").ts["0-0M"].t
            canc[pop] = cal.get_par("cl_hcc").ts["0-0M"].t
        else:
            acute[pop] = cal.get_par("cl_acu").ts[pop].vals
            cirr[pop] = cal.get_par("cl_cir").ts[pop].vals
            canc[pop] = cal.get_par("cl_hcc").ts[pop].vals

    acute["total"] = np.sum(acute.iloc[:, 1:], axis=1)
    cirr["total"] = np.sum(cirr.iloc[:, 1:], axis=1)
    canc["total"] = np.sum(canc.iloc[:, 1:], axis=1)

    # Total mortality
    hbvm = fig4.add_subplot(2,2,1)
    hbvm.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, pops="total", outputs=":dd_hbv", pop_aggregation="sum",
              t_bins=1).series[0].vals, color="red", label="Model Output")
    #hbvm.scatter(acute["year"], acute["total"]+cirr["total"]+canc["total"], alpha=0.2, label="WHO GHE data")
    #hbvm.legend(loc="best")
    hbvm.set_ylabel("Annual Deaths")
    hbvm.set_title("Total HBV Mortality")


    # Acute HBV mortality
    acum = fig4.add_subplot(2,2,2)
    acum.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, pops="total", outputs="cl_acu", pop_aggregation="sum",
              t_bins=1).series[0].vals, color="red", label="Model Output")
    acum.scatter(acute["year"], acute["total"], alpha=0.2, label="WHO GHE data")
    acum.legend(loc="best")
    acum.set_title("Acute Mortality")
    acum.set_ylabel("Annual Deaths")


    # Cirrhosis mortality
    cirm = fig4.add_subplot(2,2,3)
    cirm.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, pops="total", outputs="cl_cir", pop_aggregation="sum",
              t_bins=1).series[0].vals, color="red", label="Model Output")
    cirm.scatter(cirr["year"], cirr["total"], alpha=0.2, label="WHO GHE data")
    cirm.set_title("Cirrhosis Mortality")
    cirm.set_ylabel("Annual Deaths")


    # HCC mortality
    hccm = fig4.add_subplot(2,2,4)
    hccm.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, pops="total", outputs="cl_hcc", pop_aggregation="sum",
              t_bins=1).series[0].vals, color="red", label="Model Output")
    hccm.scatter(canc["year"], canc["total"], alpha=0.2, label="WHO GHE data")
    hccm.set_title("HCC Mortality")
    hccm.set_ylabel("Annual Deaths")

    fig4.suptitle(f"{country} HBV mortality", fontsize=25, style="normal", color="blue")

    pp.savefig(fig4)

    # Rate analysis

    fig5 = plt.figure(figsize=(20,20))
    pops = list(cons_pars.columns)

    for idx,pop in enumerate(pops):
        f5=fig5.add_subplot(int(np.ceil(len(pops)/4)), 4, idx+1)
        f5.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs={"hcc_inc":"(flw_hcc/alive)*1e5"}, pops=pop,
                t_bins=1).series[0].vals, color="red", label="HCC incidence: Model")
        f5.plot(np.arange(s_yr, e_yr, 1), at.PlotData(res_cal, outputs={"hcc_inc": "(cl_hcc/alive)*1e5"}, pops=pop,
                t_bins=1).series[0].vals, color="blue", label="HCC mortality: Model")
        f5.set_title(pop)
        f5.set_ylabel("HCC Incidence Rate (per 100k)")
        f5.set_ylim(bottom=0)

        if idx==0:
            f5.legend(loc="best")

    fig5.suptitle(f"{country} HCC incidence rate: Model Projections", fontsize=25, style="normal", color="blue")
    fig5.tight_layout(pad=2)
    fig5.subplots_adjust(top=0.95)
    pp.savefig(fig5)

    pp.close()


def compare_outcomes(countries):  

    """Visualization of outcomes from modelled scenarios compared with old estimates"""  

    scenarios = ['Baseline', 'No Vaccination'] 
    ages = ['all', 'under5'] 
    outcomes = ['deaths', 'dalys'] 

    # Load the data 
    old_estimates_all = pd.read_csv(hbv.root/"inputs"/"old_estimates_data_by_cross_section.csv") 
    old_estimates_all = old_estimates_all[old_estimates_all['disease'] == 'HepB'] 

    for scen_name in scenarios:  
        # Output name and PDF initiation 
        filename = hbv.root/"outputs"/"comparisons"/f"{scen_name} outcomes comparisons.pdf"
        pp = PdfPages(filename) 

        # Process each country 
        for country_batch in [countries[i:i + 10] for i in range(0, len(countries), 10)]: 
            fig, axes = plt.subplots(len(country_batch), len(outcomes)*len(ages), figsize=(20, 5*len(country_batch))) 
            fig.suptitle(f"{scen_name}: Death and DALYs estimates comparison", fontsize=25) 

            for row, country in enumerate(country_batch): 
                # Filter old estimates for each country 
                old_estimate_country = old_estimates_all[old_estimates_all['country'] == country] 
                
                # Load and process current estimates for each country 
                current_estimates_all = pd.read_csv(hbv.root/"outputs"/"vimc_outs"/f"{country}_{scen_name}_stochastic.csv") 
                for col, (outcome, age) in enumerate(itertools.product(outcomes, ages)): 
                    ax = axes[row, col] if len(country_batch) > 1 else axes[col] 
                    
                    # Filter data by age group 
                    old_estimates = old_estimate_country[old_estimate_country['age_group'] == age] 
                    current_estimates = current_estimates_all.copy() 
                    if age == 'under5': 
                        current_estimates = current_estimates[current_estimates['age'] < 5] 
                  
                    # Plot old estimates 
                    if scen_name == 'Baseline':
                        ax.fill_between(old_estimates['year'], old_estimates[f'{outcome}_lo'], old_estimates[f'{outcome}_hi'], color='blue', alpha=0.3) 
                        ax.plot(old_estimates['year'], old_estimates[f'{outcome}']) 
                    else:
                        ax.fill_between(old_estimates['year'], old_estimates[f'{outcome}_nv_lo'], old_estimates[f'{outcome}_nv_hi'], color='blue', alpha=0.3) 
                        ax.plot(old_estimates['year'], old_estimates[f'{outcome}_no_vac'], label='Imperial') 
                
                    # Define bounds for current estimates
                    current_estimates_summed = current_estimates.groupby(['year','run_id'])[f'{outcome}'].sum().reset_index()
                    current_estimates_summary = current_estimates_summed.groupby('year')[f'{outcome}'].agg(['min', 'max'])
                    current_estimates_point = current_estimates_summed.groupby('year')[f'{outcome}'].median()
        
                    # Plot current estimates 
                    ax.fill_between(current_estimates_summary.index, current_estimates_summary['min'], current_estimates_summary['max'], color='red', alpha=0.3) 
                    ax.plot(current_estimates_summary.index, current_estimates_point, label='Current') 
                
                    # Set labels and titles 
                    ax.set_xlim([2000, 2030]) 
                    ax.tick_params(axis='both', which='major', labelsize=16)
                    if col == 0: 
                        ax.set_ylabel(f"{country}", fontsize=20) 
                    if row == 0: 
                        ax.set_title(f"{outcome} for {age}", fontsize=20) 
                    if row == 0 and col == len(outcomes)*len(ages)-1: 
                        ax.legend(loc='best', fontsize=16) 

            fig.tight_layout(pad=2) 
            fig.subplots_adjust(top=.95) 
            pp.savefig(fig) 
            plt.close(fig)  # Close the figure after saving 
        pp.close() 
    
        print(f"PDF saved: {filename}")

def compare_scenarios(countries): 
    """Visualization of scenario comparisons for each country""" 

    scenarios = ['Baseline', 'No Vaccination'] 
    models = ['Imperial', 'Current']
    ages = ['under5', 'all'] 
    outcomes = ['deaths'] 

    # Load the data 
    old_estimates_all = pd.read_csv(hbv.root/"inputs"/"old_estimates_data_by_cross_section.csv") 
    old_estimates_all = old_estimates_all[old_estimates_all['disease'] == 'HepB'] 

    # Output name and PDF initiation 
    filename = hbv.root/"outputs"/"comparisons"/"Scenario deaths differences (trends).pdf" 
    pp = PdfPages(filename) 

    # Process countries in batches of 10 
    for country_batch in [countries[i:i + 10] for i in range(0, len(countries), 10)]:         
        fig, axes = plt.subplots(len(country_batch), len(ages) * 2, figsize=(20, 5 * len(country_batch)))  # 4 columns per country 
        fig.suptitle("Deaths Comparisons per Country", fontsize=25) 

        for row, country in enumerate(country_batch): 
            # Global min and max for each country & age group 
            global_min_max = {} 
            for age in ages: 
                # Initialize with large and small numbers 
                old_estimates = old_estimates_all[old_estimates_all['country'] == country] 
                old_estimates = old_estimates[old_estimates['age_group'] == age] 
                old_estimates = old_estimates[(old_estimates['year'] >= 2000) & (old_estimates['year'] <= 2030)] 
                current_estimates_baseline = pd.read_csv(hbv.root/"outputs"/"vimc_outs"/f"{country}_{scenarios[0]}_stochastic.csv") 
                current_estimates_novax = pd.read_csv(hbv.root/"outputs"/"vimc_outs"/f"{country}_{scenarios[1]}_stochastic.csv") 
        
                # Calculate min and max for old estimates 
                min_val_old = old_estimates[[f'{outcomes[0]}_lo', f'{outcomes[0]}_nv_lo']].min().min() 
                max_val_old = old_estimates[[f'{outcomes[0]}_hi', f'{outcomes[0]}_nv_hi']].max().max() 
                
                # Current estimates calculations 
                current_estimates = current_estimates_baseline[current_estimates_baseline['age'] < 5] if age == 'under5' else current_estimates_baseline 
                current_estimates = current_estimates[(current_estimates['year'] >= 2000) & (current_estimates['year'] <= 2030)] 
                current_estimates_summed = current_estimates.groupby(['year','run_id'])[outcomes[0]].sum().reset_index() 
                min_val_current_baseline = current_estimates_summed[f'{outcomes[0]}'].min() 
                max_val_current_baseline = current_estimates_summed[f'{outcomes[0]}'].max() 
                current_estimates = current_estimates_novax[current_estimates_novax['age'] < 5] if age == 'under5' else current_estimates_novax 
                current_estimates = current_estimates[(current_estimates['year'] >= 2000) & (current_estimates['year'] <= 2030)] 
                current_estimates_summed = current_estimates.groupby(['year','run_id'])[outcomes[0]].sum().reset_index() 
                min_val_current_novax = current_estimates_summed[f'{outcomes[0]}'].min() 
                max_val_current_novax = current_estimates_summed[f'{outcomes[0]}'].max() 

                global_min = min(min_val_old, min_val_current_baseline, min_val_current_novax) 
                global_max = max(max_val_old, max_val_current_baseline, max_val_current_novax) 
                global_min_max[age] = (global_min, global_max) 
            
            # Filter old estimates for the country 
            old_estimate_country = old_estimates_all[old_estimates_all['country'] == country] 

            for col, age in enumerate(ages): 
                for j, model in enumerate(models): 
                    old_estimates = old_estimate_country[old_estimate_country['age_group'] == age] 
                    current_estimates_baseline = pd.read_csv(hbv.root/"outputs"/"vimc_outs"/f"{country}_{scenarios[0]}_stochastic.csv") 
                    current_estimates_novax = pd.read_csv(hbv.root/"outputs"/"vimc_outs"/f"{country}_{scenarios[1]}_stochastic.csv") 
                    ax = axes[row, col * 2 + j] 
                    if model == 'Imperial':
                        # Plot old estimates 
                        ax.fill_between(old_estimates['year'], old_estimates[f'{outcomes[0]}_lo'], old_estimates[f'{outcomes[0]}_hi'], color='blue', alpha=0.3) 
                        ax.plot(old_estimates['year'], old_estimates[f'{outcomes[0]}'], label='Baseline') 
                        ax.fill_between(old_estimates['year'], old_estimates[f'{outcomes[0]}_nv_lo'], old_estimates[f'{outcomes[0]}_nv_hi'], color='blue', alpha=0.3) 
                        ax.plot(old_estimates['year'], old_estimates[f'{outcomes[0]}_no_vac'], label='No Vaccination') 
                    elif model == 'Current':
                        # Plot current estimates 
                        current_estimates = current_estimates_baseline[current_estimates_baseline['age'] < 5] if age == 'under5' else current_estimates_baseline 
                        current_estimates_summed = current_estimates.groupby(['year','run_id'])[outcomes[0]].sum() .reset_index()
                        current_estimates_summary = current_estimates_summed.groupby('year')[f'{outcomes[0]}'].agg(['min', 'max'])
                        current_estimates_point = current_estimates_summed.groupby('year')[f'{outcomes[0]}'].median()
            
                        ax.fill_between(current_estimates_summary.index, current_estimates_summary['min'], current_estimates_summary['max'], color='red', alpha=0.3) 
                        ax.plot(current_estimates_summary.index, current_estimates_point, label='Baseline') 
                        
                        current_estimates = current_estimates_novax[current_estimates_novax['age'] < 5] if age == 'under5' else current_estimates_novax 
                        current_estimates_summed = current_estimates.groupby(['year','run_id'])[outcomes[0]].sum() .reset_index()
                        current_estimates_summary = current_estimates_summed.groupby('year')[f'{outcomes[0]}'].agg(['min', 'max'])
                        current_estimates_point = current_estimates_summed.groupby('year')[f'{outcomes[0]}'].median()
            
                        ax.fill_between(current_estimates_summary.index, current_estimates_summary['min'], current_estimates_summary['max'], color='red', alpha=0.3) 
                        ax.plot(current_estimates_summary.index, current_estimates_point, label='No Vaccination') 
                    
                    # Set labels and titles 
                    ax.set_xlim([2000, 2030]) 
                    ax.set_ylim(global_min_max[age])
                    ax.tick_params(axis='both', which='major', labelsize=16)
                    if col == 0 and j == 0: 
                        ax.set_ylabel(f"{country}", fontsize=20) 
                    if row == 0: 
                        ax.set_title(f"{age.capitalize()} - {model}", fontsize=20) 
                    if row == 0 and col == 1 and j == 1: 
                        ax.legend(loc='best', fontsize=16) 
        
        fig.tight_layout(pad=2) 
        fig.subplots_adjust(top=.95) 
        pp.savefig(fig) 
        plt.close(fig) 

    pp.close()
    print(f"PDF saved: {filename}")

def difference_scenarios(countries):
    """Visualization of difference in deaths for each country""" 
    ages = ['under5', 'all'] 
    outcomes = ['deaths'] 

    # Load the data 
    old_estimates_all = pd.read_csv(hbv.root/"inputs"/"old_estimates_data_by_cross_section.csv") 
    old_estimates_all = old_estimates_all[old_estimates_all['disease'] == 'HepB'] 

    # Output name and PDF initiation 
    filename = hbv.root/"outputs"/"comparisons"/"Scenario deaths differences (scatter, table).pdf" 
    pp = PdfPages(filename) 
    
    outlier_threshold = 0.2*100 # highlight any point with an absolute difference that is greater than 20%
    
    old_estimates_diff = {age: {} for age in ages}
    current_estimates_diff = {age: {} for age in ages}
    difference = {age: {} for age in ages}
    data = {age: [] for age in ages}
    
    for country in countries:
        for age in ages:
            old_estimate_country = old_estimates_all[old_estimates_all['country'] == country] 
            old_estimates = old_estimate_country[(old_estimate_country['age_group'] == age) & (old_estimate_country['year'] == 2030)] 
            old_estimates_baseline = old_estimates[f'{outcomes[0]}']
            old_estimates_novax = old_estimates[f'{outcomes[0]}_no_vac']
            old_estimates_diff[age][country] = (old_estimates_novax.values-old_estimates_baseline.values)/old_estimates_novax.values*100
            
            current_estimates_baseline = pd.read_csv(hbv.root/"outputs"/"vimc_outs"/f"{country}_Baseline_stochastic.csv") 
            current_estimates_novax = pd.read_csv(hbv.root/"outputs"/"vimc_outs"/f"{country}_No Vaccination_stochastic.csv") 
            current_estimates_baseline = current_estimates_baseline[current_estimates_baseline['age'] < 5] if age == 'under5' else current_estimates_baseline 
            current_estimates_novax = current_estimates_novax[current_estimates_novax['age'] < 5] if age == 'under5' else current_estimates_novax 
            current_estimates_baseline = current_estimates_baseline[current_estimates_baseline['year'] == 2030] 
            current_estimates_novax = current_estimates_novax[current_estimates_novax['year'] == 2030] 
            
            current_estimates_baseline_summed = current_estimates_baseline.groupby(['year','run_id'])[outcomes[0]].sum().reset_index()
            current_estimates_baseline = current_estimates_baseline_summed.groupby('year')[f'{outcomes[0]}'].median()
            current_estimates_novax_summed = current_estimates_novax.groupby(['year','run_id'])[outcomes[0]].sum().reset_index()
            current_estimates_novax = current_estimates_novax_summed.groupby('year')[f'{outcomes[0]}'].median()
    
            current_estimates_diff[age][country] = (current_estimates_novax.values-current_estimates_baseline.values)/current_estimates_novax.values*100
            
            difference[age][country] = abs(current_estimates_diff[age][country] - old_estimates_diff[age][country])
           
            data[age].append([country, age, f"{current_estimates_diff[age][country][0]:.0f}%", f"{old_estimates_diff[age][country][0]:.0f}%", f"{difference[age][country][0]:.0f}%"])
        
    fig, ax = plt.subplots(figsize=(20,20))
    colors = ['m', 'k']
    for i, age in enumerate(ages):
        ax.scatter(150, 150, color=colors[i], label=age)
        for country in countries:
            x = current_estimates_diff[age][country]
            y = old_estimates_diff[age][country]
            ax.scatter(x, y, c=colors[i])
            if difference[age][country] > outlier_threshold:
                ax.annotate(country, (x, y), textcoords="offset points", xytext=(0,10), ha='center', color=colors[i])
        ax.set_xlim([0,100])
        ax.set_ylim([0,100])
        ax.legend(loc='best', fontsize=16) 
        ax.xaxis.set_major_formatter(mticker.PercentFormatter(xmax=100))
        ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=100))
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel("Current estimates", fontsize=20)
        ax.set_ylabel("Imperial estimates", fontsize=20)
        add_identity(ax, ls="--", c=".3")
    
    fig.suptitle("Reduction in HBV-related deaths in 2030 due to vaccine uptake", fontsize=25, style="normal", color="blue")
    pp.savefig(fig)
    plt.close(fig)

    add_table_to_pdf(pp, data, ages)
    
    pp.close()
    
    print(f"PDF saved: {filename}")

def add_table_to_pdf(pdf, table_data, ages, max_rows_per_page=30): 
    for age in ages: 
        sorted_data = sorted(table_data[age], key=lambda x: abs(float(x[-1].strip('%'))), reverse=True)
        num_pages = math.ceil(len(sorted_data) / max_rows_per_page)
        
        for page in range(num_pages): 
            start_row = page * max_rows_per_page 
            end_row = min(start_row + max_rows_per_page, len(sorted_data)) 
            table_data_page = sorted_data[start_row:end_row]

            if table_data_page: 
                fig, ax = plt.subplots(figsize=(20,10))
                ax.axis('tight') 
                ax.axis('off') 
                table = ax.table(cellText=table_data_page, colLabels=['Country', 'Age groups', 'Current reduction in deaths', 'Imperial reduction in deaths', 'Absolute difference between'], cellLoc='center', loc='center') 
                table.auto_set_font_size(False) 
                table.set_fontsize(14) 
                table.scale(1.5, 5) 

                plt.subplots_adjust(left=0.2, right=0.8, top=0.6, bottom=0.4) 
                pdf.savefig(fig) 
                plt.close(fig) 
     

def add_identity(axes, *line_args, **line_kwargs):
    identity, = axes.plot([], [], *line_args, **line_kwargs)
    def callback(axes):
        low_x, high_x = axes.get_xlim()
        low_y, high_y = axes.get_ylim()
        low = max(low_x, low_y)
        high = min(high_x, high_y)
        identity.set_data([low, high], [low, high])
    callback(axes)
    axes.callbacks.connect('xlim_changed', callback)
    axes.callbacks.connect('ylim_changed', callback)
    return axes


def comparison_plots(countries):
    compare_outcomes(countries)
    compare_scenarios(countries)
    difference_scenarios(countries)
    



    
