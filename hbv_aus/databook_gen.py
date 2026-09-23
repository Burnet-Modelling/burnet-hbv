import atomica as at
import numpy as np
import pandas as pd

from hbv_aus.claude_databook_fill import hepbd_cov_vals
from hbv_aus.utils import _get_github_folder,_get_sharepoint_folder

# Version 1.0 (basic population dynamics, no behaviour populations)
def make_age_bands(max_age=84, band_width=5):
    bands = [f"{i}-{i+band_width-1}" for i in range(0, max_age+1, band_width)]
    bands.append(f"{max_age+1}+")
    return bands


def gen_db_poptest():
    """
    Generate a databook only for population testing - easier than trying to bug solve with a million different
    parameters! For now can just fill manually.
    """
    FW_PATH = _get_github_folder() + f"/framework/fw_popsizes.xlsx"
    DB_PATH = _get_github_folder() + f"/databook/popsize_test_210926.xlsx"

    F = at.ProjectFramework(FW_PATH)

    age_bins = ["0-14", "15-64", "65+"]# only for use in testing, will be hard coded in practice
    sexes = ["_M", "_F"]
    demog = ["atsi_", "ausb_", "lros_", "hros_"] # Aboriginal Torres Strait Islander, Aus Born, Low Risk Overseas, High Risk Overseas
    pops = [f"{demo}{age}{sex}" for demo in demog for sex in sexes for age in age_bins]

    D = at.ProjectData.new(framework=F, pops=1, transfers=0, tvec=np.arange(1980, 2072, 1))
    for idx, val in enumerate (pops):
        if idx == 0:
            D.rename_pop("pop_0", new_code_name=val, new_full_name=val)
        else:
            D.add_pop(val, val)
    D.add_transfer("age", "aging")

    D.save(DB_PATH)


def gen_db_hbv_v2_1():
    """
    Produces databook using hbv_fw_v2.1.xlsx (HepAus model informed by HepB-BD in SA model)
    """
    # Import Framework to produce databook
    FW_PATH = _get_github_folder() + f"/framework/hbv_fw_v2.1_autosave.xlsx"
    DB_PATH  = _get_github_folder() + f"/databook/hbv_db_hepaus_220926.xlsx"
    DATA_PATH  = _get_github_folder() + f"input data/"

    F = at.ProjectFramework(FW_PATH)

    age_bins = ["0-14", "15-64", "65+"]  # only for use in testing, will be hard coded in practice
    sexes = ["_M", "_F"]
    demog = ["atsi_", "ausb_", "lros_",
             "hros_"]  # Aboriginal Torres Strait Islander, Aus Born, Low Risk Overseas, High Risk Overseas
    pops = [f"{demo}{age}{sex}" for demo in demog for sex in sexes for age in age_bins]

    D = at.ProjectData.new(framework=F, pops=1, transfers=0, tvec=np.arange(1980, 2072, 1))
    for idx, val in enumerate(pops):
        if idx == 0:
            D.rename_pop("pop_0", new_code_name=val, new_full_name=val)
        else:
            D.add_pop(val, val)

    D.add_transfer("age", "aging")

    # Population Data (total_pop, mig_rate, acm, and j_init)
    pd.ExcelFile(DATA_PATH+"populations_in.xlsx").sheet_names
    pop_size = pd.read_excel(DATA_PATH+"populations_in.xlsx", sheet_name="pop_size")
    acm = pd.read_excel(DATA_PATH+"populations_in.xlsx", sheet_name="acm")
    migration = pd.read_excel(DATA_PATH+"populations_in.xlsx", sheet_name="migration")
    for pop in pops:
        D.tdve["total_pop"].ts[pop] = at.TimeSeries(t=pop_size.year, vals=pop_size[pop], units="Number")
        D.tdve["acm"].ts[pop] = at.TimeSeries(t=acm.year, vals=acm[pop], units="Probability (per year)")
        D.tdve["mig_rate"].ts[pop] = at.TimeSeries(t=migration.year, vals=migration[pop], units="N.A.")
        D.tdve["j_init"].ts[pop] = at.TimeSeries(t=pop_size["year"].iloc[0], vals = pop_size[pop].iloc[0], units = "Number")

    # Births, Pregnancies and MTCT interaction matrix
    b_rate = pd.read_excel(DATA_PATH+"populations_in.xlsx", sheet_name="births")
    for pop in pops:
        if pop in b_rate.columns:
            D.tdve["b_rate"].ts[pop] = at.TimeSeries(t=b_rate.year, vals=b_rate[pop], units="N.A.")
        else:
            D.tdve["b_rate"].ts[pop].assumption = 0

    # atsi births are only divided among 15-64 ATSI females, ausb are among all other 15-64 female pops (pregs)
    atsi_preg = (b_rate["atsi_0-14_M"] + b_rate["atsi_0-14_F"])/pop_size["atsi_15-64_F"]
    ausb_preg = (b_rate["ausb_0-14_M"] + b_rate["ausb_0-14_F"])/(pop_size["ausb_15-64_F"]+pop_size["hros_15-64_F"]+pop_size["lros_15-64_F"])

    for pop in pops:
        if pop == "atsi_15-64_F":
            D.tdve["pregs"].ts[pop] = at.TimeSeries(t=np.arange(1980, 2071, 1), vals = atsi_preg, units = "N.A.")
        elif pop == "ausb_15-64_F" or pop == "hros_15-64_F" or pop == "lros_15-64_F":
            D.tdve["pregs"].ts[pop] = at.TimeSeries(t=np.arange(1980, 2071, 1), vals=ausb_preg, units="N.A.")
        else:
            D.tdve["pregs"].ts[pop].assumption =  0

    # Interaction matrices for atsi births are 1, for ausb will be proportional to 15-64 pop size at time (t)
    interaction_mtct = next(ip for ip in D.interpops if ip.code_name=="mx_mtct")
    years = list(np.arange(1980, 2071, 1))
    b_prop_ausb = list(pop_size["ausb_15-64_F"] / (pop_size["ausb_15-64_F"]+pop_size["hros_15-64_F"]+pop_size["lros_15-64_F"]))
    b_prop_hros = list(pop_size["hros_15-64_F"]/ (pop_size["ausb_15-64_F"]+pop_size["hros_15-64_F"]+pop_size["lros_15-64_F"]))
    b_prop_lros = list(pop_size["lros_15-64_F"]/ (pop_size["ausb_15-64_F"]+pop_size["hros_15-64_F"]+pop_size["lros_15-64_F"]))

    mtct_edges = [
        ("atsi_0-14_M", "atsi_0-14_M", 0),
        ("atsi_15-64_F", "atsi_0-14_M", 1.0),
        ("atsi_15-64_F", "atsi_0-14_F", 1.0),

        ("ausb_15-64_F", "ausb_0-14_M", (years, b_prop_ausb)),
        ("ausb_15-64_F", "ausb_0-14_F", (years, b_prop_ausb)),
        ("hros_15-64_F", "ausb_0-14_M", (years, b_prop_hros)),
        ("hros_15-64_F", "ausb_0-14_F", (years, b_prop_hros)),
        ("lros_15-64_F", "ausb_0-14_M", (years, b_prop_lros)),
        ("lros_15-64_F", "ausb_0-14_F", (years, b_prop_lros))]

    for from_pop, to_pop, weight in mtct_edges:
        if isinstance(weight, tuple):
            y, v = weight
            ts = at.TimeSeries(t=y, vals=v, units="N.A.")
        else:
            ts = at.TimeSeries(assumption=weight, units="N.A.")
        interaction_mtct.ts[(from_pop, to_pop)] = ts

    # Interaction Matrix: Horizontal Transmission (assume within population only)
    interaction_horiz = next(ip for ip in D.interpops if ip.code_name=="mx_horiz")

    for age_f in age_bins:
        for age_t in age_bins:
            # ATSI
            interaction_horiz.ts[(f"atsi_{age_f}_M", f"atsi_{age_t}_M")] = at.TimeSeries(assumption = 1.0, units = "N.A.")
            interaction_horiz.ts[(f"atsi_{age_f}_F", f"atsi_{age_t}_F")] = at.TimeSeries(assumption = 1.0, units = "N.A.")
            interaction_horiz.ts[(f"atsi_{age_f}_F", f"atsi_{age_t}_M")] = at.TimeSeries(assumption = 1.0, units = "N.A.")
            interaction_horiz.ts[(f"atsi_{age_f}_M", f"atsi_{age_t}_F")] = at.TimeSeries(assumption = 1.0, units = "N.A.")
            # Aus Born
            interaction_horiz.ts[(f"ausb_{age_f}_M", f"ausb_{age_t}_M")] = at.TimeSeries(assumption = 1.0, units = "N.A.")
            interaction_horiz.ts[(f"ausb_{age_f}_F", f"ausb_{age_t}_F")] = at.TimeSeries(assumption = 1.0, units = "N.A.")
            interaction_horiz.ts[(f"ausb_{age_f}_F", f"ausb_{age_t}_M")] = at.TimeSeries(assumption = 1.0, units = "N.A.")
            interaction_horiz.ts[(f"ausb_{age_f}_M", f"ausb_{age_t}_F")] = at.TimeSeries(assumption = 1.0, units = "N.A.")
            # HROS
            interaction_horiz.ts[(f"hros_{age_f}_M", f"hros_{age_t}_M")] = at.TimeSeries(assumption = 1.0, units = "N.A.")
            interaction_horiz.ts[(f"hros_{age_f}_F", f"hros_{age_t}_F")] = at.TimeSeries(assumption = 1.0, units = "N.A.")
            interaction_horiz.ts[(f"hros_{age_f}_F", f"hros_{age_t}_M")] = at.TimeSeries(assumption = 1.0, units = "N.A.")
            interaction_horiz.ts[(f"hros_{age_f}_M", f"hros_{age_t}_F")] = at.TimeSeries(assumption = 1.0, units = "N.A.")
            #LROS
            interaction_horiz.ts[(f"lros_{age_f}_M", f"lros_{age_t}_M")] = at.TimeSeries(assumption=1.0, units="N.A.")
            interaction_horiz.ts[(f"lros_{age_f}_F", f"lros_{age_t}_F")] = at.TimeSeries(assumption=1.0, units="N.A.")
            interaction_horiz.ts[(f"lros_{age_f}_F", f"lros_{age_t}_M")] = at.TimeSeries(assumption=1.0, units="N.A.")
            interaction_horiz.ts[(f"lros_{age_f}_M", f"lros_{age_t}_F")] = at.TimeSeries(assumption=1.0, units="N.A.")

    # Population Transfer (Aging)
    age_from, age_to = age_bins[:-1], age_bins[1:]
    aging = pd.read_excel(DATA_PATH+"populations_in.xlsx", sheet_name="aging")

    for i, age in enumerate(age_from):
        for dem in demog:
            for sex in sexes:
                D.transfers[0].ts.append((f"{dem}{age}{sex}", f"{dem}{age_to[i]}{sex}"),
                                         at.TimeSeries(t = aging.year,
                                                       vals = aging[f"{dem}{age}{sex}"],
                                                       units="Rate (per year)"))

    # Disease Progression & Treatment Effectiveness Parameters

    # Disease Progression Parameters
    pd.ExcelFile(DATA_PATH+"epidemiology_in.xlsx").sheet_names
    dis_prog = pd.read_excel(DATA_PATH+"epidemiology_in.xlsx", sheet_name="dispars_in")
    dis_prog = dis_prog[dis_prog['est']=="p.e"]
    dis_pars = list(dis_prog.par)

    for par in dis_pars:
        dis_temp = dis_prog[dis_prog["par"]==par]
        for pop in pops:
            D.tdve[par].ts[pop].assumption = dis_temp[pop]

    # Treatment Effectiveness Parameters
    treat_eff = pd.read_excel(DATA_PATH+"epidemiology_in.xlsx", sheet_name="treat eff_in")
    te_pars = list(treat_eff.par)

    for par in te_pars:
        te_temp = treat_eff[treat_eff["par"]==par]
        for pop in pops:
            D.tdve[par].ts[pop].assumption=te_temp["p.e."]

    # MTCT parameters (fixed)
    mtct_tiv = pd.read_excel(DATA_PATH+"epidemiology_in.xlsx", sheet_name="mtct_in")
    bir_pops = ["atsi_0-14_M", "atsi_0-14_F", "ausb_0-14_M", "ausb_0-14_F"]
    mtct_pars = list(mtct_tiv.par)
    moth_pars = ["nullipar", "preg_ltc_rate"]
    moth_pops = ["atsi_15-64_F", "ausb_15-64_F", "hros_15-64_F", "lros_15-64_F"]

    set_moth = set(moth_pars)
    mtct_pars = [item for item in mtct_pars if item not in set_moth]

    for par in mtct_pars:
        mtct_temp = mtct_tiv[mtct_tiv["par"]==par]
        for pop in pops:
            if pop in bir_pops:
                D.tdve[par].ts[pop].assumption = mtct_temp["p.e."]
            else:
                D.tdve[par].ts[pop].assumption = 0

    for par in moth_pars:
        moth_temp = mtct_tiv[mtct_tiv["par"]==par]
        for pop in pops:
            if pop in moth_pops:
                D.tdve[par].ts[pop].assumption = moth_temp["p.e."]
            else:
                D.tdve[par].ts[pop].assumption = 0

    # MTCT Pars - Time Varying (hepbd coverage, anc screening, hepb3 coverage)
    mtct_tv = pd.read_excel(DATA_PATH+"epidemiology_in.xlsx", sheet_name="mtct_tvar")

    # Antenatal Screening Coverage
    for pop in pops:
        anc_scr = mtct_tv[mtct_tv["par"]=="anc_scr"]
        if pop in moth_pops:
            D.tdve["anc_scr"].ts[pop] = at.TimeSeries(t=[1980,1999,2000,2025], vals = anc_scr.iloc[:,1:].values.flatten().tolist(),
                                                      units = "N.A.")
        else:
            D.tdve["anc_scr"].ts[pop].assumption = 0

    # HepB-BD and HepB3 coverage
    for pop in pops:
        hepbd_cov = mtct_tv[mtct_tv["par"]=="hepbd_cov"]
        hepb3_cov = mtct_tv[mtct_tv["par"]=="hepb3_cov"]

        if pop in bir_pops:
            D.tdve["hepbd_cov"].ts[pop] = at.TimeSeries(t=[1980, 1999, 2000, 2025],
                                                      vals=hepbd_cov.iloc[:, 1:].values.flatten().tolist(),
                                                      units="N.A.")
            D.tdve["hepb3_cov"].ts[pop] = at.TimeSeries(t=[1980, 1999, 2000, 2025],
                                                      vals=hepb3_cov.iloc[:, 1:].values.flatten().tolist(),
                                                      units="N.A.")
        else:
            D.tdve["hepbd_cov"].ts[pop].assumption = 0
            D.tdve["hepb3_cov"].ts[pop].assumption = 0

    # Care Cascade

    # Assuming no loss-to-follow up or treatment stopping (captured in net flow rates)
    for pop in pops:
        D.tdve["ts_rate"].ts[pop].assumption = 0
        D.tdve["ltf_rate"].ts[pop].assumption = 0

    care_data = pd.read_excel(DATA_PATH+"care cascade_in.xlsx")
    no_trt_pops = ["atsi_0-14_M", "atsi_0-14_F", "ausb_0-14_M", "ausb_0-14_F", "hros_0-14_M", "hros_0-14_F", "lros_0-14_M", 'lros_0-14_F']
    set_notrt = set(no_trt_pops)
    trt_pops = [item for item in pops if item not in set_notrt]

    for pop in pops:
        D.tdve["diag_measure"].ts[pop] = at.TimeSeries(t=[1980, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024],
                                                       vals = care_data[care_data.par=="diag_measure"].iloc[:,1:].values.flatten().tolist(),
                                                       units = "N.A.")
        D.tdve["ltc_measure"].ts[pop]= at.TimeSeries(t=[1980, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024],
                                                       vals = care_data[care_data.par=="ltc_measure"].iloc[:,1:].values.flatten().tolist(),
                                                       units = "N.A.")
        D.tdve["diagnosed_meas"].ts[pop] = at.TimeSeries(t=[1980, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024],
                                                       vals = care_data[care_data.par=="diagnosed_meas"].iloc[:,1:].values.flatten().tolist(),
                                                       units = "Fraction")
        D.tdve["linked_meas"].ts[pop] = at.TimeSeries(t=[1980, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024],
                                                       vals = care_data[care_data.par=="linked_meas"].iloc[:,1:].values.flatten().tolist(),
                                                       units = "Fraction")
        if pop in trt_pops:
            D.tdve["treat_measure"].ts[pop] = at.TimeSeries(t=[1980, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024],
                                                            vals=care_data[care_data.par == "treat_measure"].iloc[:, 1:].values.flatten().tolist(),
                                                            units="N.A.")
            D.tdve["treat_cov"].ts[pop] = at.TimeSeries(t=[1980, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024],
                                                            vals=care_data[care_data.par == "treat_cov"].iloc[:, 1:].values.flatten().tolist(),
                                                            units="Fraction")
        else:
            D.tdve["treat_measure"].ts[pop] = at.TimeSeries(
                t=[1980, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024],
                vals=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                units="N.A.")
            D.tdve["treat_cov"].ts[pop] = at.TimeSeries(
                t=[1980, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024],
                vals=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                units="Fraction")

    # Initial Conditions & Calibration Values (informed by HepB Mapping Model)
    init_data = pd.read_excel(DATA_PATH+"init_calib_in.xlsx", sheet_name = "init_conds")
    calib_data = pd.read_excel(DATA_PATH+"init_calib_in.xlsx", sheet_name = "calib_pars")

    # Initial Conditions
    init_pars = list(init_data.par)
    for par in init_pars:
        temp = init_data[init_data.par==par]
        if par == "init_clr":
            for pop in pops:
                 D.tdve[par].ts[pop] = at.TimeSeries(t=[1980], vals=temp[pop].iloc[0], units = "Proportion")
        else:
            for pop in pops:
                 D.tdve[par].ts[pop] = at.TimeSeries(t=[1980], vals=temp[pop].iloc[0], units = "N.A.")


    # Calibration Data
    cal_pars = list(pd.unique(calib_data.par))

    for par in cal_pars:
        temp = calib_data[calib_data.par==par]
        if par == "hep_dth":
            for pop in pops:
                D.tdve[par].ts[pop] = at.TimeSeries(t=temp["year"], vals=temp[pop], units="N.A.")
        else:
            for pop in pops:
                D.tdve[par].ts[pop] = at.TimeSeries(t=temp["year"], vals=temp[pop], units="Fraction")

    D.save(DB_PATH)





    #D.validate(F)











    F = at.ProjectFramework(_get_github_folder() + f"framework/hbv_fw_v2.1_autosave.xlsx")

    age_bins = ["0-14", "15-64", "65+"]# only for use in testing, will be hard coded in practice
    sexes = ["_M", "_F"]
    demog = ["atsi_", "ausb_", "lros_", "hros_"] # Aboriginal Torres Strait Islander, Aus Born, Low Risk Overseas, High Risk Overseas
    pops = [f"{demo}{age}{sex}" for demo in demog for sex in sexes for age in age_bins]

    # Generate Databook
    D = at.ProjectData.new(framework=F, pops=1, transfers=0, tvec=np.arange(1980, 2072, 1))
    for idx, val in enumerate (pops):
        if idx == 0:
            D.rename_pop("pop_0", new_code_name=val, new_full_name=val)
        else:
            D.add_pop(val, val)
    D.add_transfer("age", "aging")

    ### Population Demographics ###

    #pop_data = # Add excel sheet as reference (save in SharePoint)

    # Total Population Size

    # Births (also used for estimating preg_rate for MTCT)

    # All Cause Mortality

    # Emigration Rate (0 for ATSI)

    ### Initial Conditions ###

    ### Immigration Model (note: numbers handled in population demographics) ###

    ### MTCT Model (inc. interaction matrix) ###

    ## Horizontal Transmission Matrix ###

    ### Natural History and Effectiveness Values ###

    ### Care Cascade ###

    ### Calibration ###


    D.save(_get_github_folder() + f"databook/hbv_hepaus_db.xlsx")






    # Import Data from Excel
    # Aboriginal Torres Strait Islander Population Data
    atsi_pop_in = pd.read_excel(_get_sharepoint_folder()+f"atsi_pop_proj.xlsx", sheet_name="Total Population")
    atsi_death_in = pd.read_excel(_get_sharepoint_folder()+f"atsi_pop_proj.xlsx", sheet_name="Deaths")
    atsi_birth_in = pd.read_excel(_get_sharepoint_folder()+f"atsi_pop_proj.xlsx", sheet_name="Births")
    atsi_age_in = pd.read_excel(_get_sharepoint_folder()+f"atsi_pop_proj.xlsx", sheet_name="Aging")

    # Migrants

    # Other Australian Pop


    # Total Population (temp_alive)
    for age in age_bins:
        D.tdve["temp_alive"].ts[f"atsi_{age}_M"] = at.TimeSeries(t=atsi_pop_in[atsi_pop_in["pop"]==age]["Year"],
                                                                 vals = atsi_pop_in[atsi_pop_in["pop"]==age]["Male"],
                                                                 units="Number")

        D.tdve["temp_alive"].ts[f"atsi_{age}_F"] = at.TimeSeries(t=atsi_pop_in[atsi_pop_in["pop"] == age]["Year"],
                                                                 vals=atsi_pop_in[atsi_pop_in["pop"] == age]["Female"],
                                                                 units="Number")
    # All Cause Mortality (death_rate):
    for age in age_bins:
        D.tdve["death_rate"].ts[f"atsi_{age}_M"] = at.TimeSeries(t=atsi_death_in[atsi_death_in["pop"]==age]["Year"],
                                                                 vals = atsi_death_in[atsi_death_in["pop"]==age]["Male"],
                                                                 units="Rate (per year)")

        D.tdve["death_rate"].ts[f"atsi_{age}_F"] = at.TimeSeries(t=atsi_death_in[atsi_death_in["pop"]==age]["Year"],
                                                                 vals = atsi_death_in[atsi_death_in["pop"]==age]["Female"],
                                                                 units="Rate (per year)")
    # Births (birth_rate):
    #atsi_births_sum = atsi_birth_in.iloc[:, np.r_[0, 4,5]]
    #atsi_births_sum = atsi_births_sum.groupby('Year', as_index=False).sum()

    for age in age_bins:
        if age == "0-4":
            D.tdve["birth_rate"].ts[f"atsi_{age}_M"] = at.TimeSeries(t=atsi_birth_in["Year"], vals = atsi_birth_in["Male"], units = "Number (per year)")
            D.tdve["birth_rate"].ts[f"atsi_{age}_F"] = at.TimeSeries(t=atsi_birth_in["Year"], vals = atsi_birth_in["Female"], units = "Number (per year)")
        else:
            D.tdve["birth_rate"].ts[f"atsi_{age}_M"].assumption = 0
            D.tdve["birth_rate"].ts[f"atsi_{age}_F"].assumption = 0

    # Migration (emig_rate, imig_rate)
    for age in age_bins:
        D.tdve["imig_rate"].ts[f"atsi_{age}_M"].assumption = 0
        D.tdve["emig_rate"].ts[f"atsi_{age}_M"].assumption = 0
        D.tdve["imig_rate"].ts[f"atsi_{age}_F"].assumption = 0
        D.tdve["emig_rate"].ts[f"atsi_{age}_F"].assumption = 0

    # Aging
    #atsi_age_in = atsi_age_in[atsi_age_in["Age group"]!="85+"]
    pop_from, pop_to = age_bins[:-1], age_bins[1:]

    for idx, age in enumerate(pop_from):
        D.transfers[0].ts.append((f"atsi_{age}_M", f"atsi_{pop_to[idx]}_M"), at.TimeSeries(atsi_age_in[atsi_age_in["pop_from"]==age]["Year"],
                                                                                    atsi_age_in[atsi_age_in["pop_from"]==age]["Male"],
                                                                                    units="Rate (per year)"))
        D.transfers[0].ts.append((f"atsi_{age}_F", f"atsi_{pop_to[idx]}_F"), at.TimeSeries(atsi_age_in[atsi_age_in["pop_from"]==age]["Year"],
                                                                                    atsi_age_in[atsi_age_in["pop_from"]==age]["Female"],
                                                                                    units="Rate (per year)"))


    D.save(_get_github_folder() + f"databook/hbv_db_v2.0_test.xlsx")










def gen_databook():

    # TODO: Add ability to change reference points for getting/saving files, save names etc.

    # Populations and age bins (as they will be in databook)
    age_bins = ["0-4", "5-14", "15-29", "30-49", "50-64", "65+"]
    sexes = ["M", "F"]
    demog = ["_aus", "_oth", "_fns"]  # Born in AUS, Born OS, First Nations (Aboriginal and/or Torres Strait Islander)
    pops = [f"{age}{sex}{demo}" for demo in demog for sex in sexes for age in age_bins]

    # Import data from Sharepoint Folder
    aus_tabs = ["auspop", "ausbirth", "austrans", "ausdeath", "ausarrive", "ausdepart"]
    oth_tabs = ["othpop", "othtrans", "othdeath", "otharrive", "othdepart"]
    fns_tabs = ["atsipop", "atsibirth", "atsitrans", "atsideath"]
    tabs = aus_tabs + oth_tabs + fns_tabs

    pop_data = {}
    for tab in tabs:
        # Loops through table names and saves each as a pandas dataframe
        pop_data[tab] = read_table(
            _get_sharepoint_folder() + f"Framework Development/Populations/Population Parameters_v1.0.xlsx",
            tab)
        pop_data[tab] = pop_data[tab][pop_data[tab]["year"] >= 1980]

    # Generate Databook
    F = at.ProjectFramework(_get_github_folder() + f"framework/fw_popsizes.xlsx")
    D = at.ProjectData.new(framework=F, pops=1, transfers=0, tvec=np.arange(1980, 2072, 1))
    for idx, val in enumerate(pops):
        if idx == 0:
            D.rename_pop("pop_0", new_code_name=val, new_full_name=val)
        else:
            D.add_pop(val, val)
    D.add_transfer("age", "aging")

    # Fill the databook and save to GitHub repository

    # Australian-Born Input Parameters
    for age in age_bins:
        for sex in sexes:
            # Popsize
            D.tdve["aus_pop"].ts[f"{age}{sex}_aus"] = at.TimeSeries(t=pop_data["auspop"]["year"],
                                                                    vals=pop_data["auspop"][f"{age}{sex}"],
                                                                    units="Number")
            # Mortality
            D.tdve["acm_rate"].ts[f"{age}{sex}_aus"] = at.TimeSeries(t=pop_data["ausdeath"]["year"],
                                                                     vals=pop_data["ausdeath"][f"{age}{sex}"],
                                                                     units="Probability (per year)")
            # Arrivals
            D.tdve["pop_arrive"].ts[f"{age}{sex}_aus"] = at.TimeSeries(t=pop_data["ausarrive"]["year"],
                                                                       vals=pop_data["ausarrive"][f"{age}{sex}"],
                                                                       units="Number (per year)")
            # Departures
            D.tdve["pop_depart"].ts[f"{age}{sex}_aus"] = at.TimeSeries(t=pop_data["ausdepart"]["year"],
                                                                       vals=pop_data["ausdepart"][f"{age}{sex}"],
                                                                       units="Number (per year)")
            # Births
            if age == "0-4":
                D.tdve["b_rate"].ts[f"{age}{sex}_aus"] = at.TimeSeries(t=pop_data["ausbirth"]["year"],
                                                                       vals=pop_data["ausbirth"][f"{age}{sex}"],
                                                                       units="Number (per year)")
            else:
                D.tdve["b_rate"].ts[f"{age}{sex}_aus"].assumption = 0

    # Australian-Born Transition (Aging) Matrix
    pop_from, pop_to = age_bins[:-1], age_bins[1:]

    for i, pf in enumerate(pop_from):
        D.transfers[0].ts.append((f"{pf}M_aus", f"{pop_to[i]}M_aus"), at.TimeSeries(pop_data["austrans"]["year"],
                                                                                    pop_data["austrans"][f"{pf}M"],
                                                                                    units="Rate (per year)"))
        D.transfers[0].ts.append((f"{pf}F_aus", f"{pop_to[i]}F_aus"), at.TimeSeries(pop_data["austrans"]["year"],
                                                                                    pop_data["austrans"][f"{pf}F"],
                                                                                    units="Rate (per year)"))
    # Overseas-Born Input Parameters
    for age in age_bins:
        for sex in sexes:
            # Popsize
            D.tdve["aus_pop"].ts[f"{age}{sex}_oth"] = at.TimeSeries(t=pop_data["othpop"]["year"],
                                                                    vals=pop_data["othpop"][f"{age}{sex}"],
                                                                    units="Number")
            # Mortality
            D.tdve["acm_rate"].ts[f"{age}{sex}_oth"] = at.TimeSeries(t=pop_data["othdeath"]["year"],
                                                                     vals=pop_data["othdeath"][f"{age}{sex}"],
                                                                     units="Probability (per year)")
            # Arrivals
            D.tdve["pop_arrive"].ts[f"{age}{sex}_oth"] = at.TimeSeries(t=pop_data["otharrive"]["year"],
                                                                       vals=pop_data["otharrive"][f"{age}{sex}"],
                                                                       units="Number (per year)")
            # Departures
            D.tdve["pop_depart"].ts[f"{age}{sex}_oth"] = at.TimeSeries(t=pop_data["othdepart"]["year"],
                                                                       vals=pop_data["othdepart"][f"{age}{sex}"],
                                                                       units="Number (per year)")
            # Births
            D.tdve["b_rate"].ts[f"{age}{sex}_oth"].assumption = 0

    # Overseas-Born Transition (Aging) Matrix
    for i, pf in enumerate(pop_from):
        D.transfers[0].ts.append((f"{pf}M_oth", f"{pop_to[i]}M_oth"), at.TimeSeries(pop_data["othtrans"]["year"],
                                                                                    pop_data["othtrans"][f"{pf}M"],
                                                                                    units="Rate (per year)"))
        D.transfers[0].ts.append((f"{pf}F_oth", f"{pop_to[i]}F_oth"), at.TimeSeries(pop_data["othtrans"]["year"],
                                                                                    pop_data["othtrans"][f"{pf}F"],
                                                                                    units="Rate (per year)"))

    # Aboriginal Torres Strait Islanders Input Parameters
    for age in age_bins:
        for sex in sexes:
            # Popsize
            D.tdve["aus_pop"].ts[f"{age}{sex}_fns"] = at.TimeSeries(t=pop_data["atsipop"]["year"],
                                                                    vals=pop_data["atsipop"][f"{age}{sex}"],
                                                                    units="Number")
            # Mortality
            D.tdve["acm_rate"].ts[f"{age}{sex}_fns"] = at.TimeSeries(t=pop_data["atsideath"]["year"],
                                                                     vals=pop_data["atsideath"][f"{age}{sex}"],
                                                                     units="Probability (per year)")
            # Arrivals
            D.tdve["pop_arrive"].ts[f"{age}{sex}_fns"].assumption = 0
            # Departures
            D.tdve["pop_depart"].ts[f"{age}{sex}_fns"].assumption = 0
            # Births
            if age == "0-4":
                D.tdve["b_rate"].ts[f"{age}{sex}_fns"] = at.TimeSeries(t=pop_data["atsibirth"]["year"],
                                                                       vals=pop_data["atsibirth"][f"{age}{sex}"],
                                                                       units="Number (per year)")
            else:
                D.tdve["b_rate"].ts[f"{age}{sex}_fns"].assumption = 0

    # Aboriginal Torres Strait Islanders (Aging) Matrix
    for i, pf in enumerate(pop_from):
        D.transfers[0].ts.append((f"{pf}M_fns", f"{pop_to[i]}M_fns"), at.TimeSeries(pop_data["atsitrans"]["year"],
                                                                                    pop_data["atsitrans"][f"{pf}M"],
                                                                                    units="Rate (per year)"))
        D.transfers[0].ts.append((f"{pf}F_fns", f"{pop_to[i]}F_fns"), at.TimeSeries(pop_data["atsitrans"]["year"],
                                                                                    pop_data["atsitrans"][f"{pf}F"],
                                                                                    units="Rate (per year)"))
    D.save(_get_github_folder() + f"databook/db_demographics.xlsx")


