import atomica as at
from hbv_aus.utils import _get_github_folder,extract_hbv_effects_by_measure
import pandas as pd
import sciris as sc
import numpy as np
from matplotlib import pyplot as plt


def poptest_model_run():
    FW_PATH = _get_github_folder() + f"/framework/fw_popsizes.xlsx"
    DB_PATH = _get_github_folder() + f"/databook/popsize_test_20926.xlsx"
    YAML_PATH =  _get_github_folder()+f"/calibrations/"

    F = at.ProjectFramework(FW_PATH)

    P = at.Project(framework=F, databook=DB_PATH, do_run=False,
                   sim_start=1980, sim_end=2071, sim_dt=1)

    res = P.run_sim(parset="default", result_name="Uncalibrated")

    cal = P.parsets[0].copy()
    cal = P.calibrate(parset=cal, yaml=YAML_PATH + "YAML/calibrate_populations.yaml")

    cal_res = P.run_sim(parset=cal, result_name="Calibrated Test")

    born_os = {"BOS":["lros_0-14_M", "lros_15-64_M", "lros_65+_M", "lros_0-14_F", "lros_15-64_F", "lros_65+_F",
                                                         "hros_0-14_M", "hros_15-64_M", "hros_65+_M","hros_0-14_F", "hros_15-64_F", "hros_65+_F"]}
    atsi = {"ATSI":["atsi_0-14_M", "atsi_15-64_M", "atsi_65+_M", "atsi_0-14_F", "atsi_15-64_F", "atsi_65+_F"]}


    a = at.PlotData([res,cal_res], outputs= "aus_pop", pops=atsi, pop_aggregation="sum",t_bins=1).series[0].vals
    c = at.PlotData([res,cal_res], outputs= "aus_pop", pops=born_os, pop_aggregation="sum",t_bins=1).series[0].vals
    b = at.PlotData([res,cal_res], outputs= "aus_pop", pops="total", pop_aggregation="sum",t_bins=1).series[0].vals

    prop_atsi = a/b
    prop_bos = c/b


    plt.plot(np.linspace(1980,2071, 91), prop_bos)
    plt.plot(np.linspace(1980,2071, 91), prop_atsi)
    plt.ylim(0,0.50)


    at.plot_series(d, data=P.data, axis="results")


def hepaus_model_cal():

    """ Run the calibration using YAML files and save results to a spreadsheet"""

    # Set Paths
    FW_PATH = _get_github_folder()+f"/framework/hbv_fw_v2.1_autosave.xlsx"
    DB_PATH =  _get_github_folder()+f"/databook/hbv_db_hepaus_220926.xlsx"
    YAML_PATH =  _get_github_folder()+f"/calibrations/"

    # Get Items
    F = at.ProjectFramework(FW_PATH)
    P = at.Project(framework=F, databook = DB_PATH, do_run = False,
                            sim_start = 1980, sim_end = 2071, sim_dt=0.25)



    # Run calibrations and save y_factors
    cal = P.parsets[0].copy()

    # Currently following same approach as SA HepB-BD analysis but using deaths instead of HCC incidence.
    cal = P.calibrate(parset = cal, yaml = YAML_PATH+"YAML/hepaus_calibrate_populations.yaml")
    cal = P.calibrate(parset = cal, yaml = YAML_PATH+"YAML/hbv_prevalence_calibrate.yaml")
    cal = P.calibrate(parset = cal, yaml = YAML_PATH+"YAML/hbv_burden_calibrate.yaml")
    cal.save_calibration(YAML_PATH+"Y-factors/hbv_hepaus_calibrations.xlsx")

    res = P.run_sim(parset=cal, result_name="Calibration")
    d = at.PlotData(res, "hep_dth", t_bins=1, pop_aggregation="sum")
    at.plot_series(d, data=P.data)

    cal = P.calibrate(parset = cal, yaml = YAML_PATH+"YAML/hbv_prevalence_calibrate.yaml")
    cal = P.calibrate(parset = cal, yaml = YAML_PATH+"YAML/hbv_burden_calibrate.yaml")

    cal.save_calibration(YAML_PATH+"Y-factors/hbv_hepaus_calibrations.xlsx")

    # cal = P.calibrate(parset = cal, yaml = _get_github_folder()+f"calibrations/YAML/calibrate_care.yaml") # only if needed



def run_hepaus_uncalibrated():

    """ For rapid testing of changes to Framework etc, can run an uncalibrated model and
    returns a dictionary of project & results set for analysis"""
    FW_PATH = _get_github_folder() + f"/framework/hbv_fw_v2.1_autosave.xlsx"
    DB_PATH = _get_github_folder() + f"/databook/claude_hbv_hepaus_db.xlsx"

    F = at.ProjectFramework(FW_PATH)
    P = at.Project(framework=F, databook=DB_PATH, do_run=False,
                   sim_start=1980, sim_end=2071, sim_dt=1)

    res = P.run_sim(parset="default", result_name="Uncalibrated Test")

    return {"Project": P, "Results": res}


def export_hepaus_baselines():

    """ Export data needed for the HepAus spreadsheet used for running scenarios in a .xlsx sheet"""
    FW_PATH = _get_github_folder()+f"/framework/hbv_fw_v2.1_autosave.xlsx"
    DB_PATH =  _get_github_folder()+f"/databook/claude_hbv_hepaus_db.xlsx"
    CAL_PATH =  _get_github_folder()+f"/calibrations/"

    # Run the model and extract: testing, treatment, and linkage to care rates (for now)
    F = at.ProjectFramework(FW_PATH)
    P = at.Project(framework=F, databook=DB_PATH, do_run=False,
                   sim_start=1980, sim_end=2071, sim_dt=1)

    cal = P.parsets[0].load_calibration(CAL_PATH+"Y-factors/hbv_hepaus_calibrations.xlsx")
    res = P.run_sim(parset=cal, result_name="Calibrated")

    # Export data to excel (2026 rate for each model population)
    populations = P.parsets[0].pop_names
    pars = ["diag_rate", "treat_rate", "ltc_rate"]
    export_pars = pd.DataFrame(columns = ["par"]+populations)
    export_pars.par = pars
    export_pars = export_pars.set_index('par', drop=True)

    for pop in populations:
        for par in pars:
            store = pd.DataFrame(columns = ["year", "val"])
            store.year = at.PlotData(res, outputs = par, t_bins =1, pops=pop).series[0].tvec
            store.val = at.PlotData(res, outputs=par, t_bins=1, pops=pop).series[0].vals

            for i in range(len(store)):
                if store.year[i] == 2026.5:
                    export_pars.at[par, pop] = store.val[i]

    # Save to excel
    export_pars.to_excel(_get_github_folder()+f"baseline_hepaus_pars.xlsx")
    return export_pars



def run_hepaus_scenarios(local_ref = "C:/Users/chris.seaman/Burnet Institute/WG-Modelling - Documents/Viral hep modelling/Strategy implementation/"):
    # TODO: Incorporate uncertainty analysis

    """
    Runs the full analysis and saves results as a pickle (.pkl) file for post-processing (plotting, export to excel, etc)
    :return:
    """

    FW_PATH = _get_github_folder() + f"/framework/hbv_fw_v2.1_autosave.xlsx"
    DB_PATH = _get_github_folder() + f"/databook/claude_hbv_hepaus_db.xlsx"
    CAL_PATH = _get_github_folder() + f"/calibrations/"

    # Import scenario scale-up values
    scenario_inputs = extract_hbv_effects_by_measure(local_ref+"Hepatitis_Stategy_Cost_Worksheet v0.2.xlsx")

    col_names = {'HBV-Born Overseas (high risk)-0-14 Female': 'hros_0-14_F', 'HBV-Born Overseas (high risk)-15-64 Female': 'hros_15-64_F',
                 'HBV-Born Overseas (high risk)-65+ Female': 'hros_65+_F',  'HBV-Born Overseas (high risk)-0-14 Male': 'hros_0-14_M',
                 'HBV-Born Overseas (high risk)-15-64 Male': 'hros_15-64_M',  'HBV-Born Overseas (high risk)-65+ Male': 'hros_65+_M',
                 'HBV-Born Overseas (low risk)-0-14 Female': 'lros_0-14_F',  'HBV-Born Overseas (low risk)-15-64 Female': 'lros_15-64_F',
                 'HBV-Born Overseas (low risk)-65+ Female': 'lros_65+_F',  'HBV-Born Overseas (low risk)-0-14 Male': 'lros_0-14_M',
                 'HBV-Born Overseas (low risk)-15-64 Male': 'lros_15-64_M', 'HBV-Born Overseas (low risk)-65+ Male': 'lros_65+_M',
                 'HBV-ATSI-0-14 Female': 'atsi_0-14_F',  'HBV-ATSI-15-64 Female': 'atsi_15-64_F',  'HBV-ATSI-65+ Female': 'atsi_65+_F',
                 'HBV-ATSI-0-14 Male': 'atsi_0-14_M',  'HBV-ATSI-15-64 Male': 'atsi_15-64_M',  'HBV-ATSI-65+ Male': 'atsi_65+_M',
                 'HBV-AusBorn-0-14 Female': 'ausb_0-14_F',  'HBV-AusBorn-15-64 Female': 'ausb_15-64_F', 'HBV-AusBorn-65+ Female': 'ausb_65+_F',
                 'HBV-AusBorn-0-14 Male': 'ausb_0-14_M', 'HBV-AusBorn-15-64 Male': 'ausb_15-64_M', 'HBV-AusBorn-65+ Male': 'ausb_65+_M'}

    scenario_inputs['HBV Prevention'] =  scenario_inputs['HBV Prevention'].rename(col_names, axis=1)
    scenario_inputs['HBV Vaccine uptake'] =  scenario_inputs['HBV Vaccine uptake'].rename(col_names, axis=1)
    scenario_inputs['HBV Testing'] =  scenario_inputs['HBV Testing'].rename(col_names, axis=1)
    scenario_inputs['HBV Treatment initiation'] =  scenario_inputs['HBV Treatment initiation'].rename(col_names, axis=1)
    scenario_inputs['HBV Linkage to care'] =  scenario_inputs['HBV Linkage to care'].rename(col_names, axis=1)
    scenario_inputs['HBV Re-engagement with care'] =  scenario_inputs['HBV Re-engagement with care'].rename(col_names, axis=1)

    # Set up and run scenarios

    ### Baseline (can be turned into a loop for uncertainty) ###
    F = at.ProjectFramework(FW_PATH)
    P = at.Project(framework=F, databook=DB_PATH, do_run=False,
                   sim_start=1980, sim_end=2071, sim_dt=1)
    pop_names = list(P.data.pops.keys())

    bl_parset = P.parsets[0].copy()
    bl_parset.load_calibration(CAL_PATH+"Y-factors/hbv_hepaus_calibrations.xlsx")

    # Import Current (Baseline) Rates of Coverage in 2026
    baseline_pars = export_hepaus_baselines()
    for pop in pop_names:
        bl_parset.pars["diag_rate"].ts[pop].insert(2026.5, baseline_pars.at['diag_rate', pop])
        bl_parset.pars["ltc_rate"].ts[pop].insert(2026.5,baseline_pars.at['ltc_rate', pop])
        bl_parset.pars["treat_rate"].ts[pop].insert(2026.5,baseline_pars.at['treat_rate', pop])

    for pop in pop_names:
        bl_parset.pars["diag_rate"].skip_function[pop] = (2026.5, np.inf)
        bl_parset.pars["ltc_rate"].skip_function[pop] = (2026.5, np.inf)
        bl_parset.pars["treat_rate"].skip_function[pop] = (2026.5, np.inf)

    res_bl = P.run_sim(bl_parset, result_name = "Current Rates Continued")

    ### Intervention Scenario 1 (can add more, turned into a loop for uncertainty) ###
    s1_parset =  P.parsets[0].copy()
    s1_parset.load_calibration(CAL_PATH + "Y-factors/hbv_hepaus_calibrations.xlsx")

    for pop in pop_names:
        s1_parset.pars["diag_rate"].ts[pop].insert([2027.5, 2028.5,2029.5,2030.5], list(scenario_inputs["HBV Testing"][pop].values))
        s1_parset.pars["ltc_rate"].ts[pop].insert([2027.5, 2028.5,2029.5,2030.5], list(scenario_inputs["HBV Linkage to care"][pop].values))
        s1_parset.pars["treat_rate"].ts[pop].insert([2027.5, 2028.5,2029.5,2030.5], list(scenario_inputs["HBV Treatment initiation"][pop].values))

    for pop in pop_names:
        s1_parset.pars["diag_rate"].skip_function[pop] = (2026.5, np.inf)
        s1_parset.pars["ltc_rate"].skip_function[pop] = (2026.5, np.inf)
        s1_parset.pars["treat_rate"].skip_function[pop] = (2026.5, np.inf)

    res_s1 = P.run_sim(s1_parset, result_name = "Intervention S1")

    # Temporary (while no uncertainty) - just return P and res items for data extraction

    return {"Baseline": [P, bl_parset, res_bl], "Scenario 1": [P, s1_parset, res_s1]}















