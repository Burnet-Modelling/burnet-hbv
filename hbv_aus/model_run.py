import atomica as at
from hbv_aus.utils import _get_github_folder,extract_hbv_effects_by_measure
import pandas as pd


def hepaus_model_cal():

    """ Run the calibration using YAML files and save results to a spreadsheet"""

    # Set Paths
    FW_PATH = _get_github_folder()+f"/framework/hbv_fw_v2.1_autosave.xlsx"
    DB_PATH =  _get_github_folder()+f"/databook/claude_hbv_hepaus_db.xlsx"
    YAML_PATH =  _get_github_folder()+f"/calibrations/"

    # Get Items
    F = at.ProjectFramework(FW_PATH)
    P = at.Project(framework=F, databook = DB_PATH, do_run = False,
                            sim_start = 1980, sim_end = 2071, sim_dt=1)

    # Run calibrations and save y_factors
    cal = P.parsets[0].copy()

    cal = P.calibrate(parset = cal, yaml = YAML_PATH+"YAML/calibrate_populations.yaml")
    #cal = P.calibrate(parset = cal, yaml = YAML_PATH+"YAML/calibrate_epidemiology.yaml")
    # cal = P.calibrate(parset = cal, yaml = _get_github_folder()+f"calibrations/YAML/calibrate_care.yaml") # only if needed
    cal.save_calibration(YAML_PATH+"Y-factors/hbv_hepaus_calibrations.xlsx")


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














