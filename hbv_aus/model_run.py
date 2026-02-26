import atomica as at
import numpy as np
from hbv_aus.utils import _get_github_folder,_get_sharepoint_folder

def calibrate_model():
    # TODO: Update to include ABC approach once ready
    # TODO: Allow external definition of FW and DB to be used

    #Set up Atomica project from latest update
    F = at.ProjectFramework(_get_github_folder()+f"framework/fw_popsizes.xlsx")
    P = at.Project(framework=F, databook=_get_github_folder() + f"databook/db_demographics.xlsx",
                                do_run = False, sim_start = 1980, sim_end = 2071, sim_dt = 1)
    cal = P.parsets[0].copy()

    # Population Calibrations
    cal = P.calibrate(parset = cal, yaml = _get_github_folder()+f"calibrations/YAML/calibrate_populations.yaml")
    # TODO: Disease/Prevalence Calibrations
    # TODO: Care Calibrations

    # Save Calibration to GitHub
    cal.save_calibration(_get_github_folder()+f"calibrations/Y-factors/calibrate_populations.xlsx")










def extract_data()