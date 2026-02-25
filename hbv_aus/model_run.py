import atomica as at
import numpy as np
from hbv_aus.utils import get_github_folder,_get_sharepoint_folder

def calibrate_model():
    # TODO: Update to include ABC approach once ready
    # TODO: Allow external definition of FW and DB to be used
    # TODO: Add functionality (script called utils_model) to assess calibration fits

    #Set up Atomica project from latest update
    F = at.ProjectFramework(get_github_folder()+f"framework/fw_popsizes.xlsx")
    P = at.Project(framework=F, databook=_get_github_folder() + f"databook/db_demographics.xlsx",
                                do_run = False, sim_start = 1980, sim_end = 2071, sim_dt = 1)
    cal = P.parsets[0].copy()
    cal = P.calibrate(max_time=300, parset=cal, adjustables=["acm_rate"],
                      measurables=["aus_pop"],
                      default_min_scale=0.5, default_max_scale=2)


def extract_data()