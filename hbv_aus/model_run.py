import atomica as at
from hbv_aus.utils import _get_github_folder


def run_model(db_name = "hbv_db_v2.0_test"):
    """
    Runs model without calibration, to be used for testing purposes (e.g., population dynamics)
    :return:
    """
    F = at.ProjectFramework(_get_github_folder()+f"framework/hbv_fw_v2.0.xlsx") # import framework
    P = at.Project(framework=F, databook = _get_github_folder()+f"databook/{db_name}.xlsx", do_run=False,
                   sim_start = 1980, sim_end = 2071, sim_dt = 1)

    res = P.run_sim(parset="default", result_name="Uncalibrated Test")

    d = at.PlotData(res, "temp_alive", t_bins=1, pops="total")
    at.plot_series(d, data=P.data)


    return res



def calibrate_model():
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










