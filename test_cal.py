## YAML calibration for HBV model

# %% Import Required Packages
import hbv_vimc as hbv
import atomica as at
import sciris as sc
import numpy as np
import at_tools


def plot_calibration(country: str, calibration_file):
    P = hbv.project(country, calibrated=True)
    P.parsets[0].load_calibration(calibration_file)
    res_cal = P.run_sim(result_name="Auto_Cal")
    hbv.calib_diag_plots(country, P.parsets[0], res_cal)  # returns pdf of calibration plots
    hbv.mape_calib(country, P.parsets[0], res_cal)
    hbv.pop_calibs(country, P.parsets[0], res_cal)


def manual_calibration(country : str):
    P = hbv.project(country, calibrated=False)
    calibration_file = hbv.root / 'calibrations' / f'{country}_cal b.xlsx'
    plot_calibration(country, calibration_file)


if __name__ == '__main__':
    sc.parallelize(manual_calibration, hbv.develop.keys(), die=False)