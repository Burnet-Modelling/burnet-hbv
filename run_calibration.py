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


def run_calibrate_country(country: str, n_iterations=6, max_time_per_step=1000, cal_prev=True):
    np.random.seed(20231002)  # yyyymmdd
    P = hbv.project(country, calibrated=False)

    if cal_prev:

        P.settings.update_time_vector(end=2100)
        cal = P.parsets[0]
        for _ in range(n_iterations):
            cal_a = P.settings.run_calibration(P, parset=cal,
                                               max_time_per_step=max_time_per_step)  # change time per step for longer/shorter runs
        calibration_file = hbv.root / 'calibrations' / f'{country}_cal b.xlsx'  # altered to only calibrate population in missing countries
        cal_a.save_calibration(calibration_file)

        P.settings_b.update_time_vector(end=2100)
        for _ in range(n_iterations):
            cal = P.settings_b.run_calibration(P, parset=cal_a,
                                               max_time_per_step=max_time_per_step)  # change time per step for longer/shorter runs

        calibration_file = hbv.root / 'calibrations' / f'{country}_cal.xlsx'
        cal.save_calibration(calibration_file)  # returns y_factors from YAML calibration in .xlsx spreadsheet

        hbv.calibrate_burden(country)  # scale for burden (ussd the above saved calibration) and overrides with scaled
        calibration_file = hbv.root / 'calibrations' / f'{country}_cal b.xlsx'
        plot_calibration(country, calibration_file)

    if not cal_prev:

        calibration_file = hbv.root / 'calibrations' / f'{country}_cal_prev.xlsx'
        cal_a = P.parsets[0].load_calibration(calibration_file)

        P.settings_b.update_time_vector(end=2100)
        for _ in range(n_iterations):
            cal = P.settings_b.run_calibration(P, parset=cal_a,
                                               max_time_per_step=max_time_per_step)  # change time per step for longer/shorter runs

        calibration_file = hbv.root / 'calibrations' / f'{country}_cal.xlsx'
        cal.save_calibration(calibration_file)  # returns y_factors from YAML calibration in .xlsx spreadsheet

        hbv.calibrate_burden(country)  # scale for burden (ussd the above saved calibration) and overrides with scaled
        calibration_file = hbv.root / 'calibrations' / f'{country}_cal b.xlsx'
        plot_calibration(country, calibration_file)


if __name__ == '__main__':
    kwargs = {'n_iterations': 1, 'max_time_per_step': 400, 'cal_prev': True}  # Debug mode
    sc.parallelize(run_calibrate_country, hbv.vimc.keys(), die=False, kwargs=kwargs)

# hbv.calib_summary(hbv.develop)
