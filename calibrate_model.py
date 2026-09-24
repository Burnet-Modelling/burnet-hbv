"""
Calibrate the HepAus model against historical population size, HBsAg prevalence, and
HBV-attributable mortality data, and save the resulting Y-factors.

Runs two calibration passes (each an ASD search over per-population Y-factors):
    1. foi_cal vs hepb_prev   - horizontal transmission force-of-infection vs HBsAg prevalence
    2. y_hcc, y_cc, init_cc, init_dc vs hep_dth - disease-progression scaling vs HBV deaths

foi_cal is bounded at (0, 0.05), not the framework's wider (0.01, 10): with wide bounds the
optimizer can fit the prevalence *curve* using implausibly high local transmission, since
prevalence alone can't distinguish "sustained by local transmission" from "sustained by
immigration of already-infected arrivals" (a separate, untouched pathway via
imig_ei/eh/si/sh). The tight bound keeps local transmission near-negligible and lets
immigration + initial conditions + natural disease progression carry prevalence instead -
check the resulting incidence against run_scenarios.py's incidence panel.

Note: hbv_aus.model_run.hepaus_model_cal() calls P.calibrate(parset=cal, yaml=...), but the
installed atomica version has no yaml= support on Project.calibrate() - it silently falls
back to whichever parameters are flagged "Calibrate" in the framework (currently none), so
that path fails immediately with "ValueError: ASD: input vector cannot be zero". This script
calls P.calibrate() directly with explicit adjustables/measurables instead.

Does not touch diag_rate/ltc_rate/treat_rate (the testing/linkage/treatment cascade) or
total population size (acm/emig_rate aren't currently flagged calibratable in the
framework - a separate, unresolved framework/YAML mismatch).

Usage:
    python calibrate_model.py
"""
import os
import time
import sciris as sc
import atomica as at
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from hbv_aus.utils import _get_github_folder

FW_PATH = "framework/hbv_fw_v2.1_autosave.xlsx"
DB_PATH = "databook/hbv_db_hepaus_240926.xlsx"
CAL_PATH = "calibrations/Y-factors/hbv_hepaus_calibrations.xlsx"

SIM_START = 1980
FULL_SIM_END = 2071
CALIBRATION_SIM_END = 2026  # a small buffer past the last data point (2024) - see docstring
CALIBRATION_RELTOL = 1e-3

CALIBRATION_STEPS = [
    ("foi_cal", "hepb_prev", (0, 0.05), 800),
    (["y_hcc", "y_cc", "init_cc", "init_dc"], "hep_dth", None, 800),
]

# Characteristics to sanity-check pre- vs post-calibration fit against data
FIT_CHECKS = [
    ("total_pop", "Total Population"),
    ("hepb_prev", "HBsAg Prevalence"),
    ("hep_dth", "HBV-Attributable Deaths"),
]


def calibrate(P, cal):
    for par_names, measurable, bounds, max_time in CALIBRATION_STEPS:
        par_names = par_names if isinstance(par_names, list) else [par_names]
        bounds = bounds or (0, 5)
        print(f"Calibrating {par_names} against {measurable} (maxtime={max_time}s)...")
        t0 = time.time()
        cal = P.calibrate(
            parset=cal,
            adjustables=[(name, None, *bounds) for name in par_names],
            measurables=[(measurable, None, 1.0, "fractional")],
            max_time=max_time,
            reltol=CALIBRATION_RELTOL,
        )
        print(f"  done in {time.time() - t0:.0f}s")
    return cal


def plot_fit(P, uncal_res, cal_res, out_dir):
    """One figure per fit-check outcome: uncalibrated vs calibrated model line, data as scatter."""
    saved = []
    for code_name, label in FIT_CHECKS:
        try:
            d = at.PlotData([uncal_res, cal_res], outputs=code_name, pops="total", t_bins=1)
        except Exception as e:
            print(f"  Could not plot {code_name}: {e}")
            continue
        figs = at.plot_series(d, data=P.data, axis="results")
        fig = figs[0]
        fig.suptitle(label)
        fig.tight_layout()
        path = os.path.join(out_dir, f"calibration_fit_{code_name}.png")
        fig.savefig(path)
        plt.close(fig)
        saved.append(path)
    return saved


def main():
    out_dir = os.path.join(_get_github_folder(), "outputs")
    os.makedirs(out_dir, exist_ok=True)

    # Run the ASD search itself on a truncated project - every calibration data point is at
    # or before 2024, so this is faster per-iteration with no fit cost (see module docstring).
    # Separate ProjectFramework instances per Project, in case Project mutates its framework.
    P_cal = at.Project(framework=at.ProjectFramework(FW_PATH), databook=DB_PATH, do_run=False, sim_start=SIM_START, sim_end=CALIBRATION_SIM_END, sim_dt=1)
    cal = P_cal.parsets[0].copy()
    if os.path.exists(CAL_PATH):
        cal.load_calibration(CAL_PATH)
    cal = calibrate(P_cal, cal)
    cal.save_calibration(CAL_PATH)
    print(f"Saved calibrated Y-factors to {CAL_PATH}")

    # Y-factors are parset properties, independent of sim_end - reapply to the full-length
    # project for the fit-check plots and saved results.
    P = at.Project(framework=at.ProjectFramework(FW_PATH), databook=DB_PATH, do_run=False, sim_start=SIM_START, sim_end=FULL_SIM_END, sim_dt=1)
    cal_full = P.parsets[0].copy()
    cal_full.load_calibration(CAL_PATH)

    uncal_res = P.run_sim(parset="default", result_name="Uncalibrated")
    cal_res = P.run_sim(parset=cal_full, result_name="Calibrated")
    cal = cal_full

    saved_plots = plot_fit(P, uncal_res, cal_res, out_dir)
    print(f"Saved calibration fit plots: {saved_plots}")

    pkl_path = os.path.join(out_dir, "calibration_result.pkl")
    sc.save(pkl_path, {"Project": P, "Calibrated Parset": cal, "Uncalibrated Result": uncal_res, "Calibrated Result": cal_res})
    print(f"Saved calibration result to {pkl_path}")


if __name__ == "__main__":
    main()
