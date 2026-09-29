import atomica as at
import numpy as np

FW_PATH = "framework/hbv_fw_v2.1_autosave.xlsx"
DB_PATH = "databook/hbv_db_hepaus_240926.xlsx"
CAL_PATH = "calibrations/Y-factors/hbv_hepaus_calibrations.xlsx"

F = at.ProjectFramework(FW_PATH)
P = at.Project(framework=F, databook=DB_PATH, do_run=False, sim_start=1980, sim_end=2071, sim_dt=1)
parset = P.parsets[0].copy()
parset.load_calibration(CAL_PATH)
res = P.run_sim(parset=parset, result_name="Calibrated")

t_list = np.array(res.model.t)
post2024_idx = np.where(t_list >= 2024)[0]

for quantity in ("diag_rate", "ltc_rate", "treat_rate"):
    zero_frac_by_pop = {}
    for pop in parset.pop_names:
        vals = np.array(res.model.get_pop(pop).get_par(quantity).vals)[post2024_idx]
        zero_frac_by_pop[pop] = np.mean(vals == 0)
    overall = np.mean(list(zero_frac_by_pop.values()))
    print(f"{quantity}: overall avg fraction of post-2024 years at exactly zero = {overall:.2%}")
    worst = sorted(zero_frac_by_pop.items(), key=lambda x: -x[1])[:5]
    for pop, frac in worst:
        print(f"    {pop:16s} {frac:.1%} of years are exactly zero")
    print()
