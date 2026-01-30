import hbv_vimc as hbv
import atomica as at
import sciris as sc
import pandas as pd


def initial_update(country: str, update_yr=2000):

    # Run model using existing project set up (cal_b)

    s_yr = hbv.s_yr  # for reference to update_year
    P = at.Project(framework=hbv.fw, databook = hbv.root/'databooks'/f"{country}_db.xlsx",
                   sim_start=hbv.s_yr, sim_end=hbv.e_yr, sim_dt=hbv.dt, do_run=False)

    calibration_file = hbv.root / 'ref_cals' / f"{country}_cal b.xlsx"
    cal = P.parsets[0]
    cal.load_calibration(calibration_file)

    res_cal = P.run_sim(cal, result_name="init_cond_run")

    # Extract initial conditions using update_yr
    pops = ["0-0M", "0-0F", "1-4M", "1-4F", "5-9M", "5-9F", "10-19M", "10-19F",
            "20-29M", "20-29F", "30-39M", "30-39F", "40-49M", "40-49F", "50-59M",
            "50-59F", "60-69M", "60-69F", "70-79M", "70-79F", "80-89M", "80-89F",
            "90+M", "90+F"]

    sdis = pd.DataFrame(columns=pops)
    edis = pd.DataFrame(columns=pops)
    ccp = pd.DataFrame(columns=pops)
    dcp = pd.DataFrame(columns=pops)
    hcp = pd.DataFrame(columns=pops)


    for pop in pops:
        sdis[pop] = at.PlotData(res_cal, pops=pop, outputs={"i_sdis":"(ict+ict_dx+ict_tx)/(ict+ict_dx+ict_tx+ie+ie_dx+ie_tx)"},
                                t_bins=1). series[0].vals
        edis[pop] = at.PlotData(res_cal, pops=pop, outputs={"i_edis":"(it+it_dx)/(it+it_dx+icl+icl_dx+icl_tx)"},
                                t_bins=1). series[0].vals
        ccp[pop] = at.PlotData(res_cal, pops=pop, outputs={"i_ccp":"comp/chb_pop"},
                               t_bins=1). series[0].vals
        dcp[pop] = at.PlotData(res_cal, pops=pop, outputs={"i_dcp":"decomp/chb_pop"},
                               t_bins=1). series[0].vals
        hcp[pop] = at.PlotData(res_cal, pops=pop, outputs={"i_hcp":"hepcc/chb_pop"},
                               t_bins=1). series[0].vals

    # Insert into databoook

    D = at.ProjectData.from_spreadsheet(hbv.root/'databooks'/f"{country}_db.xlsx", framework= at.ProjectFramework(hbv.fw))

    for i, pop in enumerate(pops):
        D.tdve["i_sdis"].ts[i].assumption = sdis.at[update_yr-s_yr, pop]
        D.tdve["i_edis"].ts[i].assumption = edis.at[update_yr-s_yr, pop]
        D.tdve["i_ccp"].ts[i].assumption = ccp.at[update_yr-s_yr, pop]
        D.tdve["i_dcp"].ts[i].assumption = dcp.at[update_yr-s_yr, pop]
        D.tdve["i_hcp"].ts[i].assumption = hcp.at[update_yr-s_yr, pop]

    # Save databook for use

    D.save(hbv.root / 'databooks' / f"{country}_db.xlsx")


if __name__ == '__main__':
    kwargs = {"update_yr":2000}  # Debug mode
    sc.parallelize(initial_update, hbv.vimc.keys(), die=False, kwargs=kwargs)
