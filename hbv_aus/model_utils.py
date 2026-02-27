import atomica as at
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from hbv_aus.utils import _get_github_folder, _get_sharepoint_folder

def pop_cal_check():

    # Run the model
    F = at.ProjectFramework(_get_github_folder()+f"framework/fw_popsizes.xlsx")
    P = at.Project(framework=F, databook=_get_github_folder() + f"databook/db_demographics.xlsx",
                   do_run=False, sim_start=1980, sim_end=2071, sim_dt=1)
    cal = P.parsets[0].load_calibration(_get_github_folder()+f"calibrations/Y-factors/calibrate_populations.xlsx")
    res = P.run_sim(parset=cal, result_name="Calibration: Populations")

    # Aggregate Indices (total population, % COB, % Aboriginal Torres Strait Islander)
    pop_data = pd.read_excel(_get_sharepoint_folder()+f"Data/Calibration Data/Population Calibrations.xlsx")

    # Total Australian Population Size
    total_pop = pd.DataFrame(columns = ["year", "model"])
    total_pop.year = np.arange(1980, 2071, 1)
    total_pop.model = at.PlotData(res, "aus_pop", pops="total", pop_aggregation="sum", t_bins=1).series[0].vals

    sns.set_theme(font_scale=1.5)
    sns.set_style("darkgrid")


    fig1 = plt.figure(figsize=(12,8))
    ax = fig1.add_subplot(1,1,1)
    ax.scatter(pop_data["year"], pop_data["total_pop"]/1e6, label="Data", color="green", alpha=0.5)
    ax.plot(total_pop["year"], total_pop["model"]/1e6, label="Model", color="black")
    ax.set_ylim(bottom=0)
    ax.legend(loc="best")
    ax.set_ylabel("Total Australian Population \n (millions, Medium Series Projections [2022-2071])")
    fig1.tight_layout()
    fig1.savefig(_get_sharepoint_folder()+f"Outputs/Calibrations/Total_aggregates.png")
    plt.close()


    # % country of birth (inc total pop size)
    age_bins = ["0-4", "5-14", "15-29", "30-49", "50-64", "65+"]
    sex = ["M", "F"]
    os_pops = [f"{a}{s}_oth" for a in age_bins for s in sex]

    cob_prop = pd.DataFrame(columns = ["year", "mod_prop", "mod_pop"])
    cob_prop.year = np.arange(1980, 2071, 1)
    cob_prop.mod_prop = at.PlotData(res, "aus_pop", pops={"os":os_pops}, pop_aggregation="sum", t_bins=1).series[0].vals/at.PlotData(res, "aus_pop", pops="total", pop_aggregation="sum", t_bins=1).series[0].vals
    cob_prop.mod_pop = at.PlotData(res, "aus_pop", pops={"os":os_pops}, pop_aggregation="sum", t_bins=1).series[0].vals

    fig2 = plt.figure(figsize=(12,8))
    # Proportion Born Overseas
    prp = fig2.add_subplot(1,2,1)
    prp.scatter(pop_data["year"], pop_data["born_os"]*1e2, label="Data", color="green", alpha=0.5)
    prp.plot(cob_prop.year, cob_prop.mod_prop*1e2, label="Model", color="black" )
    prp.legend(loc="best")
    prp.set_ylabel("Proportion Born Overseas (%)")
    prp.set_ylim(bottom=0, top=100)
    # Number Born Overseas
    pop = fig2.add_subplot(1,2,2)
    pop.plot(cob_prop.year, cob_prop.mod_pop/1e6 , label="Model", color="black")
    pop.set_ylabel("Number Born Overseas (millions)")
    pop.set_ylim(bottom=0)
    fig2.tight_layout()
    fig2.savefig(_get_sharepoint_folder()+f"Outputs/Calibrations/BOS_aggregates.png")
    plt.close()



    # Total Aboriginal and/or Torres Strait Islander (inc as % of pop)
    fn_pops = [f"{a}{s}_fns" for a in age_bins for s in sex]

    fn_prop = pd.DataFrame(columns = ["year", "mod_prop", "mod_pop"])
    fn_prop.year = np.arange(1980, 2071, 1)
    fn_prop.mod_prop = at.PlotData(res, "aus_pop", pops={"os":fn_pops}, pop_aggregation="sum", t_bins=1).series[0].vals/at.PlotData(res, "aus_pop", pops="total", pop_aggregation="sum", t_bins=1).series[0].vals
    fn_prop.mod_pop = at.PlotData(res, "aus_pop", pops={"os":fn_pops}, pop_aggregation="sum", t_bins=1).series[0].vals

    fig3 = plt.figure(figsize=(12,8))
    # Proportion Born Overseas
    prp = fig3.add_subplot(1,2,1)
    prp.plot(fn_prop.year, fn_prop.mod_prop*1e2, label="Model", color="black" )
    prp.set_ylabel("Proportion Aboriginal and/or Torres Strait Islander (%)")
    prp.set_ylim(bottom=0, top=100)
    # Number Born Overseas
    pop = fig3.add_subplot(1,2,2)
    pop.scatter(pop_data["year"], pop_data["first_nations"]/1e6, label="Data", color="green", alpha=0.5)
    pop.plot(fn_prop.year, fn_prop.mod_pop/1e6 , label="Model", color="black")
    pop.legend(loc="best")
    pop.set_ylabel("Number Aboriginal and/or Torres Strait Islander \n(millions)")
    pop.set_ylim(bottom=0)
    fig3.tight_layout()
    fig3.savefig(_get_sharepoint_folder()+f"Outputs/Calibrations/FN_aggregates.png")
    plt.close()


    # Individual Population Indices (each population saved as own sheet)
    pop_bins = ["_aus", "_oth", "_fns"]
    pop_names = ["Australian Born", "Born Overseas", "Aboriginal and Torres Strait Islander"]
    age_pop = [f"{a}{s}" for a in age_bins for s in sex]

    for  n,p in enumerate(pop_bins):
        fig = plt.figure(figsize=(24,16))
        for idx, ab in enumerate(age_pop):
            ax = fig.add_subplot(6,2,idx+1)
            ax.plot(at.PlotData(res, "aus_pop", pops=f"{ab}{p}", t_bins=1).series[0].tvec,
                        at.PlotData(res, "aus_pop", pops=f"{ab}{p}", t_bins=1).series[0].vals,
                        color="black", label="Model")
            ax.scatter(P.parsets[0].get_par("aus_pop").ts[f"{ab}{p}"].t, P.parsets[0].get_par("aus_pop").ts[f"{ab}{p}"].vals,
                       color="green", alpha=0.5, label="Data")
            ax.set_ylabel("Population")
            ax.set_title(f"{ab}")
            ax.set_ylim(bottom=0)
            if idx == 0:
                ax.legend(loc="best")
        fig.suptitle(pop_names[n])
        fig.tight_layout()
        fig.savefig(_get_sharepoint_folder()+f"Outputs/Calibrations/{pop_names[n]}_individual.png")
        plt.close()


