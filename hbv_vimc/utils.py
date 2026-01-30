import atomica as at
import hbv_vimc as hbv
import pandas as pd
import sciris as sc
import numpy as np

def project(country:str, calibrated = True, eff='high', scr='norm') -> at.Project:

    # Return a pre-configured project instance
    # Argument could be added here to apply changes associated with scenarios (or otherwise
    # scenario projects could be loaded via a function that calls this one and applies the
    # scenario afterwards)

    #TODO: Remember to intergrate changes into non-dev FW file (add argument for dev)
    framework = hbv.fw
    databook = hbv.root/"databooks"/f'{country}_db.xlsx'

    P = at.Project(framework=framework, databook=databook, do_run=False)
    P.settings = sc.dcp(hbv.settings)
    P.settings_b = sc.dcp(hbv.settings_b)

    if calibrated:
        calibration_file = hbv.root/'calibrations'/f'{country}_cal b.xlsx'
        if not calibration_file.exists():
            at.logger.warn(f'Project {country} was requested to be loaded with calibration, but the calibration file is not present. Proceeding without calibration')
        else:
            P.parsets[0].load_calibration(calibration_file)

    return P


def age_bins(bins):
    """
    Function to generate relevant structures for the population age-bins to be used
    in a HBV modelling application. Allows for easy changes of model structure as
    age-dependent variables are recalculated using outputs of this function.

    Parameters
    ----------
    bins :
        "5y" = 5y age bins (0-4, 5-9 etc)
        "GIC" = Global Investment Case (0-4, 5-14, 15-49, 50-69, 70+)
        "GHE" = Global Health Estimates (0-4, 5-14, 15-29, 30-49, 50-59, 60-69, 70+)
        "VIMC_test" = Test bins for VIMC work (0-4, 5-9, 10-19, 20-29 etc)

    Returns
    -------
    pop_bins (dict) : Used for aggregating data
    birth_bins (dict): Used for aggregating birth data
    cons_pars (pd.dataframe): Constant (setting independent) parameters
    var_pars (pd.dataframe): Constant (setting independent) parameter variance
    """

    ##Get Data
    # Note - column changes will not impact function. Row changes may cause issues.
    vals= pd.read_excel(hbv.root / "inputs" / "vals.xlsx", sheet_name="val", index_col=0)
    error= pd.read_excel(hbv.root / "inputs" / "vals.xlsx", sheet_name="uncert", index_col=0)

    # vals = pd.read_excel(at.parent_dir() / "vals.xlsx", sheet_name="val", index_col=0)
    # error = pd.read_excel(at.parent_dir() / "vals.xlsx", sheet_name="uncert", index_col=0)

    # 5y age bins
    if bins == "5y":
        pop_bins = {"0-4": [0, 4], "5-9": [5, 9], "10-14": [10, 14], "15-19": [15, 19], "20-24": [20, 24],
                    "25-29": [25, 29], "30-34": [30, 34], "35-39": [35, 39], "40-44": [40, 44], "45-49": [45, 49],
                    "50-54": [50, 54], "55-59": [55, 59], "60-64": [60, 64], "65-69": [65, 69], "70-74": [70, 74],
                    "75-79": [75, 79], "80-84": [80, 84], "85-89": [85, 89], "90-94": [90, 94], "95+": [95, 120]}
        birth_bins = {"15-19": [15, 19], "20-24": [20, 24],
                      "25-29": [25, 29], "30-34": [30, 34], "35-39": [35, 39], "40-44": [40, 44], "45-49": [45, 49]}
        cons_pars = vals.T
        var_pars = error.T

    if bins == "GIC":
        pop_bins = {"0-4": [0, 4], "5-14": [5, 14], "15-49": [15, 49], "50-69": [50, 69],
                    "70+": [70, 120]}
        birth_bins = {"15-49": [15, 49]}

        # Constants for GIC age bins (format kept constant too)
        cons_pars = pd.DataFrame()
        cons_pars["0-4M"] = vals.loc[["0-4M"], :].mean(axis=0)
        cons_pars["0-4F"] = vals.loc[["0-4F"], :].mean(axis=0)
        cons_pars["5-14M"] = vals.loc[["5-9M", "10-14M"], :].mean(axis=0)
        cons_pars["5-14F"] = vals.loc[["5-9F", "10-14F"], :].mean(axis=0)
        cons_pars["15-49M"] = vals.loc[["15-19M", "20-24M", "25-29M", "30-34M", "35-39M", "40-44M", "45-49M"], :].mean(
            axis=0)
        cons_pars["15-49F"] = vals.loc[["15-19F", "20-24F", "25-29F", "30-34F", "35-39F", "40-44F", "45-49F"], :].mean(
            axis=0)
        cons_pars["50-69M"] = vals.loc[["50-54M", "55-59M", "60-64M", "65-69M"], :].mean(axis=0)
        cons_pars["50-69F"] = vals.loc[["50-54F", "55-59F", "60-64F", "65-69F"], :].mean(axis=0)
        cons_pars["70+M"] = vals.loc[["70-74M", "75-79M", "80-84M", "85-89M", "90-94M", "95+M"], :].mean(axis=0)
        cons_pars["70+F"] = vals.loc[["70-74F", "75-79F", "80-84F", "85-89F", "90-94F", "95+F"], :].mean(axis=0)

        # Errors for GIC age bins (format kept constant too)
        var_pars = pd.DataFrame()
        var_pars["0-4M"] = error.loc[["0-4M"], :].mean(axis=0)
        var_pars["0-4F"] = error.loc[["0-4F"], :].mean(axis=0)
        var_pars["5-14M"] = error.loc[["5-9M", "10-14M"], :].mean(axis=0)
        var_pars["5-14F"] = error.loc[["5-9F", "10-14F"], :].mean(axis=0)
        var_pars["15-49M"] = error.loc[["15-19M", "20-24M", "25-29M", "30-34M", "35-39M", "40-44M", "45-49M"], :].mean(
            axis=0)
        var_pars["15-49F"] = error.loc[["15-19F", "20-24F", "25-29F", "30-34F", "35-39F", "40-44F", "45-49F"], :].mean(
            axis=0)
        var_pars["50-69M"] = error.loc[["50-54M", "55-59M", "60-64M", "65-69M"], :].mean(axis=0)
        var_pars["50-69F"] = error.loc[["50-54F", "55-59F", "60-64F", "65-69F"], :].mean(axis=0)
        var_pars["70+M"] = error.loc[["70-74M", "75-79M", "80-84M", "85-89M", "90-94M", "95+M"], :].mean(axis=0)
        var_pars["70+F"] = error.loc[["70-74F", "75-79F", "80-84F", "85-89F", "90-94F", "95+F"], :].mean(axis=0)

    if bins == "GHE":
        pop_bins = {"0-4": [0, 4], "5-14": [5, 14], "15-29": [15, 29], "30-49": [30, 49],
                    "50-59": [50, 59], "60-69": [60, 69], "70+": [70, 100]}
        birth_bins = {"15-29": [15, 29], "30-49": [30, 49]}

        # Constants for GHE age bins (format kept constant too)
        cons_pars = pd.DataFrame()
        cons_pars["0-4M"] = vals.loc[["0-4M"], :].mean(axis=0)
        cons_pars["0-4F"] = vals.loc[["0-4F"], :].mean(axis=0)
        cons_pars["5-14M"] = vals.loc[["5-9M", "10-14M"], :].mean(axis=0)
        cons_pars["5-14F"] = vals.loc[["5-9F", "10-14F"], :].mean(axis=0)
        cons_pars["15-29M"] = vals.loc[["15-19M", "20-24M", "25-29M"], :].mean(axis=0)
        cons_pars["15-29F"] = vals.loc[["15-19F", "20-24F", "25-29F"], :].mean(axis=0)
        cons_pars["30-49M"] = vals.loc[["30-34M", "35-39M", "40-44M", "45-49M"], :].mean(axis=0)
        cons_pars["30-49F"] = vals.loc[["30-34F", "35-39F", "40-44F", "45-49F"], :].mean(axis=0)
        cons_pars["50-69M"] = vals.loc[["50-54M", "55-59M", "60-64M", "65-69M"], :].mean(axis=0)
        cons_pars["50-69F"] = vals.loc[["50-54F", "55-59F", "60-64F", "65-69F"], :].mean(axis=0)
        cons_pars["70+M"] = vals.loc[["70-74M", "75-79M", "80-84M", "85-89M", "90-94M", "95+M"], :].mean(axis=0)
        cons_pars["70+F"] = vals.loc[["70-74F", "75-79F", "80-84F", "85-89F", "90-94F", "95+F"], :].mean(axis=0)

        # Errors for GHE age bins (format kept constant too)
        var_pars = pd.DataFrame()
        var_pars["0-4M"] = error.loc[["0-4M"], :].mean(axis=0)
        var_pars["0-4F"] = error.loc[["0-4F"], :].mean(axis=0)
        var_pars["5-14M"] = error.loc[["5-9M", "10-14M"], :].mean(axis=0)
        var_pars["5-14F"] = error.loc[["5-9F", "10-14F"], :].mean(axis=0)
        var_pars["15-29M"] = error.loc[["15-19M", "20-24M", "25-29M"], :].mean(axis=0)
        var_pars["15-29F"] = error.loc[["15-19F", "20-24F", "25-29F"], :].mean(axis=0)
        var_pars["30-49M"] = error.loc[["30-34M", "35-39M", "40-44M", "45-49M"], :].mean(axis=0)
        var_pars["30-49F"] = error.loc[["30-34F", "35-39F", "40-44F", "45-49F"], :].mean(axis=0)
        var_pars["50-69M"] = error.loc[["50-54M", "55-59M", "60-64M", "65-69M"], :].mean(axis=0)
        var_pars["50-69F"] = error.loc[["50-54F", "55-59F", "60-64F", "65-69F"], :].mean(axis=0)
        var_pars["70+M"] = error.loc[["70-74M", "75-79M", "80-84M", "85-89M", "90-94M", "95+M"], :].mean(axis=0)
        var_pars["70+F"] = error.loc[["70-74F", "75-79F", "80-84F", "85-89F", "90-94F", "95+F"], :].mean(axis=0)

    if bins == "vimc_test":
        pop_bins = {"0-0": [0, 0], "1-4": [1, 4], "5-9": [5, 9], "10-19": [10, 19], "20-29": [20, 29],
                    "30-39": [30, 39], "40-49": [40, 49], "50-59": [50, 59],
                    "60-69": [60, 69], "70-79": [70, 79], "80-89": [80, 89],
                    "90+": [90, 120]}
        birth_bins = {"10-19": [10, 19], "20-29": [20, 29], "30-39": [30, 39],
                      "40-49": [40, 49], "50-59": [50, 59]}

        # Constants for VIMC_test age bins (format kept constant too)
        cons_pars = pd.DataFrame()
        cons_pars["0-0M"] = vals.loc[["0-0M"], :].mean(axis=0)
        cons_pars["0-0F"] = vals.loc[["0-0F"], :].mean(axis=0)
        cons_pars["1-4M"] = vals.loc[["1-4M"], :].mean(axis=0)
        cons_pars["1-4F"] = vals.loc[["1-4F"], :].mean(axis=0)
        cons_pars["5-9M"] = vals.loc[["5-9M"], :].mean(axis=0)
        cons_pars["5-9F"] = vals.loc[["5-9F"], :].mean(axis=0)
        cons_pars["10-19M"] = vals.loc[["10-14M", "15-19M"], :].mean(axis=0)
        cons_pars["10-19F"] = vals.loc[["10-14F", "15-19F"], :].mean(axis=0)
        cons_pars["20-29M"] = vals.loc[["20-24M", "25-29M"], :].mean(axis=0)
        cons_pars["20-29F"] = vals.loc[["20-24F", "25-29F"], :].mean(axis=0)
        cons_pars["30-39M"] = vals.loc[["30-34M", "35-39M"], :].mean(axis=0)
        cons_pars["30-39F"] = vals.loc[["30-34F", "35-39F"], :].mean(axis=0)
        cons_pars["40-49M"] = vals.loc[["40-44M", "45-49M"], :].mean(axis=0)
        cons_pars["40-49F"] = vals.loc[["40-44F", "45-49F"], :].mean(axis=0)
        cons_pars["50-59M"] = vals.loc[["50-54M", "55-59M"], :].mean(axis=0)
        cons_pars["50-59F"] = vals.loc[["50-54F", "55-59F"], :].mean(axis=0)
        cons_pars["60-69M"] = vals.loc[["60-64M", "65-69M"], :].mean(axis=0)
        cons_pars["60-69F"] = vals.loc[["60-64F", "65-69F"], :].mean(axis=0)
        cons_pars["70-79M"] = vals.loc[["70-74M", "75-79M"], :].mean(axis=0)
        cons_pars["70-79F"] = vals.loc[["70-74F", "75-79F"], :].mean(axis=0)
        cons_pars["80-89M"] = vals.loc[["80-84M", "85-89M"], :].mean(axis=0)
        cons_pars["80-89F"] = vals.loc[["80-84F", "85-89F"], :].mean(axis=0)
        cons_pars["90+M"] = vals.loc[["90-94M", "95+M"], :].mean(axis=0)
        cons_pars["90+F"] = vals.loc[["90-94F", "95+F"], :].mean(axis=0)

        # Errors for VIMC_test age bins (format kept constant too)
        var_pars = pd.DataFrame()
        var_pars["0-0M"] = error.loc[["0-0M"], :].mean(axis=0)
        var_pars["0-0F"] = error.loc[["0-0F"], :].mean(axis=0)
        var_pars["1-4M"] = error.loc[["1-4M"], :].mean(axis=0)
        var_pars["1-4F"] = error.loc[["1-4F"], :].mean(axis=0)
        var_pars["5-9M"] = error.loc[["5-9M"], :].mean(axis=0)
        var_pars["5-9F"] = error.loc[["5-9F"], :].mean(axis=0)
        var_pars["10-19M"] = error.loc[["10-14M", "15-19M"], :].mean(axis=0)
        var_pars["10-19F"] = error.loc[["10-14F", "15-19F"], :].mean(axis=0)
        var_pars["20-29M"] = error.loc[["20-24M", "25-29M"], :].mean(axis=0)
        var_pars["20-29F"] = error.loc[["20-24F", "25-29F"], :].mean(axis=0)
        var_pars["30-39M"] = error.loc[["30-34M", "35-39M"], :].mean(axis=0)
        var_pars["30-39F"] = error.loc[["30-34F", "35-39F"], :].mean(axis=0)
        var_pars["40-49M"] = error.loc[["40-44M", "45-49M"], :].mean(axis=0)
        var_pars["40-49F"] = error.loc[["40-44F", "45-49F"], :].mean(axis=0)
        var_pars["50-59M"] = error.loc[["50-54M", "55-59M"], :].mean(axis=0)
        var_pars["50-59F"] = error.loc[["50-54F", "55-59F"], :].mean(axis=0)
        var_pars["60-69M"] = error.loc[["60-64M", "65-69M"], :].mean(axis=0)
        var_pars["60-69F"] = error.loc[["60-64F", "65-69F"], :].mean(axis=0)
        var_pars["70-79M"] = error.loc[["70-74M", "75-79M"], :].mean(axis=0)
        var_pars["70-79F"] = error.loc[["70-74F", "75-79F"], :].mean(axis=0)
        var_pars["80-89M"] = error.loc[["80-84M", "85-89M"], :].mean(axis=0)
        var_pars["80-89F"] = error.loc[["80-84F", "85-89F"], :].mean(axis=0)
        var_pars["90+M"] = error.loc[["90-94M", "95+M"], :].mean(axis=0)
        var_pars["90+F"] = error.loc[["90-94F", "95+F"], :].mean(axis=0)

    return pop_bins, birth_bins, cons_pars, var_pars


def databook_gen(country: str):
    """
    Atomica HBV databook generation script
    Chris Seaman, 11 September 2023

    Produces databooks consistent with framework file hbv_v15.xlsx
    """

    # TODO: Remove the need to double code for this - can lead to errors when developing FW
    F = at.ProjectFramework(hbv.fw)  # import framework file
    pop_bins, birth_bins, cons_pars, var_pars = hbv.age_bins("vimc_test")
    pops = list(cons_pars.columns)
    s_yr = hbv.s_yr
    e_yr = hbv.e_yr

    ## Import data for filling
    raw = pd.read_csv(hbv.root / "inputs" / country / f"{country}_dbins.csv")

    # Initiate Databook
    D = at.ProjectData.new(framework=F, pops=1, transfers=0, tvec=np.arange(s_yr, e_yr, 1))

    ###Set up populations
    for idx, val in enumerate(pop_bins):
        if idx == 0:
            D.rename_pop('pop_0', new_code_name=val + "M", new_full_name=val + "M")
            D.add_pop(val + "F", val + "F")
        else:
            D.add_pop(val + "M", val + "M")
            D.add_pop(val + "F", val + "F")

    ## Add Population Transfers
    t_rates = raw[raw["par_name"] == "pop_t_rate"]
    p_temp = list(pop_bins.keys())

    pop_trans = {}

    for idx in range(len(p_temp) - 1):
        pop_trans[p_temp[idx] + "M_" + p_temp[idx + 1] + "M"] = [p_temp[idx] + "M", p_temp[idx + 1] + "M", idx,
                                                                 list(t_rates["year"][s_yr - 1950:e_yr - 1950]), list(
                t_rates[p_temp[idx] + "M"][s_yr - 1950:e_yr - 1950])]
        pop_trans[p_temp[idx] + "F_" + p_temp[idx + 1] + "F"] = [p_temp[idx] + "F", p_temp[idx + 1] + "F", idx,
                                                                 list(t_rates["year"][s_yr - 1950:e_yr - 1950]), list(
                t_rates[p_temp[idx] + "F"][s_yr - 1950:e_yr - 1950])]

    D.add_transfer("p_trans", "population movement")
    keys = list(pop_trans.keys())

    for idx, val in enumerate(keys):
        D.transfers[0].ts.append((pop_trans[val][0], pop_trans[val][1]), at.TimeSeries())
        D.transfers[0].ts[idx] = at.TimeSeries(pop_trans[val][3], pop_trans[val][4], units="Rate (per year)")

    ## Add MTCT interaction matrix
    mtct_rates = raw[raw["par_name"] == "mtct_prop"]
    mtct_rates = mtct_rates.dropna(axis=1)
    b_pops = list(mtct_rates.iloc[:, 2:].columns.values)

    mtct_mx = {}

    for idx in range(len(b_pops)):
        mtct_mx[b_pops[idx] + "_" + p_temp[0] + "_M"] = [b_pops[idx], p_temp[0] + "M", idx + 1,
                                                         list(mtct_rates["year"][s_yr - 1950:e_yr - 1950]),
                                                         list(mtct_rates[b_pops[idx]][s_yr - 1950:e_yr - 1950])]
        mtct_mx[b_pops[idx] + "_" + p_temp[0] + "_F"] = [b_pops[idx], p_temp[0] + "F", idx + 1,
                                                         list(mtct_rates["year"][s_yr - 1950:e_yr - 1950]),
                                                         list(mtct_rates[b_pops[idx]][s_yr - 1950:e_yr - 1950])]

    keys = list(mtct_mx.keys())

    for idx, val in enumerate(keys):
        D.interpops[0].ts.append((mtct_mx[val][0], mtct_mx[val][1]), at.TimeSeries())

    for idx, val in enumerate(keys):
        D.interpops[0].ts[idx + 1] = at.TimeSeries(mtct_mx[val][3], mtct_mx[val][4], units="N.A.")

    ## Horizontal Transmission matrices
    ec_pops = []
    at_pops_a = []
    at_pops_b = []

    for key, val in pop_bins.items():
        if val[1] <= 5:
            ec_pops.append(key)
        else:
            at_pops_a.append(key)
            at_pops_b.append(key)

    echt_mx = {}
    for i in range(len(ec_pops)):
        for j in range(len(ec_pops)):
            echt_mx[ec_pops[j] + "M_" + ec_pops[i] + "M"] = [ec_pops[i] + "M", ec_pops[j] + "M", j + 1]
            echt_mx[ec_pops[j]+ "F_" + ec_pops[i] + "M"] = [ec_pops[i] + "F", ec_pops[j] + "M", len(p_temp) + j + 1]
            echt_mx[ec_pops[j]+ "M_" + ec_pops[i] + "F"] = [ec_pops[i] + "M", ec_pops[j] + "F", 2 * len(p_temp) + j + 1]
            echt_mx[ec_pops[j] + "F_" + ec_pops[i] + "F"] = [ec_pops[i] + "F", ec_pops[j] + "F", 3 * len(p_temp) + j + 1]

    popt_mx = {}

    for i, x in enumerate(at_pops_a):
        for j, y in enumerate(at_pops_b):
            popt_mx[x + "M_" + y + "M"] = [x + "M", y + "M", j + 1]
            popt_mx[x + "F_" + y + "M"] = [x + "F", y + "M", len(at_pops_a) * len(at_pops_b) + j + 1]
            popt_mx[x + "M_" + y + "F"] = [x + "M", y + "F", 2 * len(at_pops_a) * len(at_pops_b) + j + 1]
            popt_mx[x + "F_" + y + "F"] = [x + "F", y + "F", 3 * len(at_pops_a) * len(at_pops_b) + j + 1]

    keys = list(popt_mx.keys())

    for idx, val in enumerate(keys):
        D.interpops[2].ts.append((popt_mx[val][1], popt_mx[val][0]), at.TimeSeries())
        D.interpops[2].ts[idx + 1].assumption = 1
        D.interpops[2].ts[idx + 1].units = "N.A."

    keys = list(echt_mx.keys())

    for idx, val in enumerate(keys):
        D.interpops[1].ts.append((echt_mx[val][1], echt_mx[val][0]), at.TimeSeries())
        D.interpops[1].ts[idx].assumption = 1
        D.interpops[1].ts[idx].units = "N.A."

    ## Add constants
    cons = cons_pars.index.values.tolist()
    pops = cons_pars.columns.tolist()

    for i in cons:  # loop through vars
        for j in pops:  # loop through pops
            D.tdve[i].ts[j].assumption = cons_pars.at[i, j]
            D.tdve[i].ts[j].sigma = var_pars.at[i, j]

    ## Add initialization
    init_pop = raw[((raw["par_name"] == "alive") & (raw["year"] == 1990))]
    init_pop = init_pop.reset_index()
    init_chb = raw[((raw["par_name"] == "prev") & (raw["year"] == 1990))]
    init_chb = init_chb.reset_index()
    init_hbe = raw[((raw["par_name"] == "eag_ott") & (raw["year"] == 1990))]
    init_hbe = init_hbe.reset_index()
    init_cc = raw[((raw["par_name"] == "i_cc") & (raw["year"] == 1990))]
    init_cc = init_cc.reset_index()
    init_dc = raw[((raw["par_name"] == "i_dc") & (raw["year"] == 1990))]
    init_dc = init_dc.reset_index()
    init_hcc = raw[((raw["par_name"] == "i_hcc") & (raw["year"] == 1990))]
    init_hcc = init_hcc.reset_index()

    for i, pop in enumerate(pops):
        D.tdve["j_init"].ts[i].assumption = init_pop.at[0, pop]
        D.tdve["i_uchb"].ts[i].assumption = init_chb.at[0, pop]
        D.tdve["i_epos"].ts[i].assumption = init_hbe.at[0, pop]
        D.tdve["i_ccp"].ts[i].assumption = init_cc.at[0, pop]
        D.tdve["i_dcp"].ts[i].assumption = init_dc.at[0, pop]
        D.tdve["i_hcp"].ts[i].assumption = init_hcc.at[0, pop]

    # Add demography
    alive = raw[((raw["par_name"] == "alive") & (raw["year"] >= s_yr))]
    births = raw[((raw["par_name"] == "b_rate") & (raw["year"] >= s_yr))]
    births = births.fillna(0)
    preg = raw[((raw["par_name"] == "preg") & (raw["year"] >= s_yr))]
    preg = preg.fillna(0)
    acm = raw[((raw["par_name"] == "acm") & (raw["year"] >= s_yr))]
    lexp = raw[((raw["par_name"] == "life_exp") & (raw["year"] >= s_yr))]
    mig_rate = raw[((raw["par_name"] == "mig_rate") & (raw["year"] >= s_yr))]


    for i, pop in enumerate(pops):
        D.tdve["alive"].ts[i] = at.TimeSeries(t=alive["year"], vals=alive[pop], units="Number")  # total population
        D.tdve["acm"].ts[i] = at.TimeSeries(t=acm["year"], vals=acm[pop], units="Probability (per year)")  # acm
        D.tdve["b_rate"].ts[i] = at.TimeSeries(t=births["year"], vals=births[pop], units="Number")  # births
        D.tdve["pregs"].ts[i] = at.TimeSeries(t=preg["year"], vals=preg[pop])
        D.tdve["life_exp"].ts[i] = at.TimeSeries(t=lexp["year"], vals=lexp[pop], units="Number")  # life expectancy
        D.tdve["mig_rate"].ts[i] = at.TimeSeries(t=mig_rate["year"], vals=mig_rate[pop])

    # Add calibration data
    prev = raw[((raw["par_name"] == "prev") & (raw["year"] >= s_yr))]
    chb_pop = raw[((raw["par_name"] == "chb_pop") & (raw["year"] >= s_yr))]
    eag = raw[((raw["par_name"] == "eag_ott") & (raw["year"] >= s_yr))]
    hcc_inc = raw[((raw["par_name"] == "flw_hcc") & (raw["year"] >= s_yr))]
    acu_d = raw[((raw["par_name"] == "cl_acu") & (raw["year"] >= s_yr))]
    cir_d = raw[((raw["par_name"] == "cl_cir") & (raw["year"] >= s_yr))]
    cir_d = cir_d.fillna(0)
    hcc_d = raw[((raw["par_name"] == "cl_hcc") & (raw["year"] >= s_yr))]
    hcc_d = hcc_d.fillna(0)

    # Expand out dc and hcc

    # # Decompensated Cirrhosis
    #
    # cir_yr = max(cir_d["year"])
    # cir_rep = cir_d[cir_d["year"] == cir_yr]
    # cir_rep = pd.DataFrame(np.repeat(cir_rep.values, 10, axis=0))
    # cir_rep.columns = cir_d.columns
    #
    # for idx, val in enumerate(cir_rep["year"]):
    #     cir_rep["year"][idx] = cir_rep["year"][idx] + (idx + 1)
    #
    # cir_d = pd.concat([cir_d, cir_rep])
    # cir_d = cir_d.reset_index()
    # cir_d = cir_d.iloc[:, 1:]
    #
    # # Hepatocellular Carcinoma
    #
    # hcc_yr = max(hcc_d["year"])
    # hcc_rep = hcc_d[hcc_d["year"] == hcc_yr]
    # hcc_rep = pd.DataFrame(np.repeat(hcc_rep.values, 10, axis=0))
    # hcc_rep.columns = hcc_d.columns
    #
    # for idx, val in enumerate(hcc_rep["year"]):
    #     hcc_rep["year"][idx] = hcc_rep["year"][idx] + (idx + 1)
    #
    # hcc_d = pd.concat([hcc_d, hcc_rep])
    # hcc_d = hcc_d.reset_index()
    # hcc_d = hcc_d.iloc[:, 1:]

    # Impute 1995 and 2015 coverage
    alive = raw[((raw["par_name"] == "alive") & (raw["year"] >= s_yr))]
    mort_yrs = list(cir_d["year"])

    # Years for denomionator in getting rate
    denom = alive[alive["year"].isin(mort_yrs)]
    denom = denom.set_index("year")
    denom = denom.iloc[:,1:]

    # Populations to impute deaths
    imp_years = [1995, 2024]
    mult_pop = alive[alive["year"].isin(imp_years)]
    mult_pop = mult_pop.set_index("year")
    mult_pop = mult_pop.iloc[:,1:]


    cir_imp = cir_d.set_index("year")
    cir_imp = cir_imp.iloc[:,1:]

    hcc_imp = hcc_d.set_index("year")
    hcc_imp = hcc_imp.iloc[:,1:]

    for pop in pops:
        for yrs in mort_yrs:
            cir_imp.at[yrs, pop] = np.nan_to_num(cir_imp.at[yrs, pop] / denom.at[yrs, pop])
            hcc_imp.at[yrs, pop] = np.nan_to_num(hcc_imp.at[yrs, pop] / denom.at[yrs, pop])

    cir_d = cir_d.set_index("year")
    cir_d = cir_d.iloc[:, 1:]

    hcc_d = hcc_d.set_index("year")
    hcc_d = hcc_d.iloc[:, 1:]

    for pop in pops:
        cir_d.at[1995, pop] = cir_imp.at[2000, pop] * mult_pop.at[1995, pop]
        hcc_d.at[1995, pop] = hcc_imp.at[2000, pop] * mult_pop.at[1995, pop]

        cir_d.at[2024, pop] = cir_imp.at[2019, pop] * mult_pop.at[2024, pop]
        hcc_d.at[2024, pop] = hcc_imp.at[2019, pop] * mult_pop.at[2024, pop]

    cir_d = cir_d.reset_index()
    hcc_d = hcc_d.reset_index()

    # Fill Databook

    for i, pop in enumerate(pops):
        D.tdve["prev"].ts[i] = at.TimeSeries(t=prev["year"], vals=prev[pop], units="Fraction")  # HBsAg prevalence
        D.tdve["chb_pop"].ts[i] = at.TimeSeries(t=chb_pop["year"], vals=chb_pop[pop],
                                                units="Number")  # CHB positive population
        D.tdve["eag_ott"].ts[i] = at.TimeSeries(t=eag["year"], vals=eag[pop],
                                                units="Fraction")  # HBeAg positive (among HBsAg)
        D.tdve["flw_hcc"].ts[i] = at.TimeSeries(t=hcc_inc["year"], vals=hcc_inc[pop], units="Number")  # HCC incidence
        D.tdve["cl_acu"].ts[i] = at.TimeSeries(t=acu_d["year"], vals=acu_d[pop], units="Number")  # Acute deaths
        D.tdve["cl_cir"].ts[i] = at.TimeSeries(t=cir_d["year"], vals=cir_d[pop], units="Number")  # Cirrhosis deaths
        D.tdve["cl_hcc"].ts[i] = at.TimeSeries(t=hcc_d["year"], vals=hcc_d[pop], units="Number")  # HCC deaths

    ## Add care coverage data

    diag = raw[((raw["par_name"] == "dx_cov") & (raw["year"] >= s_yr))]
    treat = raw[((raw["par_name"] == "tx_cov") & (raw["year"] >= s_yr))]
    add_row = ["hepb_bd", 1990, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    raw.loc[len(raw)] = add_row
    add_row = ["hepb_bd", 1999, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    raw.loc[len(raw)] = add_row
    bd = raw[((raw["par_name"] == "hepb_bd") & (raw["year"] >= s_yr))]
    bd = bd.fillna(0)
    hb3 = raw[((raw["par_name"] == "hepb") & (raw["year"] >= s_yr))]
    hb3 = hb3.fillna(0)

    for i, pop in enumerate(pops):
        D.tdve["dx_cov_inp"].ts[i] = at.TimeSeries(t=diag["year"], vals=diag[pop])  # HBsAg prevalence
        D.tdve["tx_cov_inp_elig"].ts[i] = at.TimeSeries(t=treat["year"], vals=treat[pop])  # CHB positive population
        D.tdve["bd"].ts[i] = at.TimeSeries(t=bd["year"], vals=bd[pop])  # HBeAg positive (among HBsAg)
        D.tdve["hb3"].ts[i] = at.TimeSeries(t=hb3["year"], vals=hb3[pop])  # HCC incidence

    D.save(hbv.root / "databooks" / f"{country}_db.xlsx")


def calibrate_burden(country: str):
    # Run model using population and prevalence/transmission calibration
    P = at.Project(framework=hbv.fw, databook=hbv.root / 'databooks' / f"{country}_db.xlsx",
                   sim_start=hbv.s_yr, sim_end=hbv.e_yr, sim_dt=hbv.dt, do_run=False)
    calibration_file = hbv.root / 'calibrations' / f'{country}_cal.xlsx'
    cal_sheet = pd.read_excel(calibration_file)
    cal_sheet = cal_sheet.set_index("par")  # allows use of pd.at functionality
    cal = P.parsets[0]
    cal.load_calibration(calibration_file)
    res = P.run_sim(cal, result_name="Pre-Cal")

    # Set-up to extract model data
    model = {}
    data = {}
    ratio = {}
    outputs = ["flw_hcc", "cl_cir", "cl_acu"]
    pops_in = ["year", "0-0M", "0-0F", "1-4M", "1-4F", "5-9M", "5-9F", "10-19M", "10-19F", "20-29M", "20-29F", "30-39M",
               "30-39F", "40-49M", "40-49F", "50-59M", "50-59F", "60-69M", "60-69F", "70-79M", "70-79F", "80-89M",
               "80-89F", "90+M", "90+F"]

    for out in outputs:
        model[out] = pd.DataFrame(columns=pops_in)
        data[out] = pd.DataFrame(columns=pops_in)
        ratio[out] = pd.DataFrame(columns=pops_in)

    # Extract the model data (model)

    for out in outputs:
        for pop in pops_in:
            if pop == "year":
                model[out][pop] = np.floor(at.PlotData(res, out, "total", t_bins=1).series[0].tvec)
            else:
                model[out][pop] = at.PlotData(res, out, pop, t_bins=1).series[0].vals

    # Extract the input data

    for out in outputs:
        for pop in pops_in:
            if pop == "year":
                data[out][pop] = cal.get_par(out).ts["0-0M"].t
            else:
                data[out][pop] = cal.get_par(out).ts[pop].vals

    # Keep only matching years

    for out in outputs:
        for pop in pops_in:
            years = list(data[out]["year"])
        if pop != "year":
            model[out] = model[out][model[out]["year"].isin(years)]

    for out in outputs:
        model[out] = model[out].reset_index()

    # Calculate ratios and get average

    for out in outputs:
        for pop in pops_in:
            if pop != "year":
                ratio[out][pop] = np.nan_to_num(data[out][pop] / model[out][pop], nan=0)
        ratio[out] = ratio[out].drop(["year"], axis=1)
        # ratio[out] = ratio[out].mean()

    # Set y_factor for relevant parameters as ratio
    for pop in pops_in:
        if pop != "year":
            # Acute Mortality
            cal_sheet.at["m_acu", pop] = ratio["cl_acu"][pop].loc[2:].mean()
            # Cirrhosis
            # cal_sheet.at["i_dcp", pop] = ratio["cl_cir"][pop].loc[0]
            # cal_sheet.at["i_ccp", pop] = ratio["cl_cir"][pop].loc[0]
            # cal_sheet.at["ie_cc", pop] = ratio["cl_cir"][pop].mean()
            # # HCC
            # cal_sheet.at["i_hcp", pop] = ratio["flw_hcc"][pop].loc[0]
            # cal_sheet.at["dc_hcc", pop] = ratio["flw_hcc"][pop].mean()

    # Save calibration file
    cal_sheet.to_excel(hbv.root / 'calibrations' / f"{country}_cal b.xlsx")