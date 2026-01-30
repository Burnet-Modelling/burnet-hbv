import pandas as pd
import atomica as at
import numpy as np
import hbv_vimc as hbv
import hbv_vimc.constants as cts
import hbv_vimc.sampling as smp
import sciris as sc

def get_scenarios(country:str):
    scenarios = {}

    # Baselines (with and without HepB-BD)

    baseline = pd.read_csv(hbv.root/"baseline.csv")
    baseline_nbd = pd.read_csv(hbv.root / "baseline_nbd.csv")

    bl, bl_nbd = baseline[baseline["country_code"] == country], baseline_nbd[baseline_nbd["country_code"] == country]
    bl, bl_nbd = bl.set_index("vaccine"), bl_nbd.set_index("vaccine")

    scen =  at.ParameterScenario(name="Baseline")
    scen.add("bd", "0-0M", list(np.arange(1990, 2101, 1)), list(bl.iloc[1, 6:]))
    scen.add("bd", "0-0F", list(np.arange(1990, 2101, 1)), list(bl.iloc[1, 6:]))
    scen.add("hb3", "0-0M", list(np.arange(1990, 2101, 1)), list(bl.iloc[0, 6:]))
    scen.add("hb3", "0-0F", list(np.arange(1990, 2101, 1)), list(bl.iloc[0, 6:]))
    scenarios[scen.name] = scen

    scen =  at.ParameterScenario(name="Baseline No BD")
    scen.add("bd", "0-0M", list(np.arange(1990, 2101, 1)), list(bl_nbd.iloc[1, 6:]))
    scen.add("bd", "0-0F", list(np.arange(1990, 2101, 1)), list(bl_nbd.iloc[1, 6:]))
    scen.add("hb3", "0-0M", list(np.arange(1990, 2101, 1)), list(bl_nbd.iloc[0, 6:]))
    scen.add("hb3", "0-0F", list(np.arange(1990, 2101, 1)), list(bl_nbd.iloc[0, 6:]))
    scenarios[scen.name] = scen

    # No Vaccine Counterfactual

    no_vax = pd.read_csv(hbv.root/"no_vax.csv")
    nv = no_vax[no_vax["country_code"] == country]
    nv = nv.set_index("vaccine")

    scen = at.ParameterScenario(name="No Vaccination")
    scen.add("bd", "0-0M", list(np.arange(1990, 2101, 1)), list(nv.iloc[1, 6:]))
    scen.add("bd", "0-0F", list(np.arange(1990, 2101, 1)), list(nv.iloc[1, 6:]))
    scen.add("hb3", "0-0M", list(np.arange(1990, 2101, 1)), list(nv.iloc[0, 6:]))
    scen.add("hb3", "0-0F", list(np.arange(1990, 2101, 1)), list(nv.iloc[0, 6:]))
    scenarios[scen.name] = scen

    # IA

    ia = pd.read_csv(hbv.root/"ia.csv")
    ia_nbd = pd.read_csv(hbv.root/"ia_nbd.csv")

    ia, ia_nbd = ia[ia["country_code"] == country], ia_nbd[ia_nbd["country_code"] == country]
    ia, ia_nbd = ia.set_index("vaccine"), ia_nbd.set_index("vaccine")

    scen = at.ParameterScenario(name="IA")
    scen.add("bd", "0-0M", list(np.arange(1990, 2101, 1)), list(ia.iloc[1, 6:]))
    scen.add("bd", "0-0F", list(np.arange(1990, 2101, 1)), list(ia.iloc[1, 6:]))
    scen.add("hb3", "0-0M", list(np.arange(1990, 2101, 1)), list(ia.iloc[0, 6:]))
    scen.add("hb3", "0-0F", list(np.arange(1990, 2101, 1)), list(ia.iloc[0, 6:]))
    scenarios[scen.name] = scen

    scen = at.ParameterScenario(name="IA No BD")
    scen.add("bd", "0-0M", list(np.arange(1990, 2101, 1)), list(ia_nbd.iloc[1, 6:]))
    scen.add("bd", "0-0F", list(np.arange(1990, 2101, 1)), list(ia_nbd.iloc[1, 6:]))
    scen.add("hb3", "0-0M", list(np.arange(1990, 2101, 1)), list(ia_nbd.iloc[0, 6:]))
    scen.add("hb3", "0-0F", list(np.arange(1990, 2101, 1)), list(ia_nbd.iloc[0, 6:]))
    scenarios[scen.name] = scen

    # Bluesky

    bsky = pd.read_csv(hbv.root / "bsky.csv")
    bsky_nbd = pd.read_csv(hbv.root / "bsky_nbd.csv")

    bsky, bsky_nbd = bsky[bsky["country_code"] == country], bsky_nbd[bsky_nbd["country_code"] == country]
    bsky, bsky_nbd = bsky.set_index("vaccine"), bsky_nbd.set_index("vaccine")

    scen = at.ParameterScenario(name="Bluesky")
    scen.add("bd", "0-0M", list(np.arange(1990, 2101, 1)), list(bsky.iloc[1, 6:]))
    scen.add("bd", "0-0F", list(np.arange(1990, 2101, 1)), list(bsky.iloc[1, 6:]))
    scen.add("hb3", "0-0M", list(np.arange(1990, 2101, 1)), list(bsky.iloc[0, 6:]))
    scen.add("hb3", "0-0F", list(np.arange(1990, 2101, 1)), list(bsky.iloc[0, 6:]))
    scenarios[scen.name] = scen

    scen = at.ParameterScenario(name="Bluesky No BD")
    scen.add("bd", "0-0M", list(np.arange(1990, 2101, 1)), list(bsky_nbd.iloc[1, 6:]))
    scen.add("bd", "0-0F", list(np.arange(1990, 2101, 1)), list(bsky_nbd.iloc[1, 6:]))
    scen.add("hb3", "0-0M", list(np.arange(1990, 2101, 1)), list(bsky_nbd.iloc[0, 6:]))
    scen.add("hb3", "0-0F", list(np.arange(1990, 2101, 1)), list(bsky_nbd.iloc[0, 6:]))
    scenarios[scen.name] = scen

    return scenarios

def run_central_scenarios(country:str):
    
    """ Produce central scenario estimates for each country.
    Basic plotting functionality for these outputs (without )
    """
    # Import Project and Calibration
    P=hbv.project(country, calibrated=True)
    P.settings.update_time_vector(dt=1)
    
    # calibration_file = hbv.root/'calibrations'/f'{country}_cal.xlsx'
    # P.parsets[0].load_calibration(calibration_file)
    
    # Run Scenarios
    scens=hbv.get_scenarios(country)
    central={}
    
    for scen in scens.values():
        central[scen.name]=scen.run(P, P.parsets[0])
    
    return central 


def run_stochastic_scenarios(country:str):
    
    """ Runs the stochastic scenarios for each country and produces
    outputs for conversion into required documentation and plotting"""
        
    # Import Project and Calibration
    P=hbv.project(country, calibrated=True)
    P.settings.update_time_vector(dt=1)
    pops = list(P.data.pops.keys())
    parset = P.parsets[0]
    s_yr=hbv.s_yr
    e_yr=hbv.e_yr
    
    
    # Sample parameters and plot outcomes
    #np.random.seed(20230906)
    runs= cts.runs #define the number of runs in the constants.py script.
    samples, point_estimates = smp.pars_stochastic(country)
    smp.sample_plots(country, samples, point_estimates, smp.parameters) #return plots to show that sampling has worked
    
    
    # MTCT parameters
    ch_pop = ["0-0M", "0-0F"]
    ch_var = ["ci_p", "eag_hvl", "eag_ve", "hvl_trisk", "sag_ve", "sag_hvl" ]
    
    # Calibrated parameters
    m_020_pop = ["0-0M", "1-4M", "5-9M", "10-19M"]
    m_020_var =  ["cc_020_m", "eag_020_m", "hcc_020_m", "sag_020_m"]
    
    f_020_pop = ["0-0F", "1-4F", "5-9F", "10-19F"]
    f_020_var =  ["cc_020_f", "eag_020_f", "hcc_020_f", "sag_020_f"]
    
    m_2040_pop = ["20-29M", "30-39M"]
    m_2040_var =  ["cc_2040_m", "eag_2040_m", "hcc_2040_m", "sag_2040_m"]
    
    f_2040_pop = ["20-29F", "30-39F"]
    f_2040_var =  ["cc_2040_f", "eag_2040_f", "hcc_2040_f", "sag_2040_f"]

    m_4060_pop = ["40-49M", "50-59M"]
    m_4060_var =  ["cc_4060_m", "eag_4060_m", "hcc_4060_m", "sag_4060_m"]
    
    f_4060_pop = ["40-49F", "50-59F"]
    f_4060_var =  ["cc_4060_f", "eag_4060_f", "hcc_4060_f", "sag_4060_f"]
    
    m_60_pop = ["60-69M", "70-79M", "80-89M", "90+M"]
    m_60_var =  ["cc_60_m", "eag_60_m", "hcc_60_m", "sag_60_m"]
    
    f_60_pop = ["60-69F", "70-79F", "80-89F", "90+F"]
    f_60_var =  ["cc_60_f", "eag_60_f", "hcc_60_f", "sag_60_f"]
    
    
    # All population parameters
    all_var = ["m_dc", "m_hcc", "te_cc_dc", "te_cc_hcc", "te_dc_cc", "te_dc_hcc",
               "te_icl_cc", "te_icl_hcc", "te_ict_hcc", "te_ie_cc", "te_ie_hcc",
               "te_m_dc", "te_m_hcc"]
    
    

    vimc_outputs=["alive", "tot_inc", ":dd_hbv", "dalys", "yll"]
    vimc_name=["cohort_size", "cases", "deaths", "dalys", "yll"]
    age_bins = ["year", "0-0", "1-4", "5-9", "10-19", "20-29", "30-39", 
                "40-49", "50-59", "60-69", "70-79", "80-89", "90+"]
    
    # FOR LOOP TO WORK, NEEDS TO HAVE CONSTANT INDEXATION
    plot_names = ["t_prev", "c_prev", "i_hcc", "d_hbv"]
    plot_outs = ["prev", "prev", "flw_hcc", ":dd_hbv"]
    plot_pops = ["total", {"0-4":["0-0F", "0-0M","1-4F", "1-4M"]}, "total", "total"]
    plot_weight = ["weighted", "weighted", "sum", "sum"]
    
    scens=hbv.get_scenarios(country)
    stochastic={}
    stoch_out = {}
    stoch_plot={}
    st_plot={}
    
    # Initialize VIMC output dictionaries
    for scen in scens.values():
        stochastic[scen.name]={}
        stoch_out[scen.name]={}
        for name in vimc_name:
            stoch_out[scen.name][name]=pd.DataFrame()
            for run in range(1,runs+1):
                stochastic[scen.name][f"{name}_{run}"]=pd.DataFrame(columns=age_bins)
    
    # Initialize plotting dictionaries
    for scen in scens.values():
        stoch_plot[scen.name] = {}
        st_plot[scen.name]={}
        for name in plot_names:
            st_plot[scen.name][name]=pd.DataFrame()
            for run in range(1, runs+1):
                stoch_plot[scen.name][f"{name}_{run}"]=pd.DataFrame(columns=[f"run_{run}"])
    
    # Fill dictionaries with data
    for scen in scens.values():
        for run in range(1, runs+1):
            
            #TODO: set constant population in samples - assures constant values are kept constant

            # 0-0 pops only
            for cp in ch_pop:
                for cv in ch_var:
                    parset.get_par(cv).ts[cp].assumption = samples[cv][cp][run-1]
            
            # 0-20M calibration
            for ap in m_020_pop:
                parset.get_par("cc_020_m").ts[ap].assumption = samples["cc_020_m"]["0-0M"][run-1]
                parset.get_par("eag_020_m").ts[ap].assumption = samples["eag_020_m"]["0-0M"][run-1]
                parset.get_par("hcc_020_m").ts[ap].assumption = samples["hcc_020_m"]["0-0M"][run-1]
                parset.get_par("sag_020_m").ts[ap].assumption = samples["sag_020_m"]["0-0M"][run-1]
                #for av in m_020_var:
                    #parset.get_par(av).ts[ap].assumption = samples[av][ap][run-1]
            
            # 0-20M calibration
            for bp in f_020_pop:
                parset.get_par("cc_020_f").ts[bp].assumption = samples["cc_020_m"]["0-0M"][run-1]
                parset.get_par("eag_020_f").ts[bp].assumption = samples["eag_020_m"]["0-0M"][run-1]
                parset.get_par("hcc_020_f").ts[bp].assumption = samples["hcc_020_m"]["0-0M"][run-1]
                parset.get_par("sag_020_f").ts[bp].assumption = samples["sag_020_m"]["0-0M"][run-1]
                #for bv in f_020_var:
                    #parset.get_par(bv).ts[bp].assumption = samples[bv][bp][run-1]
                    
            
            # 20-40M calibration
            for dp in m_2040_pop:
                parset.get_par("cc_2040_m").ts[dp].assumption = samples["cc_020_m"]["0-0M"][run-1]
                parset.get_par("eag_2040_m").ts[dp].assumption = samples["eag_020_m"]["0-0M"][run-1]
                parset.get_par("hcc_2040_m").ts[dp].assumption = samples["hcc_020_m"]["0-0M"][run-1]
                parset.get_par("sag_2040_m").ts[dp].assumption = samples["sag_020_m"]["0-0M"][run-1]
                #for dv in m_2040_var:
                    #parset.get_par(dv).ts[dp].assumption = samples[dv][dp][run-1]
            
            # 20-40F calibration
            for ep in f_2040_pop:
                parset.get_par("cc_2040_f").ts[ep].assumption = samples["cc_020_m"]["0-0M"][run-1]
                parset.get_par("eag_2040_f").ts[ep].assumption = samples["eag_020_m"]["0-0M"][run-1]
                parset.get_par("hcc_2040_f").ts[ep].assumption = samples["hcc_020_m"]["0-0M"][run-1]
                parset.get_par("sag_2040_f").ts[ep].assumption = samples["sag_020_m"]["0-0M"][run-1]
                #for ev in f_2040_var:
                    #parset.get_par(ev).ts[ep].assumption = samples[ev][ep][run-1]
            
            # 40-60M calibration
            for fp in m_4060_pop:
                parset.get_par("cc_4060_m").ts[fp].assumption = samples["cc_020_m"]["0-0M"][run-1]
                parset.get_par("eag_4060_m").ts[fp].assumption = samples["eag_020_m"]["0-0M"][run-1]
                parset.get_par("hcc_4060_m").ts[fp].assumption = samples["hcc_020_m"]["0-0M"][run-1]
                parset.get_par("sag_4060_m").ts[fp].assumption = samples["sag_020_m"]["0-0M"][run-1]
                #for fv in m_4060_var:
                    #parset.get_par(fv).ts[fp].assumption = samples[fv][fp][run-1]
            
            # 40-60F calibration
            for gp in f_4060_pop:
                parset.get_par("cc_4060_f").ts[gp].assumption = samples["cc_020_m"]["0-0M"][run-1]
                parset.get_par("eag_4060_f").ts[gp].assumption = samples["eag_020_m"]["0-0M"][run-1]
                parset.get_par("hcc_4060_f").ts[gp].assumption = samples["hcc_020_m"]["0-0M"][run-1]
                parset.get_par("sag_4060_f").ts[gp].assumption = samples["sag_020_m"]["0-0M"][run-1]
                #for gv in f_4060_var:
                    #parset.get_par(fv).ts[fp].assumption = samples[fv][fp][run-1]
            
            # 60+M calibration
            for hp in m_60_pop:
                parset.get_par("cc_60_m").ts[hp].assumption = samples["cc_020_m"]["0-0M"][run-1]
                parset.get_par("eag_60_m").ts[hp].assumption = samples["eag_020_m"]["0-0M"][run-1]
                parset.get_par("hcc_60_m").ts[hp].assumption = samples["hcc_020_m"]["0-0M"][run-1]
                parset.get_par("sag_60_m").ts[hp].assumption = samples["sag_020_m"]["0-0M"][run-1]
                #for hv in m_60_var:
                    #parset.get_par(hv).ts[hp].assumption = samples[hv][hp][run-1]
            
            # 60+F calibration
            for ip in f_60_pop:
                parset.get_par("cc_60_f").ts[ip].assumption = samples["cc_020_m"]["0-0M"][run-1]
                parset.get_par("eag_60_f").ts[ip].assumption = samples["eag_020_m"]["0-0M"][run-1]
                parset.get_par("hcc_60_f").ts[ip].assumption = samples["hcc_020_m"]["0-0M"][run-1]
                parset.get_par("sag_60_f").ts[ip].assumption = samples["sag_020_m"]["0-0M"][run-1]
                #for iv in f_60_var:
                    #parset.get_par(iv).ts[ip].assumption = samples[iv][ip][run-1]
            
            # All pops
            for tp in pops:
                for tv in all_var:
                    parset.get_par(tv).ts[tp].assumption = samples[tv]["0-0M"][run-1] 
                    
            # Calibration factors (remember - multiply, not direct replacement!)
            
            # Run model
            stoc_res = scen.run(P, parset)
            
            # Extract results for VIMC output (alive, chb_cases, deaths, dalys, yll)
            for idx, out in enumerate(vimc_outputs):
                for age in age_bins:
                    if age == "year":
                        stochastic[scen.name][f"{vimc_name[idx]}_{run}"][age] = np.arange(s_yr, e_yr, 1)
                    else:
                        stochastic[scen.name][f"{vimc_name[idx]}_{run}"][age] = at.PlotData(stoc_res, out, {"pop": [f"{age}M", f"{age}F"]}, t_bins=1).series[0].vals
                    
                    stochastic[scen.name][f"{vimc_name[idx]}_{run}"]["run_id"] = run
            
            # Extract results for plotting 
            for i, pname in enumerate(plot_names):
                    stoch_plot[scen.name][f"{pname}_{run}"][f"run_{run}"] = at.PlotData(stoc_res, plot_outs[i], plot_pops[i], pop_aggregation=plot_weight[i], t_bins=1).series[0].vals
    
    # Output to use for VIMC upload format
    for scen in scens.values():
        for name in vimc_name:
            for run in range(1, runs+1):
                stoch_out[scen.name][name]=pd.concat([stoch_out[scen.name][name], stochastic[scen.name][f"{name}_{run}"]])
    
    # Output for plotting stochastic outcomes
    for scen in scens.values():
        for name in plot_names:
            for run in range(1, runs+1):
                st_plot[scen.name][name]=pd.concat([st_plot[scen.name][name], stoch_plot[scen.name][f"{name}_{run}"]], axis=1)
    
    return samples, stoch_out, stoch_plot
