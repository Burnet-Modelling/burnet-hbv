import atomica as at
import numpy as np
import pandas as pd
from hbv_aus.utils import read_table,_get_github_folder,_get_sharepoint_folder

# Version 1.0 (basic population dynamics, no behaviour populations)

def gen_databook():

    # TODO: Add ability to change reference points for getting/saving files, save names etc.

    # Populations and age bins (as they will be in databook)
    age_bins = ["0-4", "5-14", "15-29", "30-49", "50-64", "65+"]
    sexes = ["M", "F"]
    demog = ["_aus", "_oth", "_fns"] #Born in AUS, Born OS, First Nations (Aboriginal and/or Torres Strait Islander)
    pops = [f"{age}{sex}{demo}" for demo in demog for sex in sexes for age in age_bins]

    # Import data from Sharepoint Folder
    aus_tabs = ["auspop", "ausbirth", "austrans", "ausdeath", "ausarrive", "ausdepart"]
    oth_tabs = ["othpop", "othtrans", "othdeath", "otharrive", "othdepart"]
    fns_tabs = ["atsipop", "atsibirth", "atsitrans", "atsideath"]
    tabs = aus_tabs + oth_tabs + fns_tabs

    pop_data = {}
    for tab in tabs:
        # Loops through table names and saves each as a pandas dataframe
        pop_data[tab] = read_table(_get_sharepoint_folder()+f"Framework Development/Populations/Population Parameters_v1.0.xlsx",
                                   tab)
        pop_data[tab] = pop_data[tab][pop_data[tab]["year"]>=1980]

    # Generate Databook
    F = at.ProjectFramework(_get_github_folder()+f"framework/fw_popsizes.xlsx")
    D = at.ProjectData.new(framework=F, pops=1, transfers=0, tvec=np.arange(1980,2072,1))
    for idx, val in enumerate(pops):
        if idx==0:
            D.rename_pop("pop_0", new_code_name = val, new_full_name = val)
        else:
            D.add_pop(val, val)
    D.add_transfer("age", "aging")

    # Fill the databook and save to GitHub repository

    # Australian-Born Input Parameters
    for age in age_bins:
        for sex in sexes:
            # Popsize
            D.tdve["aus_pop"].ts[f"{age}{sex}_aus"] = at.TimeSeries(t=pop_data["auspop"]["year"],
                                                                    vals = pop_data["auspop"][f"{age}{sex}"],
                                                                    units = "Number")
            # Mortality
            D.tdve["acm_rate"].ts[f"{age}{sex}_aus"] = at.TimeSeries(t=pop_data["ausdeath"]["year"],
                                                                     vals = pop_data["ausdeath"][f"{age}{sex}"],
                                                                     units = "Probability (per year)")
            # Arrivals
            D.tdve["pop_arrive"].ts[f"{age}{sex}_aus"] = at.TimeSeries(t=pop_data["ausarrive"]["year"],
                                                                       vals = pop_data["ausarrive"][f"{age}{sex}"],
                                                                       units = "Number (per year)")
            # Departures
            D.tdve["pop_depart"].ts[f"{age}{sex}_aus"] = at.TimeSeries(t=pop_data["ausdepart"]["year"],
                                                                       vals = pop_data["ausdepart"][f"{age}{sex}"],
                                                                       units = "Number (per year)")
            # Births
            if age == "0-4":
                D.tdve["b_rate"].ts[f"{age}{sex}_aus"] = at.TimeSeries(t=pop_data["ausbirth"]["year"],
                                                                           vals=pop_data["ausbirth"][f"{age}{sex}"],
                                                                           units="Number (per year)")
            else:
                D.tdve["b_rate"].ts[f"{age}{sex}_aus"].assumption = 0

    # Australian-Born Transition (Aging) Matrix
    pop_from, pop_to = age_bins[:-1], age_bins[1:]

    for i,pf in enumerate(pop_from):
        D.transfers[0].ts.append((f"{pf}M_aus", f"{pop_to[i]}M_aus"), at.TimeSeries(pop_data["austrans"]["year"],
                                                                                    pop_data["austrans"][f"{pf}M"],
                                                                                    units = "Rate (per year)"))
        D.transfers[0].ts.append((f"{pf}F_aus", f"{pop_to[i]}F_aus"), at.TimeSeries(pop_data["austrans"]["year"],
                                                                                    pop_data["austrans"][f"{pf}F"],
                                                                                    units = "Rate (per year)"))
    # Overseas-Born Input Parameters
    for age in age_bins:
        for sex in sexes:
            # Popsize
            D.tdve["aus_pop"].ts[f"{age}{sex}_oth"] = at.TimeSeries(t=pop_data["othpop"]["year"],
                                                                    vals = pop_data["othpop"][f"{age}{sex}"],
                                                                    units = "Number")
            # Mortality
            D.tdve["acm_rate"].ts[f"{age}{sex}_oth"] = at.TimeSeries(t=pop_data["othdeath"]["year"],
                                                                     vals = pop_data["othdeath"][f"{age}{sex}"],
                                                                     units = "Probability (per year)")
            # Arrivals
            D.tdve["pop_arrive"].ts[f"{age}{sex}_oth"] = at.TimeSeries(t=pop_data["otharrive"]["year"],
                                                                       vals = pop_data["otharrive"][f"{age}{sex}"],
                                                                       units = "Number (per year)")
            # Departures
            D.tdve["pop_depart"].ts[f"{age}{sex}_oth"] = at.TimeSeries(t=pop_data["othdepart"]["year"],
                                                                       vals = pop_data["othdepart"][f"{age}{sex}"],
                                                                       units = "Number (per year)")
            # Births
            D.tdve["b_rate"].ts[f"{age}{sex}_oth"].assumption = 0

    # Overseas-Born Transition (Aging) Matrix
    for i,pf in enumerate(pop_from):
        D.transfers[0].ts.append((f"{pf}M_oth", f"{pop_to[i]}M_oth"), at.TimeSeries(pop_data["othtrans"]["year"],
                                                                                    pop_data["othtrans"][f"{pf}M"],
                                                                                    units = "Rate (per year)"))
        D.transfers[0].ts.append((f"{pf}F_oth", f"{pop_to[i]}F_oth"), at.TimeSeries(pop_data["othtrans"]["year"],
                                                                                    pop_data["othtrans"][f"{pf}F"],
                                                                                    units = "Rate (per year)"))

    # Aboriginal Torres Strait Islanders Input Parameters
    for age in age_bins:
        for sex in sexes:
            # Popsize
            D.tdve["aus_pop"].ts[f"{age}{sex}_fns"] = at.TimeSeries(t=pop_data["atsipop"]["year"],
                                                                    vals = pop_data["atsipop"][f"{age}{sex}"],
                                                                    units = "Number")
            # Mortality
            D.tdve["acm_rate"].ts[f"{age}{sex}_fns"] = at.TimeSeries(t=pop_data["atsideath"]["year"],
                                                                     vals = pop_data["atsideath"][f"{age}{sex}"],
                                                                     units = "Probability (per year)")
            # Arrivals
            D.tdve["pop_arrive"].ts[f"{age}{sex}_fns"].assumption = 0
            # Departures
            D.tdve["pop_depart"].ts[f"{age}{sex}_fns"].assumption = 0
            # Births
            if age == "0-4":
                D.tdve["b_rate"].ts[f"{age}{sex}_fns"] = at.TimeSeries(t=pop_data["atsibirth"]["year"],
                                                                           vals=pop_data["atsibirth"][f"{age}{sex}"],
                                                                           units="Number (per year)")
            else:
                D.tdve["b_rate"].ts[f"{age}{sex}_fns"].assumption = 0

    # Aboriginal Torres Strait Islanders (Aging) Matrix
    for i,pf in enumerate(pop_from):
        D.transfers[0].ts.append((f"{pf}M_fns", f"{pop_to[i]}M_fns"), at.TimeSeries(pop_data["atsitrans"]["year"],
                                                                                    pop_data["atsitrans"][f"{pf}M"],
                                                                                    units = "Rate (per year)"))
        D.transfers[0].ts.append((f"{pf}F_fns", f"{pop_to[i]}F_fns"), at.TimeSeries(pop_data["atsitrans"]["year"],
                                                                                    pop_data["atsitrans"][f"{pf}F"],
                                                                                    units = "Rate (per year)"))
    D.save(_get_github_folder() + f"databook/db_demographics.xlsx")
