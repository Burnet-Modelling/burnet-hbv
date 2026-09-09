import atomica as at
import numpy as np
import pandas as pd
from hbv_aus.utils import _get_github_folder,_get_sharepoint_folder

# Version 1.0 (basic population dynamics, no behaviour populations)
def make_age_bands(max_age=84, band_width=5):
    bands = [f"{i}-{i+band_width-1}" for i in range(0, max_age+1, band_width)]
    bands.append(f"{max_age+1}+")
    return bands


def gen_db_hbv_v2_0():
    """
    Produces databook using hbv_fw_v2.0.xlsx (version control for archiving)
    """
    # Import Framework (hbv_fw_v2.0.xlsx) to produce databook
    F = at.ProjectFramework(_get_github_folder() + f"framework/hbv_fw_v2.0.xlsx")

    age_bins = ["0-4", "5-14", "15-29", "30-49", "50-64", "65+"]# only for use in testing, will be hard coded in practice
    sexes = ["_M", "_F"]
    demog = ["atsi_"]
    pops = [f"{demo}{age}{sex}" for demo in demog for sex in sexes for age in age_bins]

    # Generate Databook
    D = at.ProjectData.new(framework=F, pops=1, transfers=0, tvec=np.arange(1980, 2072, 1))
    for idx, val in enumerate (pops):
        if idx == 0:
            D.rename_pop("pop_0", new_code_name=val, new_full_name=val)
        else:
            D.add_pop(val, val)
    D.add_transfer("age", "aging")

    # Import Data from Excel
    # Aboriginal Torres Strait Islander Population Data
    atsi_pop_in = pd.read_excel(_get_sharepoint_folder()+f"atsi_pop_proj.xlsx", sheet_name="Total Population")
    atsi_death_in = pd.read_excel(_get_sharepoint_folder()+f"atsi_pop_proj.xlsx", sheet_name="Deaths")
    atsi_birth_in = pd.read_excel(_get_sharepoint_folder()+f"atsi_pop_proj.xlsx", sheet_name="Births")
    atsi_age_in = pd.read_excel(_get_sharepoint_folder()+f"atsi_pop_proj.xlsx", sheet_name="Aging")

    # Migrants

    # Other Australian Pop


    # Total Population (temp_alive)
    for age in age_bins:
        D.tdve["temp_alive"].ts[f"atsi_{age}_M"] = at.TimeSeries(t=atsi_pop_in[atsi_pop_in["pop"]==age]["Year"],
                                                                 vals = atsi_pop_in[atsi_pop_in["pop"]==age]["Male"],
                                                                 units="Number")

        D.tdve["temp_alive"].ts[f"atsi_{age}_F"] = at.TimeSeries(t=atsi_pop_in[atsi_pop_in["pop"] == age]["Year"],
                                                                 vals=atsi_pop_in[atsi_pop_in["pop"] == age]["Female"],
                                                                 units="Number")
    # All Cause Mortality (death_rate):
    for age in age_bins:
        D.tdve["death_rate"].ts[f"atsi_{age}_M"] = at.TimeSeries(t=atsi_death_in[atsi_death_in["pop"]==age]["Year"],
                                                                 vals = atsi_death_in[atsi_death_in["pop"]==age]["Male"],
                                                                 units="Rate (per year)")

        D.tdve["death_rate"].ts[f"atsi_{age}_F"] = at.TimeSeries(t=atsi_death_in[atsi_death_in["pop"]==age]["Year"],
                                                                 vals = atsi_death_in[atsi_death_in["pop"]==age]["Female"],
                                                                 units="Rate (per year)")
    # Births (birth_rate):
    #atsi_births_sum = atsi_birth_in.iloc[:, np.r_[0, 4,5]]
    #atsi_births_sum = atsi_births_sum.groupby('Year', as_index=False).sum()

    for age in age_bins:
        if age == "0-4":
            D.tdve["birth_rate"].ts[f"atsi_{age}_M"] = at.TimeSeries(t=atsi_birth_in["Year"], vals = atsi_birth_in["Male"], units = "Number (per year)")
            D.tdve["birth_rate"].ts[f"atsi_{age}_F"] = at.TimeSeries(t=atsi_birth_in["Year"], vals = atsi_birth_in["Female"], units = "Number (per year)")
        else:
            D.tdve["birth_rate"].ts[f"atsi_{age}_M"].assumption = 0
            D.tdve["birth_rate"].ts[f"atsi_{age}_F"].assumption = 0

    # Migration (emig_rate, imig_rate)
    for age in age_bins:
        D.tdve["imig_rate"].ts[f"atsi_{age}_M"].assumption = 0
        D.tdve["emig_rate"].ts[f"atsi_{age}_M"].assumption = 0
        D.tdve["imig_rate"].ts[f"atsi_{age}_F"].assumption = 0
        D.tdve["emig_rate"].ts[f"atsi_{age}_F"].assumption = 0

    # Aging
    #atsi_age_in = atsi_age_in[atsi_age_in["Age group"]!="85+"]
    pop_from, pop_to = age_bins[:-1], age_bins[1:]

    for idx, age in enumerate(pop_from):
        D.transfers[0].ts.append((f"atsi_{age}_M", f"atsi_{pop_to[idx]}_M"), at.TimeSeries(atsi_age_in[atsi_age_in["pop_from"]==age]["Year"],
                                                                                    atsi_age_in[atsi_age_in["pop_from"]==age]["Male"],
                                                                                    units="Rate (per year)"))
        D.transfers[0].ts.append((f"atsi_{age}_F", f"atsi_{pop_to[idx]}_F"), at.TimeSeries(atsi_age_in[atsi_age_in["pop_from"]==age]["Year"],
                                                                                    atsi_age_in[atsi_age_in["pop_from"]==age]["Female"],
                                                                                    units="Rate (per year)"))


    D.save(_get_github_folder() + f"databook/hbv_db_v2.0_test.xlsx")










def gen_databook():

    # TODO: Add ability to change reference points for getting/saving files, save names etc.

    # Populations and age bins (as they will be in databook)
    age_bins = ["0-4", "5-14", "15-29", "30-49", "50-64", "65+"]
    sexes = ["M", "F"]
    demog = ["_aus", "_oth", "_fns"]  # Born in AUS, Born OS, First Nations (Aboriginal and/or Torres Strait Islander)
    pops = [f"{age}{sex}{demo}" for demo in demog for sex in sexes for age in age_bins]

    # Import data from Sharepoint Folder
    aus_tabs = ["auspop", "ausbirth", "austrans", "ausdeath", "ausarrive", "ausdepart"]
    oth_tabs = ["othpop", "othtrans", "othdeath", "otharrive", "othdepart"]
    fns_tabs = ["atsipop", "atsibirth", "atsitrans", "atsideath"]
    tabs = aus_tabs + oth_tabs + fns_tabs

    pop_data = {}
    for tab in tabs:
        # Loops through table names and saves each as a pandas dataframe
        pop_data[tab] = read_table(
            _get_sharepoint_folder() + f"Framework Development/Populations/Population Parameters_v1.0.xlsx",
            tab)
        pop_data[tab] = pop_data[tab][pop_data[tab]["year"] >= 1980]

    # Generate Databook
    F = at.ProjectFramework(_get_github_folder() + f"framework/fw_popsizes.xlsx")
    D = at.ProjectData.new(framework=F, pops=1, transfers=0, tvec=np.arange(1980, 2072, 1))
    for idx, val in enumerate(pops):
        if idx == 0:
            D.rename_pop("pop_0", new_code_name=val, new_full_name=val)
        else:
            D.add_pop(val, val)
    D.add_transfer("age", "aging")

    # Fill the databook and save to GitHub repository

    # Australian-Born Input Parameters
    for age in age_bins:
        for sex in sexes:
            # Popsize
            D.tdve["aus_pop"].ts[f"{age}{sex}_aus"] = at.TimeSeries(t=pop_data["auspop"]["year"],
                                                                    vals=pop_data["auspop"][f"{age}{sex}"],
                                                                    units="Number")
            # Mortality
            D.tdve["acm_rate"].ts[f"{age}{sex}_aus"] = at.TimeSeries(t=pop_data["ausdeath"]["year"],
                                                                     vals=pop_data["ausdeath"][f"{age}{sex}"],
                                                                     units="Probability (per year)")
            # Arrivals
            D.tdve["pop_arrive"].ts[f"{age}{sex}_aus"] = at.TimeSeries(t=pop_data["ausarrive"]["year"],
                                                                       vals=pop_data["ausarrive"][f"{age}{sex}"],
                                                                       units="Number (per year)")
            # Departures
            D.tdve["pop_depart"].ts[f"{age}{sex}_aus"] = at.TimeSeries(t=pop_data["ausdepart"]["year"],
                                                                       vals=pop_data["ausdepart"][f"{age}{sex}"],
                                                                       units="Number (per year)")
            # Births
            if age == "0-4":
                D.tdve["b_rate"].ts[f"{age}{sex}_aus"] = at.TimeSeries(t=pop_data["ausbirth"]["year"],
                                                                       vals=pop_data["ausbirth"][f"{age}{sex}"],
                                                                       units="Number (per year)")
            else:
                D.tdve["b_rate"].ts[f"{age}{sex}_aus"].assumption = 0

    # Australian-Born Transition (Aging) Matrix
    pop_from, pop_to = age_bins[:-1], age_bins[1:]

    for i, pf in enumerate(pop_from):
        D.transfers[0].ts.append((f"{pf}M_aus", f"{pop_to[i]}M_aus"), at.TimeSeries(pop_data["austrans"]["year"],
                                                                                    pop_data["austrans"][f"{pf}M"],
                                                                                    units="Rate (per year)"))
        D.transfers[0].ts.append((f"{pf}F_aus", f"{pop_to[i]}F_aus"), at.TimeSeries(pop_data["austrans"]["year"],
                                                                                    pop_data["austrans"][f"{pf}F"],
                                                                                    units="Rate (per year)"))
    # Overseas-Born Input Parameters
    for age in age_bins:
        for sex in sexes:
            # Popsize
            D.tdve["aus_pop"].ts[f"{age}{sex}_oth"] = at.TimeSeries(t=pop_data["othpop"]["year"],
                                                                    vals=pop_data["othpop"][f"{age}{sex}"],
                                                                    units="Number")
            # Mortality
            D.tdve["acm_rate"].ts[f"{age}{sex}_oth"] = at.TimeSeries(t=pop_data["othdeath"]["year"],
                                                                     vals=pop_data["othdeath"][f"{age}{sex}"],
                                                                     units="Probability (per year)")
            # Arrivals
            D.tdve["pop_arrive"].ts[f"{age}{sex}_oth"] = at.TimeSeries(t=pop_data["otharrive"]["year"],
                                                                       vals=pop_data["otharrive"][f"{age}{sex}"],
                                                                       units="Number (per year)")
            # Departures
            D.tdve["pop_depart"].ts[f"{age}{sex}_oth"] = at.TimeSeries(t=pop_data["othdepart"]["year"],
                                                                       vals=pop_data["othdepart"][f"{age}{sex}"],
                                                                       units="Number (per year)")
            # Births
            D.tdve["b_rate"].ts[f"{age}{sex}_oth"].assumption = 0

    # Overseas-Born Transition (Aging) Matrix
    for i, pf in enumerate(pop_from):
        D.transfers[0].ts.append((f"{pf}M_oth", f"{pop_to[i]}M_oth"), at.TimeSeries(pop_data["othtrans"]["year"],
                                                                                    pop_data["othtrans"][f"{pf}M"],
                                                                                    units="Rate (per year)"))
        D.transfers[0].ts.append((f"{pf}F_oth", f"{pop_to[i]}F_oth"), at.TimeSeries(pop_data["othtrans"]["year"],
                                                                                    pop_data["othtrans"][f"{pf}F"],
                                                                                    units="Rate (per year)"))

    # Aboriginal Torres Strait Islanders Input Parameters
    for age in age_bins:
        for sex in sexes:
            # Popsize
            D.tdve["aus_pop"].ts[f"{age}{sex}_fns"] = at.TimeSeries(t=pop_data["atsipop"]["year"],
                                                                    vals=pop_data["atsipop"][f"{age}{sex}"],
                                                                    units="Number")
            # Mortality
            D.tdve["acm_rate"].ts[f"{age}{sex}_fns"] = at.TimeSeries(t=pop_data["atsideath"]["year"],
                                                                     vals=pop_data["atsideath"][f"{age}{sex}"],
                                                                     units="Probability (per year)")
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
    for i, pf in enumerate(pop_from):
        D.transfers[0].ts.append((f"{pf}M_fns", f"{pop_to[i]}M_fns"), at.TimeSeries(pop_data["atsitrans"]["year"],
                                                                                    pop_data["atsitrans"][f"{pf}M"],
                                                                                    units="Rate (per year)"))
        D.transfers[0].ts.append((f"{pf}F_fns", f"{pop_to[i]}F_fns"), at.TimeSeries(pop_data["atsitrans"]["year"],
                                                                                    pop_data["atsitrans"][f"{pf}F"],
                                                                                    units="Rate (per year)"))
    D.save(_get_github_folder() + f"databook/db_demographics.xlsx")


