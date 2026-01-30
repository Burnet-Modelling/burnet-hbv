import atomica as at
from at_tools import *
import sciris as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import hbv_vimc.scenarios as hbv
import hbv_vimc.constants as cts
import hbv_vimc.vimc_outputs as vimc
import hbv_vimc.plotting as plot

# %% Run vaccination scenarios (central only currently; uncertainity and VIMC outputs will be done in here)
uncert = True


def central_model_run(country:str):

    # Run Scenarios (central estimates)
    central = hbv.run_central_scenarios(country)

    # Return Results in VIMC format
    vimc.central_results(country, central)
    # TODO: Add basic plotting functionality to visualize central results (for faster debugging)


def stochastic_model_run(country:str):

    # Run Scenarios (central estimates)
    central = hbv.run_central_scenarios(country)

    # Return Results in VIMC format
    vimc.central_results(country, central)

    # Run Scenarios (Stochastic)
    samples, stoch_out, st_plot = hbv.run_stochastic_scenarios(country)

    # Return Sampled Paramaters in VIMC format
    vimc.input_results(country, samples)

    # Return Results in VIMC format
    vimc.stochastic_results(country, stoch_out)

    # Produce basic plot with 95%CrI for outcomes and central measures
    plot.stochastic_outcome_plots(country, central, st_plot)



if __name__ == '__main__':

    if uncert is False:
       sc.parallelize(central_model_run, cts.vimc.keys(), die=False)
   
    if uncert is True:
       sc.parallelize(stochastic_model_run, cts.vimc.keys(), die=False)  # cts.countries.keys()

