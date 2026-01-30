from run_scenarios import stochastic_model_run
from hbv_vimc.constants import vimc
from hbv_vimc.plotting import comparison_plots
import sciris as sc

countries = vimc.keys()

if __name__ == '__main__':
   comparison_plots(countries)
