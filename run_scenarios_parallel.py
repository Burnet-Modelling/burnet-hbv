from run_scenarios import central_model_run
from hbv_vimc.constants import vimc
import sciris as sc
from hbv_vimc.vimc_outputs import central_results_combined

countries = vimc.keys() 

if __name__ == '__main__':
    sc.parallelize(central_model_run, countries, ncpus=40)
    
central_results_combined(countries)
