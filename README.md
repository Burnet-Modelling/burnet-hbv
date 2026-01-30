# HBV VIMC analysis repository

### Setup

1. Clone repository
2. Go into the `hbv-vimc` folder (the one containing `setup.py`)
3. Run `pip install -e .`

### Usage

```python
import hbv_vimc as hbv
P = hbv.project('AFR')
res = P.run_sim()
```

Some useful module-level variables are 

- `hbv.project(country)` - Load a project, automatically loading the calibration as well
- `hbv.countries` - Dictionary of country/region codes and full names

### Key folders

- `frameworks` contains versios of the framework and YAML configuration
- `databooks` contains databooks for each country/region and corresponding calibrations
  
### Key scripts

- `run_calibration.py` to run the auto calibration routine. Files are saved to `outputs/calibrations` and once finalised, should be copied into the `databooks` folder for use in scenarios/other scripts


## Note: To change databooks (add settings, add/remove populations, or make framework based changes) needs to be done locally. 

Requires: Atomica, Atomica-Tools.

Calibrations should run automatically from the vimc_cal_yaml.py script. Settings run are determined from a list within the script (second last line).

To change calibrations, edits can be made to the model_config_hbv.yaml file. While not currently specified, upper and lower scaling bounds can be defined for each parameter using a list ["par", y_lb, y_ub] format. 

Outputs are saved to the output folder. For each setting a pdf containing calibration plots and spreadsheet of y_factors are saved. 