import numpy as np
import hbv_vimc as hbv
import hbv_vimc.constants as cts
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scipy.stats as stats

# Define sampling parameters and their format
# Optionally specify a lower bound & upper bound for probability and proportion parameters (otherwise enter None)
# Optionally specify a standard deviation for other parameters (otherwise enter None)

parameters = { 
    # Natural history
    'nathis': { 
        'ci_p': {'format': 'proportion', 'lower_bound': 0.84, 'upper_bound': 0.93, 'std_dev': None},
        #'m_acu': {'format': 'proportion', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'm_dc': {'format': 'probability', 'lower_bound': 0.16, 'upper_bound': 0.5, 'std_dev': None},
        'm_hcc': {'format': 'probability', 'lower_bound': 0.3, 'upper_bound': 0.7, 'std_dev': None}
    }, 
    # Treatment efficacy
    'trtinv': { 
        'te_dc_cc': {'format': 'probability', 'lower_bound': 0.1, 'upper_bound': 0.25, 'std_dev': None},
        'te_icl_cc': {'format': 'probability', 'lower_bound': 0.1, 'upper_bound': 0.4, 'std_dev': None},
        'te_cc_dc': {'format': 'probability', 'lower_bound': 0.35, 'upper_bound': 0.75, 'std_dev': None},
        'te_cc_hcc': {'format': 'probability', 'lower_bound': 0.5, 'upper_bound': 0.85, 'std_dev': None},
        'te_ie_cc': {'format': 'probability', 'lower_bound': 0.1, 'upper_bound': 0.4, 'std_dev': None},
        'te_m_dc': {'format': 'probability', 'lower_bound': 0.25, 'upper_bound': 0.7, 'std_dev': None},
        'te_m_hcc': {'format': 'probability', 'lower_bound': 0.5, 'upper_bound': 1, 'std_dev': None}, # point estimate was outside bounds
        'te_ict_hcc': {'format': 'probability', 'lower_bound': 0.05, 'upper_bound': 0.5, 'std_dev': None},
        'te_ie_hcc': {'format': 'probability', 'lower_bound': 0.05, 'upper_bound': 0.5, 'std_dev': None},
        'te_icl_hcc': {'format': 'probability', 'lower_bound': 0.05, 'upper_bound': 0.5, 'std_dev': None},
        'te_dc_hcc': {'format': 'probability', 'lower_bound': 0.5, 'upper_bound': 0.85, 'std_dev': None} # point estimate was outside bounds
    }, 
    # Vaccination efficacy
    'vacinv': {
        'eag_ve': {'format': 'probability', 'lower_bound': 0.6, 'upper_bound': 0.9, 'std_dev': None},
        'sag_ve': {'format': 'probability', 'lower_bound': 0.85, 'upper_bound': 1, 'std_dev': None},
        'eag_hvl': {'format': 'proportion', 'lower_bound': 0.84, 'upper_bound': 0.96, 'std_dev': None},
        'sag_hvl': {'format': 'proportion', 'lower_bound': 0.05, 'upper_bound': 0.16, 'std_dev': None},
        'hvl_trisk': {'format': 'proportion', 'lower_bound': 0.7, 'upper_bound': 1, 'std_dev': None}
        #,'hb3_ve': {'format': 'probability', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        #'mav_ve': {'format': 'probability', 'lower_bound': None, 'upper_bound': None, 'std_dev': None}
    }
    # # Calibration parameters
    , 'calibration': {
    #     'ch_pop_sus': {'format': 'proportion', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
    #     'ad_pop_sus': {'format': 'probability', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
    #     # 'hvl_cal': {'format': 'proportion', 'lower_bound': None, 'upper_bound': None, 'std_dev': None}, # not in databook, use 'hvl_trisk' instead?
        'eag_020_m': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'eag_020_f': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'eag_2040_m': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'eag_2040_f': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'eag_4060_m': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'eag_4060_f': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'eag_60_m': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'eag_60_f': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'sag_020_m': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'sag_020_f': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'sag_2040_m': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'sag_2040_f': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'sag_4060_m': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'sag_4060_f': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'sag_60_m': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'sag_60_f': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'cc_020_m': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'cc_020_f': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'cc_2040_m': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'cc_2040_f': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'cc_4060_m': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'cc_4060_f': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'cc_60_m': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'cc_60_f': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'hcc_020_m': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'hcc_020_f': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'hcc_2040_m': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'hcc_2040_f': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'hcc_4060_m': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'hcc_4060_f': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'hcc_60_m': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
        'hcc_60_f': {'format': 'number', 'lower_bound': None, 'upper_bound': None, 'std_dev': None},
     }
} 


def pars_stochastic(country):
    np.random.seed(20230906)
    sample_size = cts.runs
    # sample_size = 100 #define the sample size
    sampled_values = {}
    point_estimates = {}
    
    # Import Project and Calibration
    P = hbv.project(country, calibrated=True)
    pop_groups = list(P.data.pops.keys())
    
    for category, cat_parameters in parameters.items():
        for par, specs in cat_parameters.items():
            unit = specs['format']
            sampled_values[par] = {}
            point_estimates[par] = {}
            for pop in pop_groups:
                # Point estimate of parameter per pop group
                point_estimate = P.data.tdve[par].ts[pop].assumption 
                
                # Skip if point estimate is zero
                if point_estimate == 0:
                    continue
                
                point_estimates[par][pop] = point_estimate
                
                # Calculate lb, ub, std_dev
                lower_bound = specs['lower_bound'] if specs['lower_bound'] is not None else 0.7*point_estimate
                upper_bound = specs['upper_bound'] if specs['upper_bound'] is not None else (min(1, 1.3*point_estimate) if unit in ['probability', 'proportion'] else 1.3*point_estimate)
                std_dev = specs['std_dev'] if specs['std_dev'] is not None else (upper_bound-lower_bound)/4
                
                # Stochastic sampling based on the distribution
                if unit in ['probability', 'proportion']: 
                    # Truncated normal distribution
                    sampled_values[par][pop] = stats.truncnorm.rvs((lower_bound-point_estimate)/std_dev, (upper_bound-point_estimate)/std_dev, loc=point_estimate, scale=std_dev, size=sample_size)
                elif unit == 'rate': 
                    # Normal distribution
                    sampled_values[par][pop] = np.random.normal(point_estimate, std_dev, size=sample_size)
                else:
                    # Uniform distribution for numbers and others
                    sampled_values[par][pop] = np.random.uniform(lower_bound, upper_bound, size=sample_size)
    
    # Duplicate treatment and vaccine efficacy parameters for male onto females
    ages = ['0-0', '1-4', '5-9','10-19', '20-29','30-39','40-49','50-59','60-69','70-79','80-89','90+'] # Age groups
    tot_pars = parameters['trtinv'].copy()
    tot_pars.update(parameters['vacinv'])
    for age in ages:
        for par in tot_pars:
            if par in sampled_values and f'{age}M' in sampled_values[par]:
                sampled_values[par][f'{age}F'] = sampled_values[par][f'{age}M'] 
                    
    return sampled_values, point_estimates

def dict_flat(test_dict, rem_keys):
    if not isinstance(test_dict, dict):
        return test_dict
    res = {}
 
    for key, val in test_dict.items():
        rem = dict_flat(val, rem_keys)
 
        # performing removal
        if key not in rem_keys:
            res[key] = rem
        else:
            if isinstance(rem, dict):
                res.update(rem)
    return res


def sample_plots(country, samples, point_estimates, parameters):
    
    filename = hbv.root/"outputs"/"sampling"/f"{country}_stochastic samples.pdf"
    pp = PdfPages(filename)
    
    # Run stochastic sampling method
    #sampled = pars_stochastic(country)
    
    # Populations to iterate through
    P = hbv.project(country, calibrated=True)
    pop_groups = list(P.data.pops.keys())
    
    # List of parameters sampled to iterate through
    pars = []
    for category, cat_parameters in parameters.items():
        pars.append(list(cat_parameters.keys()))
    pars = [item for sublist in pars for item in sublist]
    
    parameters = dict_flat(parameters, ['nathis', 'trtinv', 'vacinv', 'calibration'])


# Plot outcomes
    for pop in pop_groups:
        
        fig1 = plt.figure(figsize=(20,20))
        
        for idx, par in enumerate(pars):
            f1 = fig1.add_subplot(int(np.ceil(len(pars)/5)),5, idx+1)
            try: 
                f1.hist(samples[par][pop], alpha=0.4)
                ymin, ymax = f1.get_ylim()
                f1.vlines(x=parameters[par]['lower_bound'], ymin=ymin, ymax=ymax, linestyle="dashed", color="red")
                f1.vlines(x=parameters[par]['upper_bound'], ymin=ymin, ymax=ymax, linestyle="dashed", color="green")
                f1.vlines(x=point_estimates[par][pop], ymin=ymin, ymax=ymax, linestyle="dashed", color="black")
                f1.set_title(par)
            except KeyError:
                pass
        
        fig1.suptitle(f"{pop} Sampling Distributions", fontsize=14, style="normal", color="blue")
        fig1.tight_layout(pad=2)
        fig1.subplots_adjust(top=0.95)
        pp.savefig(fig1)
        plt.close(fig1)# Attaches plot to PDF output
        
    pp.close()


    