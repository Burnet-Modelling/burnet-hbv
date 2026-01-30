## VIMC outputs
import pandas as pd
import atomica as at
import numpy as np
import hbv_vimc as hbv
import hbv_vimc.constants as cts 
from openpyxl import load_workbook 
from openpyxl.worksheet.table import Table, TableStyleInfo 

def input_results(country:str, samples):
    '''
    Updated: 17-01-2024
    
    Uses new stochastic sampling approach (via dictionary) to fill input
    parameter sheet
    
    More streamlined version for upload to Montagu
    '''
    
    # Generate databook to fill
    columns=[ #0-0 populations only
             "ve_eag", "ve_sag", "eag_hvl", "sag_hvl", "hvl_trans", "mtct_chronic",
             #Equal in all populations
             "dc_death", "te_cc_dc", "te_cc_hcc", "te_dc_cc", "te_dc_hcc",
             "te_icl_cc", "te_icl_hcc", "te_ict_hcc", "te_ie_cc", "te_ie_hcc","te_m_dc", "te_m_hcc",
             # Population Specific
             "0_4_it_icl","5_14_it_icl","15_49_it_icl","50+_it_icl", #IT to ICL
             "5_14_it_hcc", "15_49_it_hcc", "50+_it_hcc", #IT to HCC
             "0_4_icl_ict", "5_14_icl_ict", "15_49_icl_ict", "50+_icl_ict", #ICL to ICT
             "0_14M_icl_cc", "0_14F_icl_cc", "15_49M_icl_cc", "15_49F_icl_cc", "50+M_icl_cc", "50+F_icl_cc", #ICL to CC
             "5_14_icl_hcc", "15_49M_icl_hcc", "15_49F_icl_hcc", "50+M_icl_hcc", "50+F_icl_hcc", #ICL to HCC
             "0_14_ict_ie",  "15_49_ict_ie", "50+_ict_ie", #ICT to IE
             "0_4_ict_icl", "5_14_ict_icl", "15_49_ict_icl", "50+_ict_icl",#ICT to ICL
             "5_14_ict_hcc", "15_49M_ict_hcc", "15_49F_ict_hcc", "50_69M_ict_hcc", "50_69F_ict_hcc", "70+M_ict_hcc", "70+F_ict_hcc",#ICT to HCC
             "0_14M_ie_cc", "0_14F_ie_cc", "15_49M_ie_cc", "15_49F_ie_cc","50_69M_ie_cc", "50_69F_ie_cc","70+M_ie_cc", "70+F_ie_cc", #IE to CC
             "5_14_ie_hcc", "15_49M_ie_hcc","15_49F_ie_hcc","50_69M_ie_hcc","50_69F_ie_hcc","70+M_ie_hcc","70+F_ie_hcc", #IE to HCC
             "5_14_cc_hcc", "15_49M_cc_hcc","15_49F_cc_hcc","50_69M_cc_hcc","50_69F_cc_hcc","70+M_cc_hcc","70+F_cc_hcc",   # CC to HCC
             "5+_dc_hcc", # DC to HCC
             "5+_hcc_death"] #HCC death
    
    in_df = pd.DataFrame(columns =columns)
    runs=cts.runs
    in_df["run_id"] = np.arange(1,runs+1, 1)
    first_column = in_df.pop('run_id') 
    in_df.insert(0, 'run_id', first_column) 
    
    # Import project for values
    P = hbv.project(country, calibrated=True)

    # Fill databook 
    for run in range (runs):
        # Mother to Child Transmission
        in_df["ve_eag"][run] = samples["eag_ve"]["0-0M"][run]
        in_df["ve_sag"][run] = samples["sag_ve"]["0-0M"][run]
        in_df["eag_hvl"][run] = samples["eag_hvl"]["0-0M"][run]
        in_df["sag_hvl"][run] = samples["sag_hvl"]["0-0M"][run]
        in_df["hvl_trans"][run] = samples["hvl_trisk"]["0-0M"][run]
        in_df["mtct_chronic"][run] = samples["ci_p"]["0-0M"][run]
        # Equal in all populations (note: all effectiveness now expressed in same way)
        in_df["dc_death"][run] = samples["m_dc"]["0-0M"][run]
        in_df["te_cc_dc"][run] = 1- samples["te_cc_dc"]["0-0M"][run]
        in_df["te_cc_hcc"][run] =1- samples["te_cc_hcc"]["0-0M"][run]
        in_df["te_dc_cc"][run] = 1-samples["te_dc_cc"]["0-0M"][run]
        in_df["te_dc_hcc"][run] = 1-samples["te_dc_hcc"]["0-0M"][run]
        in_df["te_icl_cc"][run] = 1-samples["te_icl_cc"]["0-0M"][run]
        in_df["te_icl_hcc"][run] = 1-samples["te_icl_hcc"]["0-0M"][run]
        in_df["te_ict_hcc"][run] = 1-samples["te_ict_hcc"]["0-0M"][run]
        in_df["te_ie_cc"][run] = 1-samples["te_ie_cc"]["0-0M"][run]
        in_df["te_ie_hcc"][run] = 1-samples["te_ie_hcc"]["0-0M"][run]
        in_df["te_m_dc"][run] = 1-samples["te_m_dc"]["0-0M"][run]
        in_df["te_m_hcc"][run] = 1-samples["te_m_hcc"]["0-0M"][run]
        # Population Specific (should all be +/- 30% point estimate except hcc_death)
        in_df["5+_dc_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption
        in_df["5+_hcc_death"][run] = samples["m_hcc"]["0-0M"][run]
        
        in_df["0_4_it_icl"][run] = samples["eag_020_m"]["0-0M"][run]*P.data.tdve["it_icl_in"].ts["0-0M"].assumption
        in_df["5_14_it_icl"][run] = samples["eag_020_m"]["0-0M"][run]*P.data.tdve["it_icl_in"].ts["5-9M"].assumption
        in_df["15_49_it_icl"][run] = samples["eag_020_m"]["0-0M"][run]*P.data.tdve["it_icl_in"].ts["20-29M"].assumption
        in_df["50+_it_icl"][run] = samples["eag_020_m"]["0-0M"][run]*P.data.tdve["it_icl_in"].ts["50-59M"].assumption
        
        in_df["5_14_it_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_it_hcc"].ts["5-9M"].assumption
        in_df["15_49_it_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_it_hcc"].ts["20-29M"].assumption
        in_df["50+_it_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_it_hcc"].ts["50-59M"].assumption
        
        in_df["0_4_icl_ict"][run] = samples["eag_020_m"]["0-0M"][run]*P.data.tdve["icl_ict_in"].ts["0-0M"].assumption
        in_df["5_14_icl_ict"][run] = samples["eag_020_m"]["0-0M"][run]*P.data.tdve["icl_ict_in"].ts["5-9M"].assumption
        in_df["15_49_icl_ict"][run] = samples["eag_020_m"]["0-0M"][run]*P.data.tdve["icl_ict_in"].ts["20-29M"].assumption
        in_df["50+_icl_ict"][run] = samples["eag_020_m"]["0-0M"][run]*P.data.tdve["icl_ict_in"].ts["50-59M"].assumption
        
        in_df["0_14M_icl_cc"][run] = samples["cc_020_m"]["0-0M"][run]*P.data.tdve["ie_cc_in"].ts["0-0M"].assumption*P.data.tdve["rr_icl_cc"].ts["0-0M"].assumption
        in_df["0_14F_icl_cc"][run] = samples["cc_020_m"]["0-0M"][run]*P.data.tdve["ie_cc_in"].ts["0-0F"].assumption*P.data.tdve["rr_icl_cc"].ts["0-0F"].assumption
        in_df["15_49M_icl_cc"][run] = samples["cc_020_m"]["0-0M"][run]*P.data.tdve["ie_cc_in"].ts["20-29M"].assumption*P.data.tdve["rr_icl_cc"].ts["20-29M"].assumption
        in_df["15_49F_icl_cc"][run] = samples["cc_020_m"]["0-0M"][run]*P.data.tdve["ie_cc_in"].ts["20-29F"].assumption*P.data.tdve["rr_icl_cc"].ts["20-29F"].assumption
        in_df["50+M_icl_cc"][run] = samples["cc_020_m"]["0-0M"][run]*P.data.tdve["ie_cc_in"].ts["50-59M"].assumption*P.data.tdve["rr_icl_cc"].ts["50-59M"].assumption
        in_df["50+F_icl_cc"][run] = samples["cc_020_m"]["0-0M"][run]*P.data.tdve["ie_cc_in"].ts["50-59F"].assumption*P.data.tdve["rr_icl_cc"].ts["50-59F"].assumption
        
        in_df["5_14_icl_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_icl_hcc"].ts["5-9M"].assumption
        in_df["15_49M_icl_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_icl_hcc"].ts["20-29M"].assumption
        in_df["15_49F_icl_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_icl_hcc"].ts["20-29F"].assumption
        in_df["50+M_icl_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_icl_hcc"].ts["50-59M"].assumption
        in_df["50+F_icl_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_icl_hcc"].ts["50-59F"].assumption
        
        in_df["0_14_ict_ie"][run] = samples["sag_020_m"]["0-0M"][run]*P.data.tdve["ict_ie_in"].ts["0-0M"].assumption
        in_df["15_49_ict_ie"][run] = samples["sag_020_m"]["0-0M"][run]*P.data.tdve["ict_ie_in"].ts["20-29M"].assumption
        in_df["50+_ict_ie"][run] = samples["sag_020_m"]["0-0M"][run]*P.data.tdve["ict_ie_in"].ts["50-59M"].assumption
        
        in_df["0_4_ict_icl"][run] = samples["sag_020_m"]["0-0M"][run]*P.data.tdve["it_icl_in"].ts["0-0M"].assumption*P.data.tdve["rr_ict_icl"].ts["0-0M"].assumption
        in_df["5_14_ict_icl"][run] = samples["sag_020_m"]["0-0M"][run]*P.data.tdve["it_icl_in"].ts["5-9M"].assumption*P.data.tdve["rr_ict_icl"].ts["5-9M"].assumption
        in_df["15_49_ict_icl"][run] = samples["sag_020_m"]["0-0M"][run]*P.data.tdve["it_icl_in"].ts["20-29M"].assumption*P.data.tdve["rr_ict_icl"].ts["20-29M"].assumption
        in_df["50+_ict_icl"][run] = samples["sag_020_m"]["0-0M"][run]*P.data.tdve["it_icl_in"].ts["50-59M"].assumption*P.data.tdve["rr_ict_icl"].ts["50-59M"].assumption
        
        in_df["5_14_ict_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_ict_hcc"].ts["5-9M"].assumption
        in_df["15_49M_ict_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_ict_hcc"].ts["20-29M"].assumption
        in_df["15_49F_ict_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_ict_hcc"].ts["20-29F"].assumption
        in_df["50_69M_ict_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_ict_hcc"].ts["50-59M"].assumption
        in_df["50_69F_ict_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_ict_hcc"].ts["50-59F"].assumption
        in_df["70+M_ict_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_ict_hcc"].ts["70-79M"].assumption
        in_df["70+F_ict_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_ict_hcc"].ts["70-79F"].assumption
        
        in_df["0_14M_ie_cc"][run] = samples["cc_020_m"]["0-0M"][run]*P.data.tdve["ie_cc_in"].ts["0-0M"].assumption
        in_df["0_14F_ie_cc"][run] = samples["cc_020_m"]["0-0M"][run]*P.data.tdve["ie_cc_in"].ts["0-0F"].assumption
        in_df["15_49M_ie_cc"][run] = samples["cc_020_m"]["0-0M"][run]*P.data.tdve["ie_cc_in"].ts["20-29M"].assumption
        in_df[ "15_49F_ie_cc"][run] = samples["cc_020_m"]["0-0M"][run]*P.data.tdve["ie_cc_in"].ts["20-29F"].assumption
        in_df["50_69M_ie_cc"][run] = samples["cc_020_m"]["0-0M"][run]*P.data.tdve["ie_cc_in"].ts["50-59M"].assumption
        in_df["50_69F_ie_cc"][run] = samples["cc_020_m"]["0-0M"][run]*P.data.tdve["ie_cc_in"].ts["50-59F"].assumption
        in_df["70+M_ie_cc"][run] = samples["cc_020_m"]["0-0M"][run]*P.data.tdve["ie_cc_in"].ts["70-79M"].assumption
        in_df["70+F_ie_cc"][run] = samples["cc_020_m"]["0-0M"][run]*P.data.tdve["ie_cc_in"].ts["70-79F"].assumption
        
        in_df["5_14_ie_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_ie_hcc"].ts["5-9M"].assumption
        in_df["15_49M_ie_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_ie_hcc"].ts["20-29M"].assumption
        in_df["15_49F_ie_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_ie_hcc"].ts["20-29F"].assumption
        in_df["50_69M_ie_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_ie_hcc"].ts["50-59M"].assumption
        in_df["50_69F_ie_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_ie_hcc"].ts["50-59F"].assumption
        in_df["70+M_ie_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_ie_hcc"].ts["70-79M"].assumption
        in_df["70+F_ie_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_ie_hcc"].ts["70-79F"].assumption
        
        in_df["5_14_cc_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_cc_hcc"].ts["5-9M"].assumption
        in_df["15_49M_cc_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_cc_hcc"].ts["20-29M"].assumption
        in_df["15_49F_cc_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_cc_hcc"].ts["20-29F"].assumption
        in_df["50_69M_cc_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_cc_hcc"].ts["50-59M"].assumption
        in_df["50_69F_cc_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_cc_hcc"].ts["50-59F"].assumption
        in_df["70+M_cc_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_cc_hcc"].ts["70-79M"].assumption
        in_df["70+F_cc_hcc"][run] = samples["hcc_020_m"]["0-0M"][run]*P.data.tdve["dc_hcc_in"].ts["50-59M"].assumption*P.data.tdve["rr_cc_hcc"].ts["70-79F"].assumption
        
    
    in_df.to_csv(hbv.root/"outputs"/"vimc_outs"/f"{country}_inputs.csv", index=False)


def central_results(country:str, central):
    '''
    Created by: Phillip Luong (https://github.com/phillipluong)
    Last Updated: 31/01/23
    Generates a DataFrame listing all central output results.

    Inputs:
    - fw:       Excel sheet of the Atomica Framework of the model (str)
    - db:       Excel sheet of the Atomica Databook of the model (str)
    - cl:       Excel sheet of the Model Calibration of the model - specific to a region estimate (str)

    Outputs:
    - cen_df:   DataFrame listing the central estimates of 1-year age groups between 2000-2100 (DataFrame)
    '''
    s_yr = hbv.s_yr
    e_yr = hbv.e_yr
    scens = hbv.get_scenarios(country)

    ## Load and process the input data
    df1 = pd.read_csv(hbv.root/"inputs"/f'{country}'/f'{country}_pop.csv') # This will likely change by country
    df1 = df1[(df1.year >= 1990) & (df1.year <=2100)]
    df1 = df1.groupby(by=["country_code_numeric", "country_code", "country", "age_from", "age_to", "year"], as_index=False).sum()
    df1 = df1.drop('gender',axis=1)

    #add ability to put in full country name
    age_groups = ['0-0', '1-4', '5-9','10-19', '20-29','30-39','40-49','50-59','60-69','70-79','80-89','90+']

    ## Assort data by age group
    for idx,row in df1.iterrows():
        if row.age_from < 1:
            df1.loc[idx,'age_group'] = '0-0'
        elif row.age_from < 5:
            df1.loc[idx,'age_group'] = '1-4'
        elif row.age_from < 10:
            df1.loc[idx,'age_group'] = '5-9'
        elif row.age_from < 20:
            df1.loc[idx,'age_group'] = '10-19'
        elif row.age_from < 30:
            df1.loc[idx,'age_group'] = '20-29'
        elif row.age_from < 40:
            df1.loc[idx,'age_group'] = '30-39'
        elif row.age_from < 50:
            df1.loc[idx,'age_group'] = '40-49'
        elif row.age_from < 60:
            df1.loc[idx,'age_group'] = '50-59'
        elif row.age_from < 70:
            df1.loc[idx,'age_group'] = '60-69'
        elif row.age_from < 80:
            df1.loc[idx,'age_group'] = '70-79'
        elif row.age_from < 90:
             df1.loc[idx,'age_group'] = '80-89'
        else:
            df1.loc[idx,'age_group'] = '90+'

    ## Generate Population Distribution Weights (this could be done externally and imported if needed)
    df2 = df1.groupby(by=["country_code_numeric", "country_code", "country", "age_group", "year"],as_index=False).sum()
    df2 = df2.drop(columns=["age_from", "age_to"], axis=1)
    df2 = df2.rename(columns={"value":"total"})
        
    df1 = pd.merge(df1, df2)
    df1["wt"] = df1["value"]/df1["total"]
    
    ## Extract modelled data
    col_names = ["year"]+age_groups
    outputs = ["alive", "tot_inc", ":dd_hbv", "dalys", "yll"]
    out_name = ["cohort_size", "cases", "deaths", "dalys", "yll"]
    
    data = {}
    for scen in scens.values():
        data[scen.name] = {}
        #Initiate each data frame
        for out in outputs:
            data[scen.name][out]=pd.DataFrame(columns=col_names)
            data[scen.name][out].year=np.arange(s_yr, e_yr, 1)
            # Fill with data from the model results
            for idx,age in enumerate(age_groups):
                data[scen.name][out].iloc[:,idx+1]=at.PlotData(central[scen.name], out, pops={"pop":[age+"M", age+"F"]}, pop_aggregation="sum", t_bins=1).series[0].vals

    ## Reshape data 
    for scen in scens.values():
        for idx, out in enumerate(outputs):
            data[scen.name][out]=pd.melt(data[scen.name][out],id_vars="year", value_vars=age_groups)
            data[scen.name][out]=data[scen.name][out].rename(columns={"variable":"age_group", "value":out_name[idx]})

    ## Merge Data, distribute across age groups, and format for csv output
    for scen in scens.values():
        data[scen.name]["df1"] = df1
        for idx,out in enumerate(outputs):
            data[scen.name]["df1"] = pd.merge(data[scen.name]["df1"],data[scen.name][out])
            data[scen.name]["df1"][out_name[idx]] = data[scen.name]["df1"][out_name[idx]]*data[scen.name]["df1"]["wt"]
            
    for scen in scens.values():
        data[scen.name]["df1"]["disease"] = "HepB"
        data[scen.name]["df1"] = data[scen.name]["df1"].rename(columns={"age_from":"age", "country":"country_name", "country_code":"country"})
        # data[scen.name]["df1"]=data[scen.name]["df1"].iloc[:,np.r_[[14,5,3,1,2,10,11,13,12 ]]]
        data[scen.name]["df1"] = data[scen.name]["df1"][['disease','year','age','country','country_name','cohort_size','cases','dalys','deaths','yll']]
        data[scen.name]["df1"] = data[scen.name]["df1"][(data[scen.name]["df1"]["year"]>=2000)&(data[scen.name]["df1"]["year"]<=2100)]
        data[scen.name]["df1"] = data[scen.name]["df1"].sort_values(by=["age", "year"])
            
    
    ## Export to csv
    for scen in scens.values():
        data[scen.name]["df1"].to_csv(hbv.root/"outputs"/"vimc_outs"/f"{country}_{scen.name}_central.csv", index=False)



def stochastic_results(country:str, stoch_out):
    '''
    Created by: Phillip Luong (https://github.com/phillipluong)
    Last Updated: 31/01/23
    Generates a DataFrame listing all stochastic output results. Each stochastic output is indexed by a 'run_id' which can also refer to the stochastic input parameter set.

    Inputs:
    - n_samples:    Number of stochastic input samples (int)
    - fw:           Excel sheet of the Atomica Framework of the model (str)
    - db:           Excel sheet of the Atomica Databook of the model (str)
    - cl:           Excel sheet of the Model Calibration of the model - specific to a region estimate (str)
    - seed:         Random seed (int)

    Outputs:
    - final_df:     DataFrame listing all stochastic estimates of 1-year age groups between 2000-2100 (DataFrame)
    '''
    
    s_yr = hbv.s_yr
    e_yr = hbv.e_yr
    runs = hbv.runs
    scens = hbv.get_scenarios(country)
    out_name = ["cohort_size", "cases", "deaths", "dalys", "yll"]


    ## Load and process the input data
    df1 = pd.read_csv(hbv.root/"inputs"/f'{country}'/f'{country}_pop.csv') # This will likely change by country
    df1 = df1[(df1.year >= 1990) & (df1.year <=2100)]
    df1=df1.groupby(by=["country_code_numeric", "country_code", "country", "age_from", "age_to", "year"], as_index=False).sum()
    df1 = df1.drop('gender',axis=1)
    
    #add ability to put in full country name
    age_groups = ['0-0', '1-4', '5-9','10-19', '20-29','30-39','40-49','50-59','60-69','70-79','80-89','90+']

    ## Assort data by age group
    for idx,row in df1.iterrows():
        if row.age_from < 1:
            df1.loc[idx,'age_group'] = '0-0'
        elif row.age_from < 5:
            df1.loc[idx,'age_group'] = '1-4'
        elif row.age_from < 10:
            df1.loc[idx,'age_group'] = '5-9'
        elif row.age_from < 20:
            df1.loc[idx,'age_group'] = '10-19'
        elif row.age_from < 30:
            df1.loc[idx,'age_group'] = '20-29'
        elif row.age_from < 40:
            df1.loc[idx,'age_group'] = '30-39'
        elif row.age_from < 50:
            df1.loc[idx,'age_group'] = '40-49'
        elif row.age_from < 60:
            df1.loc[idx,'age_group'] = '50-59'
        elif row.age_from < 70:
            df1.loc[idx,'age_group'] = '60-69'
        elif row.age_from < 80:
            df1.loc[idx,'age_group'] = '70-79'
        elif row.age_from < 90:
             df1.loc[idx,'age_group'] = '80-89'
        else:
            df1.loc[idx,'age_group'] = '90+'

    ## Generate Population Distribution Weights (this could be done externally and imported if needed)
    df2=df1.groupby(by=["country_code_numeric", "country_code", "country", "age_group", "year"],as_index=False).sum()
    df2=df2.drop(columns=["age_from", "age_to"], axis=1)
    df2=df2.rename(columns={"value":"total"})
        
    df1=pd.merge(df1, df2)
    df1["wt"]=df1["value"]/df1["total"]
     
    ## Reshape Data
    for scen in scens.values():
        for idx, out in enumerate(out_name):
            stoch_out[scen.name][out] = pd.melt(stoch_out[scen.name][out],id_vars=["year", "run_id"], value_vars=age_groups)
            stoch_out[scen.name][out] = stoch_out[scen.name][out].rename(columns={"variable":"age_group", "value":out_name[idx]})
    
    ## Merge Data, distribute across age groups, and format for csv output
    for scen in scens.values():
        stoch_out[scen.name]["df1"] = df1
        for idx,out in enumerate(out_name):
            stoch_out[scen.name]["df1"] = pd.merge(stoch_out[scen.name]["df1"],stoch_out[scen.name][out])
            stoch_out[scen.name]["df1"][out_name[idx]] = stoch_out[scen.name]["df1"][out_name[idx]]*stoch_out[scen.name]["df1"]["wt"]
    
    for scen in scens.values():
        stoch_out[scen.name]["df1"]["disease"] = "HepB"
        stoch_out[scen.name]["df1"] = stoch_out[scen.name]["df1"].iloc[:,np.r_[[16,10,5,3,1,2,11,12,13,14,15]]]
        #stoch_out[scen.name]["df1"] = stoch_out[scen.name]["df1"][['disease','year','age','country','country_name','cohort_size','cases','dalys','deaths','yll']]
        stoch_out[scen.name]["df1"] = stoch_out[scen.name]["df1"][(stoch_out[scen.name]["df1"]["year"] >= 2000) & (stoch_out[scen.name]["df1"]["year"] <= 2100)]
        stoch_out[scen.name]["df1"] = stoch_out[scen.name]["df1"].rename(columns={"age_from":"age", "country":"country_name", "country_code":"country"})
        stoch_out[scen.name]["df1"] = stoch_out[scen.name]["df1"].sort_values(by=["run_id", "age", "year"])

        # Round values to 0 d.p.
        stoch_out[scen.name]["df1"]["cohort_size"] = stoch_out[scen.name]["df1"]["cohort_size"].astype(float).round(0)
        stoch_out[scen.name]["df1"]["cases"] = stoch_out[scen.name]["df1"]["cases"].astype(float).round(0)
        stoch_out[scen.name]["df1"]["dalys"] = stoch_out[scen.name]["df1"]["dalys"].astype(float).round(0)
        stoch_out[scen.name]["df1"]["deaths"] = stoch_out[scen.name]["df1"]["deaths"].astype(float).round(0)
        stoch_out[scen.name]["df1"]["yll"] = stoch_out[scen.name]["df1"]["yll"].astype(float).round(0)

    ## Export to csv
    for scen in scens.values():
        stoch_out[scen.name]["df1"].to_csv(hbv.root/"outputs"/"vimc_outs"/f"{country}_{scen.name}_stochastic.csv", index=False)


def combined_results(countries):  

    scenarios = ['Baseline', 'Baseline No BD', 'IA', 'IA No BD', 'Bluesky', 'Bluesky No BD', 'No Vaccination'] 
    country_to_region = cts.country_to_region
    
    df = pd.DataFrame(columns=['country','WHO region','cases','deaths','scenario','year'])
    for scenario in scenarios:  
        for country in countries: 
            current_estimates_all = pd.read_csv(hbv.root/"outputs"/"vimc_outs"/f"{country}_{scenario}_central.csv") 
            current_estimates = current_estimates_all.groupby(['year'])[['cases','deaths']].sum().reset_index()
            current_estimates['scenario'] = scenario
            current_estimates['WHO region'] = country_to_region[country]
            current_estimates['country'] = country
            df = pd.concat([df, current_estimates], ignore_index=True)
    
    filename = "combined_results.xlsx"
    filepath = hbv.root/"outputs"/"vimc_outs"/f"{filename}"
    df.to_excel(filepath, index=False) 
    
    wb = load_workbook(filepath) 
    ws = wb.active 
    
    tab = Table(displayName="Table1", ref=ws.dimensions) 
    style = TableStyleInfo(name="TableStyleLight1", showFirstColumn=False, 
                           showLastColumn=False, showRowStripes=True, showColumnStripes=False) 
    tab.tableStyleInfo = style 
    
    ws.add_table(tab) 
    wb.save(filepath) 
    
    print(f"Spreadsheet saved: {filepath}")
    

def central_results_combined(countries):  

    scenarios =  {'Baseline': 'hepb-hepb3-bd-default', 
                 'Baseline No BD': 'hepb-hepb3-default', 
                 'IA': 'hepb-hepb3-bd-ia2030', 
                 'IA No BD': 'hepb-hepb3-ia2030', 
                 'Bluesky': 'hepb-hepb3-bd-bluesky', 
                 'Bluesky No BD': 'hepb-hepb3-bluesky', 
                 'No Vaccination': 'hepb-no-vaccination'}
    
    for scen, scen_id in scenarios.items():  
        df = pd.DataFrame(columns=['disease','year','age','country','country_name','cohort_size','cases','dalys','deaths','yll'])
        for country in countries: 
            central_burden = pd.read_csv(hbv.root/"outputs"/"vimc_outs"/f"{country}_{scen}_central.csv") 
            # central_burden = pd.read_csv(f"C:/Users/farah.houdroge/Burnet Institute/WG-Modelling-Hepatitis B - Documents/Calibrations/2023-12-01/vimc_outs/{country}_{scen}_central.csv") 
            df = pd.concat([df, central_burden], ignore_index=True)
        df=df.sort_values(by=["age","country","year"])
        df['cohort_size'] =  df['cohort_size'].round(2)
        df['cases'] =  df['cases'].round(2)
        df['dalys'] =  df['dalys'].round(2)
        df['deaths'] =  df['deaths'].round(2)
        df['yll'] =  df['yll'].round(2)
        filename = f"central-burden-{scen_id}.csv"
        filepath = hbv.root/"outputs"/"vimc_outs"/f"{filename}"
        df.to_csv(filepath, index=False) 
        
        print(f"Spreadsheet saved: {filepath}")
        
        
def add_missing_countries(countries):  
    
    scen_ids = ['hepb-hepb3-bd-default', 'hepb-hepb3-default',  'hepb-hepb3-ia2030',  'hepb-hepb3-bd-ia2030', 
                 'hepb-hepb3-bd-bluesky', 'hepb-hepb3-bluesky', 'hepb-no-vaccination']
    
    countries_missing = {'DMA': 'Dominica',
                         'GRD': 'Grenada',
                         'LCA' : 'Saint Lucia',
                         'MDV': 'Maldives',
                         'MHL': 'Marshall Islands',
                         'NAM': 'Namibia',
                         'PSE': 'Palestine, State of',
                         'TUV': 'Tuvalu',
                         'VCT': 'Saint Vincent and the Grenadines',
                         'XK': 'Kosovo'}
    columns_to_zero = ['cohort_size','cases','dalys','deaths','yll']
    
    for scen_id in scen_ids: 
        filename = f"central-burden-{scen_id}.csv"
        filepath = hbv.root/"outputs"/"vimc_outs"/f"{filename}"
        df = pd.read_csv(filepath)
        for country, country_name in countries_missing.items():
            rows = df[df['country']=='EGY'].copy()
            rows['country'] = country
            rows['country_name'] = country_name
            rows[columns_to_zero] = 0
            df = pd.concat([df, rows], ignore_index=True)
        df=df.sort_values(by=["age","country","year"])
        filename = f"central-burden-{scen_id}_2.csv"
        filepath = hbv.root/"outputs"/"vimc_outs"/f"{filename}"
        df.to_csv(filepath, index=False) 
        
        print(f"Spreadsheet saved: {filepath}")