# Calibration diagnostics

import atomica as at
import numpy as np
import pandas as pd
import hbv_vimc as hbv
import os


def mape_calib(country: str, cal, res_cal):
    # Initiate dictionaries

    model = {}
    data = {}
    ape = {}
    mape = {}

    outputs = ["alive", "chb_pop", "flw_hcc", "cl_acu", "cl_cir", "cl_hcc"]

    pops_in = ["year", "0-0M", "0-0F", "1-4M", "1-4F", "5-9M", "5-9F", "10-19M", "10-19F", "20-29M", "20-29F", "30-39M",
               "30-39F", "40-49M", "40-49F", "50-59M", "50-59F", "60-69M", "60-69F", "70-79M", "70-79F", "80-89M",
               "80-89F", "90+M", "90+F"]

    pops_out = ["0-0M", "0-0F", "1-4M", "1-4F", "5-9M", "5-9F", "10-19M", "10-19F", "20-29M", "20-29F", "30-39M",
                "30-39F", "40-49M", "40-49F", "50-59M", "50-59F", "60-69M", "60-69F", "70-79M", "70-79F", "80-89M",
                "80-89F", "90+M", "90+F", "total"]

    for out in outputs:
        model[out] = pd.DataFrame(columns=pops_in)
        data[out] = pd.DataFrame(columns=pops_in)

    # Extract the model data (model)

    for out in outputs:
        for pop in pops_in:
            if pop == "year":
                model[out][pop] = np.floor(at.PlotData(res_cal, out, "total", t_bins=1).series[0].tvec)
            else:
                model[out][pop] = at.PlotData(res_cal, out, pop, t_bins=1).series[0].vals

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

    # Get total population values
    for out in outputs:
        data[out]["total"] = np.sum(data[out].iloc[:, 1:], axis=1)
        model[out]["total"] = np.sum(model[out].iloc[:, 1:], axis=1)

    # Extract results to an array
    for out in outputs:
        data[out] = np.array(data[out].iloc[:, 1:])
        model[out] = np.array(model[out].iloc[:, 1:])
        ape[out] = abs((data[out] - model[out]) / (data[out])) * 1e2
        ape[out] = np.nan_to_num(ape[out], nan=0)  # absolute percentage error

    # Calculate Mean Absolute Percentage Error (w/ allowance for zero errors)
    for out in outputs:
        mape[out] = np.sum(ape[out], axis=0) * (1 / len(ape[out]))

    # Export as detailed country data to excel and store pop summary across countries to produce
    # quick glance summary sheet

    out_detailed = pd.DataFrame(index=outputs, columns=pops_out)

    for out in outputs:
        for idx, pop in enumerate(pops_out):
            out_detailed.at[out, pop] = mape[out][idx]

    out_detailed.to_excel(hbv.root / "outputs" / "calibrations" / f"{country}_cal error.xlsx")  # save

    # out_summary = [country]
    #
    # for out in outputs:
    #     out_summary.append(mape[out][-1])


def calib_qual(val):
    if val <= 10:
        return ['background-color: green']
    if val > 10 and val <= 30:
        return ['background-color: yellow']
    else:
        return ['background-color: red']


def calib_summary(settings):
    calib_sum = []
    col_names = ["country", "country_code", "population", "prevalence", "hcc_incidence",
                 "acute_deaths", "cirrhosis_deaths", "hcc_deaths"]

    for iso, cnt in settings.items():
        if os.path.exists(hbv.root / "outputs" / "calibrations" / f"{iso}_cal error.xlsx") == False:
            pass
        else:
            cnt_cal = pd.read_excel(hbv.root / "outputs" / "calibrations" / f"{iso}_cal error.xlsx")
            calib_sum.append([cnt, iso, cnt_cal.iloc[0, -1], cnt_cal.iloc[1, -1], cnt_cal.iloc[2, -1],
                              cnt_cal.iloc[3, -1], cnt_cal.iloc[4, -1], cnt_cal.iloc[5, -1]])

    calib_sum = pd.DataFrame(calib_sum, columns=col_names)

    writer = pd.ExcelWriter(hbv.root / "outputs" / "calibrations" / "cal summary.xlsx", engine='xlsxwriter')
    calib_sum.to_excel(writer, sheet_name='Sheet1')
    workbook = writer.book
    worksheet = writer.sheets['Sheet1']

    format1 = workbook.add_format({'bg_color': 'green'})
    format2 = workbook.add_format({'bg_color': 'yellow'})
    format3 = workbook.add_format({'bg_color': 'red'})

    start_row = 1
    end_row = len(calib_sum)

    start_col = 3
    end_col = 8

    worksheet.conditional_format(start_row, start_col, end_row, end_col,
                                 {'type': 'cell',
                                  'criteria': '<',
                                  'value': 10,
                                  'format': format1})

    worksheet.conditional_format(start_row, start_col, end_row, end_col,
                                 {'type': 'cell',
                                  'criteria': 'between',
                                  'minimum': 10,
                                  'maximum': 30,
                                  'format': format2})

    worksheet.conditional_format(start_row, start_col, end_row, end_col,
                                 {'type': 'cell',
                                  'criteria': '>',
                                  'value': 30,
                                  'format': format3})
    writer.save()


def prog_trends_visual(country: str):
    """ Visual check that calibrated trends are age-sex consistent"""
