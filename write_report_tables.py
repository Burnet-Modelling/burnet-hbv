"""
Build the HBV rows of Tables 3-5 from "Cost and impact of VH actions_DRAFT report 0.1.docx"
(Viral hep modelling/Strategy implementation), as an Excel workbook, from the Status Quo vs
Action Scenario results already saved by run_scenarios.py.

Table 4 (direct/societal costs, QALYs, cost-effectiveness) needs cost, wage and
health-utility input data that doesn't exist for HBV in this repo yet - see
hcv-atomica/utils.py's write_societal_costs and write_premature_death_costs for a reusable
human-capital-approach methodology to adapt once that data is sourced. Its rows are left
"N/A" with a note rather than guessed.

The "Number of new hepatitis B infections" row (and the incidence-reduction row in Table 5)
depends on how well-constrained foi_cal is - see run_scenarios.py's incidence panel, which
exists to sanity-check this figure against plausible Australian HBV incidence before it's
used in any report.

Usage:
    python write_report_tables.py
"""
import os
import numpy as np
import sciris as sc
import atomica as at
from isort.wrap_modes import vertical_prefix_from_module_import

from hbv_aus.utils import _get_github_folder
from run_scenarios import _series_vals, _cascade_vals, SCENARIOS

BASELINE_YEAR = 2015
YEAR_START = 2027
YEAR_END = 2036
TARGET_YEAR = 2030

CSQ, ACT = "Status Quo", "Action Scenario"
NA = "N/A"

COST_NOTE = "N/A - requires HBV cost/wage/utility input data not yet available in this repository"


def _flow_sum(res, param_name, years):
    """Sum a rate/probability-formatted parameter's actual flow (people/year, not the raw
    rate) across every population and every transition it drives, for the given years."""
    total = 0.0
    for pop in res.model.pops:
        for link in pop.links:
            if link.parameter.name != param_name:
                continue
            t = list(link.t)
            for y in years:
                if y in t:
                    total += link.vals[t.index(y)]
    return total


def _value_at_year(tvec, vals, year):
    tvec = list(tvec)
    return vals[tvec.index(year + 0.5)]


def _new_infections(res, pop_names, years):
    """Horizontal (community) + mother-to-child transmission new infections, summed over
    `years`. See module docstring - sanity-check against run_scenarios.py's incidence
    panel before use."""
    horiz = _flow_sum(res, "horiz_chb", years)
    b_chb_total = 0.0
    for p in pop_names:
        s = at.PlotData(res, outputs="b_chb", pops=p, t_bins=1).series[0]
        tvec = list(s.tvec)
        for y in years:
            if y + 0.5 in tvec:
                b_chb_total += s.vals[tvec.index(y + 0.5)]
    return horiz + b_chb_total


def _new_infections_at_year(res, pop_names, year):
    return _new_infections(res, pop_names, [year])


def build_tables(results):
    project, _, _ = results[SCENARIOS[0][0]]
    pop_names = list(project.data.pops.keys())

    res = {label: results[key][2] for key, label in SCENARIOS}

    # ---- Table 3a: Service delivery and epi, 2027-2030 -------------------------------
    years_range = list(range(YEAR_START, 2030 + 1))
    treatments = {lbl: _flow_sum(res[lbl], "treat_rate", years_range) for _, lbl in SCENARIOS}
    new_infections = {lbl: _new_infections(res[lbl], pop_names, years_range) for _, lbl in SCENARIOS}
    deaths_c_2030 = {}
    hcc_c_2030 = {}
    hepb_pop_2030 = {}
    prev_2030 = {}
    death_i_2030 = {}
    hcc_i_2030 = {}
    for _, lbl in SCENARIOS:
        t_dth, v_dth = _series_vals(res[lbl], pop_names, None, "hep_dth")
        deaths_c_2030[lbl] = sum(v_dth[list(t_dth).index(y + 0.5)] for y in years_range)
        t_hcc, v_hcc = _series_vals(res[lbl], pop_names, None, "inc_hcc")
        hcc_c_2030[lbl] = sum(v_hcc[list(t_hcc).index(y + 0.5)] for y in years_range)
        t_hepb, v_hepb = _series_vals(res[lbl], pop_names, None, "hepb_pop")
        hepb_pop_2030[lbl] = _value_at_year(t_hepb, v_hepb, 2030)
        t_tot, v_tot = _series_vals(res[lbl], pop_names, None, "total_pop")
        prev_2030[lbl] = hepb_pop_2030[lbl] / _value_at_year(t_tot, v_tot, 2030)
        t_hc, v_hc = _series_vals(res[lbl], pop_names, None, "inc_hcc")
        hcc_i_2030[lbl] = _value_at_year(t_hc, v_hc, 2030)
        t_dt, v_dt = _series_vals(res[lbl], pop_names, None, "hep_dth")
        death_i_2030[lbl] = _value_at_year(t_dt, v_dt, 2030)

    def d(dct):
        return dct[ACT] - dct[CSQ]

    table3a = [
        ("Service delivery", None, None, None, "header"),
        (f"Total number of hepatitis B tests", NA, NA, NA, "text"),
        (f"Total number of hepatitis B treatments", treatments[CSQ], treatments[ACT], d(treatments), "count"),
        ("Epidemiology", None, None, None, "header"),
        (f"People with hepatitis B in {2030}", hepb_pop_2030[CSQ], hepb_pop_2030[ACT], d(hepb_pop_2030), "count"),
        (f"Number of new hepatitis B infections, {YEAR_START}-{2030} (see note below table)",
         new_infections[CSQ], new_infections[ACT], d(new_infections), "count"),
        (f"Number of hepatitis B-related deaths, {YEAR_START}-{2030}", deaths_c_2030[CSQ], deaths_c_2030[ACT], d(deaths_c_2030), "count"),
        (f"Number of hepatitis B-related liver cancers, {YEAR_START}-{2030}", hcc_c_2030[CSQ], hcc_c_2030[ACT], d(hcc_c_2030),
         "count"),
        (f"Hepatitis B prevalence among the whole population in {2030} (%)",
         prev_2030[CSQ], prev_2030[ACT], d(prev_2030), "pct"),
        (f"HCC incidence in single year of {2030} ",
         hcc_i_2030[CSQ], hcc_i_2030[ACT], d(hcc_i_2030), "count"),
        (f"HBV deaths in single year of {2030} ",
         death_i_2030[CSQ], death_i_2030[ACT], d(death_i_2030), "count"),
    ]

    # ---- Table 3b: Service delivery and epi, 2027-2036 -------------------------------
    years_range = list(range(YEAR_START, 2036 + 1))
    treatments = {lbl: _flow_sum(res[lbl], "treat_rate", years_range) for _, lbl in SCENARIOS}
    new_infections = {lbl: _new_infections(res[lbl], pop_names, years_range) for _, lbl in SCENARIOS}
    deaths_c_2036 = {}
    hcc_c_2036 = {}
    hepb_pop_2036 = {}
    prev_2036 = {}
    death_i_2036 = {}
    hcc_i_2036 = {}
    for _, lbl in SCENARIOS:
        t_dth, v_dth = _series_vals(res[lbl], pop_names, None, "hep_dth")
        deaths_c_2036[lbl] = sum(v_dth[list(t_dth).index(y + 0.5)] for y in years_range)
        t_hcc, v_hcc = _series_vals(res[lbl], pop_names, None, "inc_hcc")
        hcc_c_2036[lbl] = sum(v_hcc[list(t_hcc).index(y + 0.5)] for y in years_range)
        t_hepb, v_hepb = _series_vals(res[lbl], pop_names, None, "hepb_pop")
        hepb_pop_2036[lbl] = _value_at_year(t_hepb, v_hepb, 2036)
        t_tot, v_tot = _series_vals(res[lbl], pop_names, None, "total_pop")
        prev_2036[lbl] = hepb_pop_2036[lbl] / _value_at_year(t_tot, v_tot, 2036)
        t_hc, v_hc = _series_vals(res[lbl], pop_names, None, "inc_hcc")
        hcc_i_2036[lbl] = _value_at_year(t_hc, v_hc, 2036)
        t_dt, v_dt = _series_vals(res[lbl], pop_names, None, "hep_dth")
        death_i_2036[lbl] = _value_at_year(t_dt, v_dt, 2036)

    table3b = [
        ("Service delivery", None, None, None, "header"),
        (f"Total number of hepatitis B tests", NA, NA, NA, "text"),
        (f"Total number of hepatitis B treatments", treatments[CSQ], treatments[ACT], d(treatments), "count"),
        ("Epidemiology", None, None, None, "header"),
        (f"People with hepatitis B in {2036}", hepb_pop_2036[CSQ], hepb_pop_2036[ACT], d(hepb_pop_2036), "count"),
        (f"Number of new hepatitis B infections, {YEAR_START}-{2036} (see note below table)",
         new_infections[CSQ], new_infections[ACT], d(new_infections), "count"),
        (f"Number of hepatitis B-related deaths, {YEAR_START}-{2036}", deaths_c_2036[CSQ], deaths_c_2036[ACT], d(deaths_c_2036), "count"),
        (f"Number of hepatitis B-related liver cancers, {YEAR_START}-{2036}", hcc_c_2036[CSQ], hcc_c_2036[ACT], d(hcc_c_2036),
         "count"),
        (f"Hepatitis B prevalence among the whole population in {2036} (%)",
         prev_2036[CSQ], prev_2036[ACT], d(prev_2036), "pct"),
        (f"HCC incidence in single year of {2036} ",
         hcc_i_2036[CSQ], hcc_i_2036[ACT], d(hcc_i_2036), "count"),
        (f"HBV deaths in single year of {2036} ",
         death_i_2036[CSQ], death_i_2036[ACT], d(death_i_2036), "count"),
    ]

    # ---- Table 3c: Service delivery and epi, 2027-2040 -------------------------------
    years_range = list(range(YEAR_START, 2040 + 1))
    treatments = {lbl: _flow_sum(res[lbl], "treat_rate", years_range) for _, lbl in SCENARIOS}
    new_infections = {lbl: _new_infections(res[lbl], pop_names, years_range) for _, lbl in SCENARIOS}
    deaths_c_2040 = {}
    hcc_c_2040 = {}
    hepb_pop_2040 = {}
    prev_2040 = {}
    death_i_2040 = {}
    hcc_i_2040 = {}
    for _, lbl in SCENARIOS:
        t_dth, v_dth = _series_vals(res[lbl], pop_names, None, "hep_dth")
        deaths_c_2040[lbl] = sum(v_dth[list(t_dth).index(y + 0.5)] for y in years_range)
        t_hcc, v_hcc = _series_vals(res[lbl], pop_names, None, "inc_hcc")
        hcc_c_2040[lbl] = sum(v_hcc[list(t_hcc).index(y + 0.5)] for y in years_range)
        t_hepb, v_hepb = _series_vals(res[lbl], pop_names, None, "hepb_pop")
        hepb_pop_2040[lbl] = _value_at_year(t_hepb, v_hepb, 2040)
        t_tot, v_tot = _series_vals(res[lbl], pop_names, None, "total_pop")
        prev_2040[lbl] = hepb_pop_2036[lbl] / _value_at_year(t_tot, v_tot, 2040)
        t_hc, v_hc = _series_vals(res[lbl], pop_names, None, "inc_hcc")
        hcc_i_2040[lbl] = _value_at_year(t_hc, v_hc, 2040)
        t_dt, v_dt = _series_vals(res[lbl], pop_names, None, "hep_dth")
        death_i_2040[lbl] = _value_at_year(t_dt, v_dt, 2040)

    table3c = [
        ("Service delivery", None, None, None, "header"),
        (f"Total number of hepatitis B tests", NA, NA, NA, "text"),
        (f"Total number of hepatitis B treatments", treatments[CSQ], treatments[ACT], d(treatments), "count"),
        ("Epidemiology", None, None, None, "header"),
        (f"People with hepatitis B in {2040}", hepb_pop_2040[CSQ], hepb_pop_2040[ACT], d(hepb_pop_2040), "count"),
        (f"Number of new hepatitis B infections, {YEAR_START}-{2040} (see note below table)",
         new_infections[CSQ], new_infections[ACT], d(new_infections), "count"),
        (f"Number of hepatitis B-related deaths, {YEAR_START}-{2040}", deaths_c_2040[CSQ], deaths_c_2040[ACT], d(deaths_c_2040), "count"),
        (f"Number of hepatitis B-related liver cancers, {YEAR_START}-{2040}", hcc_c_2040[CSQ], hcc_c_2040[ACT], d(hcc_c_2040),
         "count"),
        (f"Hepatitis B prevalence among the whole population in {2040} (%)",
         prev_2040[CSQ], prev_2040[ACT], d(prev_2040), "pct"),
        (f"HCC incidence in single year of {2040} ",
         hcc_i_2040[CSQ], hcc_i_2040[ACT], d(hcc_i_2040), "count"),
        (f"HBV deaths in single year of {2040} ",
         death_i_2040[CSQ], death_i_2040[ACT], d(death_i_2040), "count"),
    ]

    # ---- Table 3d: Service delivery and epi, 2027-2050 -------------------------------
    years_range = list(range(YEAR_START, 2050 + 1))
    treatments = {lbl: _flow_sum(res[lbl], "treat_rate", years_range) for _, lbl in SCENARIOS}
    new_infections = {lbl: _new_infections(res[lbl], pop_names, years_range) for _, lbl in SCENARIOS}
    deaths_c_2050 = {}
    hcc_c_2050 = {}
    hepb_pop_2050 = {}
    prev_2050 = {}
    death_i_2050 = {}
    hcc_i_2050 = {}
    for _, lbl in SCENARIOS:
        t_dth, v_dth = _series_vals(res[lbl], pop_names, None, "hep_dth")
        deaths_c_2050[lbl] = sum(v_dth[list(t_dth).index(y + 0.5)] for y in years_range)
        t_hcc, v_hcc = _series_vals(res[lbl], pop_names, None, "inc_hcc")
        hcc_c_2050[lbl] = sum(v_hcc[list(t_hcc).index(y + 0.5)] for y in years_range)
        t_hepb, v_hepb = _series_vals(res[lbl], pop_names, None, "hepb_pop")
        hepb_pop_2050[lbl] = _value_at_year(t_hepb, v_hepb, 2050)
        t_tot, v_tot = _series_vals(res[lbl], pop_names, None, "total_pop")
        prev_2050[lbl] = hepb_pop_2036[lbl] / _value_at_year(t_tot, v_tot, 2050)
        t_hc, v_hc = _series_vals(res[lbl], pop_names, None, "inc_hcc")
        hcc_i_2050[lbl] = _value_at_year(t_hc, v_hc, 2050)
        t_dt, v_dt = _series_vals(res[lbl], pop_names, None, "hep_dth")
        death_i_2050[lbl] = _value_at_year(t_dt, v_dt, 2050)

    table3d = [
        ("Service delivery", None, None, None, "header"),
        (f"Total number of hepatitis B tests", NA, NA, NA, "text"),
        (f"Total number of hepatitis B treatments", treatments[CSQ], treatments[ACT], d(treatments), "count"),
        ("Epidemiology", None, None, None, "header"),
        (f"People with hepatitis B in {2050}", hepb_pop_2050[CSQ], hepb_pop_2050[ACT], d(hepb_pop_2050), "count"),
        (f"Number of new hepatitis B infections, {YEAR_START}-{2050} (see note below table)",
         new_infections[CSQ], new_infections[ACT], d(new_infections), "count"),
        (f"Number of hepatitis B-related deaths, {YEAR_START}-{2050}", deaths_c_2050[CSQ], deaths_c_2050[ACT], d(deaths_c_2050), "count"),
        (f"Number of hepatitis B-related liver cancers, {YEAR_START}-{2050}", hcc_c_2050[CSQ], hcc_c_2050[ACT], d(hcc_c_2050),
         "count"),
        (f"Hepatitis B prevalence among the whole population in {2050} (%)",
         prev_2050[CSQ], prev_2050[ACT], d(prev_2050), "pct"),
        (f"HCC incidence in single year of {2050} ",
         hcc_i_2050[CSQ], hcc_i_2050[ACT], d(hcc_i_2050), "count"),
        (f"HBV deaths in single year of {2050} ",
         death_i_2050[CSQ], death_i_2050[ACT], d(death_i_2050), "count"),
    ]



    # ---- Table 4a: costs and economic outcomes 2027-2030 ---------------------------------------------
    dir_c_2030, test_c_2030, treat_c_2030, dism_c_2030, prem_c_2030, qaly_2030 = {},{},{},{},{},{}
    years_range = list(range(YEAR_START, 2030 + 1))

    for _, lbl in SCENARIOS:
        t_dir, v_dir = _series_vals(res[lbl], pop_names, None, "direct_costs")
        dir_c_2030[lbl] = sum(v_dir[list(t_dir).index(y + 0.5)] for y in years_range)
        t_test, v_test = _series_vals(res[lbl], pop_names, None, "diag_cost")
        test_c_2030[lbl] = sum(v_test[list(t_test).index(y + 0.5)] for y in years_range)
        t_treat, v_treat = _series_vals(res[lbl], pop_names, None, "treat_cost")
        treat_c_2030[lbl] = sum(v_treat[list(t_treat).index(y + 0.5)] for y in years_range)
        t_dism, v_dism = _series_vals(res[lbl], pop_names, None, "dm_cost")
        dism_c_2030[lbl] = sum(v_dism[list(t_dism).index(y + 0.5)] for y in years_range)
        t_prem, v_prem = _series_vals(res[lbl], pop_names, None, "prod_loss")
        prem_c_2030[lbl] = sum(v_prem[list(t_prem).index(y + 0.5)] for y in years_range)
        t_qaly, v_qaly = _series_vals(res[lbl], pop_names, None, "qalys_total")
        qaly_2030[lbl] = sum(v_qaly[list(t_qaly).index(y + 0.5)] for y in years_range)

    table4a = [
        (f"Costs (million A$), {YEAR_START}-{2030}", None, None, None, "header"),
        ("Total direct costs", dir_c_2030[CSQ], dir_c_2030[ACT], d(dir_c_2030), "count"),
        ("HBV testing", test_c_2030[CSQ], test_c_2030[ACT], d(test_c_2030), "count"),
        ("HBV treatment", treat_c_2030[CSQ], treat_c_2030[ACT], d(treat_c_2030), "count"),
        ("HBV disease management", dism_c_2030[CSQ], dism_c_2030[ACT], d(dism_c_2030), "count"),
        ("Societal costs", None, None, None, "header"),
        ("Absenteeism + presenteeism", COST_NOTE, COST_NOTE, COST_NOTE, "text"),
        ("Premature deaths", prem_c_2030[CSQ], prem_c_2030[ACT], d(prem_c_2030), "count"),
        ("Cost-effectiveness", None, None, None, "header"),
        ("Total QALYs", qaly_2030[CSQ], qaly_2030[ACT], d(qaly_2030), "count"),
        (f"Direct costs per QALY gained at {2030}", "-", COST_NOTE, None, "text"),
        ("Net economic benefit", NA, NA, NA, "text"),
        (f"At {YEAR_END} (millions A$)", "-", NA, None, "text"),
    ]

    # ---- Table 4b: costs and economic outcomes 2027-2036 ---------------------------------------------
    dir_c_2036, test_c_2036, treat_c_2036, dism_c_2036, prem_c_2036, qaly_2036 = {}, {}, {}, {}, {}, {}
    years_range = list(range(YEAR_START, 2036 + 1))

    for _, lbl in SCENARIOS:
        t_dir, v_dir = _series_vals(res[lbl], pop_names, None, "direct_costs")
        dir_c_2036[lbl] = sum(v_dir[list(t_dir).index(y + 0.5)] for y in years_range)
        t_test, v_test = _series_vals(res[lbl], pop_names, None, "diag_cost")
        test_c_2036[lbl] = sum(v_test[list(t_test).index(y + 0.5)] for y in years_range)
        t_treat, v_treat = _series_vals(res[lbl], pop_names, None, "treat_cost")
        treat_c_2036[lbl] = sum(v_treat[list(t_treat).index(y + 0.5)] for y in years_range)
        t_dism, v_dism = _series_vals(res[lbl], pop_names, None, "dm_cost")
        dism_c_2036[lbl] = sum(v_dism[list(t_dism).index(y + 0.5)] for y in years_range)
        t_prem, v_prem = _series_vals(res[lbl], pop_names, None, "prod_loss")
        prem_c_2036[lbl] = sum(v_prem[list(t_prem).index(y + 0.5)] for y in years_range)
        t_qaly, v_qaly = _series_vals(res[lbl], pop_names, None, "qalys_total")
        qaly_2036[lbl] = sum(v_qaly[list(t_qaly).index(y + 0.5)] for y in years_range)

    table4b = [
        (f"Costs (million A$), {YEAR_START}-{2036}", None, None, None, "header"),
        ("Total direct costs", dir_c_2036[CSQ], dir_c_2036[ACT], d(dir_c_2036), "count"),
        ("HBV testing", test_c_2036[CSQ], test_c_2036[ACT], d(test_c_2036), "count"),
        ("HBV treatment", treat_c_2036[CSQ], treat_c_2036[ACT], d(treat_c_2036), "count"),
        ("HBV disease management", dism_c_2036[CSQ], dism_c_2036[ACT], d(dism_c_2036), "count"),
        ("Societal costs", None, None, None, "header"),
        ("Absenteeism + presenteeism", COST_NOTE, COST_NOTE, COST_NOTE, "text"),
        ("Premature deaths", prem_c_2036[CSQ], prem_c_2036[ACT], d(prem_c_2036), "count"),
        ("Cost-effectiveness", None, None, None, "header"),
        ("Total QALYs", qaly_2036[CSQ], qaly_2036[ACT], d(qaly_2036), "count"),
        (f"Direct costs per QALY gained at {2036}", "-", COST_NOTE, None, "text"),
        ("Net economic benefit", NA, NA, NA, "text"),
        (f"At {YEAR_END} (millions A$)", "-", NA, None, "text"),
    ]

    # ---- Table 4c: costs and economic outcomes 2027-2040 ---------------------------------------------
    dir_c_2040, test_c_2040, treat_c_2040, dism_c_2040, prem_c_2040, qaly_2040 = {}, {}, {}, {}, {}, {}
    years_range = list(range(YEAR_START, 2040 + 1))

    for _, lbl in SCENARIOS:
        t_dir, v_dir = _series_vals(res[lbl], pop_names, None, "direct_costs")
        dir_c_2040[lbl] = sum(v_dir[list(t_dir).index(y + 0.5)] for y in years_range)
        t_test, v_test = _series_vals(res[lbl], pop_names, None, "diag_cost")
        test_c_2040[lbl] = sum(v_test[list(t_test).index(y + 0.5)] for y in years_range)
        t_treat, v_treat = _series_vals(res[lbl], pop_names, None, "treat_cost")
        treat_c_2040[lbl] = sum(v_treat[list(t_treat).index(y + 0.5)] for y in years_range)
        t_dism, v_dism = _series_vals(res[lbl], pop_names, None, "dm_cost")
        dism_c_2040[lbl] = sum(v_dism[list(t_dism).index(y + 0.5)] for y in years_range)
        t_prem, v_prem = _series_vals(res[lbl], pop_names, None, "prod_loss")
        prem_c_2040[lbl] = sum(v_prem[list(t_prem).index(y + 0.5)] for y in years_range)
        t_qaly, v_qaly = _series_vals(res[lbl], pop_names, None, "qalys_total")
        qaly_2040[lbl] = sum(v_qaly[list(t_qaly).index(y + 0.5)] for y in years_range)

    table4c = [
        (f"Costs (million A$), {YEAR_START}-{2040}", None, None, None, "header"),
        ("Total direct costs", dir_c_2040[CSQ], dir_c_2040[ACT], d(dir_c_2040), "count"),
        ("HBV testing", test_c_2040[CSQ], test_c_2040[ACT], d(test_c_2040), "count"),
        ("HBV treatment", treat_c_2040[CSQ], treat_c_2040[ACT], d(treat_c_2040), "count"),
        ("HBV disease management", dism_c_2040[CSQ], dism_c_2040[ACT], d(dism_c_2040), "count"),
        ("Societal costs", None, None, None, "header"),
        ("Absenteeism + presenteeism", COST_NOTE, COST_NOTE, COST_NOTE, "text"),
        ("Premature deaths", prem_c_2040[CSQ], prem_c_2040[ACT], d(prem_c_2040), "count"),
        ("Cost-effectiveness", None, None, None, "header"),
        ("Total QALYs", qaly_2040[CSQ], qaly_2040[ACT], d(qaly_2040), "count"),
        (f"Direct costs per QALY gained at {2040}", "-", COST_NOTE, None, "text"),
        ("Net economic benefit", NA, NA, NA, "text"),
        (f"At {YEAR_END} (millions A$)", "-", NA, None, "text"),
    ]

    # ---- Table 4d: costs and economic outcomes 2027-2050 ---------------------------------------------
    dir_c_2050, test_c_2050, treat_c_2050, dism_c_2050, prem_c_2050, qaly_2050 = {}, {}, {}, {}, {}, {}
    years_range = list(range(YEAR_START, 2050 + 1))

    for _, lbl in SCENARIOS:
        t_dir, v_dir = _series_vals(res[lbl], pop_names, None, "direct_costs")
        dir_c_2050[lbl] = sum(v_dir[list(t_dir).index(y + 0.5)] for y in years_range)
        t_test, v_test = _series_vals(res[lbl], pop_names, None, "diag_cost")
        test_c_2050[lbl] = sum(v_test[list(t_test).index(y + 0.5)] for y in years_range)
        t_treat, v_treat = _series_vals(res[lbl], pop_names, None, "treat_cost")
        treat_c_2050[lbl] = sum(v_treat[list(t_treat).index(y + 0.5)] for y in years_range)
        t_dism, v_dism = _series_vals(res[lbl], pop_names, None, "dm_cost")
        dism_c_2050[lbl] = sum(v_dism[list(t_dism).index(y + 0.5)] for y in years_range)
        t_prem, v_prem = _series_vals(res[lbl], pop_names, None, "prod_loss")
        prem_c_2050[lbl] = sum(v_prem[list(t_prem).index(y + 0.5)] for y in years_range)
        t_qaly, v_qaly = _series_vals(res[lbl], pop_names, None, "qalys_total")
        qaly_2050[lbl] = sum(v_qaly[list(t_qaly).index(y + 0.5)] for y in years_range)

    table4d = [
        (f"Costs (million A$), {YEAR_START}-{2050}", None, None, None, "header"),
        ("Total direct costs", dir_c_2050[CSQ], dir_c_2050[ACT], d(dir_c_2050), "count"),
        ("HBV testing", test_c_2050[CSQ], test_c_2050[ACT], d(test_c_2050), "count"),
        ("HBV treatment", treat_c_2050[CSQ], treat_c_2050[ACT], d(treat_c_2050), "count"),
        ("HBV disease management", dism_c_2050[CSQ], dism_c_2050[ACT], d(dism_c_2050), "count"),
        ("Societal costs", None, None, None, "header"),
        ("Absenteeism + presenteeism", COST_NOTE, COST_NOTE, COST_NOTE, "text"),
        ("Premature deaths", prem_c_2050[CSQ], prem_c_2050[ACT], d(prem_c_2050), "count"),
        ("Cost-effectiveness", None, None, None, "header"),
        ("Total QALYs", qaly_2050[CSQ], qaly_2050[ACT], d(qaly_2050), "count"),
        (f"Direct costs per QALY gained at {2050}", "-", COST_NOTE, None, "text"),
        ("Net economic benefit", NA, NA, NA, "text"),
        (f"At {YEAR_END} (millions A$)", "-", NA, None, "text"),
    ]


    # ---- Table 5a: progress towards targets, 2030 ------------------------------------------
    diag_2030, ltc_2030, treat_2030 = {}, {}, {}
    inc_2030, inc_2015, mort_2030, mort_2015 = {}, {}, {}, {}
    for _, lbl in SCENARIOS:
        t_d, v_d = _cascade_vals(res[lbl], pop_names, None, "diagnosed")
        diag_2030[lbl] = _value_at_year(t_d, v_d, TARGET_YEAR)
        t_l, v_l = _cascade_vals(res[lbl], pop_names, None, "linked")
        ltc_2030[lbl] = _value_at_year(t_l, v_l, TARGET_YEAR)
        t_tr, v_tr = _cascade_vals(res[lbl], pop_names, None, "treated")
        treat_2030[lbl] = _value_at_year(t_tr, v_tr, TARGET_YEAR)

        inc_2030[lbl] = _new_infections_at_year(res[lbl], pop_names, TARGET_YEAR)
        inc_2015[lbl] = _new_infections_at_year(res[lbl], pop_names, BASELINE_YEAR)
        t_dth, v_dth = _series_vals(res[lbl], pop_names, None, "hep_dth")
        mort_2030[lbl] = _value_at_year(t_dth, v_dth, TARGET_YEAR)
        mort_2015[lbl] = _value_at_year(t_dth, v_dth, BASELINE_YEAR)

    def reduction(cur, base):
        return {lbl: 1 - (cur[lbl] / base[lbl]) for _, lbl in SCENARIOS}

    inc_reduction = reduction(inc_2030, inc_2015)
    mort_reduction = reduction(mort_2030, mort_2015)

    table5a = [
        ("Hepatitis B", None, None, None, "header"),
        ("% Diagnosed hepatitis B", "90%", diag_2030[CSQ], diag_2030[ACT], ("text", "pct", "pct")),
        ("% Linked to hepatitis B care", "80%", ltc_2030[CSQ], ltc_2030[ACT], ("text", "pct", "pct")),
        ("% On treatment for hepatitis B", "27%", treat_2030[CSQ], treat_2030[ACT], ("text", "pct", "pct")),
        (f"Reduction in hepatitis B incidence by {TARGET_YEAR} (vs {BASELINE_YEAR}) (see note below table)",
         "95%", inc_reduction[CSQ], inc_reduction[ACT], ("text", "pct", "pct")),
        (f"Reduction in hepatitis B mortality by {TARGET_YEAR} (vs {BASELINE_YEAR})",
         "30%", mort_reduction[CSQ], mort_reduction[ACT], ("text", "pct", "pct")),
    ]

    # ---- Table 5b: progress towards targets, 2036 ------------------------------------------
    diag_2036, ltc_2036, treat_2036 = {}, {}, {}
    inc_2036, mort_2036 = {}, {}
    for _, lbl in SCENARIOS:
        t_d, v_d = _cascade_vals(res[lbl], pop_names, None, "diagnosed")
        diag_2036[lbl] = _value_at_year(t_d, v_d, 2036)
        t_l, v_l = _cascade_vals(res[lbl], pop_names, None, "linked")
        ltc_2036[lbl] = _value_at_year(t_l, v_l, 2036)
        t_tr, v_tr = _cascade_vals(res[lbl], pop_names, None, "treated")
        treat_2036[lbl] = _value_at_year(t_tr, v_tr, 2036)

        inc_2036[lbl] = _new_infections_at_year(res[lbl], pop_names, 2036)
        t_dth, v_dth = _series_vals(res[lbl], pop_names, None, "hep_dth")
        mort_2036[lbl] = _value_at_year(t_dth, v_dth, 2036)


    inc_reduction = reduction(inc_2036, inc_2015)
    mort_reduction = reduction(mort_2036, mort_2015)

    table5b = [
        ("Hepatitis B", None, None, None, "header"),
        ("% Diagnosed hepatitis B", "90%", diag_2036[CSQ], diag_2036[ACT], ("text", "pct", "pct")),
        ("% Linked to hepatitis B care", "80%", ltc_2036[CSQ], ltc_2036[ACT], ("text", "pct", "pct")),
        ("% On treatment for hepatitis B", "27%", treat_2036[CSQ], treat_2036[ACT], ("text", "pct", "pct")),
        (f"Reduction in hepatitis B incidence by {2036} (vs {BASELINE_YEAR}) (see note below table)",
         "95%", inc_reduction[CSQ], inc_reduction[ACT], ("text", "pct", "pct")),
        (f"Reduction in hepatitis B mortality by {2036} (vs {BASELINE_YEAR})",
         "30%", mort_reduction[CSQ], mort_reduction[ACT], ("text", "pct", "pct")),
    ]

    # ---- Table 5c: progress towards targets, 2040 ------------------------------------------
    diag_2040, ltc_2040, treat_2040 = {}, {}, {}
    inc_2040, mort_2040 = {}, {}
    for _, lbl in SCENARIOS:
        t_d, v_d = _cascade_vals(res[lbl], pop_names, None, "diagnosed")
        diag_2040[lbl] = _value_at_year(t_d, v_d, 2040)
        t_l, v_l = _cascade_vals(res[lbl], pop_names, None, "linked")
        ltc_2040[lbl] = _value_at_year(t_l, v_l, 2040)
        t_tr, v_tr = _cascade_vals(res[lbl], pop_names, None, "treated")
        treat_2040[lbl] = _value_at_year(t_tr, v_tr, 2040)

        inc_2040[lbl] = _new_infections_at_year(res[lbl], pop_names, 2040)
        t_dth, v_dth = _series_vals(res[lbl], pop_names, None, "hep_dth")
        mort_2040[lbl] = _value_at_year(t_dth, v_dth, 2040)


    inc_reduction = reduction(inc_2040, inc_2015)
    mort_reduction = reduction(mort_2040, mort_2015)

    table5c = [
        ("Hepatitis B", None, None, None, "header"),
        ("% Diagnosed hepatitis B", "90%", diag_2040[CSQ], diag_2040[ACT], ("text", "pct", "pct")),
        ("% Linked to hepatitis B care", "80%", ltc_2040[CSQ], ltc_2040[ACT], ("text", "pct", "pct")),
        ("% On treatment for hepatitis B", "27%", treat_2040[CSQ], treat_2040[ACT], ("text", "pct", "pct")),
        (f"Reduction in hepatitis B incidence by {2040} (vs {BASELINE_YEAR}) (see note below table)",
         "95%", inc_reduction[CSQ], inc_reduction[ACT], ("text", "pct", "pct")),
        (f"Reduction in hepatitis B mortality by {2040} (vs {BASELINE_YEAR})",
         "30%", mort_reduction[CSQ], mort_reduction[ACT], ("text", "pct", "pct")),
    ]

    # ---- Table 5d: progress towards targets, 2050 ------------------------------------------
    diag_2050, ltc_2050, treat_2050 = {}, {}, {}
    inc_2050, mort_2050 = {}, {}
    for _, lbl in SCENARIOS:
        t_d, v_d = _cascade_vals(res[lbl], pop_names, None, "diagnosed")
        diag_2050[lbl] = _value_at_year(t_d, v_d, 2050)
        t_l, v_l = _cascade_vals(res[lbl], pop_names, None, "linked")
        ltc_2050[lbl] = _value_at_year(t_l, v_l, 2050)
        t_tr, v_tr = _cascade_vals(res[lbl], pop_names, None, "treated")
        treat_2050[lbl] = _value_at_year(t_tr, v_tr, 2050)

        inc_2050[lbl] = _new_infections_at_year(res[lbl], pop_names, 2050)
        t_dth, v_dth = _series_vals(res[lbl], pop_names, None, "hep_dth")
        mort_2050[lbl] = _value_at_year(t_dth, v_dth, 2050)


    inc_reduction = reduction(inc_2050, inc_2015)
    mort_reduction = reduction(mort_2050, mort_2015)

    table5d = [
        ("Hepatitis B", None, None, None, "header"),
        ("% Diagnosed hepatitis B", "90%", diag_2050[CSQ], diag_2050[ACT], ("text", "pct", "pct")),
        ("% Linked to hepatitis B care", "80%", ltc_2050[CSQ], ltc_2050[ACT], ("text", "pct", "pct")),
        ("% On treatment for hepatitis B", "27%", treat_2050[CSQ], treat_2050[ACT], ("text", "pct", "pct")),
        (f"Reduction in hepatitis B incidence by {2050} (vs {BASELINE_YEAR}) (see note below table)",
         "95%", inc_reduction[CSQ], inc_reduction[ACT], ("text", "pct", "pct")),
        (f"Reduction in hepatitis B mortality by {2050} (vs {BASELINE_YEAR})",
         "30%", mort_reduction[CSQ], mort_reduction[ACT], ("text", "pct", "pct")),
    ]


    return table3a, table3b,table3c, table3d, table4a, table4b, table4c, table4d, table5a, table5b, table5c, table5d


def write_excel(table3a, table3b, table3c, table3d, table4a, table4b, table4c, table4d, table5a, table5b, table5c, table5d, save_path):
    import pandas as pd
    writer = pd.ExcelWriter(save_path, engine="xlsxwriter")
    workbook = writer.book
    fmt_header = workbook.add_format({"bold": True})
    fmt_count = workbook.add_format({"num_format": "#,##0"})
    fmt_pct = workbook.add_format({"num_format": "0.0%"})
    fmt_note = workbook.add_format({"italic": True, "font_color": "#7F7F7F", "text_wrap": True})

    def kind_format(kind):
        return {"count": fmt_count, "pct": fmt_pct}.get(kind)

    def write_cell(worksheet, r, c, val, kind):
        if val is None:
            return
        if kind in ("count", "pct") and isinstance(val, (int, float, np.floating, np.integer)):
            worksheet.write_number(r, c, val, kind_format(kind))
        else:
            worksheet.write_string(r, c, str(val))

    def write_sheet(sheet_name, columns, rows, note=None):
        worksheet = workbook.add_worksheet(sheet_name)
        worksheet.set_column(0, 0, 55)
        worksheet.set_column(1, len(columns) - 1, 20)
        for c, col_name in enumerate(columns):
            worksheet.write(0, c, col_name, fmt_header)
        r = 1
        for row in rows:
            indicator, v1, v2, v3, kind = row
            is_header = kind == "header"
            worksheet.write(r, 0, indicator, fmt_header if is_header else None)
            kinds = kind if isinstance(kind, tuple) else (kind, kind, kind)
            for c, (val, val_kind) in enumerate(zip((v1, v2, v3), kinds), start=1):
                write_cell(worksheet, r, c, val, val_kind)
            r += 1
        if note:
            worksheet.write(r + 1, 0, note, fmt_note)
        return r

    # Table 3a (2027-2030)
    write_sheet("Table 3a_2030", ["Indicator", "Continued status quo (Total)", "Action scenario (Total)", "Difference from status quo"], table3a,
        note=("Note: 'Number of new hepatitis B infections' depends on the calibrated force-of-infection "
              "(foi_cal) - check it against run_scenarios.py's incidence panel before citing this figure."),)

    # Table 3b (2027-2036)
    write_sheet("Table 3b_2036", ["Indicator", "Continued status quo (Total)", "Action scenario (Total)", "Difference from status quo"], table3b,
        note=("Note: 'Number of new hepatitis B infections' depends on the calibrated force-of-infection "
              "(foi_cal) - check it against run_scenarios.py's incidence panel before citing this figure."),)
    # Table 3c (2027-2040)
    write_sheet("Table 3c_2040", ["Indicator", "Continued status quo (Total)", "Action scenario (Total)", "Difference from status quo"], table3c,
        note=("Note: 'Number of new hepatitis B infections' depends on the calibrated force-of-infection "
              "(foi_cal) - check it against run_scenarios.py's incidence panel before citing this figure."),)
    # Table 3d (2027-2050)
    write_sheet("Table 3d_2050", ["Indicator", "Continued status quo (Total)", "Action scenario (Total)", "Difference from status quo"], table3d,
        note=("Note: 'Number of new hepatitis B infections' depends on the calibrated force-of-infection "
              "(foi_cal) - check it against run_scenarios.py's incidence panel before citing this figure."),)

    write_sheet("Table 4a_2030", ["Indicator", "Continued status quo (Total)", "Action scenario (Total)", "Difference from status quo"], table4a)
    write_sheet("Table 4b_2036", ["Indicator", "Continued status quo (Total)", "Action scenario (Total)", "Difference from status quo"], table4b)
    write_sheet("Table 4c_2040", ["Indicator", "Continued status quo (Total)", "Action scenario (Total)", "Difference from status quo"], table4c)
    write_sheet("Table 4d_2050", ["Indicator", "Continued status quo (Total)", "Action scenario (Total)", "Difference from status quo"], table4d)



    write_sheet(
        "Table 5a_2030", ["Indicator", "Target", "Continued status quo", "Action scenario"], table5a,
        note="Note: incidence-reduction row inherits the same new-infections caveat as Table 3.",
    )

    write_sheet(
        "Table 5b_2036", ["Indicator", "Target", "Continued status quo", "Action scenario"], table5b,
        note="Note: incidence-reduction row inherits the same new-infections caveat as Table 3.",
    )

    write_sheet(
        "Table 5c_2040", ["Indicator", "Target", "Continued status quo", "Action scenario"], table5c,
        note="Note: incidence-reduction row inherits the same new-infections caveat as Table 3.",
    )
    write_sheet(
        "Table 5d_2050", ["Indicator", "Target", "Continued status quo", "Action scenario"], table5d,
        note="Note: incidence-reduction row inherits the same new-infections caveat as Table 3.",
    )

    writer.close()
    print(f"Saved: {save_path}")


def main():
    out_dir = os.path.join(_get_github_folder(), "outputs")
    results = sc.load(os.path.join(out_dir, "hepaus_scenarios.pkl"))
    table3a,table3b,table3c, table3d, table4a, table4b, table4c, table4d, table5a, table5b, table5c, table5d = build_tables(results)
    save_path = os.path.join(out_dir, "HBV_report_tables.xlsx")
    write_excel(table3a, table3b, table3c, table3d, table4a, table4b, table4c, table4d,table5a, table5b, table5c, table5d, save_path)


if __name__ == "__main__":
    main()
