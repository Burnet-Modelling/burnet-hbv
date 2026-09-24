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

YEAR_START, YEAR_END = 2027, 2036
TARGET_YEAR = 2030
BASELINE_YEAR = 2015

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
    years_range = list(range(YEAR_START, YEAR_END + 1))

    res = {label: results[key][2] for key, label in SCENARIOS}

    # ---- Table 3: service delivery + epidemiology, 2027-2036 ----------------------------
    treatments = {lbl: _flow_sum(res[lbl], "treat_rate", years_range) for _, lbl in SCENARIOS}
    new_infections = {lbl: _new_infections(res[lbl], pop_names, years_range) for _, lbl in SCENARIOS}
    deaths = {}
    hepb_pop_2036 = {}
    prev_2036 = {}
    for _, lbl in SCENARIOS:
        t_dth, v_dth = _series_vals(res[lbl], pop_names, None, "hep_dth")
        deaths[lbl] = sum(v_dth[list(t_dth).index(y + 0.5)] for y in years_range)
        t_hepb, v_hepb = _series_vals(res[lbl], pop_names, None, "hepb_pop")
        hepb_pop_2036[lbl] = _value_at_year(t_hepb, v_hepb, YEAR_END)
        t_tot, v_tot = _series_vals(res[lbl], pop_names, None, "total_pop")
        prev_2036[lbl] = hepb_pop_2036[lbl] / _value_at_year(t_tot, v_tot, YEAR_END)

    def d(dct):
        return dct[ACT] - dct[CSQ]

    table3 = [
        ("Service delivery", None, None, None, "header"),
        (f"Total number of hepatitis B tests", NA, NA, NA, "text"),
        (f"Total number of hepatitis B treatments", treatments[CSQ], treatments[ACT], d(treatments), "count"),
        ("Epidemiology", None, None, None, "header"),
        (f"People with hepatitis B in {YEAR_END}", hepb_pop_2036[CSQ], hepb_pop_2036[ACT], d(hepb_pop_2036), "count"),
        (f"Number of new hepatitis B infections, {YEAR_START}-{YEAR_END} (see note below table)",
         new_infections[CSQ], new_infections[ACT], d(new_infections), "count"),
        (f"Number of hepatitis B-related deaths, {YEAR_START}-{YEAR_END}", deaths[CSQ], deaths[ACT], d(deaths), "count"),
        (f"Hepatitis B prevalence among the whole population in {YEAR_END} (%)",
         prev_2036[CSQ], prev_2036[ACT], d(prev_2036), "pct"),
    ]

    # ---- Table 4: costs and economic outcomes ---------------------------------------------
    dir_c_2036, test_c_2036, treat_c_2036, dism_c_2036, prem_c_2036, qaly_2036 = {},{},{},{},{},{}
    for _, lbl in SCENARIOS:
        t_dir, v_dir = _series_vals(res[lbl], pop_names, None, "direct_costs")
        dir_c_2036[lbl] = sum(v_dir[list(t_dir).index(y + 0.5)] for y in years_range)
        t_test, v_test = _series_vals(res[lbl], pop_names, None, "diag_cost")
        test_c_2036[lbl] = sum(v_test[list(t_test).index(y + 0.5)] for y in years_range)
        t_treat, v_treat = _series_vals(res[lbl], pop_names, None, "treat_cost")
        treat_c_2036[lbl] = sum(v_treat[list(t_treat).index(y + 0.5)] for y in years_range)
        t_dism, v_dism = _series_vals(res[lbl], pop_names, None, "dm_cost")
        dism_c_2036[lbl] = sum(v_dism[list(t_dism).index(y + 0.5)] for y in years_range)
        t_prem, v_prem = _series_vals(res[lbl], pop_names, None, "mort_cost")
        prem_c_2036[lbl] = sum(v_prem[list(t_prem).index(y + 0.5)] for y in years_range)
        t_qaly, v_qaly = _series_vals(res[lbl], pop_names, None, "qalys_total")
        qaly_2036[lbl] = sum(v_qaly[list(t_qaly).index(y + 0.5)] for y in years_range)
    table4 = [
        (f"Costs (million A$), {YEAR_START}-{YEAR_END}", None, None, None, "header"),
        ("Total direct costs", dir_c_2036[CSQ], dir_c_2036[ACT], d(dir_c_2036), "count"),
        ("HBV testing", test_c_2036[CSQ], test_c_2036[ACT], d(test_c_2036), "count"),
        ("HBV treatment", treat_c_2036[CSQ], treat_c_2036[ACT], d(treat_c_2036), "count"),
        ("HBV disease management", dism_c_2036[CSQ], dism_c_2036[ACT], d(dism_c_2036), "count"),
        ("Societal costs", None, None, None, "header"),
        ("Absenteeism + presenteeism", COST_NOTE, COST_NOTE, COST_NOTE, "text"),
        ("Premature deaths", prem_c_2036[CSQ], prem_c_2036[ACT], d(prem_c_2036), "count"),
        ("Cost-effectiveness", None, None, None, "header"),
        ("Total QALYs", qaly_2036[CSQ], qaly_2036[ACT], d(qaly_2036), "count"),
        (f"Direct costs per QALY gained at {YEAR_END}", "-", COST_NOTE, None, "text"),
        ("Net economic benefit", NA, NA, NA, "text"),
        (f"At {YEAR_END} (millions A$)", "-", NA, None, "text"),
    ]

    # ---- Table 5: progress towards targets, 2030 ------------------------------------------
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

    table5 = [
        ("Hepatitis B", None, None, None, "header"),
        ("% Diagnosed hepatitis B", "90%", diag_2030[CSQ], diag_2030[ACT], ("text", "pct", "pct")),
        ("% Linked to hepatitis B care", "80%", ltc_2030[CSQ], ltc_2030[ACT], ("text", "pct", "pct")),
        ("% On treatment for hepatitis B", "27%", treat_2030[CSQ], treat_2030[ACT], ("text", "pct", "pct")),
        (f"Reduction in hepatitis B incidence by {TARGET_YEAR} (vs {BASELINE_YEAR}) (see note below table)",
         "95%", inc_reduction[CSQ], inc_reduction[ACT], ("text", "pct", "pct")),
        (f"Reduction in hepatitis B mortality by {TARGET_YEAR} (vs {BASELINE_YEAR})",
         "30%", mort_reduction[CSQ], mort_reduction[ACT], ("text", "pct", "pct")),
    ]

    return table3, table4, table5


def write_excel(table3, table4, table5, save_path):
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

    write_sheet(
        "Table 3", ["Indicator", "Continued status quo (Total)", "Action scenario (Total)", "Difference from status quo"],
        table3,
        note=("Note: 'Number of new hepatitis B infections' depends on the calibrated force-of-infection "
              "(foi_cal) - check it against run_scenarios.py's incidence panel before citing this figure."),
    )
    write_sheet("Table 4", ["Indicator", "Continued status quo (Total)", "Action scenario (Total)", "Difference from status quo"], table4)
    write_sheet(
        "Table 5", ["Indicator", "Target", "Continued status quo", "Action scenario"], table5,
        note="Note: incidence-reduction row inherits the same new-infections caveat as Table 3.",
    )

    writer.close()
    print(f"Saved: {save_path}")


def main():
    out_dir = os.path.join(_get_github_folder(), "outputs")
    results = sc.load(os.path.join(out_dir, "hepaus_scenarios.pkl"))
    table3, table4, table5 = build_tables(results)
    save_path = os.path.join(out_dir, "HBV_report_tables.xlsx")
    write_excel(table3, table4, table5, save_path)


if __name__ == "__main__":
    main()
