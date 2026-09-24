"""
Run the HepAus model for two scenarios - a "status quo" scenario (current testing,
linkage-to-care and treatment rates held constant) and an "action" scenario (rates
scaled up per the Hepatitis Strategy Cost Worksheet) - and save the results for
comparison.

Produces two figures in outputs/figures/:
  - hepaus_scenarios_comparison.png - national aggregate outcomes (publication-formatted)
  - hepaus_population_breakdown.png - the same care-cascade metrics broken down by
    population group (one row per population, plus a "Total" aggregate row)

Usage:
    python run_scenarios.py
"""
import os
import string
import numpy as np
import sciris as sc
import atomica as at
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator

from hbv_aus.model_run import run_hepaus_scenarios
from hbv_aus.utils import _get_github_folder

# Publication-style figure defaults
COLOR_STATUS_QUO = "#2A78D6"
COLOR_ACTION = "#EB6834"
COLOR_DATA = "#3B3B3B"

# Origin-group prefixes and colors for the population-composition stacked area chart
ORIGIN_GROUPS = [
    ("atsi", "Aboriginal and Torres Strait Islander", "#2A78D6"),
    ("ausb", "Australian-born", "#1BAF7A"),
    ("lros", "Overseas-born (low risk)", "#EDA100"),
    ("hros", "Overseas-born (high risk)", "#EB6834"),
]

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10,
    "axes.titlesize": 10.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 9.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": "#4D4D4D",
    "axes.linewidth": 0.8,
    "axes.grid": True,
    "grid.color": "#DDDDDD",
    "grid.linewidth": 0.6,
    "axes.axisbelow": True,
    "savefig.dpi": 300,
    "figure.dpi": 150,
})

# Outcomes to compare between scenarios (panel key, framework code name, display name,
# kind). "key" is the dict key used to store/retrieve this panel's series - needed
# because e.g. "diagnosed" appears twice below (once as a raw count, once as a % of all
# HBV), and both can't share one key. kind "count" = a number (people/deaths); summed
# across populations and, where the databook has real data for it, overlaid as scatter.
# kind "proportion" = a fraction (0-1), data (if any) aggregated as a total_pop-weighted
# average since summing proportions across populations is meaningless. kind "cascade" =
# code_name's national count as a fraction of national hepb_pop (% diagnosed/linked/
# treated, of everyone with HBV) - no real data exists for these.
OUTCOMES = [
    ("total_pop", "total_pop", "Total Population", "count"),
    ("hepb_prev", "hepb_prev", "HBsAg Prevalence", "proportion"),
    ("hepb_pop", "hepb_pop", "Population Living with HBV", "count"),
    ("diagnosed", "diagnosed", "Diagnosed with HBV", "count"),
    ("treated", "treated", "On Antiviral Treatment", "count"),
    ("compensated", "compensated", "Compensated Cirrhosis", "count"),
    ("decompensated", "decompensated", "Decompensated Cirrhosis", "count"),
    ("livercancer", "livercancer", "Liver Cancer", "count"),
    ("hep_dth", "hep_dth", "HBV-Attributable Deaths (annual)", "count"),
    ("incidence", None, "New HBV Infections (annual)", "incidence"),
    ("pct_diagnosed", "diagnosed", "% Diagnosed (of all HBV)", "cascade"),
    ("pct_linked", "linked", "% Linked to Care (of all HBV)", "cascade"),
    ("pct_treated", "treated", "% On Treatment (of all HBV)", "cascade"),
]

# Columns for the per-population breakdown figure (framework code, display name, kind).
# "cascade" = code_name's count as a fraction of hepb_pop for that same row (% diagnosed/
# linked/treated, of all people with HBV) - no real data exists for these.
BREAKDOWN_COLUMNS = [
    ("total_pop", "Population Size", "count"),
    ("hepb_prev", "HBV Prevalence", "prevalence"),
    ("diagnosed", "% Diagnosed", "cascade"),
    ("linked", "% Linked to Care", "cascade"),
    ("treated", "% On Treatment", "cascade"),
]

SCENARIOS = [
    ("Baseline", "Status Quo"),
    ("Scenario 1", "Action Scenario"),
]


def _series_vals(res, pop_names, pop, code_name):
    """
    Model values for `code_name`: for a specific real population (pop="atsi_0-14_M",
    etc.), or summed across every real population if pop is None (the aggregate/
    "Total" row/panel).
    """
    if pop is not None:
        s = at.PlotData(res, outputs=code_name, pops=pop, t_bins=1).series[0]
        return s.tvec, s.vals
    tvec, total = None, None
    for p in pop_names:
        s = at.PlotData(res, outputs=code_name, pops=p, t_bins=1).series[0]
        if total is None:
            tvec, total = s.tvec, s.vals.copy()
        else:
            total = total + s.vals
    return tvec, total


def _prevalence_vals(res, pop_names, pop):
    """
    HBV prevalence: for a real population this is just hepb_prev directly. For the
    aggregate row it must be the population-weighted sum(hepb_pop)/sum(total_pop) -
    NOT atomica's own pops="total" aggregation, which (for a proportion-type
    characteristic) averages populations unweighted rather than weighting by size.
    """
    if pop is not None:
        s = at.PlotData(res, outputs="hepb_prev", pops=pop, t_bins=1).series[0]
        return s.tvec, s.vals
    tvec, hepb = _series_vals(res, pop_names, None, "hepb_pop")
    _, tot = _series_vals(res, pop_names, None, "total_pop")
    return tvec, hepb / tot


def _cascade_vals(res, pop_names, pop, numerator_code):
    """numerator_code's count as a fraction of hepb_pop, for the same row (population
    or aggregate)."""
    tvec, num = _series_vals(res, pop_names, pop, numerator_code)
    _, hepb = _series_vals(res, pop_names, pop, "hepb_pop")
    return tvec, np.divide(num, hepb, out=np.zeros_like(num), where=hepb != 0)


def _incidence_vals(res, pop_names):
    """
    New HBV infections per year, nationally: horizontal (community) transmission flow +
    mother-to-child transmission (b_chb). "horiz" is a rate/probability-formatted
    parameter, so its actual per-year flow (people, not the raw rate) has to be read off
    the underlying model links rather than via PlotData. Its own panel so incidence stays
    visible and can be sanity-checked against published Australian estimates whenever
    foi_cal is recalibrated (see calibrate_model.py).
    """
    horiz_by_year = {}
    for pop in res.model.pops:
        for link in pop.links:
            if link.parameter.name != "horiz_chb":
                continue
            for ti, y in enumerate(link.t):
                horiz_by_year[y] = horiz_by_year.get(y, 0.0) + link.vals[ti]

    years = sorted(horiz_by_year.keys())
    horiz_vals = np.array([horiz_by_year[y] for y in years])

    # b_chb is already a per-year count (mother-to-child transmission); its PlotData tvec
    # uses the midpoint convention (year+0.5) - align onto horiz's integer-year keys.
    b_tvec, b_vals = _series_vals(res, pop_names, None, "b_chb")
    b_chb_by_year = {int(np.floor(t)): v for t, v in zip(b_tvec, b_vals)}

    combined = horiz_vals + np.array([b_chb_by_year.get(y, 0.0) for y in years])
    tvec_out = np.array(years) + 0.5  # match PlotData's midpoint convention for the x-axis
    return tvec_out, combined


def aggregate_data(P, code_name, kind):
    """
    Aggregate real (databook) input data for `code_name` across all populations that
    have it, into a single national series - summed for a "count", or weighted by
    total_pop for a "proportion". Returns (years, values) as numpy arrays, or None if
    no population has data for this quantity.
    """
    pop_names = list(P.data.pops.keys())
    years = None
    val_sum = None
    weight_sum = None
    for pop in pop_names:
        ts = P.data.get_ts(code_name, pop)
        if ts is None or not ts.has_time_data:
            continue
        t, v = np.array(ts.t, dtype=float), np.array(ts.vals, dtype=float)
        if years is None:
            years = t
            val_sum = np.zeros_like(t)
            weight_sum = np.zeros_like(t) if kind == "proportion" else None
        if kind == "proportion":
            w_ts = P.data.get_ts("total_pop", pop)
            if w_ts is not None and w_ts.has_time_data:
                w = np.interp(t, np.array(w_ts.t, dtype=float), np.array(w_ts.vals, dtype=float))
            else:
                w = np.ones_like(t)
            val_sum += v * w
            weight_sum += w
        else:
            val_sum += v

    if years is None:
        return None
    return years, (val_sum / weight_sum if kind == "proportion" else val_sum)


def _data_vals(P, pop, code_name, kind):
    """Real (databook) input data, if any. pop=None aggregates across every real
    population (summed for "count", total_pop-weighted for "proportion")."""
    if pop is not None:
        ts = P.data.get_ts(code_name, pop)
        if ts is None or not ts.has_time_data:
            return None
        return np.array(ts.t, dtype=float), np.array(ts.vals, dtype=float)
    return aggregate_data(P, code_name, kind)


def extract_outcomes(results):
    """Pull the national comparison outcomes (model lines + any real data) out of a
    {scenario_key: [Project, ParameterSet, Result]} dict, for the aggregate figure."""
    data = {"model": {}, "data": {}}
    project, pop_names = None, None
    for scenario_key, _ in SCENARIOS:
        project, _, res = results[scenario_key]
        pop_names = list(project.data.pops.keys())
        data["model"][scenario_key] = {}
        for key, code_name, _, kind in OUTCOMES:
            if kind == "proportion":
                t, v = _prevalence_vals(res, pop_names, None)
            elif kind == "cascade":
                t, v = _cascade_vals(res, pop_names, None, code_name)
            elif kind == "incidence":
                t, v = _incidence_vals(res, pop_names)
            else:
                t, v = _series_vals(res, pop_names, None, code_name)
            data["model"][scenario_key][key] = {"t": t, "vals": v}

    for key, code_name, _, kind in OUTCOMES:
        data["data"][key] = _data_vals(project, None, code_name, kind) if kind not in ("cascade", "incidence") else None

    return data


def extract_population_breakdown(results):
    """Per-population-group series (model lines + any real data) for the breakdown
    figure: rows = a "Total" aggregate (pop=None) followed by every real population;
    columns = BREAKDOWN_COLUMNS."""
    project, _, _ = results[SCENARIOS[0][0]]
    pop_names = list(project.data.pops.keys())
    rows = [None] + pop_names

    out = {"rows": rows, "model": {}, "data": {}}
    for scenario_key, _ in SCENARIOS:
        _, _, res = results[scenario_key]
        out["model"][scenario_key] = {}
        for row in rows:
            out["model"][scenario_key][row] = {}
            for code_name, _, kind in BREAKDOWN_COLUMNS:
                if kind == "prevalence":
                    t, v = _prevalence_vals(res, pop_names, row)
                elif kind == "cascade":
                    t, v = _cascade_vals(res, pop_names, row, code_name)
                else:
                    t, v = _series_vals(res, pop_names, row, code_name)
                out["model"][scenario_key][row][code_name] = {"t": t, "vals": v}

    for row in rows:
        out["data"][row] = {}
        for code_name, _, kind in BREAKDOWN_COLUMNS:
            if kind == "cascade":
                out["data"][row][code_name] = None
            else:
                out["data"][row][code_name] = _data_vals(project, row, code_name, "proportion" if kind == "prevalence" else "count")

    return out


def _axis_formatter(kind):
    if kind in ("proportion", "prevalence", "cascade"):
        return FuncFormatter(lambda v, _: f"{v * 100:.1f}%")
    return FuncFormatter(lambda v, _: f"{v:,.0f}")


def plot_scenario_comparison(data, save_path):
    """Auto-plot every outcome in OUTCOMES as one multipanel figure, formatted for
    publication (Status Quo vs Action Scenario model lines, with real input data
    overlaid as scatter wherever available)."""
    ncols = 3
    nrows = -(-len(OUTCOMES) // ncols)  # ceil division
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 3.1 * nrows), squeeze=False)
    axes = axes.flatten()
    panel_letters = string.ascii_uppercase

    scenario_lines = {}
    data_handle = None
    for idx, (key, code_name, display_name, kind) in enumerate(OUTCOMES):
        ax = axes[idx]
        for scenario_key, label, color in [
            (SCENARIOS[0][0], SCENARIOS[0][1], COLOR_STATUS_QUO),
            (SCENARIOS[1][0], SCENARIOS[1][1], COLOR_ACTION),
        ]:
            (line,) = ax.plot(
                data["model"][scenario_key][key]["t"],
                data["model"][scenario_key][key]["vals"],
                label=label, color=color, linewidth=1.8, solid_capstyle="round",
            )
            scenario_lines[label] = line

        observed = data["data"].get(key)
        if observed is not None:
            obs_t, obs_v = observed
            data_handle = ax.scatter(
                obs_t, obs_v, s=22, facecolor=COLOR_DATA, edgecolor="white", linewidth=0.6,
                zorder=5, label="Data",
            )

        ax.set_title(f"{panel_letters[idx]}.  {display_name}", loc="left")
        ax.set_ylim(bottom=0)
        ax.set_xlim(1980, 2071)
        ax.yaxis.set_major_formatter(_axis_formatter(kind))
        ax.tick_params(axis="both", length=3, color="#4D4D4D")

    for ax in axes[len(OUTCOMES):]:
        fig.delaxes(ax)

    handles = [scenario_lines[label] for _, label in SCENARIOS]
    labels = [label for _, label in SCENARIOS]
    if data_handle is not None:
        handles.append(data_handle)
        labels.append("Data")
    fig.legend(handles, labels, loc="upper center", ncol=len(handles), frameon=False, bbox_to_anchor=(0.5, 1.02))

    fig.suptitle("HepAus model: Status Quo vs. Action Scenario, 1980–2071", y=1.06, fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 1))
    fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)


def plot_population_breakdown(breakdown, save_path):
    """One row per population group (plus a Total aggregate row), one column per
    BREAKDOWN_COLUMNS metric - Status Quo vs Action Scenario, with real data overlaid
    on population size and prevalence wherever the databook has it."""
    rows = breakdown["rows"]
    ncols = len(BREAKDOWN_COLUMNS)
    nrows = len(rows)
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.6 * ncols, 0.95 * nrows), squeeze=False)

    scenario_lines = {}
    data_handle = None
    for r, pop in enumerate(rows):
        is_total = pop is None
        for c, (code_name, col_label, kind) in enumerate(BREAKDOWN_COLUMNS):
            ax = axes[r][c]
            for scenario_key, label, color in [
                (SCENARIOS[0][0], SCENARIOS[0][1], COLOR_STATUS_QUO),
                (SCENARIOS[1][0], SCENARIOS[1][1], COLOR_ACTION),
            ]:
                series = breakdown["model"][scenario_key][pop][code_name]
                (line,) = ax.plot(
                    series["t"], series["vals"], color=color,
                    linewidth=1.8 if is_total else 1.1, label=label,
                )
                scenario_lines[label] = line

            observed = breakdown["data"][pop].get(code_name)
            if observed is not None:
                obs_t, obs_v = observed
                data_handle = ax.scatter(
                    obs_t, obs_v, s=8, facecolor=COLOR_DATA, edgecolor="white",
                    linewidth=0.3, zorder=5, label="Data",
                )

            ax.set_ylim(bottom=0)
            ax.set_xlim(1980, 2071)
            ax.yaxis.set_major_formatter(_axis_formatter(kind))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
            ax.tick_params(axis="both", length=2, labelsize=6.5)
            if r == 0:
                ax.set_title(col_label, fontsize=8.5, fontweight="bold")
            if r < nrows - 1:
                ax.set_xticklabels([])
            if c == 0:
                row_label = "TOTAL\n(Australia)" if is_total else pop
                ax.set_ylabel(
                    row_label, fontsize=7, rotation=0, ha="right", va="center",
                    labelpad=2, fontweight="bold" if is_total else "normal",
                )
            if is_total:
                ax.set_facecolor("#EEF3FB")
                for spine in ax.spines.values():
                    spine.set_visible(True)
                    spine.set_linewidth(1.3)
                    spine.set_color("#2A78D6")

    handles = [scenario_lines[label] for _, label in SCENARIOS]
    labels = [label for _, label in SCENARIOS]
    if data_handle is not None:
        handles.append(data_handle)
        labels.append("Data")
    fig.legend(handles, labels, loc="upper center", ncol=len(handles), frameon=False, bbox_to_anchor=(0.5, 1.012), fontsize=9)
    fig.suptitle("HepAus model by population group: Status Quo vs. Action Scenario", y=1.028, fontsize=11, fontweight="bold")

    fig.subplots_adjust(left=0.15, right=0.995, top=0.965, bottom=0.015, hspace=0.35, wspace=0.45)
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_population_composition(results, save_path, scenario_key=None):
    """Stacked area chart of hepb_pop (people living with HBV) over time, stacked by
    origin group (ATSI / Australian-born / overseas-born low-risk / overseas-born
    high-risk), for one scenario only (defaults to Status Quo)."""
    scenario_key = scenario_key or SCENARIOS[0][0]
    project, _, res = results[scenario_key]
    pop_names = list(project.data.pops.keys())

    tvec = None
    series = []
    for prefix, label, color in ORIGIN_GROUPS:
        group_pops = [p for p in pop_names if p.startswith(prefix)]
        t, v = _series_vals(res, group_pops, None, "hepb_pop")
        tvec = t
        series.append(v)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.stackplot(
        tvec, series, labels=[label for _, label, _ in ORIGIN_GROUPS],
        colors=[color for _, _, color in ORIGIN_GROUPS],
        edgecolor="white", linewidth=0.4,
    )
    ax.set_xlim(1980, 2071)
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_formatter(_axis_formatter("count"))
    ax.set_ylabel("People living with HBV")
    scenario_label = dict((k, lbl) for k, lbl in SCENARIOS)[scenario_key]
    ax.set_title(f"Population Living with HBV by Population Group ({scenario_label})", loc="left")
    ax.legend(loc="upper left", frameon=False)
    ax.tick_params(axis="both", length=3, color="#4D4D4D")

    fig.tight_layout()
    fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)


def print_summary(data, summary_year=2070.5):
    print(f"\nOutcome comparison at {summary_year}:")
    header = f"{'Outcome':<34}" + "".join(f"{label:>16}" for _, label in SCENARIOS)
    print(header)
    for key, code_name, display_name, _ in OUTCOMES:
        row = f"{display_name:<34}"
        for scenario_key, _ in SCENARIOS:
            t = data["model"][scenario_key][key]["t"]
            vals = data["model"][scenario_key][key]["vals"]
            val = vals[list(t).index(summary_year)] if summary_year in t else float("nan")
            row += f"{val:>16,.3f}" if val < 1 else f"{val:>16,.0f}"
        has_data = "yes" if data["data"].get(key) is not None else "no"
        row += f"   (data: {has_data})"
        print(row)


def main():
    print("Running Status Quo and Action scenarios...")
    results = run_hepaus_scenarios()

    out_dir = os.path.join(_get_github_folder(), "outputs")
    figures_dir = os.path.join(out_dir, "figures")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    pkl_path = os.path.join(out_dir, "hepaus_scenarios.pkl")
    sc.save(pkl_path, results)
    print(f"Saved raw scenario results to {pkl_path}")

    data = extract_outcomes(results)
    plot_path = os.path.join(figures_dir, "hepaus_scenarios_comparison.png")
    plot_scenario_comparison(data, plot_path)
    print(f"Saved comparison plot to {plot_path}")

    breakdown = extract_population_breakdown(results)
    breakdown_path = os.path.join(figures_dir, "hepaus_population_breakdown.png")
    plot_population_breakdown(breakdown, breakdown_path)
    print(f"Saved population breakdown plot to {breakdown_path}")

    composition_path = os.path.join(figures_dir, "hepaus_population_composition.png")
    plot_population_composition(results, composition_path)
    print(f"Saved population composition plot to {composition_path}")

    print_summary(data)


if __name__ == "__main__":
    main()
