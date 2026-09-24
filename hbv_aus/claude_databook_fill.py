"""
Fill hbv_hepaus_db.xlsx with plausible dummy data, consistent with the
hbv_fw_v2.1_autosave.xlsx framework structure, and save as claude_hbv_hepaus_db.xlsx.

Population naming: {origin}_{age}_{sex}
  origin: atsi (Aboriginal/Torres Strait Islander), ausb (Australian-born, non-Indigenous),
          lros (low-risk-country overseas-born), hros (high-risk-country overseas-born)
  age:    0-14, 15-64, 65+   (renamed from 15-54/65+ to match updated databook_gen.py)
  sex:    M, F

Two modelling rules applied here:
  - Births to lros/hros mothers enter the model as ausb_0-14 (born in Australia
    = Australian-born by definition, regardless of parental country of birth).
    atsi births stay atsi (Indigenous status is not birthplace-defined).
  - lros/hros populations are therefore only ever populated via immigration
    (at any age, including children), never via in-model births.
"""
import atomica as at
import numpy as np

FW_PATH = 'framework/hbv_fw_v2.1_autosave.xlsx'
DB_IN = 'databook/hbv_hepaus_db.xlsx'
DB_OUT = 'databook/claude_hbv_hepaus_db.xlsx'

F = at.ProjectFramework(FW_PATH)
D = at.ProjectData.from_spreadsheet(DB_IN, framework=F)

ORIGINS = ['atsi', 'ausb', 'lros', 'hros']
AGES = ['0-14', '15-64', '65+']
SEXES = ['M', 'F']

# Rename the old 15-54 band to 15-64 (folding the 55-64 decade in), matching
# the age_bins update in hbv_aus/databook_gen.py.
for origin in ORIGINS:
    for sex in SEXES:
        old = f'{origin}_15-54_{sex}'
        new = f'{origin}_15-64_{sex}'
        D.rename_pop(old, new, new)

pop_names = list(D.pops.keys())

def parse(pop):
    origin, age, sex = pop.split('_')
    return origin, age, sex

YEARS = [1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2030, 2040, 2050, 2060, 2071]

# ---------------------------------------------------------------------------
# Demographics - anchored to real ABS figures where available:
#   - Total ERP ~27.6M at Jun 2025 (ABS National, state and territory population)
#   - ATSI population 983,700 at Jun 2021, ~3.8% of total, growing ~2.0%/yr,
#     33.1% aged under 15, median age 24 (ABS ATSI population estimates/projections)
#   - Overseas-born 32.0% of population at 2025 (8.83M); country-of-birth split
#     bucketed here into "low HBV-prevalence" (lros, e.g. England, NZ, and other
#     Europe/Americas) vs "high HBV-prevalence" (hros, e.g. China, India,
#     Philippines, Vietnam, South Africa, Nepal, Sri Lanka, Malaysia, and other
#     Asia/Africa/Pacific) per ABS "Australia's population by country of birth"
#   - Overseas-born population skews strongly working-age: only ~4.7% are aged
#     0-14 (ABS: just 8.9% of Australia's 0-14 children were born overseas, out
#     of ~411,000 overseas-born children / 8.83M overseas-born total), median
#     age 43 vs 35 for Australia-born
#   - National age structure ~0-14:17.4%, 15-64:65.1%, 65+:17.5% (ABS)
# 1980 and 2071 anchors are NOT direct ABS figures (ABS doesn't publish this
# specific 4-group split historically or in projection) - they are informed
# extrapolations from known historical migration-pattern shifts (pre-1970s/80s
# migration was much more Europe-dominated; Asian-source migration expanded
# from the 1980s onward) and continuation of recent growth trends.
# NOTE: with lros/hros births now flowing to ausb (see b_rate below), these
# origin_share anchors describe the population BEFORE that internal transfer
# plays out - the simulated ausb/lros/hros split will drift from them over
# time as lros/hros can now only grow via immigration, not natural increase.
# ---------------------------------------------------------------------------
anchor_years = [1980, 2024, 2071]
origin_share_anchors = {
    #           1980    2024    2071
    # atsi 1980 share is backward-extrapolated from the real 2021 ABS figure
    # (983,700, ~2.0%/yr growth) rather than a historical census count, since
    # historical Indigenous census counts are heavily affected by undercounting
    # and by real increases in self-identification over time, not just births
    # and deaths - a true demographic backward-projection is more internally
    # consistent for seeding population dynamics than the recorded count.
    'atsi': [0.0295, 0.038, 0.058],
    'hros': [0.030, 0.190, 0.260],
    'lros': [0.080, 0.130, 0.150],
}
origin_share_anchors['ausb'] = [1 - sum(origin_share_anchors[o][i] for o in ('atsi', 'hros', 'lros')) for i in range(3)]

total_pop_anchors = [14_700_000.0, 27_600_000.0, 46_000_000.0]

age_share = {
    'atsi': {'0-14': 0.33, '15-64': 0.62, '65+': 0.05},
    'ausb': {'0-14': 0.228, '15-64': 0.582, '65+': 0.190},
    'lros': {'0-14': 0.047, '15-64': 0.793, '65+': 0.160},
    'hros': {'0-14': 0.047, '15-64': 0.793, '65+': 0.160},
}
# Fertility rate = births per woman per year in the childbearing pool. The
# pool is now 15-64 (50 years) rather than 15-54 (40 years), but real fertility
# is concentrated in ~15-49; scaled down by ~0.8 (35/50 vs previous 35/40
# effective coverage) so the ABSOLUTE number of births stays realistic rather
# than being inflated by including more non-fertile older women in the pool.
fertility_rate = {'atsi': 0.075, 'ausb': 0.036, 'lros': 0.036, 'hros': 0.040}

def origin_pop(origin, year):
    total = np.interp(year, anchor_years, total_pop_anchors)
    share = np.interp(year, anchor_years, origin_share_anchors[origin])
    return total * share

def pop_size(pop, year):
    origin, age, sex = parse(pop)
    return origin_pop(origin, year) * age_share[origin][age] * 0.5

def set_series(name, pop, years, vals, units=None):
    ts = D.tdve[name].ts[pop]
    if units is None:
        units = ts.units
    ts.t = []
    ts.vals = []
    ts.assumption = None
    ts.units = units
    ts.insert(years, vals)

def set_flat(name, pop, val, units=None):
    ts = D.tdve[name].ts[pop]
    if units is None:
        units = ts.units
    ts.units = units
    ts.assumption = val

# total_pop
for pop in pop_names:
    vals = [pop_size(pop, y) for y in YEARS]
    set_series('total_pop', pop, YEARS, vals)

# j_init: everyone enters the model at sim_start via this junction
for pop in pop_names:
    set_flat('j_init', pop, pop_size(pop, 1980))

# b_rate: annual births "by sex" - the transition wiring (src_births -> sus
# etc.) is *within* the same population as b_rate itself, so despite the
# name this must be attached to the *_0-14_{M,F} (newborn) populations, not
# the mother's *_15-64_F population. Sex ratio at birth ~105:100 (M:F), i.e.
# ~51.2% male / 48.8% female.
#
# Births to lros/hros mothers are credited to ausb_0-14 (born in Australia =
# Australian-born, regardless of parents' country of birth); atsi and ausb
# mothers' births stay in their own origin group. lros/hros therefore have
# b_rate = 0 everywhere and can only grow via imig_rate.
sex_at_birth = {'M': 0.512, 'F': 0.488}
for pop in pop_names:
    origin, age, sex = parse(pop)
    if age != '0-14':
        set_flat('b_rate', pop, 0.0)
        continue

    if origin == 'ausb':
        contributing_mothers = ['ausb', 'lros', 'hros']
    elif origin == 'atsi':
        contributing_mothers = ['atsi']
    else:  # lros, hros: no direct births - only entry is via immigration
        set_flat('b_rate', pop, 0.0)
        continue

    vals = [
        sum(pop_size(f'{mother_origin}_15-64_F', y) * fertility_rate[mother_origin] for mother_origin in contributing_mothers) * sex_at_birth[sex]
        for y in YEARS
    ]
    set_series('b_rate', pop, YEARS, vals)

# acm: all-cause mortality probability per year. '65+' is now a genuine 65+
# band (55-64 has moved into '15-64'), so uses a proper 65+-specific rate.
# '15-64' blends the low-mortality 15-54 range with higher-mortality 55-64.
acm_base = {'0-14': 0.001, '15-64': 0.0035, '65+': 0.040}
acm_adj = {'atsi': 1.4, 'ausb': 1.0, 'lros': 1.0, 'hros': 1.05}
for pop in pop_names:
    origin, age, sex = parse(pop)
    set_flat('acm', pop, min(acm_base[age] * acm_adj[origin], 0.95))

# emigration / immigration - rates chosen so the resulting overall population
# growth roughly tracks ABS ERP growth (~1.4%/yr 1980-2024, moderating to
# ~1.1%/yr through to 2071): crude birth rate ~1.2-1.3%, crude death rate
# ~0.6-0.7%, net overseas migration contributing the remaining ~0.5-0.9%/yr
# (net migration has historically been Australia's larger growth driver
# alongside natural increase; ABS reports 445,600 net overseas migration in
# 2024 alone, ~1.6% of population, though that is an unusually high year).
emig_rate_val = {'atsi': 0.001, 'ausb': 0.004, 'lros': 0.004, 'hros': 0.004}
for pop in pop_names:
    origin, age, sex = parse(pop)
    set_flat('emig_rate', pop, emig_rate_val[origin])

# imig_rate allows only 'N.A.' units and feeds a 'src_immig -> comp' NUMBER-type
# inflow (same architecture as b_rate/src_births), i.e. it must be an ABSOLUTE
# annual headcount of arrivals, not a per-capita rate. lros/hros now depend
# ENTIRELY on this for growth (births no longer feed them - see b_rate above),
# so a flat rate can't keep pace with origin_share_anchors targets that grow
# FASTER than total population (hros: 19%->26% of pop between 2024 and 2071).
# Both total migration intake and its origin mix are therefore time-varying:
# total volume increases over time (Australia's real migration intake has
# grown substantially - 220,000/yr was already a long-run-average
# approximation of a rising trend), and the mix shifts toward hros-origin
# countries over time (matching the real shift in Australia's migration
# source countries from Europe-dominated pre-1980s towards Asia/Africa in
# recent decades - see the country-of-birth data referenced above).
migration_years = [1980, 2000, 2024, 2050, 2071]
total_net_migration_vals = [100_000.0, 180_000.0, 260_000.0, 380_000.0, 480_000.0]
imig_origin_share_vals = {
    'atsi': [0.0, 0.0, 0.0, 0.0, 0.0],
    'ausb': [0.10, 0.07, 0.04, 0.02, 0.02],
    'lros': [0.50, 0.42, 0.36, 0.32, 0.28],
    'hros': [0.40, 0.51, 0.60, 0.66, 0.70],
}
imig_age_share = {'0-14': 0.10, '15-64': 0.85, '65+': 0.05}

def imig_rate_at(origin, age, year):
    total = np.interp(year, migration_years, total_net_migration_vals)
    share = np.interp(year, migration_years, imig_origin_share_vals[origin])
    return total * share * imig_age_share[age] * 0.5

for pop in pop_names:
    origin, age, sex = parse(pop)
    vals = [imig_rate_at(origin, age, y) for y in YEARS]
    set_series('imig_rate', pop, YEARS, vals)

imig_prev_val = {'atsi': 0.01, 'ausb': 0.005, 'lros': 0.01, 'hros': 0.07}
imig_vax_val = {'atsi': 0.3, 'ausb': 0.4, 'lros': 0.5, 'hros': 0.25}
imig_clr_val = {'atsi': 0.1, 'ausb': 0.1, 'lros': 0.08, 'hros': 0.25}
for pop in pop_names:
    origin, age, sex = parse(pop)
    set_flat('imig_prev', pop, imig_prev_val[origin])
    set_flat('imig_vax_prop', pop, imig_vax_val[origin])
    set_flat('imig_clr_prop', pop, imig_clr_val[origin])

# hepb_prev: background prevalence input used for horizontal FOI weighting - mild decline post-2000
hepb_prev_base = {'atsi': 0.015, 'ausb': 0.003, 'lros': 0.008, 'hros': 0.065}
hepb_prev_years = [1980, 2000, 2020, 2040, 2071]
hepb_prev_mult = [1.2, 1.0, 0.85, 0.70, 0.6]
for pop in pop_names:
    origin, age, sex = parse(pop)
    vals = [hepb_prev_base[origin] * m for m in hepb_prev_mult]
    set_series('hepb_prev', pop, hepb_prev_years, vals)

# epos_total: HBeAg positive fraction of all HBV, by age
epos_base = {'0-14': 0.6, '15-64': 0.22, '65+': 0.08}
for pop in pop_names:
    origin, age, sex = parse(pop)
    set_flat('epos_total', pop, epos_base[age])

# ---------------------------------------------------------------------------
# Initial disease-stage distribution (sim-start conditions)
# ---------------------------------------------------------------------------
chb_init_val = {'atsi': 0.015, 'ausb': 0.004, 'lros': 0.01, 'hros': 0.07}
init_clr_val = {'atsi': 0.2, 'ausb': 0.1, 'lros': 0.2, 'hros': 0.35}
for pop in pop_names:
    origin, age, sex = parse(pop)
    set_flat('chb_init', pop, chb_init_val[origin])
    set_flat('init_clr', pop, init_clr_val[origin])
    set_flat('eag_init', pop, 0.3)
    set_flat('cc_prop', pop, 0.05)
    set_flat('dc_prop', pop, 0.01)
    set_flat('hcc_prop', pop, 0.005)
    set_flat('ei_prop', pop, 0.7)
    set_flat('si_prop', pop, 0.6)

# ---------------------------------------------------------------------------
# Natural history / clinical parameters (mostly universal, some age-varying)
# ---------------------------------------------------------------------------
flat_universal = {
    'tx_chb_cc': 0.65, 'tx_chb_hcc': 0.6, 'tx_cir_hcc': 0.55, 'tx_cc_dc': 0.7,
    'tx_dc_death': 0.5, 'dc_cc_t': 0.15, 'sh_si_t': 0.4,
    'si_sh': 0.06, 'sh_si': 0.08, 'eh_cc_in': 0.02, 'sh_cc_in': 0.015,
    'cc_dc': 0.04, 'ei_hcc_in': 0.003, 'act_hcc_in': 0.005, 'si_hcc_in': 0.004,
    'cir_hcc_in': 0.03, 'y_hcc': 1.0, 'y_cc': 1.0,
    'nullipar': 0.42, 've_hepbd_lvl': 0.7, 've_hepbd_hvl': 0.5, 've_hepbd_pap': 0.9,
    've_hepb3': 0.95, 'mtct_chb': 0.85, 'trisk_low': 0.02, 'trisk_high': 0.2,
    'hvl_epos': 0.75, 'hvl_eneg': 0.1, 'ts_rate': 0.03, 'ltf_rate': 0.07,
}
for name, val in flat_universal.items():
    for pop in pop_names:
        set_flat(name, pop, val)

age_varying = {
    'eh_si': {'0-14': 0.02, '15-64': 0.09, '65+': 0.05},
    'acu_chb': {'0-14': 0.3, '15-64': 0.05, '65+': 0.05},
    'acu_dth': {'0-14': 0.001, '15-64': 0.001, '65+': 0.01},
    'ei_eh': {'0-14': 0.02, '15-64': 0.09, '65+': 0.1},
    'dc_dth': {'0-14': 0.15, '15-64': 0.19, '65+': 0.25},
    'hcc_dth': {'0-14': 0.5, '15-64': 0.57, '65+': 0.65},
}
for name, lookup in age_varying.items():
    for pop in pop_names:
        origin, age, sex = parse(pop)
        set_flat(name, pop, lookup[age])

foi_cal_val = {'atsi': 0.04, 'ausb': 0.01, 'lros': 0.015, 'hros': 0.02}
for pop in pop_names:
    origin, age, sex = parse(pop)
    set_flat('foi_cal', pop, foi_cal_val[origin])

# ---------------------------------------------------------------------------
# Antenatal / vertical transmission block
# ---------------------------------------------------------------------------
anc_scr_years = [1980, 1990, 2000, 2010, 2071]
anc_scr_vals = [0.0, 0.5, 0.85, 0.95, 0.95]
preg_ltc_years = [1980, 1995, 2010, 2071]
preg_ltc_vals = [0.2, 0.5, 0.85, 0.85]
vax_years = [1980, 1988, 1995, 2000, 2005, 2071]
# NOTE: hepbd_cov must not be exactly 0 while b_rate>0 - adj_hepbd_cov divides
# by b_rate*hepbd_cov and hits the same sdiv nonzero-numerator/zero-denominator
# issue found earlier for treat_rate/ltc_rate. A tiny floor avoids it.
hepbd_cov_vals = [0.001, 0.001, 0.5, 0.8, 0.95, 0.95]
hepb3_cov_vals = [0.001, 0.001, 0.6, 0.85, 0.95, 0.95]

for pop in pop_names:
    origin, age, sex = parse(pop)
    set_series('anc_scr', pop, anc_scr_years, anc_scr_vals)
    set_series('preg_ltc_rate', pop, preg_ltc_years, preg_ltc_vals)
    set_series('hepbd_cov', pop, vax_years, hepbd_cov_vals)
    set_series('hepb3_cov', pop, vax_years, hepb3_cov_vals)
    if age == '15-64' and sex == 'F':
        vals = [pop_size(pop, y) * fertility_rate[origin] for y in YEARS]
        set_series('pregs', pop, YEARS, vals)
    else:
        set_flat('pregs', pop, 0.0)

# ---------------------------------------------------------------------------
# Diagnosis / linkage / treatment cascade targets
# NOTE: onset years are staggered (diag before ltc before treat) so that
# diag_n / ltc_n have a chance to become non-zero before ltc_measure / treat_measure
# demand anything of them - this avoids triggering the divide-by-zero issue
# identified in treat_rate / ltc_rate (see conversation).
# ---------------------------------------------------------------------------
ceiling_mult = {'ausb': 1.0, 'atsi': 0.85, 'lros': 0.75, 'hros': 0.7}

diag_years = [1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2071]
diag_base =  [0.0,  0.0,  0.05, 0.15, 0.30, 0.45, 0.55, 0.65, 0.72, 0.78]

ltc_years =  [1980, 1990, 1992, 1995, 2000, 2005, 2010, 2015, 2020, 2071]
ltc_base =   [0.0,  0.0,  0.02, 0.08, 0.18, 0.28, 0.36, 0.42, 0.47, 0.52]

treat_years = [1980, 1995, 2000, 2002, 2005, 2010, 2015, 2020, 2071]
treat_base =  [0.0,  0.0,  0.0,  0.02, 0.08, 0.15, 0.22, 0.28, 0.33]

for pop in pop_names:
    origin, age, sex = parse(pop)
    m = ceiling_mult[origin]
    diag_vals = [v * m for v in diag_base]
    ltc_vals = [v * m for v in ltc_base]
    treat_vals = [v * m for v in treat_base]

    set_series('diag_measure', pop, diag_years, diag_vals)
    set_series('ltc_measure', pop, ltc_years, ltc_vals)
    set_series('treat_measure', pop, treat_years, treat_vals)

    # observed comparators mirror the fit targets
    set_series('diagnosed_meas', pop, diag_years, diag_vals)
    set_series('linked_meas', pop, ltc_years, ltc_vals)
    set_series('treat_cov', pop, treat_years, treat_vals)

# ---------------------------------------------------------------------------
# Interactions: mx_mtct (mother-to-child, within origin group only,
# from *_15-64_F to *_0-14_{M,F} of the same origin), mx_horiz (horizontal,
# within origin group only, any age/sex to any age/sex of the same origin)
# ---------------------------------------------------------------------------
mx_mtct = next(ip for ip in D.interpops if ip.code_name == 'mx_mtct')
mx_horiz = next(ip for ip in D.interpops if ip.code_name == 'mx_horiz')

mx_mtct.ts.clear()
for origin in ORIGINS:
    src = f'{origin}_15-64_F'
    for age in AGES:
        for sex in SEXES:
            if age != '0-14':
                continue
            dst = f'{origin}_{age}_{sex}'
            mx_mtct.ts[(src, dst)] = at.TimeSeries(assumption=1.0, units='N.A.')

mx_horiz.ts.clear()
for origin in ORIGINS:
    group_pops = [f'{origin}_{age}_{sex}' for age in AGES for sex in SEXES]
    for src in group_pops:
        for dst in group_pops:
            mx_horiz.ts[(src, dst)] = at.TimeSeries(assumption=1.0, units='N.A.')

# ---------------------------------------------------------------------------
# Transfers: aging within each population (origin + sex fixed), 0-14 -> 15-64 -> 65+
# ---------------------------------------------------------------------------
age_transfer = D.transfers[0]
age_transfer.ts.clear()
age_rate = {'0-14': 1 / 15, '15-64': 1 / 50}
for origin in ORIGINS:
    for sex in SEXES:
        age_transfer.ts[(f'{origin}_0-14_{sex}', f'{origin}_15-64_{sex}')] = at.TimeSeries(assumption=age_rate['0-14'], units='Rate (per year)')
        age_transfer.ts[(f'{origin}_15-64_{sex}', f'{origin}_65+_{sex}')] = at.TimeSeries(assumption=age_rate['15-64'], units='Rate (per year)')

# ---------------------------------------------------------------------------
# Validate and save
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    D.validate(F)
    print('VALIDATION PASSED')
    D.save(DB_OUT)
    print('Saved to', DB_OUT)
