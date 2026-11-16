"""
Packed-column CO2 absorption lab: analysis and findings page, in one file.

Objectives
  1. Compare measured gas-phase pressure drop with the Strigle GPDC prediction.
  2. Determine the overall gas-phase coefficient K_OG*a for CO2 into dilute NaOH.
  3. Show how NaOH concentration affects the absorption driving force.

Inputs (read only, never modified)
  <data-dir>/Gas Absorption Data.xlsx   lab workbook (Summary, run sheets, GPDC Comparison)
  <data-dir>/GPDC-1.xlsx                digitized Strigle GPDC curves (Flv, CP, Dp in in H2O/ft)

Outputs (in --out-dir)
  results.json   every computed number
  index.html     self-contained findings page (charts, tables, literature comparison)

Usage
  python gas_absorption_project.py
  python gas_absorption_project.py --data-dir "C:/Users/elmho/Downloads" --out-dir .

Requires: pandas, numpy, openpyxl   (pip install pandas numpy openpyxl)
"""
import argparse
import json
import math
import re
import statistics as st
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# =============================================================================
# Inputs and constants
# =============================================================================
D_IN, Z_IN = 4.0, 48.0                       # column ID and packed height, in (lab group)
A_FT2 = math.pi * (D_IN / 12) ** 2 / 4       # cross-section, ft^2
A_M2 = A_FT2 * 0.09290304
Z_FT, Z_M = Z_IN / 12, Z_IN * 0.0254
FP = 200.0                                   # ft^-1, 1/2-in ceramic Intalox (Strigle / Kister & Gill)
T_C = 22.0                                   # assumed lab temperature, C
P_ATM = 1.0
V_STD = 0.7302 * 529.67                      # ft^3/lbmol at 70 F, 14.7 psia (rotameter standard)
RHO_AIR = 28.96 / V_STD                      # lb/ft^3
RHO_L = 62.30                                # lb/ft^3, water at 22 C
NU_L = 0.96                                  # cSt, water at 22 C
GAL_LB = 8.33 * RHO_L / 62.37                # lb per US gallon
C_NAOH = 0.076                               # mol/L, HCl titration
HE_CO2 = 29.4e-3                             # atm m^3/mol, Henry's constant CO2-water
H_X = HE_CO2 * 55.4e3                        # atm, mole-fraction basis (~1630 atm)
D_CO2 = 1.9e-9                               # m^2/s
D_OH = 5.27e-9                               # m^2/s
K2 = 10 ** (11.895 - 2382 / (T_C + 273.15))  # L/(mol s), Pohorecki & Moniuk, CO2 + OH-
DP_FLOOD = 0.12 * FP ** 0.7                  # Kister et al. 2007 Eq. 4, in H2O/ft (off-chart for Fp=200)
DP_INCIPIENT_FLOOD = 1.5                     # in H2O/ft, Perry's / Strigle

G_NAOH_SCFM, L_NAOH_GPM = 7.5, 1.35          # NaOH run settings (10-15)
STEADY_WINDOW = (60, 680)                    # s, liquid flowing, after start-up

# Literature values (Perry's 7th ed. Table 14-8, CO2-air-1 N NaOH, 25 % conversion)
PERRY_T148 = {"25 mm ceramic Intalox": 45.1, "50 mm ceramic Intalox": 30.1,
              "25 mm ceramic Raschig": 37.0, "25 mm metal Pall": 49.6, "25 mm metal IMTP": 54.8}
DRY_AREA_HALF_IN = 623                       # m^2/m^3, 13 mm ceramic Intalox (Table 14-7b)
DRY_AREA_25MM = 256                          # m^2/m^3, 25 mm ceramic Intalox


def scfm_to_lbmol_h(q):
    return q * 60 / V_STD


# =============================================================================
# Objective 1: pressure drop vs GPDC
# =============================================================================
class GPDC:
    """Strigle GPDC chart: interpolate dP/Z (in H2O/ft) at (F_LV, CP)."""

    def __init__(self, df):
        self.df = df
        self.curves = sorted(df.Dp.unique())

    def curve_cp(self, dp, flv):
        c = self.df[self.df.Dp == dp].sort_values("Flv")
        return float(np.interp(math.log10(flv), np.log10(c.Flv), c.CP))

    def dp(self, flv, cp):
        cps = np.array([self.curve_cp(d, flv) for d in self.curves])
        logd = np.log10(self.curves)
        if cp <= cps[0]:      # below lowest curve: dP ~ CP^2 (pre-loading)
            return self.curves[0] * (cp / cps[0]) ** 2, "below 0.05 curve (extrapolated)"
        if cp >= cps[-1]:     # above 1.5 curve: flooding region
            slope = (logd[-1] - logd[-2]) / (cps[-1] - cps[-2])
            return 10 ** (logd[-1] + slope * (cp - cps[-1])), "above 1.5 curve (flood region, extrapolated)"
        return float(10 ** np.interp(cp, cps, logd)), "interpolated"

    def curves_for_chart(self, every=4):
        return {str(d): self.df[self.df.Dp == d].sort_values("Flv")[["Flv", "CP"]].iloc[::every].round(4).values.tolist()
                for d in self.curves}


def hydraulics(G_scfm, L_gpm, y_co2=0.0):
    """GPDC coordinates for a run; column assumed at ~standard conditions."""
    mw = 28.96 * (1 - y_co2) + 44.01 * y_co2
    rho_g = mw / V_STD
    G_lb = G_scfm * 60 * rho_g                # lb/h
    L_lb = L_gpm * 60 * GAL_LB
    us = G_scfm / 60 / A_FT2                  # ft/s
    cs = us * math.sqrt(rho_g / (RHO_L - rho_g))
    cp = cs * FP ** 0.5 * NU_L ** 0.05
    flv = (L_lb / G_lb) * math.sqrt(rho_g / RHO_L)
    return dict(us=us, cs=cs, cp=cp, flv=flv, G_flux=G_lb / A_FT2, L_flux=L_lb / A_FT2,
                L_gpm_ft2=L_gpm / A_FT2)


def analyze_pressure(dp_rows, sheet_cmp, gpdc):
    rows = []
    for _, r in dp_rows.iterrows():
        h = hydraulics(r.G_scfm, r.L_gpm)
        pred_ft, region = gpdc.dp(h["flv"], h["cp"])
        meas_ft = r.dP / Z_FT
        sc = sheet_cmp[(sheet_cmp.Liquid_GPM == r.L_gpm) & (sheet_cmp.Gas_SCFM == r.G_scfm)]
        sheet_pred = sc.GPDC_dP_inH2O.iloc[0] if len(sc) and pd.notna(sc.GPDC_dP_inH2O.iloc[0]) else None
        rows.append(dict(L_gpm=r.L_gpm, G_scfm=r.G_scfm, dP_meas_in=r.dP, sd_dP_in=r.sd_dP,
                         dP_meas_ft=meas_ft, sd_ft=r.sd_dP / Z_FT, dP_pred_ft=pred_ft,
                         dP_pred_in=pred_ft * Z_FT, region=region,
                         dev_pct=100 * (meas_ft - pred_ft) / pred_ft,
                         sheet_pred_in=sheet_pred, **h))
    pdf = pd.DataFrame(rows).sort_values(["L_gpm", "G_scfm"])

    flood = []                                # gas rate where GPDC reaches incipient flood
    for L in sorted(pdf.L_gpm.unique()):
        for g in np.arange(2, 40, 0.05):
            h = hydraulics(g, L)
            if gpdc.dp(h["flv"], h["cp"])[0] >= DP_INCIPIENT_FLOOD:
                flood.append(dict(L_gpm=L, G_flood_scfm=round(g, 2)))
                break
    return pdf, pd.DataFrame(flood)


# =============================================================================
# Objective 2: analyzer cross-calibration and K_OG*a
# =============================================================================
def calibrate_analyzers(dp_rows, co2_rows):
    """Fit out_reading = a + b*in_reading on runs where outlet should equal inlet."""
    zero = dp_rows[["yin", "yout"]]                                  # 9-17: no CO2 fed
    water = co2_rows[co2_rows.sheet.str.contains("Water")].copy()
    x = np.concatenate([zero.yin, water.yin])
    y = np.concatenate([zero.yout, water.yout])
    b, a = np.polyfit(x, y, 1)
    resid = y - (a + b * x)
    s_res = float(np.sqrt(np.sum(resid ** 2) / (len(x) - 2)))

    water["ratio"] = water.yout / water.yin
    A_fac = []                                # physical absorption capacity of water
    for _, r in water.iterrows():
        Lm = r.L_gpm * 3.785 / 18.02                       # kmol/min water
        Gm = scfm_to_lbmol_h(r.G_scfm) * 0.4536 / 60       # kmol/min gas
        A_fac.append(Lm / (H_X / P_ATM * Gm))
    water["A_factor"] = A_fac
    water["max_removal_pct"] = 100 * water.A_factor        # upper bound for A << 1

    points = [dict(x=float(xx), y=float(yy), kind=("no CO2 (9-17)" if i < len(zero) else "water"))
              for i, (xx, yy) in enumerate(zip(x, y))]
    return dict(a=a, b=b, s_res=s_res, n=len(x), points=points), water


class KogModel:
    """Dilute absorption with y* = 0 (excess OH-): NOG = ln(y1/y2), KOGa = GM*NOG/(Z*P)."""

    def __init__(self):
        self.Gm_kmol_h = scfm_to_lbmol_h(G_NAOH_SCFM) * 0.4536
        self.G_flux = self.Gm_kmol_h / A_M2            # kmol/(h m^2)
        self.L_L_min = L_NAOH_GPM * 3.785

    def kog(self, yin, yout):
        nog = math.log(yin / yout)
        hog = Z_M / nog
        return nog, hog, self.G_flux / (hog * P_ATM)


def monte_carlo_kog(model, cal, yin_s, yout_s, corrected, rng, n=4000):
    """Uncertainty: analyzer time variation (~30-s autocorrelation), calibration fit, +/-3 % gas flow."""
    m = len(yin_s)
    out = []
    for _ in range(n):
        yi = rng.normal(yin_s.mean(), yin_s.std(ddof=1) / math.sqrt(m / 30))
        yo = rng.normal(yout_s.mean(), yout_s.std(ddof=1) / math.sqrt(m / 30))
        if corrected:                     # back to raw reading, re-apply a perturbed calibration
            bb = rng.normal(cal["b"], 0.015 * cal["b"])
            aa = rng.normal(cal["a"], cal["s_res"] / 2)
            yo = (yo * cal["b"] + cal["a"] - aa) / bb
        g = rng.normal(1, 0.03)
        if yo <= 0 or yi <= yo:
            continue
        out.append(g * model.kog(yi, yo)[2])
    return float(np.mean(out)), float(np.std(out))


def analyze_naoh_run(naoh_ts, cal, rng):
    model = KogModel()
    t = naoh_ts["Elapsed Time (s)"].values
    yin_raw = naoh_ts["CO2 In (%)"].values
    yout_raw = naoh_ts["CO2 Out (%)"].values
    dp_ts = naoh_ts["Pressure Drop (in. H2O)"].values
    yout_cor = (yout_raw - cal["a"]) / cal["b"]
    steady = (t >= STEADY_WINDOW[0]) & (t < STEADY_WINDOW[1])

    res = {}
    for label, yo_series, corr in [("raw", yout_raw, False), ("corrected", yout_cor, True)]:
        yi_m, yo_m = yin_raw[steady].mean(), yo_series[steady].mean()
        nog, hog, k = model.kog(yi_m, yo_m)
        _, mc_sd = monte_carlo_kog(model, cal, yin_raw[steady], yo_series[steady], corr, rng)
        co2_abs = model.Gm_kmol_h * 1000 / 60 * (yi_m - yo_m) / 100     # mol/min
        res[label] = dict(yin=yi_m, yout=yo_m, NOG=nog, HOG_m=hog, KOGa=k, KOGa_sd=mc_sd,
                          removal_pct=100 * (1 - yo_m / yi_m), co2_abs_mol_min=co2_abs,
                          naoh_conv_pct=100 * 2 * co2_abs / (model.L_L_min * C_NAOH))

    win = (t // 10) * 10                      # 10-s averages for the time-series chart
    ts = (pd.DataFrame(dict(t=win, yin=yin_raw, yout=yout_raw, yout_cor=yout_cor, dp=dp_ts))
          .groupby("t").mean().reset_index())
    ts["k_raw"] = [model.kog(r.yin, r.yout)[2] if r.yin > r.yout > 0 else None for r in ts.itertuples()]
    ts["k_cor"] = [model.kog(r.yin, r.yout_cor)[2] if r.yin > r.yout_cor > 0 else None for r in ts.itertuples()]
    return model, res, ts


# =============================================================================
# Objective 3: NaOH concentration and driving force
# =============================================================================
def analyze_concentration(model, kc):
    oh_top = C_NAOH
    oh_bot = C_NAOH * (1 - kc["naoh_conv_pct"] / 100)
    oh_mean = (oh_top + oh_bot) / 2

    # reaction regime at the measured conditions
    k1 = K2 * oh_mean
    sqrt_k1d = math.sqrt(k1 * D_CO2)                       # m/s = E*kL0 (fast pseudo-first-order)
    ha = sqrt_k1d / 1.0e-4                                 # kL0 ~ 1e-4 m/s
    c_ai = (kc["yin"] / 100) * P_ATM / HE_CO2 / 1000       # mol/L at gas inlet (upper bound)
    e_inf = 1 + (D_OH * oh_mean) / (2 * D_CO2 * c_ai)
    a_eff = (kc["KOGa"] * 1000 / 3600) / (sqrt_k1d / HE_CO2)   # m^2/m^3

    kog_phys = kc["KOGa"] / math.sqrt(1 + ha ** 2)         # same packing and flows, E = 1
    co2_in = model.Gm_kmol_h * 1000 / 60 * kc["yin"] / 100  # mol/min
    A_w = model.L_L_min / 18.02 / (H_X / P_ATM * model.Gm_kmol_h / 60)

    def kog_from_oh(oh_m):
        """E = sqrt(1 + Ha^2), Ha ~ sqrt([OH-]); tends to physical absorption as [OH-] -> 0."""
        ha_c = ha * math.sqrt(max(oh_m, 0) / oh_mean)
        return kog_phys * math.sqrt(1 + ha_c ** 2)

    def predicted_kog(c_feed):
        if c_feed <= 0:
            return kog_phys
        conv = 0.3
        for _ in range(50):                                # iterate hydroxide conversion
            kog = kog_from_oh(c_feed * (1 - conv / 2))
            nog = kog * Z_M * P_ATM / model.G_flux
            co2 = co2_in * (1 - math.exp(-nog))
            conv_new = min(0.999, 2 * co2 / (model.L_L_min * c_feed))
            if abs(conv_new - conv) < 1e-6:
                break
            conv = 0.5 * conv + 0.5 * conv_new
        return kog

    def predicted_removal(c_feed):
        """Rate-limited removal capped by NaOH supply plus water's physical capacity."""
        nog = predicted_kog(c_feed) * Z_M * P_ATM / model.G_flux
        stoich = (model.L_L_min * c_feed / 2) / co2_in
        return 100 * min(1 - math.exp(-nog), stoich + A_w)

    grid = [0, 0.002, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.076, 0.09, 0.1, 0.125,
            0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.875, 1.0]
    conc_model = [dict(C=c, KOGa=round(predicted_kog(c), 2), removal=round(predicted_removal(c), 2)) for c in grid]

    # driving force y - y* along the bed: NaOH (y* = 0, measured NOG) vs water (predicted)
    z = np.linspace(0, 1, 21)
    drive_naoh = [dict(z=float(zz), y=float(v), ystar=0.0, drive=float(v))
                  for zz, v in zip(z, kc["yin"] * np.exp(-kc["NOG"] * z))]
    nog_w = kog_phys * Z_M * P_ATM / model.G_flux
    yin_w = kc["yin"]

    def water_profile(y_top, n=400):
        """Integrate dy/dz = NOG*(y - y*) from top (z=1) down; y* = (y - y_top)/A."""
        y, prof = y_top, [(1.0, y_top, 0.0)]
        for i in range(n):
            y += nog_w * (y - (y - y_top) / A_w) / n
            prof.append((1 - (i + 1) / n, y, (y - y_top) / A_w))
        return prof

    lo, hi = yin_w * (1 - A_w), yin_w
    for _ in range(80):                                    # shoot on y_top so y(bottom) = y_in
        mid = (lo + hi) / 2
        if water_profile(mid)[-1][1] > yin_w:
            hi = mid
        else:
            lo = mid
    wp = sorted(water_profile((lo + hi) / 2)[::20])
    drive_water = [dict(z=zz, y=yy, ystar=ys, drive=yy - ys) for zz, yy, ys in wp]

    return dict(
        kinetics=dict(oh_top=oh_top, oh_bot=oh_bot, oh_mean=oh_mean, k1=k1, Ha=ha, E_inf=e_inf,
                      a_eff=a_eff, a_dry=DRY_AREA_HALF_IN),
        conc_model=conc_model,
        c_stoich_full=2 * co2_in / model.L_L_min,
        perry=PERRY_T148,
        perry_scaled=PERRY_T148["25 mm ceramic Intalox"] * math.sqrt(oh_mean / 0.875),
        driving=dict(naoh=drive_naoh, water=drive_water, A_water=A_w, kog_water=kog_phys,
                     nog_water=nog_w, water_removal_pct=100 * (1 - wp[-1][1] / yin_w)),
    )


# =============================================================================
# Run the analysis
# =============================================================================
def load_summary(wb):
    summ = pd.read_excel(wb, sheet_name="Summary", header=None)
    dp_rows = summ.iloc[1:19, :9].copy()
    dp_rows.columns = ["sheet", "G_scfm", "L_gpm", "yin", "yout", "sd_out", "dP", "sd_dP", "sd_in"]
    co2_rows = summ.iloc[22:32, :9].copy()
    co2_rows.columns = ["sheet", "G_scfm", "L_gpm", "yin", "yout", "sd_out", "sd_in", "dP", "sd_dP"]
    for df in (dp_rows, co2_rows):
        for c in df.columns[1:]:
            df[c] = df[c].astype(float)
    return dp_rows, co2_rows


def run_analysis(data_dir):
    wb = data_dir / "Gas Absorption Data.xlsx"
    gpdc = GPDC(pd.read_excel(data_dir / "GPDC-1.xlsx"))
    rng = np.random.default_rng(42)

    dp_rows, co2_rows = load_summary(wb)
    sheet_cmp = pd.read_excel(wb, sheet_name="GPDC Comparison")
    pdf, flood = analyze_pressure(dp_rows, sheet_cmp, gpdc)

    cal, water = calibrate_analyzers(dp_rows, co2_rows)
    model, kog, ts = analyze_naoh_run(pd.read_excel(wb, sheet_name="10-15 NaOH CO2"), cal, rng)
    conc = analyze_concentration(model, kog["corrected"])

    return dict(
        inputs=dict(D_in=D_IN, Z_in=Z_IN, area_ft2=A_FT2, Fp=FP, T_C=T_C, C_NaOH=C_NAOH, nu_cSt=NU_L,
                    rho_L=RHO_L, rho_air=RHO_AIR, He=HE_CO2, k2=K2, dp_flood_kister=DP_FLOOD),
        pressure=pdf.round(4).to_dict(orient="records"),
        flood=flood.to_dict(orient="records"),
        gpdc_curves=gpdc.curves_for_chart(),
        calibration=cal,
        water=water[["sheet", "G_scfm", "L_gpm", "yin", "yout", "ratio", "A_factor", "max_removal_pct"]]
            .round(4).to_dict(orient="records"),
        kog=kog,
        kog_ts=ts.round(4).replace({np.nan: None}).to_dict(orient="records"),
        flows=dict(G_scfm=G_NAOH_SCFM, L_gpm=L_NAOH_GPM, Gm_kmol_h=model.Gm_kmol_h, G_flux=model.G_flux,
                   L_L_min=model.L_L_min, L_flux_kg_s_m2=model.L_L_min / 60 / A_M2,
                   G_flux_kg_s_m2=model.Gm_kmol_h * 29.4 / 3600 / A_M2),
        **conc,
    )


def clean(o):
    """numpy -> python, NaN/inf -> None, so JSON and the page see plain values."""
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(x) for x in o]
    if isinstance(o, (np.floating, np.integer)):
        o = o.item()
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return o


# =============================================================================
# Findings page
# =============================================================================
def page_values(R):
    """Numbers quoted in the page prose ({{key}} placeholders)."""
    P = R["pressure"]
    devs = [r["dev_pct"] for r in P if r["region"] == "interpolated"]
    sheet_diffs = [abs(r["dP_pred_in"] - r["sheet_pred_in"]) / r["sheet_pred_in"] * 100
                   for r in P if r["sheet_pred_in"]]
    row = lambda L, G: next(r for r in P if r["L_gpm"] == L and r["G_scfm"] == G)  # noqa: E731
    k, kin = R["kog"], R["kinetics"]
    flood15 = next(f["G_flood_scfm"] for f in R["flood"] if f["L_gpm"] == 1.5)
    return dict(
        dev_median=f"{st.median(devs):.0f}", n_interp=len(devs),
        dev_min=f"{min(devs):.0f}", dev_max=f"{max(devs):.0f}",
        sheet_diff=f"{st.mean(sheet_diffs):.0f}",
        nu=f"{R['inputs']['nu_cSt']:.2f}",
        n_flv_high=sum(r["flv"] > 0.3 for r in P),
        dp8=f"{row(1.5, 8)['dP_meas_ft']:.2f}", dp9=f"{row(1.5, 9)['dP_meas_ft']:.2f}",
        gp8=f"{row(1.5, 8)['dP_pred_ft']:.2f}", gp9=f"{row(1.5, 9)['dP_pred_ft']:.2f}",
        flood15=f"{flood15:.0f}",
        kog_cor=f"{k['corrected']['KOGa']:.1f}", kog_cor_sd=f"{k['corrected']['KOGa_sd']:.1f}",
        kog_raw=f"{k['raw']['KOGa']:.1f}", kog_raw_sd=f"{k['raw']['KOGa_sd']:.1f}",
        rem_cor=f"{k['corrected']['removal_pct']:.0f}",
        water_rem=f"{R['driving']['water_removal_pct']:.1f}",
        gm_flux=f"{R['flows']['G_flux']:.1f}",
        cal_a=f"{R['calibration']['a']:.3f}", cal_b=f"{R['calibration']['b']:.3f}",
        cal_s=f"{R['calibration']['s_res']:.3f}",
        oh_mean=f"{kin['oh_mean']:.3f}", oh_bot=f"{kin['oh_bot']:.3f}",
        perry_scaled=f"{R['perry_scaled']:.0f}",
        area_ratio=f"{DRY_AREA_HALF_IN / DRY_AREA_25MM:.1f}",
        a_eff=f"{kin['a_eff']:.0f}", a_frac=f"{100 * kin['a_eff'] / DRY_AREA_HALF_IN:.0f}",
        A_water=f"{R['driving']['A_water']:.3f}",
        k2=f"{R['inputs']['k2']:,.0f}", Ha=f"{kin['Ha']:.0f}", Einf=f"{kin['E_inf']:.0f}",
        c_stoich=f"{R['c_stoich_full']:.3f}",
        conv=f"{k['corrected']['naoh_conv_pct']:.0f}",
        nog_cor=f"{k['corrected']['NOG']:.3f}",
        kog_water=f"{R['driving']['kog_water']:.2f}",
    )


def literature_rows(R, v):
    """Comparison table: quantity | this work | reference | source."""
    return [
        ["ΔP/Z vs GPDC",
         f"Measured {v['dev_min']}–{v['dev_max']}% above prediction (median {v['dev_median']}%); {v['n_flv_high']}/18 runs at F<sub>LV</sub> &gt; 0.3",
         "GPDC is optimistic (predicts too little ΔP) for F<sub>LV</sub> &gt; 0.3 and limited at low liquid rates; deviations are systematic",
         "Kister et al. 2007, CEP pp. 29–30"],
        ["Packing factor",
         "F<sub>p</sub> = 200 ft⁻¹ reproduces the lab group's GPDC sheet within " + v["sheet_diff"] + "%",
         "Packing factors from Kister & Gill and Strigle; values should be used with caution",
         "Perry's Table 14-7b, p. 14-43"],
        ["Flooding",
         f"Loading break at 1.5 gpm between 8 and 9 SCFM (ΔP/Z {v['dp8']} → {v['dp9']} in/ft)",
         "1.5 in H2O/ft = incipient flooding; measured floods at 2.0–2.5 in H2O/ft. GPDC flood at 1.5 gpm ≈ " + v["flood15"] + " SCFM",
         "Perry's p. 14-41/42; Strigle GPDC"],
        ["K<sub>OG</sub>a, CO2–NaOH",
         f"{v['kog_cor']} ± {v['kog_cor_sd']} (corrected), {v['kog_raw']} ± {v['kog_raw_sd']} (as recorded); 0.076 M, ½-in ceramic Intalox",
         "25 mm ceramic Intalox 45.1; 50 mm 30.1 (1 N NaOH, 25% conversion, 24 °C). 38 mm Intalox 38 (4 wt%, 33% conversion)",
         "Perry's Table 14-8, p. 14-52; Table 14-3, p. 14-20"],
        ["Liquid-rate dependence",
         f"L = {R['flows']['L_flux_kg_s_m2']:.1f} kg/(s·m²), G = {R['flows']['G_flux_kg_s_m2']:.2f} kg/(s·m²)",
         "K<sub>G</sub>a rises 50–65% from 2.7 to 13.6 kg/(s·m²) liquid for 25 mm packings",
         "Perry's Table 14-2, p. 14-19"],
        ["Concentration effect",
         f"Ha ≈ {v['Ha']}, E<sub>∞</sub> ≈ {v['Einf']}: fast pseudo-first-order, K<sub>OG</sub>a ∝ √[OH⁻]; scaled Perry value ≈ {v['perry_scaled']}",
         "K<sub>G</sub>a depends on reagent concentration and conversion; liquid-phase controlled, k<sub>L</sub> = φ·k<sub>L</sub>⁰",
         "Perry's Fig. 14-12 and Eq. 14-72 to 14-78, pp. 14-18 to 14-21"],
        ["Driving force model",
         "Plug-flow gas, ln(y<sub>in</sub>/y<sub>out</sub>) with gas-phase depletion accounted for",
         "Using a no-depletion gas model distorts CO2 k<sub>L</sub>a; CO2 is 26× more soluble than O2",
         "Kordač & Linek, Ind. Eng. Chem. Res."],
        ["Solvent strength",
         "Removal rises from ~2% (water) to 31% (0.076 M NaOH)",
         "CO2 absorption rate increases with solvent concentration (1.11–3.24 mol/L) and with gas flow rate",
         "Daneshvar et al., Chem. Biochem. Eng. Q."],
        ["Analyzer accuracy",
         f"Outlet/inlet offset: out = {v['cal_a']} + {v['cal_b']}·in",
         "NDIR reads ∝ p/T: +10 °C → −3%; 980 hPa → −3%. Water vapour dilutes the reading (Table 5)",
         "Vaisala application note; Vaisala CARBOCAP note"],
        ["Flow measurement",
         "±3% assumed on G in the Monte Carlo",
         "Rotameter accuracy ±3% of full scale; standard 70 °F, 14.7 psia; correct for gas density",
         "Mosedale, rotameter note"],
    ]


def build_page(R):
    v = page_values(R)
    html = TEMPLATE.replace("__DATA__", json.dumps(R))
    html = html.replace("__LIT__", json.dumps(literature_rows(R, v), ensure_ascii=False))
    return re.sub(r"\{\{(\w+)\}\}", lambda m: str(v[m.group(1)]), html)


def print_summary(R):
    k = R["kog"]
    print(f"Column: {D_IN} in ID x {Z_IN} in, Fp {FP} ft^-1, area {A_FT2:.4f} ft2")
    print("\nPressure drop (in H2O/ft):")
    print(pd.DataFrame(R["pressure"])[["L_gpm", "G_scfm", "flv", "cp", "dP_meas_ft", "dP_pred_ft", "dev_pct", "region"]]
          .round(3).to_string(index=False))
    print("\nGPDC flood (1.5 in H2O/ft):", {f["L_gpm"]: f["G_flood_scfm"] for f in R["flood"]})
    c = R["calibration"]
    print(f"\nAnalyzer calibration: out = {c['a']:.4f} + {c['b']:.4f}*in  (s = {c['s_res']:.4f}, n = {c['n']})")
    for label in ("raw", "corrected"):
        d = k[label]
        print(f"KOGa {label:9s}: {d['KOGa']:.2f} +/- {d['KOGa_sd']:.2f} kmol/(h m3 atm)  "
              f"NOG {d['NOG']:.3f}  removal {d['removal_pct']:.1f}%  NaOH conv {d['naoh_conv_pct']:.0f}%")
    kin = R["kinetics"]
    print(f"[OH-] mean {kin['oh_mean']:.4f} M, Ha {kin['Ha']:.1f}, E_inf {kin['E_inf']:.0f}, a_eff {kin['a_eff']:.0f} m2/m3")
    print(f"Water at same G, L: removal {R['driving']['water_removal_pct']:.2f}%, "
          f"NaOH for 100% removal {R['c_stoich_full']:.3f} M")


def main(argv=None):
    sys.stdout.reconfigure(encoding="utf-8")
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data-dir", type=Path, default=Path.home() / "Downloads",
                    help="folder holding 'Gas Absorption Data.xlsx' and 'GPDC-1.xlsx'")
    ap.add_argument("--out-dir", type=Path, default=here, help="where results.json and index.html go")
    args = ap.parse_args(argv)

    R = clean(run_analysis(args.data_dir))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "results.json").write_text(json.dumps(R, indent=1), encoding="utf-8")
    (args.out_dir / "index.html").write_text(build_page(R), encoding="utf-8")
    print_summary(R)
    print(f"\nWrote {args.out_dir / 'results.json'} and {args.out_dir / 'index.html'}")


# =============================================================================
# Page template: {{key}} = prose numbers, __DATA__ = results, __LIT__ = comparison rows
# =============================================================================
TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Packed Column CO2 Absorption Findings</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap" rel="stylesheet">
<style>
/* Layout idea: a lab report read top to bottom; summary tiles first, then one section per objective, each chart in a hairline card with its data table beneath. */
:root {
  color-scheme: light;
  --page: #f6f7f8;
  --surface: #fcfcfb;
  --ink: #0b0b0b;
  --ink-2: #4c4e52;
  --muted: #84868a;
  --grid: #e3e4e2;
  --axis: #c3c4c0;
  --border: rgba(11, 11, 11, 0.10);
  --band: rgba(42, 120, 214, 0.07);
  --link: #1c5cab;
  --s1: #2a78d6;
  --s2: #eb6834;
  --s3: #1baf7a;
  --warn-ink: #8a5a00;
  --font-body: "IBM Plex Sans", system-ui, -apple-system, "Segoe UI", sans-serif;
  --font-mono: "IBM Plex Mono", ui-monospace, "Cascadia Mono", Consolas, monospace;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --page: #0e0f10;
    --surface: #1a1b1c;
    --ink: #f0efec;
    --ink-2: #c3c2b7;
    --muted: #8d8f93;
    --grid: #2c2d2e;
    --axis: #3a3b3c;
    --border: rgba(255, 255, 255, 0.10);
    --band: rgba(57, 135, 229, 0.12);
    --link: #86b6ef;
    --s1: #3987e5;
    --s2: #d95926;
    --s3: #199e70;
    --warn-ink: #fab219;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0e0f10;
  --surface: #1a1b1c;
  --ink: #f0efec;
  --ink-2: #c3c2b7;
  --muted: #8d8f93;
  --grid: #2c2d2e;
  --axis: #3a3b3c;
  --border: rgba(255, 255, 255, 0.10);
  --band: rgba(57, 135, 229, 0.12);
  --link: #86b6ef;
  --s1: #3987e5;
  --s2: #d95926;
  --s3: #199e70;
  --warn-ink: #fab219;
}
* { box-sizing: border-box; }
html, body { margin: 0; }
body {
  background: var(--page);
  color: var(--ink);
  font-family: var(--font-body);
  font-size: 16px;
  line-height: 1.55;
}
a { color: var(--link); }
a:focus-visible, summary:focus-visible, svg:focus-visible {
  outline: 2px solid var(--link);
  outline-offset: 2px;
  border-radius: 4px;
}
.wrap { max-width: 1100px; margin: 0 auto; padding: 40px 20px 72px; display: grid; gap: 56px; }
header { display: grid; gap: 12px; }
.kicker { font-family: var(--font-mono); font-size: 13px; color: var(--ink-2); letter-spacing: 0.02em; }
h1 { font-size: clamp(28px, 4vw, 40px); line-height: 1.15; margin: 0; font-weight: 600; text-wrap: balance; }
h2 { font-size: 26px; line-height: 1.2; margin: 0; font-weight: 600; text-wrap: balance; }
h3 { font-size: 17px; margin: 0; font-weight: 600; text-wrap: balance; }
p { margin: 0; max-width: 68ch; }
.lede { font-size: 18px; color: var(--ink-2); max-width: 70ch; }
.meta { display: flex; flex-wrap: wrap; gap: 8px; }
.chip {
  font-family: var(--font-mono); font-size: 12.5px; color: var(--ink-2);
  border: 1px solid var(--border); background: var(--surface);
  padding: 3px 10px; border-radius: 999px;
}
section { display: grid; gap: 20px; }
.sec-head { display: grid; gap: 8px; border-top: 1px solid var(--axis); padding-top: 20px; }
.sec-num { font-family: var(--font-mono); font-size: 13px; color: var(--muted); }
.prose { display: grid; gap: 12px; }
.tiles { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 12px; }
.tile {
  background: var(--surface); border: 1px solid var(--border); border-radius: 10px;
  padding: 16px 18px; display: grid; gap: 4px; align-content: start;
}
.tile .label { font-size: 13.5px; color: var(--ink-2); }
.tile .value { font-size: 30px; font-weight: 600; line-height: 1.15; }
.tile .unit { font-size: 14px; font-weight: 400; color: var(--ink-2); }
.tile .note { font-size: 13px; color: var(--muted); }
.grid2 { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 16px; }
.grid4 { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px; }
@media (max-width: 900px) { .grid4 { grid-template-columns: repeat(2, minmax(0, 1fr)); } }
@media (max-width: 760px) { .grid2 { grid-template-columns: 1fr; } }
@media (max-width: 460px) { .grid4 { grid-template-columns: 1fr; } }
.card {
  background: var(--surface); border: 1px solid var(--border); border-radius: 10px;
  padding: 16px 16px 12px; display: grid; gap: 10px; min-width: 0;
}
.card .cap { font-size: 13.5px; color: var(--ink-2); max-width: none; }
.chart { position: relative; min-width: 0; }
.chart svg { display: block; width: 100%; height: auto; overflow: visible; }
.chart text { font-family: var(--font-body); }
.tick { font-size: 11.5px; fill: var(--muted); font-variant-numeric: tabular-nums; }
.axis-label { font-size: 12.5px; fill: var(--ink-2); }
.ref-label { font-size: 11px; fill: var(--muted); }
.anno { font-size: 12px; fill: var(--ink-2); }
.panel-title { font-size: 13px; font-weight: 600; fill: var(--ink); }
.legend { display: flex; flex-wrap: wrap; gap: 6px 16px; font-size: 13px; color: var(--ink-2); }
.legend span { display: inline-flex; align-items: center; gap: 6px; }
.sw { position: relative; width: 18px; height: 10px; display: inline-block; }
.sw i { position: absolute; left: 0; right: 0; top: 4px; height: 2px; border-radius: 1px; }
.sw b { position: absolute; left: 5px; top: 1px; width: 8px; height: 8px; border-radius: 50%; }
.tip {
  position: absolute; pointer-events: none; z-index: 5;
  background: var(--surface); color: var(--ink); border: 1px solid var(--border);
  box-shadow: 0 4px 16px rgba(0, 0, 0, 0.12);
  border-radius: 8px; padding: 8px 10px; font-size: 12.5px; line-height: 1.45;
  min-width: 140px; max-width: 260px; display: none; font-variant-numeric: tabular-nums;
}
.tip .t-h { font-weight: 600; margin-bottom: 2px; }
.tip .t-row { display: flex; gap: 6px; align-items: center; }
.tip .t-row b { width: 8px; height: 8px; border-radius: 50%; flex: none; }
details { font-size: 13.5px; }
summary { cursor: pointer; color: var(--link); width: fit-content; }
.tbl-wrap { overflow-x: auto; margin-top: 8px; }
table { border-collapse: collapse; font-size: 13px; font-variant-numeric: tabular-nums; width: 100%; }
th, td { text-align: right; padding: 5px 10px; border-bottom: 1px solid var(--grid); white-space: nowrap; }
th { color: var(--ink-2); font-weight: 500; }
th:first-child, td:first-child { text-align: left; }
td.wrap-cell, th.wrap-cell { white-space: normal; text-align: left; min-width: 180px; }
.eq {
  font-family: var(--font-mono); font-size: 13.5px; background: var(--surface);
  border: 1px solid var(--border); border-radius: 8px; padding: 12px 14px;
  overflow-x: auto; white-space: pre; line-height: 1.7; color: var(--ink);
}
.bars { display: grid; gap: 10px; }
.bar-row { display: grid; grid-template-columns: minmax(150px, 230px) 1fr; gap: 12px; align-items: center; }
.bar-row .bl { font-size: 13px; color: var(--ink-2); line-height: 1.3; }
.bar-track { position: relative; height: 20px; }
.bar-fill { position: absolute; left: 0; top: 2px; height: 16px; border-radius: 0 4px 4px 0; }
.bar-err { position: absolute; top: 5px; height: 10px; border-left: 1.5px solid var(--ink); border-right: 1.5px solid var(--ink); }
.bar-err::after { content: ""; position: absolute; left: 0; right: 0; top: 4px; border-top: 1.5px solid var(--ink); }
.bar-val { position: absolute; top: 0; font-size: 12.5px; color: var(--ink); white-space: nowrap; font-variant-numeric: tabular-nums; line-height: 20px; }
.bar-axis { display: grid; grid-template-columns: minmax(150px, 230px) 1fr; gap: 12px; }
.bar-axis .ticks { position: relative; height: 16px; border-top: 1px solid var(--axis); }
.bar-axis .ticks span { position: absolute; top: 2px; transform: translateX(-50%); font-size: 11.5px; color: var(--muted); }
.callout {
  border-left: 3px solid var(--s1); background: var(--surface); padding: 12px 16px;
  border-radius: 0 8px 8px 0; display: grid; gap: 6px;
}
.callout strong { font-weight: 600; }
ul.clean { margin: 0; padding-left: 20px; display: grid; gap: 6px; max-width: 75ch; }
.src { font-size: 13.5px; }
.src li { color: var(--ink-2); }
.flag { color: var(--warn-ink); font-weight: 500; }
footer { font-size: 13px; color: var(--muted); }
@media (prefers-reduced-motion: no-preference) { .tip { transition: opacity 0.08s; } }
</style>
</head>
<body>
<div class="wrap">

<header>
  <div class="kicker">Unit operations lab · CO2–air absorption · runs 17 Sep – 15 Oct 2026</div>
  <h1>CO2 absorption in a 4-inch packed column: pressure drop, K<sub>OG</sub>a and the effect of NaOH</h1>
  <p class="lede">We measured gas-side pressure drop over 18 gas and liquid rate combinations and absorbed CO2 into water and into dilute NaOH. This page compares the hydraulics with the generalized pressure drop correlation (GPDC), reports the overall gas-phase mass-transfer coefficient for CO2 into 0.076 M NaOH, and explains why the hydroxide controls the driving force.</p>
  <div class="meta">
    <span class="chip">ID 4 in · packed height 48 in</span>
    <span class="chip">½-in ceramic Intalox saddles · F<sub>p</sub> = 200 ft⁻¹</span>
    <span class="chip">NaOH 0.076 M (HCl titration)</span>
    <span class="chip">P ≈ 1 atm · T ≈ 22 °C (assumed)</span>
  </div>
</header>

<section aria-labelledby="sum-h">
  <h2 id="sum-h" class="visually-plain">Key findings</h2>
  <div class="tiles">
    <div class="tile">
      <div class="label">Measured ΔP vs GPDC prediction</div>
      <div class="value">+{{dev_median}}%</div>
      <div class="note">median over {{n_interp}} runs on the chart; every run sits above the prediction</div>
    </div>
    <div class="tile">
      <div class="label">K<sub>OG</sub>a, CO2 into 0.076 M NaOH</div>
      <div class="value">{{kog_cor}} <span class="unit">± {{kog_cor_sd}}</span></div>
      <div class="note">kmol/(h·m³·atm), sensor-corrected; {{kog_raw}} ± {{kog_raw_sd}} uncorrected</div>
    </div>
    <div class="tile">
      <div class="label">CO2 removed, 7.5 SCFM gas · 1.35 gpm liquid</div>
      <div class="value">{{rem_cor}}%</div>
      <div class="note">0.076 M NaOH (corrected). Water at the same rates is limited to about {{water_rem}}%</div>
    </div>
    <div class="tile">
      <div class="label">Loading onset at 1.5 gpm</div>
      <div class="value">8–9 <span class="unit">SCFM</span></div>
      <div class="note">ΔP/Z doubled from {{dp8}} to {{dp9}} in H2O/ft; GPDC puts flooding at {{flood15}} SCFM</div>
    </div>
  </div>
</section>

<!-- ============================ 1. PRESSURE DROP ============================ -->
<section aria-labelledby="s1-h">
  <div class="sec-head">
    <div class="sec-num">Objective 1</div>
    <h2 id="s1-h">The GPDC under-predicts pressure drop in this column</h2>
  </div>
  <div class="prose">
    <p>For each run on 17 Sep we averaged the logged pressure drop, divided by the 4 ft bed height, and read the Strigle version of the GPDC (Perry's Fig. 14-48, digitized in <code>GPDC-1.xlsx</code>) at the run's flow parameter and capacity parameter. Measured ΔP/Z is higher than predicted at every one of the 18 settings. Where the point falls inside the chart, the excess runs from {{dev_min}}% to {{dev_max}}%, with a median of {{dev_median}}%. The predictions closely reproduce the values in the lab group's own "GPDC Comparison" sheet (mean difference {{sheet_diff}}%), so the gap comes from the column itself, not from the arithmetic.</p>
  </div>
  <div class="eq">C_s = U_s·[ρ_G/(ρ_L − ρ_G)]^0.5            capacity parameter  CP = C_s·F_p^0.5·ν^0.05
F_LV = (L/G)·(ρ_G/ρ_L)^0.5                  U_s in ft/s, F_p = 200 ft⁻¹, ν = {{nu}} cSt
ΔP/Z (GPDC) = curve value at (F_LV, CP), log-interpolated between the 0.05 … 1.5 in H2O/ft curves</div>

  <div class="card">
    <h3>Pressure drop per foot of packing vs gas rate, one panel per liquid rate</h3>
    <div class="legend" id="leg-dp"></div>
    <div class="grid4" id="dp-multiples"></div>
    <p class="cap">Error bars are ±1 standard deviation of the 1-s readings at each setting. Panels share the same y-axis.</p>
    <details><summary>Show data table</summary><div class="tbl-wrap" id="tbl-dp"></div></details>
  </div>

  <div class="grid2">
    <div class="card">
      <h3>Measured vs predicted ΔP/Z</h3>
      <div class="chart" id="ch-parity"></div>
      <p class="cap">Every point lies above the 1 : 1 line. The two highest relative errors are at 3 SCFM, where ΔP is below 0.1 in H2O/ft and close to the transducer's resolution.</p>
    </div>
    <div class="card">
      <h3>The runs on the Strigle GPDC chart</h3>
      <div class="chart" id="ch-gpdc"></div>
      <p class="cap">Grey curves are constant ΔP/Z (in H2O/ft). Hover a point to compare its measured ΔP/Z with the curve value at its position.</p>
    </div>
  </div>

  <div class="callout">
    <strong>Why the measured ΔP runs high</strong>
    <ul class="clean">
      <li><b>High flow parameter.</b> {{n_flv_high}} of 18 runs have F<sub>LV</sub> &gt; 0.3. Kister et al. (2007) report that the GPDC is optimistic, meaning it predicts too little ΔP, above F<sub>LV</sub> ≈ 0.3. It also loses accuracy at low liquid rates.</li>
      <li><b>Small column, small packing.</b> The column-to-packing diameter ratio is only 8. Wall channelling and a non-uniform bed are expected, and published F<sub>p</sub> values come from much larger towers. Kister et al. recommend GPDC interpolation charts built from data for the specific packing for exactly this reason.</li>
      <li><b>What the taps measure.</b> If the pressure taps span the support plate and liquid distributor as well as the bed, the measurement includes losses that the per-foot correlation does not.</li>
      <li><b>Early loading.</b> At 1.5 gpm, ΔP/Z jumped from {{dp8}} to {{dp9}} in H2O/ft between 8 and 9 SCFM, while the GPDC predicted only {{gp8}} to {{gp9}}. Perry's §14 treats 1.5 in H2O/ft as incipient flooding. The 9 SCFM point is close to that, about 2 SCFM below the GPDC flood estimate of {{flood15}} SCFM.</li>
    </ul>
  </div>
</section>

<!-- ============================ 2. KOGa ============================ -->
<section aria-labelledby="s2-h">
  <div class="sec-head">
    <div class="sec-num">Objective 2</div>
    <h2 id="s2-h">K<sub>OG</sub>a for CO2 into 0.076 M NaOH is about 20 kmol/(h·m³·atm)</h2>
  </div>
  <div class="prose">
    <p>Because the hydroxide reacts with CO2 quickly and stays in excess, the liquid holds essentially no free CO2, so y* ≈ 0 throughout the bed. With only about 3% CO2 the system is dilute, so the integrated rate equation reduces to Perry's Eq. 14-76:</p>
  </div>
  <div class="eq">N_OG = ln(y_in / y_out)        H_OG = Z / N_OG        K_OG·a = G_M / (H_OG · P)
G_M = {{gm_flux}} kmol/(h·m²)   (7.5 SCFM at 70 °F, 14.7 psia over 0.0873 ft²)      Z = 1.219 m      P = 1 atm</div>

  <div class="prose">
    <p><b>The CO2 sensors needed a correction first.</b> With no CO2 fed (17 Sep) and in every water run, the outlet analyzer read higher than the inlet analyzer, by 6–10% in the water runs. Water cannot produce that: at these flows its absorption factor is only about 0.01–0.03, so it can remove at most 1–3% of the CO2. The difference is therefore an offset between the two instruments. A straight-line fit of outlet reading against inlet reading over those 27 runs gives <code>out = {{cal_a}} + {{cal_b}}·in</code>. Inverting this line puts the NaOH outlet readings on the inlet analyzer's scale.</p>
  </div>

  <div class="grid2">
    <div class="card">
      <h3>Inlet and outlet CO2 during the NaOH run</h3>
      <div class="legend" id="leg-ts"></div>
      <div class="chart" id="ch-ts"></div>
      <p class="cap">10-s averages. The shaded window (60–680 s) is the steady period used for K<sub>OG</sub>a. Pressure drop collapsed at about 690 s (below) as liquid flow stopped, and absorption stopped with it.</p>
    </div>
    <div class="card">
      <h3>Column pressure drop during the NaOH run</h3>
      <div class="chart" id="ch-tsdp"></div>
      <p class="cap">Steady wet-bed ΔP of about 2.1 in H2O fell to about 0.6 in H2O at the same gas rate, as liquid holdup drained from the bed. This indicates that the loss of absorption after 690 s is a liquid-supply event, not a change in mass transfer.</p>
    </div>
  </div>

  <div class="grid2">
    <div class="card">
      <h3>Analyzer cross-calibration</h3>
      <div class="legend" id="leg-cal"></div>
      <div class="chart" id="ch-cal"></div>
      <p class="cap">Runs with no CO2 fed and water runs, where outlet ≈ inlet is expected. Residual standard deviation {{cal_s}} mol%.</p>
      <details><summary>Show data table</summary><div class="tbl-wrap" id="tbl-cal"></div></details>
    </div>
    <div class="card">
      <h3>K<sub>OG</sub>a compared with Perry's</h3>
      <div class="legend" id="leg-bars"></div>
      <div id="bars-kog" class="bars"></div>
      <p class="cap">Our error bars are ±1σ from a Monte Carlo over the time variation of the analyzers, the calibration fit and ±3% gas flow. The Perry's values are for 1 N NaOH at 25% conversion, the vendor standard test (Table 14-8). The scaled bar adjusts the 25 mm value to our mean [OH⁻] using K<sub>OG</sub>a ∝ √[OH⁻] (objective 3).</p>
    </div>
  </div>

  <div class="card">
    <h3>K<sub>OG</sub>a calculation, steady window 60–680 s</h3>
    <div class="tbl-wrap" id="tbl-kog"></div>
  </div>

  <div class="callout">
    <strong>How it compares</strong>
    <p>Our corrected {{kog_cor}} kmol/(h·m³·atm) is well below Perry's 1 N values (45.1 for 25 mm ceramic Intalox, 38 for 38 mm). That is expected at 0.076 M: CO2 absorption into NaOH is liquid-film controlled (Perry's p. 14-18), and the reaction-enhanced liquid coefficient grows only with √[OH⁻]. Scaling the 25 mm value to our mean [OH⁻] of {{oh_mean}} M gives about {{perry_scaled}}. Our ½-in saddles carry {{area_ratio}}× the surface area per volume of 25 mm saddles (623 vs 256 m²/m³, Table 14-7b), so a result somewhat above {{perry_scaled}} is consistent with the literature. The back-calculated effective interfacial area is {{a_eff}} m²/m³, about {{a_frac}}% of the dry packing area.</p>
  </div>
</section>

<!-- ============================ 3. NaOH concentration ============================ -->
<section aria-labelledby="s3-h">
  <div class="sec-head">
    <div class="sec-num">Objective 3</div>
    <h2 id="s3-h">NaOH keeps the driving force at its maximum, and its concentration sets how fast the column uses it</h2>
  </div>
  <div class="prose">
    <p>The driving force for absorption is y − y*, the gap between the CO2 in the gas and the CO2 that would be in equilibrium with the liquid. In water, dissolved CO2 builds up as the liquid descends: y* climbs toward y and the driving force pinches near the bottom. At our flows water has an absorption factor A = L/(mG) of only {{A_water}}, so it saturates after removing about {{water_rem}}% of the CO2. Hydroxide converts dissolved CO2 to carbonate (CO2 + 2 OH⁻ → CO3²⁻ + H2O). Free CO2 stays near zero, so y* ≈ 0 and the whole gas-phase concentration is available as driving force, from the bottom to the top of the bed.</p>
    <p>Hydroxide concentration then sets the rate. With k₂ = {{k2}} L/(mol·s) at 22 °C and a mean [OH⁻] of {{oh_mean}} M, the Hatta number is about {{Ha}}. That is well below the instantaneous-reaction limit E<sub>∞</sub> ≈ {{Einf}}, so the reaction is <i>fast pseudo-first-order</i>: absorption is enhanced by E ≈ √(1 + Ha²), which scales with √[OH⁻]. A second limit is stoichiometric. Each mole of CO2 consumes two moles of OH⁻, so complete removal at these flows needs at least {{c_stoich}} M NaOH. At 0.076 M we converted about {{conv}}% of the hydroxide, leaving the bottom of the bed at about {{oh_bot}} M.</p>
  </div>

  <div class="card">
    <h3>Driving force along the bed: water vs 0.076 M NaOH (same gas and liquid rates)</h3>
    <div class="legend" id="leg-df"></div>
    <div class="chart" id="ch-df"></div>
    <p class="cap">NaOH: measured N<sub>OG</sub> = {{nog_cor}} with y* = 0. Water: countercurrent integration with Henry's law (H ≈ 1,630 atm) and the physical coefficient K<sub>OG</sub>a = {{kog_water}} kmol/(h·m³·atm), the NaOH value divided by the enhancement factor. The water curve is a prediction, because the analyzers cannot resolve a 2% change. At the top, fresh liquid enters, so y* = 0 for both liquids. Water's curve is higher there only because it has removed almost no CO2 from the gas. Lower in the bed, its y* climbs until little driving force is left.</p>
  </div>

  <div class="grid2">
    <div class="card">
      <h3>Predicted K<sub>OG</sub>a vs feed NaOH concentration</h3>
      <div class="legend" id="leg-kc"></div>
      <div class="chart" id="ch-kc"></div>
      <p class="cap">Model: E = √(1 + Ha²), with Ha ∝ √[OH⁻], anchored at our measured point. Above about 0.5 N, ionic-strength and viscosity effects (Perry's Fig. 14-12) would flatten the real curve, so treat that range as indicative only.</p>
    </div>
    <div class="card">
      <h3>Predicted CO2 removal vs feed NaOH concentration</h3>
      <div class="legend" id="leg-rc"></div>
      <div class="chart" id="ch-rc"></div>
      <p class="cap">The grey line is the most CO2 the hydroxide supply could neutralize. The model sits below it everywhere except under about 0.01 M, so across this range the reaction rate limits removal, not the hydroxide supply. Our 0.076 M run is rate-limited.</p>
    </div>
  </div>
  <details class="card"><summary>Show model table</summary><div class="tbl-wrap" id="tbl-conc"></div></details>
</section>

<!-- ============================ Literature ============================ -->
<section aria-labelledby="s4-h">
  <div class="sec-head">
    <div class="sec-num">Comparison</div>
    <h2 id="s4-h">Results against the reference material</h2>
  </div>
  <div class="card"><div class="tbl-wrap" id="tbl-lit"></div></div>
</section>

<!-- ============================ Uncertainty ============================ -->
<section aria-labelledby="s5-h">
  <div class="sec-head">
    <div class="sec-num">Uncertainty</div>
    <h2 id="s5-h">Error sources and how large they are</h2>
  </div>
  <ul class="clean">
    <li><b>Analyzer offset, the largest effect.</b> The correction raises K<sub>OG</sub>a from {{kog_raw}} to {{kog_cor}}. Water runs remove up to 1–3% of the CO2, which slightly inflates the fitted slope, so the corrected value may be high by a similar few percent.</li>
    <li><b>Pressure and temperature at the analyzers.</b> NDIR sensors count molecules per volume, so readings scale with p/T (Vaisala). A 2.2 in H2O pressure difference between inlet and outlet changes the reading by about 0.5%. A 10 °C temperature difference changes it by about 3%.</li>
    <li><b>Humidity.</b> Gas leaving the bed is close to saturation (about 2.6% H2O at 22 °C), which dilutes the outlet CO2 reading by about 2–3% if the inlet gas is dry (Vaisala Table 5). The cross-calibration absorbs this, because the water runs were humidified the same way.</li>
    <li><b>Rotameters.</b> Variable-area meters are typically ±3% of full scale (Mosedale), so percentage errors are largest at the low end, such as 3 SCFM and 0.5 gpm. Gas containing 3% CO2 is denser than the air the meter was calibrated for, so it reads about 0.8% high (correction √(M<sub>air</sub>/M)).</li>
    <li><b>Assumed properties.</b> Taking T = 22 °C sets ν, ρ, Henry's constant and k₂. ±3 °C moves k₂ by about ±20%, but K<sub>OG</sub>a comes from the measured compositions and does not depend on it.</li>
    <li><b>Transients.</b> Start-up (0–60 s) and the liquid-supply loss after 680 s are excluded. Pump start and stop transients (Perry's §6, hydraulic transients) disturb flow briefly.</li>
    <li><b>Once-through NaOH.</b> At 1.35 gpm the steady window used about 53 L of solution. The calculation assumes 0.076 M was the feed concentration throughout. If the solution was recirculated, the hydroxide conversion figure would understate depletion.</li>
  </ul>
</section>

<section aria-labelledby="s6-h">
  <div class="sec-head">
    <div class="sec-num">Sources</div>
    <h2 id="s6-h">Data and references</h2>
  </div>
  <ul class="clean src">
    <li>Lab data: <code>Gas Absorption Data.xlsx</code> (Summary, 9-17 Pressure Drop, 9-24 and 10-1 Water CO2, 10-15 NaOH CO2, GPDC Comparison); <code>GPDC-1.xlsx</code> (digitized Strigle GPDC).</li>
    <li>Perry's Chemical Engineers' Handbook, 7th ed., §14 Gas Absorption and Gas-Liquid System Design: Eq. 14-76, Tables 14-2, 14-3, 14-7b and 14-8, Figs. 14-12 and 14-48, pp. 14-18 to 14-20, 14-41 to 14-43 and 14-52.</li>
    <li>Kister, Scherffius, Afshar and Abkar, "Realistically Predict Capacity and Pressure Drop for Packed Columns," CEP, July 2007, from p. 28: GPDC limits, Eq. 4 flood ΔP.</li>
    <li>Kordač and Linek, "Dynamic Measurement of Carbon Dioxide Volumetric Mass Transfer Coefficient…": gas-phase depletion and the driving-force model in CO2 k<sub>L</sub>a measurement.</li>
    <li>Daneshvar et al., "Absorption of Carbon Dioxide into Aqueous Solutions of N,N-dimethylethanolamine…": absorption rate rises with solvent concentration and gas rate.</li>
    <li>Vaisala, "CARBOCAP Carbon Dioxide Sensors: What is CO2?" and Häkkänen, "Atmospheric Conditions Affect CO2 Measurement Accuracy": p/T and humidity effects on NDIR readings.</li>
    <li>Mosedale, "The Surprisingly Confusing World of Rotameters": standard conditions (70 °F, 14.7 psia), gas-density correction, ±3% of full scale.</li>
    <li>Perry's 7th ed., §6 (hydraulic transients): pump start and stop effects.</li>
    <li>Rate constant: Pohorecki and Moniuk, log₁₀k₂ = 11.895 − 2382/T. Henry's constant for CO2–water and diffusivities are standard literature values at 22 °C.</li>
  </ul>
</section>

<footer>All values were calculated from the lab workbooks by <code>analyze.py</code>. Units: ΔP in in H2O, K<sub>OG</sub>a in kmol/(h·m³·atm), compositions in mol% CO2.</footer>
</div>

<script>
const D = __DATA__;

/* ------------------------------------------------------------ helpers */
const NS = "http://www.w3.org/2000/svg";
function el(tag, attrs, parent) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs || {}) e.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(e);
  return e;
}
function h(tag, attrs, html) {
  const e = document.createElement(tag);
  for (const k in attrs || {}) e.setAttribute(k, attrs[k]);
  if (html != null) e.innerHTML = html;
  return e;
}
function niceTicks(min, max, n) {
  const span = max - min, raw = span / (n || 5);
  const mag = Math.pow(10, Math.floor(Math.log10(raw))), r = raw / mag;
  const step = (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * mag;
  const out = [];
  for (let v = Math.ceil(min / step - 1e-9) * step; v <= max + span * 1e-9; v += step) out.push(+v.toFixed(10));
  return out;
}
function logTicks(min, max) {
  const out = [];
  for (let e = Math.floor(Math.log10(min)); e <= Math.ceil(Math.log10(max)); e++)
    for (const m of [1, 2, 5]) { const v = m * Math.pow(10, e); if (v >= min * 0.999 && v <= max * 1.001) out.push(v); }
  return out;
}
function tickFmt(v) {
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1000) return v.toLocaleString("en-US", { maximumFractionDigits: 0 });
  if (a >= 1) return String(+v.toFixed(2));
  return String(+v.toPrecision(2));
}
const f = (v, d) => (v == null ? "–" : Number(v).toFixed(d));
function scale(a, p0, p1) {
  if (a.log) {
    const l0 = Math.log10(a.min), l1 = Math.log10(a.max);
    return v => p0 + (Math.log10(v) - l0) / (l1 - l0) * (p1 - p0);
  }
  return v => p0 + (v - a.min) / (a.max - a.min) * (p1 - p0);
}
let clipN = 0;

function legend(id, items) {
  const host = document.getElementById(id);
  if (!host) return;
  host.innerHTML = "";
  for (const it of items) {
    const sw = it.kind === "ref"
      ? `<span class="sw"><i style="background:var(--muted);height:1px;top:5px"></i></span>`
      : `<span class="sw">${it.line ? `<i style="background:var(${it.color})"></i>` : ""}${it.points !== false ? `<b style="background:var(${it.color})"></b>` : ""}</span>`;
    host.appendChild(h("span", {}, sw + it.name));
  }
}

/* ------------------------------------------------------------ chart */
function chart(host, o) {
  host.style.position = "relative";
  const draw = () => {
    host.innerHTML = "";
    const W = Math.max(240, host.clientWidth || 300), H = o.height || 300;
    const m = { l: o.ml || 58, r: o.mr || 16, t: o.mt || 14, b: o.mb || 46 };
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H, role: "img", tabindex: "0", "aria-label": o.aria || "" }, host);
    const sx = scale(o.x, m.l, W - m.r), sy = scale(o.y, H - m.b, m.t);
    const cid = "clip" + (++clipN);
    el("rect", { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b }, el("clipPath", { id: cid }, el("defs", {}, svg)));

    if (o.panelTitle) { const t = el("text", { x: m.l, y: 11, class: "panel-title" }, svg); t.textContent = o.panelTitle; }

    // bands
    for (const b of o.bands || []) {
      el("rect", { x: sx(b.x0), y: m.t, width: sx(b.x1) - sx(b.x0), height: H - m.t - m.b, fill: "var(--band)" }, svg);
      if (b.label) { const t = el("text", { x: sx(b.x0) + 6, y: m.t + 14, class: "ref-label" }, svg); t.textContent = b.label; }
    }
    // grid + ticks
    const xt = o.x.ticks || (o.x.log ? logTicks(o.x.min, o.x.max) : niceTicks(o.x.min, o.x.max, o.x.n || 6));
    const yt = o.y.ticks || (o.y.log ? logTicks(o.y.min, o.y.max) : niceTicks(o.y.min, o.y.max, o.y.n || 5));
    for (const v of yt) {
      const y = sy(v);
      el("line", { x1: m.l, x2: W - m.r, y1: y, y2: y, stroke: "var(--grid)", "stroke-width": 1 }, svg);
      if (!o.hideY) { const t = el("text", { x: m.l - 8, y: y + 4, "text-anchor": "end", class: "tick" }, svg); t.textContent = (o.y.fmt || tickFmt)(v); }
    }
    for (const v of xt) {
      const x = sx(v);
      if (o.x.grid) el("line", { x1: x, x2: x, y1: m.t, y2: H - m.b, stroke: "var(--grid)", "stroke-width": 1 }, svg);
      el("line", { x1: x, x2: x, y1: H - m.b, y2: H - m.b + 4, stroke: "var(--axis)" }, svg);
      const t = el("text", { x: x, y: H - m.b + 17, "text-anchor": "middle", class: "tick" }, svg); t.textContent = (o.x.fmt || tickFmt)(v);
    }
    el("line", { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, stroke: "var(--axis)", "stroke-width": 1 }, svg);
    if (o.x.label) { const t = el("text", { x: (m.l + W - m.r) / 2, y: H - 8, "text-anchor": "middle", class: "axis-label" }, svg); t.textContent = o.x.label; }
    if (o.y.label && !o.hideY) { const t = el("text", { transform: `translate(14 ${(m.t + H - m.b) / 2}) rotate(-90)`, "text-anchor": "middle", class: "axis-label" }, svg); t.textContent = o.y.label; }

    const plot = el("g", { "clip-path": `url(#${cid})` }, svg);
    const path = pts => pts.filter(p => p.y != null && isFinite(p.y)).map((p, i) => (i ? "L" : "M") + sx(p.x).toFixed(1) + " " + sy(p.y).toFixed(1)).join("");

    // reference curves
    for (const r of o.refs || []) {
      el("path", { d: path(r.data), fill: "none", stroke: "var(--muted)", "stroke-width": r.width || 1 }, plot);
      if (r.label) {
        const at = r.labelAt || r.data[r.data.length - 1];
        const t = el("text", { x: sx(at.x) + (r.dx || 4), y: sy(at.y) + (r.dy || -4), "text-anchor": r.anchor || "start", class: "ref-label" }, svg);
        t.textContent = r.label;
      }
    }
    // vertical markers
    for (const a of o.vlines || []) {
      el("line", { x1: sx(a.x), x2: sx(a.x), y1: m.t, y2: H - m.b, stroke: "var(--ink-2)", "stroke-width": 1 }, svg);
      const t = el("text", { x: sx(a.x) + (a.dx || -6), y: m.t + (a.dy || 14), "text-anchor": a.anchor || "end", class: "anno" }, svg);
      t.textContent = a.text;
    }

    const hits = [];
    for (const s of o.series) {
      const col = `var(${s.color})`;
      if (s.line) el("path", { d: path(s.data), fill: "none", stroke: col, "stroke-width": 2, "stroke-linejoin": "round", "stroke-linecap": "round" }, plot);
      for (const p of s.data) {
        if (p.y == null) continue;
        const px = sx(p.x), py = sy(p.y);
        if (p.lo != null && p.hi != null) {
          el("line", { x1: px, x2: px, y1: sy(p.lo), y2: sy(p.hi), stroke: col, "stroke-width": 1.25 }, plot);
          for (const yy of [p.lo, p.hi]) el("line", { x1: px - 3, x2: px + 3, y1: sy(yy), y2: sy(yy), stroke: col, "stroke-width": 1.25 }, plot);
        }
        if (s.points !== false) el("circle", { cx: px, cy: py, r: s.r || 4, fill: col, stroke: "var(--surface)", "stroke-width": 2 }, plot);
        hits.push({ px, py, x: p.x, s, p });
      }
      if (s.endLabel) {
        const last = s.data.filter(p => p.y != null).slice(-1)[0];
        const t = el("text", { x: sx(last.x) + 6, y: sy(last.y) + 4, class: "anno" }, svg); t.textContent = s.endLabel;
      }
    }
    for (const a of o.annos || []) {
      const t = el("text", { x: sx(a.x) + (a.dx || 0), y: sy(a.y) + (a.dy || 0), "text-anchor": a.anchor || "start", class: "anno" }, svg);
      t.textContent = a.text;
    }

    // tooltip layer
    const tip = h("div", { class: "tip", role: "status" });
    host.appendChild(tip);
    const hl = el("circle", { r: 6, fill: "none", stroke: "var(--ink)", "stroke-width": 1.5, visibility: "hidden" }, svg);
    const xhair = el("line", { y1: m.t, y2: H - m.b, stroke: "var(--ink-2)", "stroke-width": 1, visibility: "hidden" }, svg);
    const show = (hit, rows) => {
      tip.innerHTML = rows;
      tip.style.display = "block";
      const k = svg.getBoundingClientRect().width / W;
      let left = hit.px * k + 14, top = hit.py * k - 10;
      if (left + tip.offsetWidth > host.clientWidth) left = hit.px * k - tip.offsetWidth - 14;
      tip.style.left = Math.max(0, left) + "px";
      tip.style.top = Math.max(0, top) + "px";
    };
    const hide = () => { tip.style.display = "none"; hl.setAttribute("visibility", "hidden"); xhair.setAttribute("visibility", "hidden"); };
    const row = (s, txt) => `<div class="t-row"><b style="background:var(${s.color})"></b>${txt}</div>`;
    const showHit = hit => {
      if (o.crosshair) {
        const same = hits.filter(q => q.x === hit.x);
        xhair.setAttribute("x1", hit.px); xhair.setAttribute("x2", hit.px); xhair.setAttribute("visibility", "visible");
        show(hit, `<div class="t-h">${o.crosshair(hit.x)}</div>` + same.map(q => row(q.s, q.s.tip ? q.s.tip(q.p) : `${q.s.name}: ${f(q.p.y, 3)}`)).join(""));
      } else {
        hl.setAttribute("cx", hit.px); hl.setAttribute("cy", hit.py); hl.setAttribute("visibility", "visible");
        show(hit, row(hit.s, `<span><span class="t-h">${hit.s.name}</span></span>`) + (hit.s.tip ? hit.s.tip(hit.p) : ""));
      }
    };
    svg.addEventListener("pointermove", ev => {
      const r = svg.getBoundingClientRect(), k = W / r.width;
      const mx = (ev.clientX - r.left) * k, my = (ev.clientY - r.top) * k;
      let best = null, bd = Infinity;
      for (const q of hits) {
        const d = o.crosshair ? Math.abs(q.px - mx) : Math.hypot(q.px - mx, q.py - my);
        if (d < bd) { bd = d; best = q; }
      }
      if (best && (o.crosshair ? mx >= m.l && mx <= W - m.r : bd < 28)) showHit(best); else hide();
    });
    svg.addEventListener("pointerleave", hide);
    const order = hits.slice().sort((a, b) => a.px - b.px || a.py - b.py);
    let ki = -1;
    svg.addEventListener("keydown", ev => {
      if (ev.key !== "ArrowRight" && ev.key !== "ArrowLeft") return;
      ev.preventDefault();
      ki = (ki + (ev.key === "ArrowRight" ? 1 : -1) + order.length) % order.length;
      showHit(order[ki]);
    });
    svg.addEventListener("blur", hide);
  };
  draw();
  let w = host.clientWidth;
  new ResizeObserver(() => { if (Math.abs(host.clientWidth - w) > 4) { w = host.clientWidth; draw(); } }).observe(host);
}

function table(id, cols, rows) {
  const host = document.getElementById(id);
  const t = h("table");
  t.appendChild(h("thead", {}, "<tr>" + cols.map(c => `<th class="${c.cls || ""}">${c.h}</th>`).join("") + "</tr>"));
  const tb = h("tbody");
  for (const r of rows) tb.appendChild(h("tr", {}, cols.map(c => `<td class="${c.cls || ""}">${c.v(r)}</td>`).join("")));
  t.appendChild(tb);
  host.appendChild(t);
}

/* ------------------------------------------------------------ objective 1 */
const P = D.pressure;
const Ls = [...new Set(P.map(r => r.L_gpm))].sort((a, b) => a - b);
legend("leg-dp", [
  { name: "Measured (±1σ)", color: "--s1", line: true },
  { name: "GPDC prediction", color: "--s2", line: true },
]);
const yMaxDp = 1.5;
const mult = document.getElementById("dp-multiples");
Ls.forEach((L, i) => {
  const rows = P.filter(r => r.L_gpm === L).sort((a, b) => a.G_scfm - b.G_scfm);
  const c = h("div", { class: "chart" });
  mult.appendChild(c);
  chart(c, {
    height: 250, ml: i === 0 ? 50 : 34, mt: 22, panelTitle: `L = ${L} gpm (${f(L / D.inputs.area_ft2, 1)} gpm/ft²)`,
    aria: `Pressure drop per foot vs gas rate at ${L} gpm`,
    x: { min: 0, max: 14, ticks: [0, 4, 8, 12], label: "Gas, SCFM" },
    y: { min: 0, max: yMaxDp, ticks: [0, 0.5, 1, 1.5], label: i === 0 ? "ΔP/Z, in H2O/ft" : "" },
    series: [
      { name: "Measured", color: "--s1", line: true,
        data: rows.map(r => ({ x: r.G_scfm, y: r.dP_meas_ft, lo: r.dP_meas_ft - r.sd_ft, hi: r.dP_meas_ft + r.sd_ft, r })),
        tip: p => `${p.r.G_scfm} SCFM · ${L} gpm<br>Measured ${f(p.y, 3)} ± ${f(p.r.sd_ft, 3)} in/ft<br>GPDC ${f(p.r.dP_pred_ft, 3)} in/ft` },
      { name: "GPDC", color: "--s2", line: true,
        data: rows.map(r => ({ x: r.G_scfm, y: r.dP_pred_ft, r })),
        tip: p => `${p.r.G_scfm} SCFM · ${L} gpm<br>Predicted ${f(p.y, 3)} in/ft<br><span style="color:var(--muted)">${p.r.region}</span>` },
    ],
  });
});
table("tbl-dp", [
  { h: "Liquid, gpm", v: r => f(r.L_gpm, 2) },
  { h: "Gas, SCFM", v: r => f(r.G_scfm, 0) },
  { h: "F<sub>LV</sub>", v: r => f(r.flv, 3) },
  { h: "CP", v: r => f(r.cp, 3) },
  { h: "ΔP meas, in H2O", v: r => f(r.dP_meas_in, 3) },
  { h: "ΔP/Z meas, in/ft", v: r => f(r.dP_meas_ft, 3) },
  { h: "ΔP/Z GPDC, in/ft", v: r => f(r.dP_pred_ft, 3) },
  { h: "Deviation", v: r => "+" + f(r.dev_pct, 0) + "%" },
  { h: "Sheet's GPDC ΔP, in H2O", v: r => r.sheet_pred_in == null ? "–" : f(r.sheet_pred_in, 3) },
  { h: "Chart region", cls: "wrap-cell", v: r => r.region },
], P);

chart(document.getElementById("ch-parity"), {
  height: 320, aria: "Measured versus predicted pressure drop, log scales",
  x: { log: true, min: 0.02, max: 2, label: "GPDC ΔP/Z, in H2O/ft" },
  y: { log: true, min: 0.02, max: 2, label: "Measured ΔP/Z, in H2O/ft" },
  refs: [
    { data: [{ x: 0.02, y: 0.02 }, { x: 2, y: 2 }], label: "1 : 1", labelAt: { x: 1.15, y: 1.15 }, dx: 6, dy: 12 },
    { data: [{ x: 0.02, y: 0.04 }, { x: 1, y: 2 }], label: "2× prediction", labelAt: { x: 0.6, y: 1.2 }, dx: -8, dy: -2, anchor: "end" },
  ],
  series: [{ name: "Run", color: "--s1",
    data: P.map(r => ({ x: r.dP_pred_ft, y: r.dP_meas_ft, r })),
    tip: p => `${p.r.L_gpm} gpm · ${p.r.G_scfm} SCFM<br>Measured ${f(p.y, 3)} in/ft<br>GPDC ${f(p.x, 3)} in/ft (+${f(p.r.dev_pct, 0)}%)` }],
});

const curveRefs = Object.entries(D.gpdc_curves).map(([dp, pts]) => {
  const data = pts.map(([x, y]) => ({ x, y })).filter(p => p.x >= 0.1 && p.x <= 4);
  const dy = { "1.5": -2, "1.0": 10, "1": 10 }[dp] ?? 4;   // 1.0 and 1.5 curves nearly touch at the left edge
  return { data, label: `${+dp}`, labelAt: data[0], dx: -4, dy, anchor: "end" };
});
chart(document.getElementById("ch-gpdc"), {
  height: 320, ml: 64, aria: "Runs plotted on the generalized pressure drop correlation",
  x: { log: true, min: 0.06, max: 4, ticks: [0.1, 0.2, 0.5, 1, 2], label: "Flow parameter F_LV" },
  y: { min: 0, max: 2.4, ticks: [0, 0.4, 0.8, 1.2, 1.6, 2, 2.4], label: "Capacity parameter CP" },
  refs: curveRefs,
  annos: [{ x: 0.075, y: 2.28, text: "ΔP/Z, in H2O/ft", anchor: "start" }],
  series: [{ name: "Run", color: "--s1",
    data: P.map(r => ({ x: r.flv, y: r.cp, r })),
    tip: p => `${p.r.L_gpm} gpm · ${p.r.G_scfm} SCFM<br>F<sub>LV</sub> ${f(p.x, 3)}, CP ${f(p.y, 3)}<br>Measured ${f(p.r.dP_meas_ft, 3)} in/ft · GPDC ${f(p.r.dP_pred_ft, 3)}` }],
});

/* ------------------------------------------------------------ objective 2 */
const TS = D.kog_ts.filter(r => r.t <= 820);
legend("leg-ts", [
  { name: "Inlet", color: "--s1", line: true, points: false },
  { name: "Outlet, as recorded", color: "--s2", line: true, points: false },
  { name: "Outlet, corrected", color: "--s3", line: true, points: false },
]);
chart(document.getElementById("ch-ts"), {
  height: 290, aria: "Inlet and outlet CO2 during the NaOH run",
  x: { min: 0, max: 830, ticks: [0, 200, 400, 600, 800], label: "Elapsed time, s" },
  y: { min: 1.6, max: 3.0, ticks: [1.6, 2, 2.4, 2.8], label: "CO2, mol%" },
  bands: [{ x0: 60, x1: 680, label: "Steady window" }],
  vlines: [{ x: 690, text: "Liquid lost", dx: 6, anchor: "start", dy: 30 }],
  crosshair: x => `t = ${x}–${x + 10} s`,
  series: [
    { name: "Inlet", color: "--s1", line: true, points: false, data: TS.map(r => ({ x: r.t, y: r.yin })), tip: p => `Inlet ${f(p.y, 3)}%` },
    { name: "Outlet, as recorded", color: "--s2", line: true, points: false, data: TS.map(r => ({ x: r.t, y: r.yout })), tip: p => `Outlet, as recorded ${f(p.y, 3)}%` },
    { name: "Outlet, corrected", color: "--s3", line: true, points: false, data: TS.map(r => ({ x: r.t, y: r.yout_cor })), tip: p => `Outlet, corrected ${f(p.y, 3)}%` },
  ],
});
chart(document.getElementById("ch-tsdp"), {
  height: 290, aria: "Pressure drop during the NaOH run",
  x: { min: 0, max: 830, ticks: [0, 200, 400, 600, 800], label: "Elapsed time, s" },
  y: { min: 0, max: 2.5, ticks: [0, 0.5, 1, 1.5, 2, 2.5], label: "ΔP, in H2O" },
  bands: [{ x0: 60, x1: 680, label: "Steady window" }],
  vlines: [{ x: 690, text: "Liquid lost", dx: 6, anchor: "start", dy: 30 }],
  crosshair: x => `t = ${x}–${x + 10} s`,
  series: [{ name: "ΔP", color: "--s1", line: true, points: false, data: TS.map(r => ({ x: r.t, y: r.dp })), tip: p => `ΔP ${f(p.y, 2)} in H2O` }],
});

const C = D.calibration;
legend("leg-cal", [
  { name: "No CO2 fed (17 Sep)", color: "--s1" },
  { name: "Water runs (24 Sep, 1 Oct)", color: "--s2" },
  { name: "1 : 1 and fitted line", kind: "ref" },
]);
chart(document.getElementById("ch-cal"), {
  height: 300, aria: "Outlet analyzer reading versus inlet analyzer reading",
  x: { min: 0, max: 3.5, ticks: [0, 1, 2, 3], label: "Inlet analyzer, mol%" },
  y: { min: 0, max: 3.8, ticks: [0, 1, 2, 3], label: "Outlet analyzer, mol%" },
  refs: [
    { data: [{ x: 0, y: 0 }, { x: 3.5, y: 3.5 }], label: "1 : 1", labelAt: { x: 3.2, y: 3.2 }, dx: 6, dy: 12 },
    { data: [{ x: 0, y: C.a }, { x: 3.5, y: C.a + C.b * 3.5 }], label: `out = ${f(C.a, 3)} + ${f(C.b, 3)}·in`, labelAt: { x: 1.2, y: C.a + C.b * 1.2 }, dx: -6, dy: -8, anchor: "end" },
  ],
  series: [
    { name: "No CO2 fed", color: "--s1", data: C.points.filter(p => p.kind !== "water").map(p => ({ x: p.x, y: p.y })), tip: p => `In ${f(p.x, 4)}% · Out ${f(p.y, 4)}%` },
    { name: "Water run", color: "--s2", data: C.points.filter(p => p.kind === "water").map(p => ({ x: p.x, y: p.y })), tip: p => `In ${f(p.x, 3)}% · Out ${f(p.y, 3)}%<br>Ratio ${f(p.y / p.x, 3)}` },
  ],
});
table("tbl-cal", [
  { h: "Run", v: r => r.sheet },
  { h: "Gas, SCFM", v: r => f(r.G_scfm, 1) },
  { h: "Liquid, gpm", v: r => f(r.L_gpm, 2) },
  { h: "In, %", v: r => f(r.yin, 3) },
  { h: "Out, %", v: r => f(r.yout, 3) },
  { h: "Out/In", v: r => f(r.ratio, 3) },
  { h: "Absorption factor A", v: r => f(r.A_factor, 4) },
  { h: "Max removal, %", v: r => f(r.max_removal_pct, 2) },
], D.water);

// KOGa bars
const K = D.kog;
const barsData = [
  { name: "This run, corrected", v: K.corrected.KOGa, sd: K.corrected.KOGa_sd, ours: true },
  { name: "This run, as recorded", v: K.raw.KOGa, sd: K.raw.KOGa_sd, ours: true },
  { name: `Perry 25 mm Intalox scaled to ${f(D.kinetics.oh_mean, 3)} M OH⁻`, v: D.perry_scaled, ours: false },
  { name: "Perry Table 14-3: 38 mm Intalox, 1 N", v: 38, ours: false },
  { name: "Perry Table 14-8: 50 mm Intalox, 1 N", v: D.perry["50 mm ceramic Intalox"], ours: false },
  { name: "Perry Table 14-8: 25 mm Intalox, 1 N", v: D.perry["25 mm ceramic Intalox"], ours: false },
];
legend("leg-bars", [{ name: "This work (±1σ)", color: "--s1", points: true }, { name: "Literature (Perry's)", color: "--muted", points: true }]);
(function () {
  const host = document.getElementById("bars-kog"), max = 50;
  for (const b of barsData) {
    const pct = v => (100 * v / max) + "%";
    const row = h("div", { class: "bar-row" });
    row.appendChild(h("div", { class: "bl" }, b.name));
    const tr = h("div", { class: "bar-track", tabindex: "0", "aria-label": `${b.name}: ${f(b.v, 1)} kmol per hour per cubic metre per atm` });
    tr.appendChild(h("div", { class: "bar-fill", style: `width:${pct(b.v)};background:var(${b.ours ? "--s1" : "--muted"})` }));
    let end = b.v;
    if (b.sd) { tr.appendChild(h("div", { class: "bar-err", style: `left:${pct(b.v - b.sd)};width:${pct(2 * b.sd)}` })); end = b.v + b.sd; }
    tr.appendChild(h("div", { class: "bar-val", style: `left:calc(${pct(end)} + 6px)` }, f(b.v, 1) + (b.sd ? ` ± ${f(b.sd, 1)}` : "")));
    row.appendChild(tr);
    host.appendChild(row);
  }
  const ax = h("div", { class: "bar-axis" });
  ax.appendChild(h("div", { class: "bl", style: "font-size:12px;color:var(--muted)" }, "kmol/(h·m³·atm)"));
  const tk = h("div", { class: "ticks" });
  for (const v of [0, 10, 20, 30, 40, 50]) tk.appendChild(h("span", { style: `left:${2 * v}%` }, v));
  ax.appendChild(tk);
  host.appendChild(ax);
})();

table("tbl-kog", [
  { h: "Basis", v: r => r.name },
  { h: "y<sub>in</sub>, %", v: r => f(r.d.yin, 3) },
  { h: "y<sub>out</sub>, %", v: r => f(r.d.yout, 3) },
  { h: "Removal, %", v: r => f(r.d.removal_pct, 1) },
  { h: "N<sub>OG</sub>", v: r => f(r.d.NOG, 3) },
  { h: "H<sub>OG</sub>, m", v: r => f(r.d.HOG_m, 2) },
  { h: "K<sub>OG</sub>a, kmol/(h·m³·atm)", v: r => `${f(r.d.KOGa, 1)} ± ${f(r.d.KOGa_sd, 1)}` },
  { h: "CO2 absorbed, mol/min", v: r => f(r.d.co2_abs_mol_min, 4) },
  { h: "NaOH converted, %", v: r => f(r.d.naoh_conv_pct, 0) },
], [{ name: "As recorded", d: K.raw }, { name: "Analyzer-corrected", d: K.corrected }]);

/* ------------------------------------------------------------ objective 3 */
const DF = D.driving;
legend("leg-df", [
  { name: "0.076 M NaOH: y − y* (y* ≈ 0), measured", color: "--s1", line: true, points: false },
  { name: "Water: y − y*, predicted", color: "--s2", line: true, points: false },
]);
const ZFT = D.inputs.Z_in / 12;
chart(document.getElementById("ch-df"), {
  height: 300, mr: 24, aria: "Driving force y minus y-star along the bed",
  x: { min: 0, max: ZFT, ticks: [0, 1, 2, 3, 4], label: "Height above bottom of bed, ft  (gas enters at 0, liquid at 4)" },
  y: { min: 0, max: 3, ticks: [0, 0.5, 1, 1.5, 2, 2.5, 3], label: "y − y*, mol% CO2" },
  crosshair: x => `${f(x, 2)} ft above bottom`,
  series: [
    { name: "NaOH", color: "--s1", line: true, points: false,
      data: DF.naoh.map(d => ({ x: +(d.z * ZFT).toFixed(4), y: d.drive, d })),
      tip: p => `NaOH: y ${f(p.d.y, 3)}%, y* 0 → ${f(p.y, 3)}` },
    { name: "Water", color: "--s2", line: true, points: false,
      data: DF.water.map(d => ({ x: +(d.z * ZFT).toFixed(4), y: d.drive, d })),
      tip: p => `Water: y ${f(p.d.y, 3)}%, y* ${f(p.d.ystar, 3)} → ${f(p.y, 3)}` },
  ],
});

const CM = D.conc_model;
legend("leg-kc", [
  { name: "Model", color: "--s1", line: true, points: false },
  { name: "This run (±1σ)", color: "--s2" },
  { name: "Perry 25 mm Intalox, 1 N", color: "--s3" },
]);
chart(document.getElementById("ch-kc"), {
  height: 300, aria: "Predicted KOGa versus NaOH concentration",
  x: { min: -0.02, max: 1.04, ticks: [0, 0.2, 0.4, 0.6, 0.8, 1], label: "Feed NaOH, mol/L" },
  y: { min: 0, max: 90, ticks: [0, 20, 40, 60, 80], label: "K_OG·a, kmol/(h·m³·atm)" },
  series: [
    { name: "Model", color: "--s1", line: true, points: false, data: CM.map(r => ({ x: r.C, y: r.KOGa })), tip: p => `${f(p.x, 3)} M → ${f(p.y, 1)}` },
    { name: "This run", color: "--s2", data: [{ x: D.inputs.C_NaOH, y: K.corrected.KOGa, lo: K.corrected.KOGa - K.corrected.KOGa_sd, hi: K.corrected.KOGa + K.corrected.KOGa_sd }],
      tip: p => `0.076 M: ${f(p.y, 1)} ± ${f(K.corrected.KOGa_sd, 1)}<br>½-in ceramic Intalox` },
    { name: "Perry Table 14-8", color: "--s3", data: [{ x: 1, y: D.perry["25 mm ceramic Intalox"] }],
      tip: p => `1 N, 25% converted: ${f(p.y, 1)}<br>25 mm saddles, higher L and G` },
  ],
});
legend("leg-rc", [
  { name: "Model", color: "--s1", line: true, points: false },
  { name: "This run, 0.076 M", color: "--s2" },
  { name: "Water, predicted", color: "--s3" },
  { name: "NaOH supply limit", kind: "ref" },
]);
const cs = D.c_stoich_full, Aw = DF.A_water * 100;
chart(document.getElementById("ch-rc"), {
  height: 300, aria: "Predicted CO2 removal versus NaOH concentration",
  x: { min: -0.008, max: 0.305, ticks: [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3], label: "Feed NaOH, mol/L" },
  y: { min: 0, max: 100, ticks: [0, 25, 50, 75, 100], label: "CO2 removed, %" },
  refs: [{ data: [{ x: 0, y: Aw }, { x: cs * (1 - Aw / 100), y: 100 }], label: "NaOH supply limit", labelAt: { x: 0.045, y: 50 }, dx: 6, dy: 0 }],
  series: [
    { name: "Model", color: "--s1", line: true, points: false, data: CM.filter(r => r.C <= 0.3).map(r => ({ x: r.C, y: r.removal })), tip: p => `${f(p.x, 3)} M → ${f(p.y, 1)}%` },
    { name: "This run", color: "--s2", data: [{ x: D.inputs.C_NaOH, y: K.corrected.removal_pct, lo: K.raw.removal_pct, hi: K.corrected.removal_pct }],
      tip: p => `${f(p.y, 1)}% corrected (${f(K.raw.removal_pct, 1)}% as recorded)` },
    { name: "Water", color: "--s3", data: [{ x: 0, y: DF.water_removal_pct }], tip: p => `${f(p.y, 2)}%: capacity-limited, A = ${f(DF.A_water, 4)}` },
  ],
});
table("tbl-conc", [
  { h: "Feed NaOH, M", v: r => f(r.C, 3) },
  { h: "K<sub>OG</sub>a, kmol/(h·m³·atm)", v: r => f(r.KOGa, 1) },
  { h: "CO2 removed, %", v: r => f(r.removal, 1) },
], CM);

/* ------------------------------------------------------------ literature table */
table("tbl-lit", [
  { h: "Quantity", cls: "wrap-cell", v: r => r[0] },
  { h: "This work", cls: "wrap-cell", v: r => r[1] },
  { h: "Reference value or statement", cls: "wrap-cell", v: r => r[2] },
  { h: "Source", cls: "wrap-cell", v: r => r[3] },
], __LIT__);
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
