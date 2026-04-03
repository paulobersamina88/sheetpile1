
import math
from pathlib import Path

import pandas as pd
import streamlit as st

st.set_page_config(page_title="Sheet Pile Design", layout="wide")

SECTION_FILE = Path(__file__).parent / "data" / "section_database.csv"

DEFAULT_SECTIONS = [
    {"Type": "W", "Section": "W14X257", "A_in2": 75.6, "d_in": 16.4, "bf_in": 16.0, "tmin_in": 1.18, "Sx_in3": 415.0, "Zx_in3": 468.0, "rx_in": 6.36, "ry_in": 4.08},
    {"Type": "W", "Section": "W14X211", "A_in2": 62.0, "d_in": 15.9, "bf_in": 15.8, "tmin_in": 0.98, "Sx_in3": 343.0, "Zx_in3": 387.0, "rx_in": 6.21, "ry_in": 4.02},
    {"Type": "HSS", "Section": "HSS12X12X5/8", "A_in2": 27.1, "d_in": 12.0, "bf_in": 12.0, "tmin_in": 0.58, "Sx_in3": 97.8, "Zx_in3": 111.0, "rx_in": 4.58, "ry_in": 4.58},
    {"Type": "PIPE", "Section": "PIPE14STD", "A_in2": 13.5, "d_in": 14.0, "bf_in": 14.0, "tmin_in": 0.44, "Sx_in3": 43.6, "Zx_in3": 54.0, "rx_in": 4.92, "ry_in": 4.92},
]

def load_sections() -> pd.DataFrame:
    if SECTION_FILE.exists():
        try:
            df = pd.read_csv(SECTION_FILE)
            required = {"Type", "Section", "A_in2", "d_in", "bf_in", "tmin_in", "Sx_in3", "Zx_in3", "rx_in", "ry_in"}
            if required.issubset(df.columns):
                return df
        except Exception:
            pass
    return pd.DataFrame(DEFAULT_SECTIONS)

def aisc_compression_asd(Ag, rx, ry, Fy, E, K, L_ft):
    """Approximate AISC 360 ASD compression allowable strength."""
    if Ag <= 0 or min(rx, ry) <= 0:
        return 0.0, 0.0, 0.0
    r = min(rx, ry)
    slenderness = K * L_ft * 12.0 / r
    Fe = math.pi**2 * E / max(slenderness**2, 1e-9)
    lambdac = math.sqrt(Fy / Fe) if Fe > 0 else 999.0
    if lambdac <= 1.5:
        Fcr = (0.658 ** (lambdac**2)) * Fy
    else:
        Fcr = 0.877 * Fe
    omega_c = 1.67
    P_allow = Fcr * Ag / omega_c
    return P_allow / 1000.0, slenderness, Fcr

def flexure_allowable_asd(Sx, Zx, Fy, use_plastic=False):
    """Simple ASD bending allowable strength, major axis only."""
    if use_plastic and Zx > 0:
        Mn = Fy * Zx
    else:
        Mn = Fy * Sx
    omega_b = 1.67
    M_allow_kip_in = Mn / omega_b
    return M_allow_kip_in / 12.0

def design_calcs(inputs, section_row):
    H = inputs["Hcant_ft"]
    ws = inputs["ws_psf"]
    Pp = inputs["Pp_psf_per_ft"]
    Pa = inputs["Pa_pcf"]
    Pe = inputs["Pe_psf_per_ft"] if H > 12 else 0.0
    gamma_b = inputs["gamma_b_pcf"]
    D_in = inputs["D_in"]
    S_ft = inputs["S_ft"]
    Fy = inputs["Fy_ksi"]
    gamma_c = inputs["gamma_c_pcf"]
    K = inputs["K"]
    Lc = inputs["Lc_ft"]

    D_ft = D_in / 12.0

    Hb = 0.5 * S_ft * Pa * H**2 / 1000.0
    Hs = ws * S_ft * Pa * H / gamma_b / 1000.0 if gamma_b > 0 else 0.0
    He = 0.5 * S_ft * Pe * H**2 / 1000.0
    P = (S_ft * D_ft * ws + 0.25 * math.pi * gamma_c * D_ft**2 * H) / 1000.0
    V = Hb + Hs + He
    M = (Hb / 3.0 + 2.0 * He / 3.0 + Hs / 2.0) * H

    bf = float(section_row["bf_in"])
    d = float(section_row["d_in"])
    tmin = float(section_row["tmin_in"])
    Ag = float(section_row["A_in2"])
    Sx = float(section_row["Sx_in3"])
    Zx = float(section_row["Zx_in3"])
    rx = float(section_row["rx_in"])
    ry = float(section_row["ry_in"])

    Pc_allow, slenderness, Fcr = aisc_compression_asd(Ag, rx, ry, Fy, 29000.0, K, Lc)
    Mc_allow = flexure_allowable_asd(Sx, Zx, Fy, use_plastic=inputs["use_plastic_flexure"])

    pr_pc = P / Pc_allow if Pc_allow > 0 else float("inf")
    mr_mc = M / Mc_allow if Mc_allow > 0 else float("inf")
    if pr_pc >= 0.2:
        h1_ratio = pr_pc + (8.0 / 9.0) * mr_mc
        h1_eqn = "Pr/Pc + 8/9 · (Mr/Mc)"
    else:
        h1_ratio = pr_pc / 2.0 + mr_mc
        h1_eqn = "Pr/(2Pc) + Mr/Mc"

    A_term = 2.34 * V / (D_ft * (2.0 * Pp * min(max(0.01, 1.0), 12.0)/1000.0)) if False else None  # placeholder removed below
    h_term = M / V if V > 0 else 0.0

    # Solve embedment by iteration because S1 depends on d/3.
    def required_depth(d_guess):
        s1_ksf = 2.0 * Pp * min(d_guess / 3.0, 12.0) / 1000.0
        A = 2.34 * V / (D_ft * s1_ksf) if s1_ksf > 0 and D_ft > 0 else 1e9
        return A / 2.0 * (1.0 + math.sqrt(1.0 + 4.36 * h_term / A))

    d_req = max(H / 2.0, 1.0)
    for _ in range(100):
        d_new = required_depth(d_req)
        if abs(d_new - d_req) < 1e-4:
            break
        d_req = d_new

    S3 = 2.0 * Pp * min(d_req, 12.0) / 1000.0
    S1 = 2.0 * Pp * min(d_req / 3.0, 12.0) / 1000.0
    A = 2.34 * V / (D_ft * S1) if S1 > 0 and D_ft > 0 else 0.0

    return {
        "Hb_k": Hb, "Hs_k": Hs, "He_k": He, "P_k": P, "V_k": V, "M_ftk": M,
        "bf_in": bf, "d_in": d, "tmin_in": tmin,
        "thickness_check_left": 14.0 * tmin, "depth_check_left": d,
        "Pc_allow_k": Pc_allow, "Mc_allow_ftk": Mc_allow,
        "slenderness": slenderness, "Fcr_ksi": Fcr,
        "pr_pc": pr_pc, "mr_mc": mr_mc, "h1_ratio": h1_ratio, "h1_eqn": h1_eqn,
        "d_req_ft": d_req, "S3_ksf": S3, "S1_ksf": S1, "A_term": A, "h_term": h_term,
    }

def result_badge(ok):
    return "✅ ADEQUATE" if ok else "❌ NOT ADEQUATE"

st.title("Steel Sheet Pile / Soldier Pile Retaining Wall Designer")
st.caption("Streamlit tool based on the design logic shown in your sample sheets. Current version focuses on ASD-style force determination, section checks, AISC H1 interaction, and IBC embedment depth iteration.")

with st.sidebar:
    st.header("Project Info")
    project = st.text_input("Project", "Sheet Pile Wall Design")
    client = st.text_input("Client", "")
    job_no = st.text_input("Job No.", "")
    design_by = st.text_input("Design By", "")
    review_by = st.text_input("Review By", "")
    st.divider()

    st.header("Section Source")
    mode = st.radio("Section input mode", ["Built-in database", "Manual properties"], index=0)
    section_df = load_sections()

col1, col2 = st.columns([1.1, 1.0])

with col1:
    st.subheader("Input Data")
    c1, c2 = st.columns(2)
    with c1:
        Hcant_ft = st.number_input("Height of cantilever, Hcant (ft)", min_value=1.0, value=10.0, step=0.5)
        ws_psf = st.number_input("Surcharge weight, ws (psf)", min_value=0.0, value=500.0, step=25.0)
        Pp_psf_per_ft = st.number_input("Allowable lateral soil-bearing in embedment, Pp (psf/ft)", min_value=1.0, value=300.0, step=25.0)
        Pa_pcf = st.number_input("Lateral soil pressure, Pa (pcf)", min_value=1.0, value=35.0, step=1.0)
        Pe_psf_per_ft = st.number_input("Seismic ground shaking, Pe (psf/ft)", min_value=0.0, value=450.0, step=25.0)
        gamma_b_pcf = st.number_input("Soil specific weight, γb (pcf)", min_value=1.0, value=110.0, step=5.0)
        Fy_ksi = st.number_input("Steel yield stress, Fy (ksi)", min_value=20.0, value=50.0, step=1.0)
    with c2:
        D_in = st.number_input("Pile effective diameter, D (in)", min_value=1.0, value=16.0, step=0.5)
        S_ft = st.number_input("Steel pile spacing, S (ft o.c.)", min_value=0.5, value=3.0, step=0.25)
        R_kip_per_ft = st.number_input("Restrained force, R (kip/ft)", min_value=0.0, value=5.0, step=0.5)
        h_rest_ft = st.number_input("Restraint elevation parameter, h (ft)", min_value=0.0, value=2.0, step=0.5)
        gamma_c_pcf = st.number_input("Concrete unit weight, γc (pcf)", min_value=100.0, value=150.0, step=5.0)
        K = st.number_input("Compression K-factor", min_value=0.5, value=1.0, step=0.1)
        Lc_ft = st.number_input("Unbraced compression length for AISC check (ft)", min_value=1.0, value=10.0, step=0.5)

    use_plastic_flexure = st.checkbox("Use plastic section modulus Zx for flexure check", value=False)

    if mode == "Built-in database":
        section = st.selectbox("Steel section", section_df["Section"].tolist(), index=0)
        row = section_df.loc[section_df["Section"] == section].iloc[0]
        st.dataframe(pd.DataFrame([row]).T.rename(columns={row.name: "Value"}), use_container_width=True)
    else:
        st.markdown("#### Manual Section Properties")
        m1, m2 = st.columns(2)
        with m1:
            sec_name = st.text_input("Section name", "Custom Section")
            sec_type = st.selectbox("Section type", ["W", "HSS", "PIPE"])
            A_in2 = st.number_input("Area, A (in²)", min_value=0.1, value=75.6, step=0.1)
            d_in = st.number_input("Overall depth, d (in)", min_value=0.1, value=16.4, step=0.1)
            bf_in = st.number_input("Effective width / bf (in)", min_value=0.1, value=16.0, step=0.1)
        with m2:
            tmin_in = st.number_input("Minimum thickness, tmin (in)", min_value=0.05, value=1.18, step=0.01)
            Sx_in3 = st.number_input("Elastic section modulus, Sx (in³)", min_value=0.1, value=415.0, step=1.0)
            Zx_in3 = st.number_input("Plastic section modulus, Zx (in³)", min_value=0.1, value=468.0, step=1.0)
            rx_in = st.number_input("rx (in)", min_value=0.1, value=6.36, step=0.01)
            ry_in = st.number_input("ry (in)", min_value=0.1, value=4.08, step=0.01)
        row = pd.Series({
            "Type": sec_type, "Section": sec_name, "A_in2": A_in2, "d_in": d_in, "bf_in": bf_in,
            "tmin_in": tmin_in, "Sx_in3": Sx_in3, "Zx_in3": Zx_in3, "rx_in": rx_in, "ry_in": ry_in
        })

inputs = {
    "Hcant_ft": Hcant_ft,
    "ws_psf": ws_psf,
    "Pp_psf_per_ft": Pp_psf_per_ft,
    "Pa_pcf": Pa_pcf,
    "Pe_psf_per_ft": Pe_psf_per_ft,
    "gamma_b_pcf": gamma_b_pcf,
    "Fy_ksi": Fy_ksi,
    "D_in": D_in,
    "S_ft": S_ft,
    "R_kip_per_ft": R_kip_per_ft,
    "h_rest_ft": h_rest_ft,
    "gamma_c_pcf": gamma_c_pcf,
    "K": K,
    "Lc_ft": Lc_ft,
    "use_plastic_flexure": use_plastic_flexure,
}

results = design_calcs(inputs, row)

with col2:
    st.subheader("Design Summary")
    summary_df = pd.DataFrame([
        ["Project", project],
        ["Client", client],
        ["Job No.", job_no],
        ["Section", row["Section"]],
        ["Steel Type", row["Type"]],
        ["Design By", design_by],
        ["Review By", review_by],
    ], columns=["Field", "Value"])
    st.table(summary_df)

    adequate = (results["h1_ratio"] <= 1.0)
    st.metric("Pile Strength Status", result_badge(adequate))
    st.metric("Required Embedment, Hembd (ft)", f'{results["d_req_ft"]:.2f}')
    st.metric("Interaction Ratio", f'{results["h1_ratio"]:.3f}')

st.divider()
st.header("Analysis")

a1, a2, a3 = st.columns(3)
a1.metric("Hb (kips)", f'{results["Hb_k"]:.3f}')
a1.metric("Hs (kips)", f'{results["Hs_k"]:.3f}')
a2.metric("He (kips)", f'{results["He_k"]:.3f}')
a2.metric("P (kips)", f'{results["P_k"]:.3f}')
a3.metric("V (kips)", f'{results["V_k"]:.3f}')
a3.metric("M (ft-kips)", f'{results["M_ftk"]:.3f}')

with st.expander("Show governing formulas", expanded=True):
    st.markdown(f"""
**Pile section forces at cantilever bottom**

- `Hb = 0.5 × S × Pa × Hcant² = {results["Hb_k"]:.3f} kips`
- `Hs = ws × S × Pa × Hcant / γb = {results["Hs_k"]:.3f} kips`
- `He = 0.5 × S × Pe × Hcant² = {results["He_k"]:.3f} kips`  
  *(Automatically set to zero when Hcant ≤ 12 ft, following the note in your sample.)*
- `P = S × D × ws + 0.25 × π × γc × D² × Hcant = {results["P_k"]:.3f} kips`
- `V = Hb + Hs + He = {results["V_k"]:.3f} kips`
- `M = (Hb/3 + 2He/3 + Hs/2) × Hcant = {results["M_ftk"]:.3f} ft-kips`
    """)

st.subheader("Section Limitation Checks")
lim_df = pd.DataFrame([
    ["14 × tmin", f'{results["thickness_check_left"]:.2f} in', f'> bf = {results["bf_in"]:.2f} in', "OK" if results["thickness_check_left"] > results["bf_in"] else "NG"],
    ["tmin", f'{results["tmin_in"]:.2f} in', '> 0.375 in', "OK" if results["tmin_in"] > 0.375 else "NG"],
    ["d", f'{results["depth_check_left"]:.2f} in', '> 8 in', "OK" if results["depth_check_left"] > 8.0 else "NG"],
], columns=["Check", "Provided", "Required", "Status"])
st.dataframe(lim_df, use_container_width=True, hide_index=True)

st.subheader("Combined Compression and Bending Capacity (AISC 360 H1 – approximate ASD implementation)")
h1_df = pd.DataFrame([
    ["Allowable axial capacity, Pc", f'{results["Pc_allow_k"]:.2f} kips'],
    ["Allowable flexural capacity, Mc", f'{results["Mc_allow_ftk"]:.2f} ft-kips'],
    ["Pr/Pc", f'{results["pr_pc"]:.3f}'],
    ["Mr/Mc", f'{results["mr_mc"]:.3f}'],
    ["Interaction equation used", results["h1_eqn"]],
    ["Interaction ratio", f'{results["h1_ratio"]:.3f}'],
    ["Status", "Satisfactory" if results["h1_ratio"] <= 1.0 else "Not satisfactory"],
], columns=["Item", "Value"])
st.table(h1_df)

with st.expander("AISC compression details"):
    st.write(f"Slenderness KL/r = {results['slenderness']:.2f}")
    st.write(f"Computed Fcr = {results['Fcr_ksi']:.2f} ksi")
    st.caption("This is a simplified major-axis-focused section capacity check suitable for a first-pass sizing tool. Final design should still be reviewed against the full code requirements and project-specific geotechnical criteria.")

st.subheader("Embedment Depth (IBC-style iterative solution)")
emb_df = pd.DataFrame([
    ["By trial / iteration, use pile depth d = Hembd", f'{results["d_req_ft"]:.2f} ft'],
    ["Lateral bearing at bottom, S3", f'{results["S3_ksf"]:.2f} ksf'],
    ["Lateral bearing at d/3, S1", f'{results["S1_ksf"]:.2f} ksf'],
    ["A = 2.34P / (D × S1)", f'{results["A_term"]:.2f}'],
    ["h = M / V", f'{results["h_term"]:.2f} ft'],
], columns=["Parameter", "Value"])
st.table(emb_df)

st.info("Current version keeps the restraint inputs R and h available in the interface, but the core force model follows the equations visible in your sample screenshots. If you want, the next version can add an explicit restrained-top analysis branch and LRFD load combinations.")

st.download_button(
    "Download current results as CSV",
    data=pd.DataFrame([{
        "project": project,
        "client": client,
        "section": row["Section"],
        "Hb_k": results["Hb_k"],
        "Hs_k": results["Hs_k"],
        "He_k": results["He_k"],
        "P_k": results["P_k"],
        "V_k": results["V_k"],
        "M_ftk": results["M_ftk"],
        "Pc_allow_k": results["Pc_allow_k"],
        "Mc_allow_ftk": results["Mc_allow_ftk"],
        "interaction_ratio": results["h1_ratio"],
        "embedment_ft": results["d_req_ft"],
    }]).to_csv(index=False).encode("utf-8"),
    file_name="sheet_pile_design_results.csv",
    mime="text/csv",
)
