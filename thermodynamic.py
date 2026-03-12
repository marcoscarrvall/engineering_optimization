import numpy as np
from dataclasses import dataclass
from scipy.optimize import brentq

@dataclass
class Params:
    # Flight & Ambient
    mach: float
    ambient_temp: float
    ambient_press: float

    # Mass flow (no bypass — single spool turbojet)
    mass_flow_rate: float

    # Fixed pressure ratios
    inlet_pressure_ratio: float
    combustor_pressure_ratio: float

    # Combustor: fixed temperature RISE (not fixed TIT)
    delta_T_combustor: float

    # Component efficiencies
    compressor_eff: float
    turbine_eff: float
    mech_eff: float
    combustor_eff: float
    nozzle_eff: float

    # Gas properties
    gas_constant_R: float
    cp_air: float
    kappa_air: float
    cp_gas: float
    kappa_gas: float
    fuel_lhv: float


engine = Params(
    mach=0.78,
    ambient_temp=216.5,
    ambient_press=22632.0,

    mass_flow_rate=16.3,        # core-only flow (no bypass)

    inlet_pressure_ratio=0.99,
    combustor_pressure_ratio=0.96,

    delta_T_combustor=900.0,    # fixed combustor temperature rise (K)

    compressor_eff=0.87,
    turbine_eff=0.90,
    mech_eff=0.995,
    combustor_eff=0.995,
    nozzle_eff=0.98,

    gas_constant_R=287.0,
    cp_air=1005.0,
    kappa_air=1.4,
    cp_gas=1150.0,
    kappa_gas=1.33,
    fuel_lhv=43_000_000.0,
)


# ---------------------------------------------------------------------------
# Thermodynamic helpers
# ---------------------------------------------------------------------------

def compress(Ta, pa, kappa, eff, pr):
    Tb = Ta * (1 + (1/eff) * (pr**((kappa-1)/kappa) - 1))
    return Tb, pr * pa

def expand(work, eff, mdot, cp, Ta, pa, kappa):
    Tb = Ta - work / (mdot * cp)
    pb = (1 - (1/eff) * (1 - Tb/Ta))**(kappa/(kappa-1)) * pa
    return Tb, pb

def nozzle_thrust(pa, Ta, p_amb, kappa, eff, cp, mdot, V0, R):
    if pa <= p_amb:
        return mdot * (-V0)
    cpr = (1 - (1/eff) * (kappa-1)/(kappa+1))**(-kappa/(kappa-1))
    if (pa / p_amb) < cpr:                  # unchoked
        Tb = Ta * (1 + eff * ((p_amb/pa)**((kappa-1)/kappa) - 1))
        dT = Ta - Tb
        if dT <= 0:
            return mdot * (-V0)
        Vb = np.sqrt(2 * cp * dT)
        return mdot * (Vb - V0)
    else:                                    # choked
        pb  = pa / cpr
        Tb  = Ta * (2 / (kappa + 1))
        Vb  = np.sqrt(kappa * R * Tb)
        rho = pb / (R * Tb)
        Ab  = mdot / (rho * Vb)
        return mdot * (Vb - V0) + (pb - p_amb) * Ab


# ---------------------------------------------------------------------------
# Full cycle given compressor PR → (thrust_N, tsfc, T3, TIT)
# ---------------------------------------------------------------------------

def run_cycle(comp_pr: float, p: Params):
    # Station 0: ambient total conditions
    T0 = p.ambient_temp * (1 + (p.kappa_air-1)/2 * p.mach**2)
    p0 = p.ambient_press * (1 + (p.kappa_air-1)/2 * p.mach**2)**(p.kappa_air/(p.kappa_air-1))
    V0 = p.mach * np.sqrt(p.kappa_air * p.gas_constant_R * p.ambient_temp)

    # Station 2: inlet
    T2 = T0
    p2 = p.inlet_pressure_ratio * p0

    # Station 3: compressor exit
    T3, p3 = compress(T2, p2, p.kappa_air, p.compressor_eff, comp_pr)
    work_comp = p.mass_flow_rate * p.cp_air * (T3 - T2)

    # Station 4: combustor exit — fixed ΔT
    T4  = T3 + p.delta_T_combustor
    p4  = p.combustor_pressure_ratio * p3
    fuel = p.mass_flow_rate * p.cp_gas * p.delta_T_combustor / (p.combustor_eff * p.fuel_lhv)
    mdot_gas = p.mass_flow_rate + fuel

    # Station 5: turbine exit — extracts exactly enough work to drive compressor
    work_turb = work_comp / p.mech_eff
    T5, p5 = expand(work_turb, p.turbine_eff, mdot_gas, p.cp_gas, T4, p4, p.kappa_gas)

    if T5 <= 0 or p5 <= 0:
        return None

    # Station 8: nozzle exit
    F = nozzle_thrust(p5, T5, p.ambient_press, p.kappa_gas,
                      p.nozzle_eff, p.cp_gas, mdot_gas, V0, p.gas_constant_R)

    tsfc = fuel / F * 1e6 if F > 0 else np.inf   # g/kNs
    return F, tsfc, T3, T4                         # T4 = TIT


# ---------------------------------------------------------------------------
# Solver: find comp_pr that hits target thrust (brentq)
# ---------------------------------------------------------------------------

def solve(target_thrust_N: float, p: Params,
          pr_min: float = 1.5, pr_max: float = 40.0) -> dict:

    def residual(pr):
        r = run_cycle(pr, p)
        return (r[0] - target_thrust_N) if r else -target_thrust_N

    # Find achievable range and bracket direction
    r_lo = run_cycle(pr_min, p)
    r_hi = run_cycle(pr_max, p)
    t_lo = r_lo[0] if r_lo else 0.0
    t_hi = r_hi[0] if r_hi else 0.0

    t_min_pr, t_max_pr = min(t_lo, t_hi), max(t_lo, t_hi)

    if not (t_min_pr <= target_thrust_N <= t_max_pr):
        raise ValueError(
            f"Target {target_thrust_N/1000:.2f} kN is outside the achievable range.\n"
            f"  PR {pr_min} → {t_lo/1000:.2f} kN\n"
            f"  PR {pr_max} → {t_hi/1000:.2f} kN\n"
            f"Achievable range: {t_min_pr/1000:.2f} – {t_max_pr/1000:.2f} kN"
        )

    pr_sol = brentq(residual, pr_min, pr_max, xtol=1e-8)
    F, tsfc, T3, T4 = run_cycle(pr_sol, p)

    return {
        'comp_pr':    pr_sol,
        'thrust_N':   F,
        'tsfc':       tsfc,
        'T3_K':       T3,
        'TIT_K':      T4,
        'delta_T':    T4 - T3,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    # Show range
    samples = []
    for pr in np.linspace(1.5, 40, 200):
        r = run_cycle(pr, engine)
        if r and r[0] > 0:
            samples.append(r[0])

    print(f"Engine: single-spool turbojet")
    print(f"  Mass flow    : {engine.mass_flow_rate:.1f} kg/s")
    print(f"  Combustor ΔT : {engine.delta_T_combustor:.0f} K (fixed)")
    print(f"  Mach         : {engine.mach}")
    print(f"\nAchievable thrust range: {min(samples)/1000:.2f} – {max(samples)/1000:.2f} kN\n")

    target_kN = float(input("Enter target thrust (kN): "))

    result = solve(target_kN * 1000.0, engine)

    print(f"\n{'='*48}")
    print(f"  Target thrust        : {target_kN:.2f} kN")
    print(f"  Achieved thrust      : {result['thrust_N']/1000:.4f} kN")
    print(f"  ── Design variable ──────────────────────")
    print(f"  Compressor PR        : {result['comp_pr']:.4f}")
    print(f"  ── Outputs ──────────────────────────────")
    print(f"  Compressor exit T3   : {result['T3_K']:.1f} K")
    print(f"  Turbine inlet T4     : {result['TIT_K']:.1f} K")
    print(f"  Combustor ΔT (check) : {result['delta_T']:.1f} K")
    print(f"  TSFC                 : {result['tsfc']:.4f} g/kNs")
    print(f"{'='*48}\n")