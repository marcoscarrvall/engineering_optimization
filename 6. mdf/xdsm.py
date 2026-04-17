import pyxdsm.XDSM as xx
from pyxdsm.XDSM import XDSM, OPT, SOLVER, FUNC, IFUNC

# Configuration for A320 Engine Optimization
main = XDSM(use_sfmath=True)

# 1. Systems
main.add_system("opt", OPT, (r"0, 8 \rightarrow 1:", r"\text{Optimizer}"))
main.add_system("MDA", SOLVER, (r"1, 5\rightarrow 2:", r"\text{MDA}", r"\text{Converger}"))
main.add_system("thermo", FUNC, (r"2:", r"\text{Thermodynamics}"))
main.add_system("engine", FUNC, (r"3:", r"\text{Engine Sizing}"))
main.add_system("aero", FUNC, (r"4:", r"\text{Aerodynamics}"))
main.add_system("obj", IFUNC, (r"6:", r"\text{Obj. Function}"))
main.add_system("con", IFUNC, (r"7:", r"\text{Constraints}"))

# 2. Inputs 
main.add_input("opt", (r"V^0, BPR^0", r"\Pi_\text{fan}^0, \Pi_\text{LPC}^0, \Pi_\text{HPC}^0"))
main.add_input("MDA", (r"MTOM_\text{ref}", r"D_\text{ref}", r"\eta_\text{fan, ref}"))
main.add_input("thermo", (r"T_0, p_0, R", r"\text{thermo. constants}", r"N_\text{eng}"))
main.add_input("engine", (r"T_0, R, \gamma", r"\rho_\text{mat}, W_\text{eng, ref}, \eta_\text{fan, ref}", r"N_\text{eng}, OEW_\text{ref}", r"W_\text{pay}, W_\text{fuel}"))
main.add_input("aero", (r"\rho_0", r"D_\text{fan, ref}, L_\text{eng, ref}", r"C_\text{D, ref}, C_\text{Dpar, ref}, C_\text{L, ref}", r"AR, e, S_\text{ref}"))
main.add_input("obj", (r"W_\text{fuel}"))
main.add_input("con", (r"h_\text{eng}", r"\text{h}_\text{min}", r"M_\text{tip, max}", r"TIT_\text{max}"))


# 3. Connections (Data Flow)
# Optimizer to Disciplines
main.connect("opt", "thermo", (r"V, BPR", r"\Pi_\text{fan}, \Pi_\text{LPC}, \Pi_\text{HPC}"))
main.connect("opt", "engine", (r"BPR, \Pi_\text{fan}", r"\Pi_\text{LPC}, \Pi_\text{HPC}"))
main.connect("opt", "aero", r"V")
main.connect("opt", "obj", r"V")

# MDA Converger Loop (Weight consistency)
main.connect("MDA", "thermo", (r"\hat{D}, \hat{\eta}_\text{fan}"))
main.connect("engine", "MDA", (r"MTOM", r"\eta_\text{fan}"))
main.connect("aero", "MDA", r"D")

# Inter-discipline Connections
main.connect("thermo", "engine", r"\dot{m}_\text{tot}")
main.connect("engine", "aero", r"D_\text{fan}, L_\text{eng}, MTOM")

main.connect("thermo", "obj", r"TSFC")
main.connect("engine", "obj", r"MTOM")
main.connect("aero", "obj", r"C_L/C_D")

main.connect("thermo", "con", r"TIT")
main.connect("engine", "con", (r"M_\text{tip}", r"D_\text{nacelle}"))

# Feedback to Optimizer
main.connect("obj", "opt", r"f")
main.connect("con", "opt", r"g")

# 4. Outputs (Final Optimized Design)
main.add_output("opt", (r"V^*, BPR^*", r"\Pi_{fan}^*, \Pi_{LPC}^*, \Pi_{HPC}^*"))
main.add_output("obj", r"f^*")
main.add_output("con", (r"g_1^*, g_2^*, g_3^*"))

main.write("MDF", quiet=True, outdir='mdf')
