from pyxdsm.XDSM import XDSM, OPT, SOLVER, FUNC, LEFT, IFUNC


# Change `use_sfmath` to False to use computer modern
x = XDSM(use_sfmath=True)

# --- ADD THIS LINE TO SET MAX WIDTH ---
# "3cm" is the width; "align=center" centers the text.
# --------------------------------------

x.add_system("opt", OPT, r"\begin{tabular}{c}  0, 8 $\rightarrow$ 1 \\ Optimizer\end{tabular}")
x.add_system("Thermodynamics", FUNC, "2: Thermodynamics")
x.add_system("Structures", FUNC, "3: Structures")
x.add_system("Emissions", FUNC, "5: Emissions")
x.add_system("Engine Weight", FUNC, "6: Engine Weight")
x.add_system("Objective", IFUNC, "7: Objective Function")

x.add_system(
    "W/S",
    IFUNC,
    r"\begin{tabular}{c} 7: Inequality 1: \\ $W/S \leq (W/S)_{max}$ \end{tabular}"
)

x.add_system(
    "vfuel",
    IFUNC,
    r"\begin{tabular}{c} 7: Inequality 2: \\ $V_{fuel} \leq V_{tank}*f_{tank}$ \end{tabular}"
)
x.add_system(
    "range",
    IFUNC,
    r"\begin{tabular}{c} 7: Equality 1: \\ $R = R_{ref}$ \end{tabular}"
)


x.connect("opt", "Thermodynamics", "1: Design Variables")
x.connect("Thermodynamics", "Engine Weight", r"4: \dot{m}_{total}, BPR, D_{fan}, N_{stages}, \rho_{mat})")
x.connect("Thermodynamics", "Emissions", r"4: P_{T3}, T_{T3}, f, \tau_{res}, \Phi")



x.write("mdf")

import fitz  # PyMuPDF
from PIL import Image

pdf_doc = fitz.open("mdf.pdf")
page = pdf_doc.load_page(0)
pix = page.get_pixmap(dpi=300)
img = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
img.show()