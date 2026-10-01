import re

with open('src/global.F', 'r') as f:
    text = f.read()

vars = [
    "dgswe", "padapt", "pflag", "gflag", "diorism", "pl", "ph", "px", "slimit", "plimit",
    "pflag2con1", "pflag2con2", "lebesgueP", "fluxtype", "rk_stage", "rk_order",
    "modal_ic", "dghot", "dghotspool", "slopeflag", "slope_weight", "sedflag", "porosity",
    "sevdm", "layers", "reaction_rate", "mnes", "artdif", "kappa", "s0", "uniform_dif",
    "tune_by_hand", "sed_equationX", "sed_equationY", "rainfall"
]

for v in vars:
    for line in text.split('\n'):
        if re.search(r'\b'+v+r'\b', line, re.IGNORECASE) and 'POINTER' in line.upper():
            print(f"{v} is POINTER in global.F: {line.strip()}")
