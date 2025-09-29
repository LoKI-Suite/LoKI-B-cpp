import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from os import path
import glob

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
})

base_path = "./default-N2";

matlab = pd.read_csv(path.join(base_path, "matlab/lookUpTableSwarm.txt"), sep=r'\s+')
cpp = pd.read_csv(path.join(base_path, "cpp/lookUpTableSwarm.txt"), sep=r'\s+')

# Normalize column names to match
column_mapping = {
    'RedDiff((ms)^-1)': 'RedDif(1/(ms))',
    'RedMob((msV)^-1)': 'RedMob(1/(msV))',
    'DriftVelocity(ms^-1)': 'DriftVelocity(m/s)',
    'RedTow(m^2)': 'RedTow(m2)',
    'RedAtt(m^2)': 'RedAtt(m2)',
    'RedDiffE(eV(ms)^-1)': 'RedDiffE(eV/(ms))',
    'RedMobE(eV(msV)^-1)': 'RedMobE(eV/(msV))',
    'MeanE(eV)': 'MeanE(eV)',
    'CharE(eV)': 'CharE(eV)',
    'EleTemp(eV)': 'EleTemp(eV)'
}
matlab.rename(columns=column_mapping, inplace=True)

legend_mapping = {
    'RedDif(1/(ms))': r'$D_e N$',
    'RedMob(1/(msV))': r'$\mu_e N$',
    'DriftVelocity(m/s)': r'$v_d$',
    'RedTow(m2)': r'$\alpha / N$',
    'RedAtt(m2)': r'$\eta/N$',
    'RedDiffE(eV/(ms))': r'$D_\epsilon N$',
    'RedMobE(eV/(msV))': r'$\mu_\epsilon N$',
    'MeanE(eV)': r'$\epsilon$',
    'CharE(eV)': r'$u_k$',
    'EleTemp(eV)': r'$T_e$',
    'EEDF(eV^-(3/2))': r'$f(u)$'
}

def relative_error(matlab, cpp):
    return np.abs((matlab - cpp) / matlab)

def get_field(file_name):
    return float(file_name.split("_")[1])

for col in matlab.columns:
    if col == 'RedField(Td)':
        eedf_label = "EEDF(eV^-(3/2))"

        matlab_eedfs = sorted(glob.glob(path.join(base_path, "matlab", "reducedField_*")), key=get_field)
        cpp_eedfs = sorted(glob.glob(path.join(base_path, "cpp", "reducedField_*")), key=get_field)

        eedf_err = []

        for E, matlab_file, cpp_file in zip(matlab[col], matlab_eedfs, cpp_eedfs):
            matlab_eedf = pd.read_csv(path.join(matlab_file, "eedf.txt"), sep=r'\s+')
            cpp_eedf = pd.read_csv(path.join(cpp_file, "eedf.txt"), sep=r'\s+')

            eedf_err.append(np.mean(relative_error(matlab_eedf[eedf_label], cpp_eedf[eedf_label])))

        plt.loglog(matlab["RedField(Td)"], np.array(eedf_err), label=legend_mapping[eedf_label], linestyle="--")
        continue

    error = relative_error(matlab[col], cpp[col])
    plt.loglog(matlab['RedField(Td)'], error, label=legend_mapping[col])

plt.xlabel("Reduced electric field (Td)")
plt.ylabel("Absolute relative error")
plt.grid(True, which='both', linestyle='-', linewidth=0.5)
plt.tick_params(direction='in', which='both')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=False)
plt.savefig(path.join(base_path, "plot.pdf"), bbox_inches="tight")
plt.show()
