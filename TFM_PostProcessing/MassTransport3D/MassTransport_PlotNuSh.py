# PLOT NUSSELT AND SHERWOOD NUMBERS
from matplotlib import pyplot as plt, font_manager as fm
import numpy as np

# Results from Paper
Results_Paper_NU = np.loadtxt(open("MoistAir/PaperData/Nusselt_Data.txt", "r"), usecols=(0,1))
Results_Paper_SH = np.loadtxt(open("MoistAir/PaperData/Sherwood_Data.txt", "r"), usecols=(0,1))

# Results obtained
Results_Obtained_NU = np.loadtxt(open("MoistAir/ResultsObtained/Channel_Nusselt.txt", "r"), usecols=(0,1), skiprows=1)
Results_Obtained_SH = np.loadtxt(open("MoistAir/ResultsObtained/Channel_Sherwood.txt", "r"), usecols=(0,1), skiprows=1)

marker_style_NU = dict(linestyle='', color='10', markersize=15,
                    markerfacecolor="k", markeredgecolor="k", mew=1.5, ms=3)

marker_style_SH = dict(linestyle='', color='10', markersize=15,
                    markerfacecolor="k", markeredgecolor="r", mew=1.5, ms=3)


# Figura Comparativa
Figura, ax = plt.subplots()

# Plot Nusselt Number
ax.plot(Results_Paper_NU[:,0], Results_Paper_NU[:,1], marker="+", **marker_style_NU, label="Paper Nusselt")
ax.plot(Results_Obtained_NU[:,0], Results_Obtained_NU[:,1], linewidth=2.0, color="b", label="Results Nusselt")

# Plot Nusselt Number
ax.plot(Results_Paper_SH[:,0], Results_Paper_SH[:,1], marker="+", **marker_style_SH, label="Paper Sherwood")
ax.plot(Results_Obtained_SH[:,0], Results_Obtained_SH[:,1], linewidth=2.0, color="g", label="Results Sherwood")


ax.set_xscale("log")
ax.set_title("Nusselt and Sherwood Comparison",size=14, fontweight="bold")
ax.set_ylabel("Nu/Sh (Average)",size=12, fontweight='bold')
ax.set_xlabel("X/(Re Pr Dh)",size=12, fontweight='bold')
ax.grid()
ax.legend()

plt.show()
