# PLOT DRIVEN CAVITY VELOCITIES
from matplotlib import pyplot as plt, font_manager as fm
import numpy as np

# RAYLEIGH 1E3

# Results from Markatos
Markatos_U_RA1E3 = np.loadtxt(open("RA_1E3/ResultsMarkatos/VelocidadesU_Data.txt", "r"), usecols=(0,1))
Markatos_V_RA1E3 = np.loadtxt(open("RA_1E3/ResultsMarkatos/VelocidadesV_Data.txt", "r"), usecols=(0,1))

# Results obtained
Results_U_RA1E3 = np.loadtxt(open("RA_1E3/ResultsObtained/U_Velocities_CDS.txt", "r"), usecols=(0,1), skiprows=1)
Results_V_RA1E3 = np.loadtxt(open("RA_1E3/ResultsObtained/V_Velocities_CDS.txt", "r"), usecols=(0,1), skiprows=1)

Results_U_RA1E3 = np.insert(Results_U_RA1E3, 0, [[0.0, 0.0]], axis=0)
Results_U_RA1E3 = np.append(Results_U_RA1E3, [[1.0, 0.0]], axis=0)

Results_V_RA1E3 = np.insert(Results_V_RA1E3, 0, [[0.0, 0.0]], axis=0)
Results_V_RA1E3 = np.append(Results_V_RA1E3, [[1.0, 0.0]], axis=0)



marker_style_RA1E3 = dict(linestyle='', color='10', markersize=15,
                    markerfacecolor="k", markeredgecolor="k", mew=1.5, ms=3)


# RAYLEIGH 1E5

# Results from Markatos
Markatos_U_RA1E5 = np.loadtxt(open("RA_1E5/ResultsMarkatos/VelocidadesU_Data.txt", "r"), usecols=(0,1))
Markatos_V_RA1E5 = np.loadtxt(open("RA_1E5/ResultsMarkatos/VelocidadesV_Data.txt", "r"), usecols=(0,1))

# Results obtained
Results_U_RA1E5 = np.loadtxt(open("RA_1E5/ResultsObtained/Differentially_U_Results.txt", "r"), usecols=(0,1), skiprows=1)
Results_V_RA1E5 = np.loadtxt(open("RA_1E5/ResultsObtained/Differentially_V_Results.txt", "r"), usecols=(0,1), skiprows=1)

Results_U_RA1E5 = np.insert(Results_U_RA1E5, 0, [[0.0, 0.0]], axis=0)
Results_U_RA1E5 = np.append(Results_U_RA1E5, [[1.0, 0.0]], axis=0)

Results_V_RA1E5 = np.insert(Results_V_RA1E5, 0, [[0.0, 0.0]], axis=0)
Results_V_RA1E5 = np.append(Results_V_RA1E5, [[1.0, 0.0]], axis=0)

marker_style_RA1E5 = dict(linestyle='', color='10', markersize=15,
                    markerfacecolor="r", markeredgecolor="r", mew=1.5, ms=3)


# RYALEIGH 1E6

# Results from Markatos
Markatos_U_RA1E6 = np.loadtxt(open("RA_1E6/ResultsMarkatos/VelocidadesU_Data.txt", "r"), usecols=(0,1))
Markatos_V_RA1E6 = np.loadtxt(open("RA_1E6/ResultsMarkatos/VelocidadesV_Data.txt", "r"), usecols=(0,1))

# Results obtained
#Results_U_RA1E6 = np.loadtxt(open("RA_1E6/ResultsObtained/U_Velocities.txt", "r"), usecols=(0,1), skiprows=1)
#Results_V_RA1E6 = np.loadtxt(open("RA_1E6/ResultsObtained/V_Velocities.txt", "r"), usecols=(0,1), skiprows=1)

marker_style_RA1E6 = dict(linestyle='', color='10', markersize=15,
                    markerfacecolor="r", markeredgecolor="m", mew=1.5, ms=3)


# Figura Velocidades U
Figura_U, ax = plt.subplots()

# Plot Rayleigh 1E3
ax.plot(Markatos_U_RA1E3[:,0], Markatos_U_RA1E3[:,1], marker="+", **marker_style_RA1E3, label="Markatos Rayleigh 1E3")
ax.plot(Results_U_RA1E3[:,0], Results_U_RA1E3[:,1], linewidth=2.0, color="b", label="Results Rayleigh 1E3")

# Plot Rayleigh 1E5
#ax.plot(Markatos_U_RA1E5[:,0], Markatos_U_RA1E5[:,1], marker="+", **marker_style_RA1E5, label="Markatos Rayleigh 1E5")
#ax.plot(Results_U_RA1E5[:,0], Results_U_RA1E5[:,1], linewidth=2.0, color="k", label="Results Rayleigh 1E5")

# Plot Rayleigh 1E6
#ax.plot(Markatos_U_RA1E6[:,0], Markatos_U_RA1E6[:,1], marker="+", **marker_style_RA1E6, label="Markatos Rayleigh 1E6")
#ax.plot(Results_U_RA1E6[:,0], Results_U_RA1E6[:,1], linewidth=2.0, color="b", label="Results Rayleigh 1E6")

ax.set_title("U Velocities Comparison",size=14, fontweight="bold")
ax.set_ylabel("UL/k",size=12, fontweight='bold')
ax.set_xlabel("Coordinate Y (m)",size=12, fontweight='bold')
ax.set_xlim(0.0, 1.0)
ax.grid()
ax.legend()


# Figura Velocidades V
Figura_V, bx = plt.subplots()

# Plot Rayleigh 1E3
bx.plot(Markatos_V_RA1E3[:,1], Markatos_V_RA1E3[:,0], marker="+", **marker_style_RA1E3, label="Markatos Rayleigh 1E3")
bx.plot(Results_V_RA1E3[:,0], Results_V_RA1E3[:,1], linewidth=2.0, color="b", label="Results Rayleigh 1E3")

# Plot Rayleigh 1E5
#bx.plot(Markatos_V_RA1E5[:,1], Markatos_V_RA1E5[:,0], marker="+", **marker_style_RA1E5, label="Markatos Rayleigh 1E5")
#bx.plot(Results_V_RA1E5[:,0], Results_V_RA1E5[:,1], linewidth=2.0, color="k", label="Results Rayleigh 1E5")

# Plot Rayleigh 1E6
#bx.plot(Markatos_V_RA1E6[:,1], Markatos_V_RA1E6[:,0], marker="+", **marker_style_RA1E6, label="Markatos Rayleigh 1E6")
#bx.plot(Results_V_RA1E6[:,1], Results_V_RA1E6[:,0], linewidth=2.0, color="b", label="Results Rayleigh 1E6")

bx.set_title("V Velocities Comparison", size=14, fontweight="bold")
bx.set_ylabel("VL/k",size=12, fontweight='bold')
bx.set_xlabel("Coordinate X (m)",size=12, fontweight='bold')
bx.set_xlim(0.0, 1.0)
bx.grid()
bx.legend()

plt.show()
