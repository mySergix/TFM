# PLOT DRIVEN CAVITY VELOCITIES
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


# REYNOLDS 100

# Results from Ghia et Al.
Ghia_U_RE100 = np.loadtxt(open("RE_100/ResultsGhia/U_Velocities.txt", "r"), usecols=(0,1))
Ghia_V_RE100 = np.loadtxt(open("RE_100/ResultsGhia/V_Velocities.txt", "r"), usecols=(0,1))

# Results obtained
Results_U_RE100 = np.loadtxt(open("RE_100/ResultsObtained/HyperbolicMesh/U_Velocities.txt", "r"), usecols=(0,1), skiprows=1)
Results_V_RE100 = np.loadtxt(open("RE_100/ResultsObtained/HyperbolicMesh/V_Velocities.txt", "r"), usecols=(0,1), skiprows=1)

Results_U_RE100 = np.append(Results_U_RE100, [[1.0, 1.0]], axis=0)

Results_V_RE100 = np.insert(Results_V_RE100, 1, [[0.0, 0.0]], axis=0)
Results_V_RE100 = np.append(Results_V_RE100, [[1.0, 0.0]], axis=0)

marker_style_RE100 = dict(linestyle='', color='10', markersize=4,
                    markerfacecolor="k", markeredgecolor="k", mew=1.5, ms=3)

# REYNOLDS 1000

# Results from Ghia et Al.
Ghia_U_RE1000 = np.loadtxt(open("RE_1000/ResultsGhia/U_Velocities.txt", "r"), usecols=(0,1))
Ghia_V_RE1000 = np.loadtxt(open("RE_1000/ResultsGhia/V_Velocities.txt", "r"), usecols=(0,1))

# Results obtained
#Results_U_RE1000 = np.loadtxt(open("RE_1000/ResultsObtained/HyperbolicMesh/U_Velocities.txt", "r"), usecols=(0,1), skiprows=1)
#Results_V_RE1000 = np.loadtxt(open("RE_1000/ResultsObtained/HyperbolicMesh/V_Velocities.txt", "r"), usecols=(0,1), skiprows=1)

#Results_U_RE1000 = np.append(Results_U_RE1000, [[1.0, 1.0]], axis=0)

#Results_V_RE1000 = np.insert(Results_V_RE1000, 1, [[0.0, 0.0]], axis=0)
#Results_V_RE1000 = np.append(Results_V_RE1000, [[1.0, 0.0]], axis=0)

marker_style_RE1000 = dict(linestyle='', color='10', markersize=15,
                    markerfacecolor="r", markeredgecolor="r", mew=1.5, ms=3)


# REYNOLDS 3200

# Results from Ghia et Al.
Ghia_U_RE3200 = np.loadtxt(open("RE_3200/ResultsGhia/U_Velocities.txt", "r"), usecols=(0,1))
Ghia_V_RE3200 = np.loadtxt(open("RE_3200/ResultsGhia/V_Velocities.txt", "r"), usecols=(0,1))

# Results obtained
#Results_U_RE3200 = np.loadtxt(open("RE_3200/ResultsObtained/HyperbolicMesh/U_Velocities.txt", "r"), usecols=(0,1), skiprows=1)
#Results_V_RE3200 = np.loadtxt(open("RE_3200/ResultsObtained/HyperbolicMesh/V_Velocities.txt", "r"), usecols=(0,1), skiprows=1)

#Results_U_RE3200 = np.append(Results_U_RE3200, [[1.0, 1.0]], axis=0)

#Results_V_RE3200 = np.insert(Results_V_RE3200, 1, [[0.0, 0.0]], axis=0)
#Results_V_RE3200 = np.append(Results_V_RE3200, [[1.0, 0.0]], axis=0)

marker_style_RE3200 = dict(linestyle='', color='10', markersize=2,
                    markerfacecolor="r", markeredgecolor="m", mew=1.5, ms=3)


sns.set_theme(style="dark")

# Figura Velocidades U
Figura_U, ax = plt.subplots()

# Plot Reynolds 100
ax.plot(Ghia_U_RE100[:,0], Ghia_U_RE100[:,1], marker="o", **marker_style_RE100, label="Ghia Reynolds 100")
ax.plot(Results_U_RE100[:,0], Results_U_RE100[:,1], linewidth=2.0, color="b",  label="Results Reynolds 100")



# Plot Reynolds 1000
#ax.plot(Ghia_U_RE1000[:,0], Ghia_U_RE1000[:,1], marker="+", **marker_style_RE1000, label="Ghia Reynolds 1000")
#ax.plot(Results_U_RE1000[:,0], Results_U_RE1000[:,1], linewidth=2.0, color="g", label="Results Reynolds 1000")

# Plot Reynolds 3200
#ax.plot(Ghia_U_RE3200[:,0], Ghia_U_RE3200[:,1], marker="+", **marker_style_RE3200, label="Ghia Reynolds 3200")
#ax.plot(Results_U_RE3200[:,0], Results_U_RE3200[:,1], linewidth=2.0, color="p", label="Results Reynolds 3200")

ax.set_title("U Velocities Comparison",size=14, fontweight="bold")
ax.set_ylabel("Velocity U (m/s)",size=12, fontweight='bold')
ax.set_xlabel("Coordinate Y (m)",size=12, fontweight='bold')
ax.set_xlim(0, 1.0)
ax.grid()
ax.legend()

#sns.set_style("whitegrid")

# Figura Velocidades V
#Figura_V, bx = plt.subplots()

# Plot Reynolds 100
#bx.plot(Ghia_V_RE100[:,0], Ghia_V_RE100[:,1], marker="+", **marker_style_RE100, label="Ghia Reynolds 100")
#bx.plot(Results_V_RE100[:,0], Results_V_RE100[:,1], linewidth=2.0, color="b", label="Results Reynolds 100")

# Plot Reynolds 1000
#bx.plot(Ghia_V_RE1000[:,0], Ghia_V_RE1000[:,1], marker="+", **marker_style_RE1000, label="Ghia Reynolds 1000")
#bx.plot(Results_V_RE1000[:,0], Results_V_RE1000[:,1], linewidth=2.0, color="b", label="Results Reynolds 1000")

# Plot Reynolds 3200
#bx.plot(Ghia_V_RE3200[:,0], Ghia_V_RE3200[:,1], marker="+", **marker_style_RE3200, label="Ghia Reynolds 3200")
#bx.plot(Results_V_RE3200[:,0], Results_V_RE3200[:,1], linewidth=2.0, color="b", label="Results Reynolds 3200")

#bx.set_title("V Velocities Comparison", size=14, fontweight="bold")
#bx.set_ylabel("Velocity V (m/s)",size=12, fontweight='bold')
#bx.set_xlabel("Coordinate X (m)",size=12, fontweight='bold')
#bx.set_xlim(0, 1.0)
#bx.grid()
#bx.legend()

plt.show()
