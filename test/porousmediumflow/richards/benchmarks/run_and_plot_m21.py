#!/usr/bin/env python3

import sys
import numpy as np
import subprocess

# Run benchmark M2.1 infiltration scenarios
#subprocess.run(["./test_richards_benchmark_tpfa", "params_infiltration_sand.input", "-Grid.Refinement", "2"])
#subprocess.run(["./test_richards_benchmark_tpfa", "params_infiltration_loam.input", "-Grid.Refinement", "2"])
#subprocess.run(["./test_richards_benchmark_tpfa", "params_infiltration_clay.input", "-Grid.Refinement", "2"])

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Matplotlib not found so no result plots will be created.")
    sys.exit(0)

# Create Fig. 4def (right column) of Vanderborght 2005
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(5, 10), dpi=72)

# set labels
axes[0].set_title(r"$\theta$")
for ax in axes:
    ax.set_ylabel(r"$\Delta\eta$ (cm)")

# set limits
axes[0].set_ylim([5.0, -10.0])
axes[0].set_xlim([0.0, 0.5])
axes[1].set_ylim([5.0, -10.0])
axes[1].set_xlim([0.0, 0.5])
axes[2].set_ylim([2.0, -5.0])
axes[2].set_xlim([0.0, 0.5])

# first plot
for i, (t, marker) in enumerate(zip([0.1, 0.2, 0.3], ["x", "D", "o"])):
    theta, eta = np.genfromtxt("theta_deltaeta_num_sand_{}.dat".format(i+1)).T
    axes[0].plot(theta, eta, linestyle='none', fillstyle='none', marker=marker, label="{:.1f} days".format(t))
theta, eta = np.genfromtxt("theta_deltaeta_ana_sand_0.dat").T
axes[0].plot(theta, eta, "k", label="Benchmark")
axes[0].legend()

# second plot
for i, (t, marker) in enumerate(zip([0.2, 0.5, 1.0], ["x", "D", "o"])):
    theta, eta = np.genfromtxt("theta_deltaeta_num_loam_{}.dat".format(i+1)).T
    axes[1].plot(theta, eta, linestyle='none', fillstyle='none', marker=marker, label="{:.1f} days".format(t))
theta, eta = np.genfromtxt("theta_deltaeta_ana_loam_0.dat").T
axes[1].plot(theta, eta, "k", label="Benchmark")
axes[1].legend()

# third plot
for i, (t, marker) in enumerate(zip([0.1, 0.2, 0.5], ["x", "D", "o"])):
    theta, eta = np.genfromtxt("theta_deltaeta_num_clay_{}.dat".format(i+1)).T
    axes[2].plot(theta, eta, linestyle='none', fillstyle='none', marker=marker, label="{:.1f} days".format(t))
theta, eta = np.genfromtxt("theta_deltaeta_ana_clay_0.dat").T
axes[2].plot(theta, eta, "k", label="Benchmark")
axes[2].legend()

fig.tight_layout()
plt.savefig("benchmark_infiltration.png")

# Create Fig. of Schnepf et al 2020
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 5), dpi=72, sharey=True)
textOutput = ""

# set labels
axes[0].set_title(r"$\theta$")
axes[0].set_ylabel(r"depth (cm)")
axes[0].set_ylim([-200.0, 0.0])

# first plot
for i, (t, marker) in enumerate(zip([0.1, 0.2, 0.3], ["x", "D", "o"])):
    theta, depth = np.genfromtxt("theta_depth_num_sand_{}.dat".format(i+1)).T
    textOutput += ",".join(f"{i:.10e}" for i in depth) + "\n"
    textOutput += ",".join(f"{i:.10e}" for i in theta) + "\n"
    axes[0].plot(theta, depth, linestyle='none', fillstyle='none', marker=marker, label="{:.1f} days".format(t))
axes[0].legend()

# second plot
for i, (t, marker) in enumerate(zip([0.2, 0.5, 1.0], ["x", "D", "o"])):
    theta, depth = np.genfromtxt("theta_depth_num_loam_{}.dat".format(i+1)).T
    textOutput += ",".join(f"{i:.10e}" for i in depth) + "\n"
    textOutput += ",".join(f"{i:.10e}" for i in theta) + "\n"
    axes[1].plot(theta, depth, linestyle='none', fillstyle='none', marker=marker, label="{:.1f} days".format(t))
axes[1].legend()

# third plot
for i, (t, marker) in enumerate(zip([0.1, 0.2, 0.5], ["x", "D", "o"])):
    theta, depth = np.genfromtxt("theta_depth_num_clay_{}.dat".format(i+1)).T
    textOutput += ",".join(f"{i:.10e}" for i in depth) + "\n"
    textOutput += ",".join(f"{i:.10e}" for i in theta) + "\n"
    axes[2].plot(theta, depth, linestyle='none', fillstyle='none', marker=marker, label="{:.1f} days".format(t))
axes[2].legend()

with open("benchmark_infiltration_depth.txt", "w") as out:
    out.write(textOutput)

fig.tight_layout()
plt.savefig("benchmark_infiltration_depth.png")
