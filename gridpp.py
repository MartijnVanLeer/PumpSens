#%%
import flopy 
import pandas as pd
import os
import matplotlib.pyplot as plt

sim = flopy.mf6.MFSimulation.load(sim_ws=os.path.join('..','ws_10_10_closed' ))
gwf = sim.get_model('pomp')
#%%
ppdf = pd.read_csv(os.path.join('..', 'pest_files','ppdf.csv'),sep=" ", names = ['n','m', 'x', 'y','value'])

fig, axs = plt.subplots(2)
fig.set_size_inches(3.5, 6)
pmva = flopy.plot.PlotMapView(gwf, ax = axs[0])
pmvb = flopy.plot.PlotMapView(gwf, ax = axs[1])

pmva.plot_grid(linewidth=  0.3)
pmvb.plot_grid(linewidth=  0.3)
ppdf.plot.scatter(x = 'x', y = 'y', ax = axs[1],s = 2, zorder = 3, label = 'Pilot point', color = 'black', legend = False)
axs[1].scatter(0, 0, color = 'red', s = 2, label = 'Pumping well', zorder = 3)

axs[1].scatter(0, 200, color = 'lime', s = 2, label = 'Observation well', zorder = 4)
axs[1].set_xlim(-350, 350)
axs[1].set_ylim(-200,400)
axs[0].set_xlim(-2000, 2000)
axs[0].set_ylim(-2000,2000)
for ax in axs:
    ax.set_aspect('equal')
    locs = ax.get_xticks()
    labels = [float(item)/200 for item in locs]
    ax.set_xticks(locs, labels)
    ax.set_xlabel('x [$r/λ$]')
    locs = ax.get_yticks()
    labels = [float(item)/200 for item in locs]
    ax.set_yticks(locs, labels)
    ax.set_ylabel('y [$r/λ$]')
fig.tight_layout()
