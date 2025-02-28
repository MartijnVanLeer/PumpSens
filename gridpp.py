#%%
import flopy 
import pandas as pd
import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

sim = flopy.mf6.MFSimulation.load(sim_ws=os.path.join('..','ws_10_10_closed' ))
gwf = sim.get_model('pomp')
ppdf = pd.read_csv(os.path.join('..', 'pest_files','ppdf.csv'),sep=" ", names = ['n','m', 'x', 'y','value'])
#%%
fig, ax = plt.subplots(constrained_layout=True)
fig.set_size_inches(8.25*0.393701, 8.25*0.393701)


pmva = flopy.plot.PlotMapView(gwf, ax = ax)
axins = inset_axes(ax, width="30%", height="30%",bbox_to_anchor=(0, 0, 1, 1),
                      bbox_transform=ax.transAxes,
                      loc="lower left")
pmvb = flopy.plot.PlotMapView(gwf, ax = axins)

pmva.plot_grid(linewidth=  0.2)
pmvb.plot_grid(linewidth=  0.2)
pp = ppdf.plot.scatter(x = 'x', y = 'y', ax = ax,s = 1, zorder = 3, label = 'Pilot point', color = 'black')
ax.scatter(0, 0, color = 'red', s = 2, label = 'Pumping well', zorder = 3)

ax.scatter(0, 200, color = 'lime', s = 2, label = 'Observation well', zorder = 4)
r = 450
ax.set_xlim(-r, r)
ax.set_ylim(-r,r)
axins.add_patch(plt.Rectangle((-r,-r), 2*r, 2*r,fc = 'none', ec = 'black', linewidth = 1, zorder = 10))
axins.set_xlim(-2000, 2000)
axins.set_ylim(-2000,2000)

ax.set_aspect('equal')
axins.set_aspect('equal')
locs = ax.get_xticks()
labels = [float(item)/200 for item in locs]
ax.set_xticks(locs, labels)
ax.set_xlabel('$x/λ$')
locs = ax.get_yticks()
labels = [float(item)/200 for item in locs]
ax.set_yticks(locs, labels)
ax.set_ylabel('$y/λ$')
ax.set_xlim(-r, r)
ax.set_ylim(-r,r)

axins.set_xticks([])
axins.set_yticks([])
# mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
ax.legend(framealpha = 1, fontsize = 7, edgecolor = 'black', loc = 'lower right', fancybox = False)
fig.savefig(os.path.join('..','Images','Grid_pilot_points.pdf'), bbox_inches = 'tight')

# %%
