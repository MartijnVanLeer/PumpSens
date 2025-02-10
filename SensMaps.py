#%%
import xarray as xr
import matplotlib.pyplot as plt
import flopy 
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

#Load jacobians ds (JcoToNetCDF.py)
ds = xr.open_dataset(os.path.join('..','Results','jacobians.nc'))

#Select timesteps
timeid = np.array([5,20,50,-1])
ds = ds.isel(time = timeid)

#%%


l = 0
tops = ['open']
PlotFor = 'k'
vectors = True

cols = 5 if vectors else 4
fig, axs = plt.subplots(len(tops*3),cols, constrained_layout=True, dpi = 600)
fig.set_size_inches(17.15/2.54,5*len(tops))
cm = matplotlib.colormaps['viridis']
cmap = 'bone_r'
handles2 = []
dc = {'2.5' : 'K1 < K2','10' : 'K1 = K2','40' : 'K1 > K2'}
no = 0
for top in tops:# 'open']:
    for k1 in [2.5,10,40]:
        ws = os.path.join('..', f'ws_{k1}_10_{top}')
        sim = flopy.mf6.MFSimulation.load(sim_ws=ws)
        gwf = sim.get_model('pomp')
        mask = np.load(os.path.join( ws, 'zonearray.npy'))
        if vectors:
            buds = gwf.output.budget().get_data(text='DATA-SPDIS')
        for ti,t in enumerate(ds.time):
                # for no, l in enumerate([0]):#(ds.layer):
                    ax = axs[no,ti]
                    mapview = flopy.plot.PlotMapView(model=gwf, ax = ax)
                    data = ds.loc[dict(time = t, layer = l, top = top, k1 = k1)]
                    if PlotFor == 'ss':
                        map = np.where(mask ==1, data.ss.values, 0)
                        vmax = data.ss.loc[{'cellid' : [224,2786]}].max()
                    else:
                        map = np.where(mask ==1, data.k.values, 0)
                        vmax = data.k.loc[{'cellid' : [224,2786]}].max()
                    mv = mapview.plot_array(map, cmap = cmap,ax = ax, vmin = 0, vmax = vmax)
                    if ti == 0:
                        ax.set_ylabel(f'{dc[str(k1)]}\n y [$r/λ$]')
                        ax.set_xlim(-150,150)
                        ax.set_ylim(-100,300)
                        locs = ax.get_yticks()
                        labels = [float(item)/200 for item in locs]
                        ax.set_yticks(locs, labels)
                    else:
                        ax.set_yticks([])
                    ax.set_xlim(-150,150)
                    ax.set_ylim(-100,300)
                    ax.set_aspect('equal')
                    if no == 0:
                        ax.set_title(f'$t_d$ = {t:.1e}', fontsize = 10)  
                    if no != len(tops)*3-1:
                        ax.set_xticks([])
                    else:
                        locs = ax.get_xticks()
                        labels = [float(item)/200 for item in locs]
                        ax.set_xticks(locs, labels)
                        ax.set_xlabel('x [$r/λ$]')
                    ax.set_xlim(-150,150)
                    ax.set_ylim(-100,300)
                    #Wells
                    ax.scatter(0,0, color = 'red', s = 2, label = 'Pumping well')
                    ax.scatter(0,200, color = 'lime', s = 2, label = 'Observation well')

                    #colorbars
                    divider = make_axes_locatable(ax)
                    cax = ax.inset_axes((0.0, 0.05, 0.05, 0.90))
                    if PlotFor =='ss':
                        cbar = plt.colorbar(mv, cax=cax, ticks=[0, data.ss.max()],orientation='vertical', shrink = 0.8)
                        cbar.ax.set_yticklabels([0, f'{vmax:.1e}'], fontsize = 7, alpha = 0.7)
                    else:
                        cbar = plt.colorbar(mv, cax=cax, ticks=[0, data.k.max()],orientation='vertical',shrink = 0.8)
                        cbar.ax.set_yticklabels([0, f'{vmax:.1e}'], fontsize = 7, alpha = 0.7)                   

                    #Vectors
                    if vectors:
                        zbud = buds[timeid[ti]]
                        props = dict(arrowstyle="-|>,head_width=0.4,head_length=0.8", color = cm(ti*0.33),
                        shrinkA=0,shrinkB=0)
                        magnitude = np.sqrt(zbud['qy'][2786+l*14342]**2 + zbud['qz'][2786+l*14342]**2)
                        arr = axs[no,-1].annotate("", xy=(zbud['qy'][2786+l*14342]/magnitude, zbud['qz'][2786+l*14342]/magnitude), xytext = (0,0),arrowprops = props, color = cm(ti*0.33))
                        handles2.append(matplotlib.lines.Line2D([0], [0], color=cm(ti*0.33), marker=r'$\leftarrow$', markersize=15, linestyle='None'))
                        axs[no,-1].set_xlim(-1.0,0.1)
                        axs[no,-1].set_ylim(-1,0.1)
                        axs[no,-1].set_aspect('equal')
                        axs[no,-1].set_xticklabels([])
                        axs[no,-1].set_yticklabels([])
                        if no == len(tops)*3-1:
                            axs[no,-1].set_xlabel('y')
                        axs[no,-1].set_ylabel('z')

        no +=1 

#Legends
handles, labels = axs[0,0].get_legend_handles_labels()
labels2 = []
for i in ds.time.values:
    labels2.append(f"$t_d$ = {float(i):.1e}")
fig.legend(handles + handles2, labels + labels2, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=3)
# fig.tight_layout()


# %%

#%%
fig,ax = plt.subplots()
mapview = flopy.plot.PlotMapView(model=gwf, ax = ax,layer = 2)
mapview.plot_array(zbud['qx'], ax = ax)
ax.set_xlim(-125,125)
ax.set_ylim(-40,250)
