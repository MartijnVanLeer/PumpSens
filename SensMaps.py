#%%
import xarray as xr
import matplotlib.pyplot as plt
import flopy 
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from ConceptProfile import * 
import os

#Load jacobians ds (JcoToNetCDF.py)
ds = xr.open_dataset(os.path.join('..','Results','jacobians.nc'))

#Select timesteps
timeid = np.array([0,15,-1])
ds = ds.isel(time = timeid)

l = 0
tops = ['closed']
PlotFor = 'ss'

vectors = False
concepts = False

cols = len (ds.time)
cols = cols+1 if vectors else cols
cols = cols+2 if concepts else cols
fig, axs = plt.subplots(len(tops*3),cols, constrained_layout=True)
if vectors or concepts:
    fig.set_size_inches(17.15*0.393701,4*len(tops))
else:
    fig.set_size_inches(8.25*0.393701,3.5)


xmin = -175
xmax = 175
ymin = -100
ymax = 300
asp = (ymax - ymin)/(xmax - xmin)


cm = matplotlib.colormaps['viridis']
cmap = 'bone_r'
handles2 = []
dc = {'2.5' : '$K_1$ < $K_2$','10' : '$K_1$ = $K_2$','40' : '$K_1$ > $K_2$'}
no = 0
for top in tops:# 'open']:
    for k1 in [2.5,10,40]:
        ws = os.path.join('..', f'ws_{k1}_10_{top}')
        sim = flopy.mf6.MFSimulation.load(sim_ws=ws)
        gwf = sim.get_model('pomp')
        wellcell = gwf.wel.stress_period_data.data[0][0][0][1]
        obscell = gwf.obs.continuous.data['head_0_200'][0][2][1]
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
                        vmax = data.ss.loc[{'cellid' : [wellcell,obscell]}].max()
                    else:
                        map = np.where(mask ==1, data.k.values, 0)
                        vmax = data.k.loc[{'cellid' : [wellcell,obscell]}].max()
                    mv = mapview.plot_array(map, cmap = cmap,ax = ax, vmin = 0, vmax = vmax)
                    if ti == 0:
                        ax.set_ylabel(f'{dc[str(k1)]}\n $y/λ$', fontsize = 8)
                        ax.set_xlim(xmin,xmax)
                        ax.set_ylim(ymin,ymax)
                        locs = ax.get_yticks()
                        labels = [float(item)/200 for item in locs]
                        ax.set_yticks(locs, labels, fontsize = 8)
                    else:
                        ax.set_yticks([])
                    ax.set_xlim(xmin,xmax)
                    ax.set_ylim(ymin,ymax)
                    ax.set_aspect('equal')
                    if no == 0:
                        ax.set_title(f'$t_d$ = {t/24:.1e}', fontsize = 8)  
                    if no != len(tops)*3-1:
                        ax.set_xticks([])
                    else:
                        locs = ax.get_xticks()
                        labels = [float(item)/200 for item in locs]
                        ax.set_xticks(locs, labels, fontsize = 8)
                        ax.set_xlabel('$x/λ$', fontsize = 8)
                    ax.set_xlim(xmin,xmax)
                    ax.set_ylim(ymin,ymax)
                    #Wells
                    ax.scatter(0,0, color = 'red', s = 2, label = 'Pumping well')
                    ax.scatter(0,200, color = 'lime', s = 2, label = 'Observation well')

                    #colorbars
                    divider = make_axes_locatable(ax)
                    cax = ax.inset_axes((0.0, 0.05, 0.05, 0.90))
                    if PlotFor =='ss':
                        cbar = plt.colorbar(mv, cax=cax, ticks=[0, data.ss.max()],orientation='vertical', shrink = 0.72)
                        cbar.ax.set_yticklabels([0, f'{vmax:.1e}'], fontsize = 6, alpha = 1)
                    else:
                        cbar = plt.colorbar(mv, cax=cax, ticks=[0, data.k.max()],orientation='vertical',shrink = 0.72)
                        cbar.ax.set_yticklabels([0, f'{vmax:.1e}'], fontsize = 6, alpha = 1)                   

                    #Vectors
                    if vectors:
                        zbud = buds[timeid[ti]]
                        props = dict(arrowstyle="->,head_width=0.4,head_length=0.8", color = cm(ti*0.5),
                        shrinkA=0,shrinkB=0)
                        magnitude = np.sqrt(zbud['qy'][obscell+l*int(len(zbud)/3)]**2 + zbud['qz'][obscell+l*int(len(zbud)/3)]**2)
                        arr = axs[no,-3].annotate("", xy=(zbud['qy'][obscell+l*int(len(zbud)/3)]/magnitude, zbud['qz'][obscell+l*int(len(zbud)/3)]/magnitude), xytext = (0,0),arrowprops = props)
                        handles2.append(matplotlib.lines.Line2D([0], [0], color=cm(ti*0.5), marker=r'$\leftarrow$', markersize=15, linestyle='None'))
                        axs[no,-3].set_xlim(-1.0,0.1)
                        axs[no,-3].set_ylim(-1,asp-0.9)
                        axs[no,-3].set_aspect('equal')
                        axs[no,-3].set_xticklabels([])
                        axs[no,-3].set_xticks([])
                        axs[no,-3].set_yticklabels([])
                        axs[no,-3].set_yticks([])
                        if no == 0:
                            axs[no,-3].set_title('Flow direction', fontsize = 8)
                        if no == len(tops)*3-1:
                            axs[no,-3].set_xlabel('y', fontsize = 8)
                        axs[no,-3].set_ylabel('z', fontsize = 8, ha = 'right', labelpad = 2)
                    #Conceptual figures
        if concepts:
            confined = True if top == 'closed' else False
            axs[no,-2] = plot_arrowpic(axs[no,-2], k1,'early',confined)
            axs[no,-1] = plot_arrowpic(axs[no,-1], k1,'late',confined)
            if no == 0:
                axs[no,-2].set_title('Signal (early)', fontsize = 8)
                axs[no,-1].set_title('Signal (late)', fontsize = 8)
            for col in [-1,-2]:
                axs[no,col].set_xticklabels([])
                axs[no,col].set_xticks([])
                axs[no,col].set_yticklabels([])
                axs[no,col].set_yticks([])
                axs[no,col].set_aspect(250/40*asp)
                if no == 2:
                    axs[no,col].set_xlabel('y', fontsize = 8)

        no +=1 

#Legends
handles, labels = axs[0,0].get_legend_handles_labels()
labels2 = []
for i in ds.time.values:
    labels2.append(f"$t_d$ = {float(i):.1e}")
fig.legend(handles + handles2, labels + labels2, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize = 8, framealpha = 1)
# fig.tight_layout()
fig.savefig(os.path.join('..','Images',f'Sensitivity_{tops[0]}_{l}_{PlotFor}.pdf'), bbox_inches = 'tight')
