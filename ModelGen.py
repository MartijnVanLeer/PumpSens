#%% Import stuff
import os
from grid_funcs import *
from sklearn.metrics import r2_score
import shutil

name = 'pomp'

#grid settings
radius = 2000.0 #radius circle
max_area = 50000 # max area at boundary
zonearea = 15
zonedist = 350

areas = [50000,10000,2500,1200, zonearea, zonearea]
dists = [1500,1000,600,zonedist,20]

#Parameter settings
k = np.array([[10,10],[2.5,10],[40,10]]) #K in m/d
D = 10 #in m 
c = 400 #in d
ss = 0.00001 
Q = 1000 # in m³/h
nlay = 3

#Observation settings
SL = np.sqrt(k*D*c) #Spreidingslengte in m
ObsR = [200]#ObsL * SL #Afstand observatieput

#temporal settings
Tlen = 4 #Time (h)
Tsteps = 50 #Timesteps 
Tmult = 1.05 #Time multiplier

VG,GI,gridprops,tri = MakeGrid(os.path.join('..', 'grid'), radius, max_area, areas, dists)
PlotGrid(VG, zonedist, ObsR)
#%% Run MF6 models

# PlotGrid(VG, zonedist, ObsR)
gwf = []
sim = []
mls = []
hmf = pd.DataFrame()
httim = pd.DataFrame()
tdf = pd.DataFrame()
ObsCells = GetObsCells(GI,ObsR)
CellCenters = np.sqrt(VG.xyzcellcenters[1][ObsCells]**2 + VG.xyzcellcenters[0][ObsCells]**2)
rw = np.sqrt(VG.geo_dataframe.area[ObsCells[0]]/np.pi)
for km in k:
    for TBC in(['closed', 'open']):
        k1 = int(km[0]) if km[0] != 2.5 else km[0]
        ws = os.path.join('..', f'ws_{k1}_{int(km[1])}_{TBC}')
        if not os.path.isdir(ws):
            os.makedirs(ws)
        else:
            shutil.rmtree(ws, ignore_errors=True)

        gwft,simt = Init_Modflow(ws, name,gridprops, radius,GI, km,D,c, Q, TBC, ObsR, ss, nlay, Tlen, Tsteps,Tmult )
        save_zonearray(zonedist,VG,ws)

        ml = run_TTIM(gwft,km,c,TBC, rw)

        hds = gwft.output.head().get_alldata()
        t = np.array(gwft.output.head().get_times())
        td = km[1]/24 * t/(ss*ObsR[0]**2)
        for i in [0, 1]:
            lay = 2 if i == 1 else 0
            hmf[f'{km[0]}_{km[1]}_{TBC}_{lay}'] = -hds[:,lay,0, ObsCells[0]]*4*np.pi*km[1]/24*10/Q
            httim[f'{km[0]}_{km[1]}_{TBC}_{lay}'] = -ml.head(CellCenters[0],0, t,i)[0]*4*np.pi*km[1]/24*10/Q
            tdf[f'{km[0]}_{km[1]}_{TBC}_{lay}'] = td

hmf.to_csv(os.path.join('..', 'Results', 'hmf.csv'))
httim.to_csv(os.path.join('..', 'Results', 'httim.csv'))
tdf.to_csv(os.path.join('..', 'Results', 'tdf.csv'))
        
# PlotHeads(gwf, t = Tsteps-1)

#%%Compare with Ttim

import matplotlib as mpl
cm = mpl.colormaps.get_cmap('viridis')
fig,ax =plt.subplots(dpi = 600)

for no,col in enumerate(hmf.columns):
        color = cm(no/len(hmf.columns))
        ax.plot(tdf[col],hmf[col], label = f'mf + {col}', color = color)
        ax.plot(tdf[col],httim[col], label = f'ttim + {col}', linestyle = '--', color = color)
        ax.set_ylabel (r'$s_d = \frac{K_2 4\pi s}{Q}$')
        ax.set_xlabel(r'$t_d= \frac{K_2t}{S_2r^2}$')
# ax.legend(loc = 'lower right', fontsize = 8)
# ax.set_xscale("log", base=10)
# ax.set_yscale("log", base=10)
# ax[0].set_title('Aquifer 1')
# ax[1].set_title('Aquifer 2')


# %%
