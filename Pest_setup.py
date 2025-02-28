#%%
import pyemu
import pypestutils.helpers as helpers
import os
import pandas as pd
import numpy as np
from pest_funcs import *
import flopy
name = 'pomp'
PestDir = os.path.join('..','pest_files')
mfexe = os.path.join('..','exe','mf6.exe')
sen = os.path.join('..','exe','pestpp-sen')
glm = os.path.join('..','exe','pestpp-glm')
scenes = ['2.5_10_closed','10_10_closed','40_10_closed',
          '2.5_10_open','10_10_open','40_10_open']
for scene in scenes:

    OrgDir = os.path.join('..', f'ws_{scene}')
    master_dir = os.path.join('..',f'Master_{scene}');

    zonearray = np.load(os.path.join(OrgDir, 'zonearray.npy'))
    sim = flopy.mf6.MFSimulation.load(sim_ws=OrgDir)
    gwf = sim.get_model(name)


    sr = Get_SR(OrgDir, name)
    pf = pyemu.utils.PstFrom(OrgDir,PestDir, remove_existing=True,spatial_reference=sr)
    fix_k(PestDir, 'pomp.npf_k_layer2.txt')

    ppdf = generate_points_within_circle(350,12.5)
    ppdf = ppdf[(ppdf.x.between(-200,200, inclusive= 'both')) & (ppdf.y >= -100)]

    pp_opt= {'try_use_ppu' : True, 'pp_space' : ppdf, 'search_radius' : 20, 'maxpts_interp' : 5, 'num_threads'  : 1 }
    v = pyemu.geostats.SphVario(contribution=1.0, a=100, anisotropy=1, bearing=0)
    gs = pyemu.geostats.GeoStruct(variograms=v,nugget=0.0,transform = 'log')
    pf.add_parameters('pomp.npf_k_layer2.txt', 'pilotpoints', transform= 'log', geostruct=gs,upper_bound = 100, lower_bound = 0.01,
                    pp_options = pp_opt, zone_array=np.array([zonearray]))
    pf.add_parameters('pomp.sto_ss_layer2.txt', 'pilotpoints', transform= 'log', geostruct=gs,upper_bound = 100, lower_bound = 0.01,
                    pp_options = pp_opt, zone_array=np.array([zonearray]))


    for subdir, dirs, files in os.walk(OrgDir):
        for f in files:
            if f.startswith('head_0') or f.startswith('head_2'):
                pf.add_observations(f,f'{f}.ins',index_cols = 'time',ofile_sep= ',', use_cols='HEAD',zone_array=np.array([1]*50))
    pf.mod_sys_cmds.append(f'{mfexe} {name}.nam')
    pst = pf.build_pst(os.path.join(PestDir,'eg.pst'))
    os.chmod(os.path.join(PestDir,'forward_run.py'), 0o755)

    pst.control_data.noptmax = -1
    pst.control_data.facparmax = 1000
    pst.write(os.path.join(PestDir,'eg.pst'))
    print(f'Expected time for PESTPP-GLM: {0.885*len(ppdf)/8/60} h')
    pyemu.os_utils.start_workers(os.path.join(PestDir), # the folder which contains the "template" PEST dataset
                                os.path.join(glm), #the PEST software version we want to run
                                'eg.pst', # the control file to use with PEST
                                num_workers=8, #how many agents to deploy
                                worker_root='..', #where to deploy the agent directories; relative to where python is running
                                master_dir= master_dir,
                                verbose = True
                                )


    # %%
