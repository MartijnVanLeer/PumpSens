import os
import pypestutils.helpers as helpers
import pandas as pd
import numpy as np
import pyemu
import xarray as xr
def Get_SR(OrgDir, name):
    grb_fname = os.path.join(OrgDir,f"{name}.disv.grb")
    grid_info = helpers.get_2d_grid_info_from_mf6_grb(grb_fname)
    df = pd.DataFrame({'x' : grid_info['x'], 'y':grid_info['y']})
    sr = df.apply(tuple,axis=1).to_dict()
    return sr

def fix_k(PestDir, f):
    p = os.path.join(PestDir, f)
    with open(p, 'r') as file:
        lines = file.readlines()
    
    with open(os.path.join(PestDir,f), 'w') as file:
        for line in lines:
            values = line.split()
            for value in values:
                file.write(value + '\n')

def generate_points_within_circle(radius, min_distance):
    points = []
    x = np.arange(-radius, radius + min_distance, min_distance)
    y = np.arange(-radius, radius + min_distance, min_distance)
    xin = []
    yin = []
    print('Calculating pilot point locations..')
    for xi in x:
        for yi in y:
            if xi**2 + yi**2 <= radius**2:  # Check if the point lies within the circle
                xin.append(xi)
                yin.append(yi)
    df = pd.DataFrame({'x' : xin, 'y' : yin, 'name':list(range(len(xin))), 'zone' : 1, 'value' : 1})
    df['name'] = 'pp' + df.name.astype(str)
    return df


def Get_Sens_fields(k1,top, PlotFor, ds,SL = [200], layers = [0,2]):
    dir = f'Master_{k1}_10_{top}'

    #Read pilot point file
    ppdf = pd.read_csv(os.path.join('..', 'pest_files', 'p_inst0pp.dat.tpl'), sep= ' ',skiprows=1, names = ['name','x', 'y'], usecols=[0,1,2])
    ppdf['index'] = ppdf.name 
    ppdf.set_index('index', inplace= True)

    #read jacobian and get values for ss or k
    jco = pyemu.Jco.from_binary(os.path.join('..', dir, 'eg.jcb'))
    jcodf = abs(jco.to_dataframe().T)
    jco = jcodf[:int(len(jcodf)/2)] if PlotFor == 'k' else jcodf[int(len(jcodf)/2):]

    #get times from jacobian file
    times = np.unique([i.split(':')[-1] for i in jcodf.columns.unique()]).astype(float)
    factors = os.path.join('..', 'pest_files', 'p_inst0pp.fac')
    for tno, t in enumerate(times):
        for no,l in enumerate(layers): 
            for o in SL: 
                tdf = ppdf.copy()
                obsname = f'oname:head_{l}_{o}_otype:lst_usecol:head_time:{t}'            
                tdf[obsname] = jco[obsname].values
                ppdfpath = os.path.join('..', 'pest_files','ppdf.csv')
                tdf.to_csv(ppdfpath, sep = ' ', header = False)
                c = pyemu.geostats.fac2real(ppdfpath, factors_file=factors,out_file=None, fill_value = 0)
                cdf  = c[0]*4*np.pi*10/24*10/1000
                if not isinstance(ds, xr.Dataset):
                    ds = xr.Dataset(coords={'cellid': np.arange(len(cdf)),
                        'k1': [2.5,10,40],
                        'time': times*10/(0.00001*SL[0]**2),
                        'top': ['open', 'closed'],
                        'layer': layers})
                    emp = np.full((len(cdf), len([2.5,10,40]), len(times), len(['open', 'closed']), len(layers)), np.nan)
                    ds['k'] = (['cellid', 'k1', 'time', 'top', 'layer'],emp)
                    ds['ss'] = (['cellid', 'k1', 'time', 'top', 'layer'],emp)
                if PlotFor == 'k':
                    ds['k'].loc[dict(k1=k1, time=t*10/(0.00001*SL[0]**2), top=top, layer=l)] = cdf
                elif PlotFor == 'ss':
                    ds['ss'].loc[dict(k1=k1, time=t*10/(0.00001*SL[0]**2), top=top, layer=l)] = cdf
    return ds