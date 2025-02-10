#%%
import pandas as pd
import matplotlib.pyplot as plt
import os 
import numpy as np
import xarray as xr
from pest_funcs import Get_Sens_fields

ds = None
for k1 in [2.5, 10, 40]:
    for top in ['open', 'closed']:
        for PlotFor in ['ss', 'k']:
            ds = Get_Sens_fields(k1,top, PlotFor,ds)

ds.to_netcdf(os.path.join('..', 'Results', 'jacobians.nc'))

# %%
