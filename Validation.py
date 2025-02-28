#%%
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
cm = mpl.colormaps.get_cmap('viridis')


moddf = pd.read_csv(os.path.join('..', 'Results','hmf.csv'), index_col = 0)
ttimdf = pd.read_csv(os.path.join('..', 'Results','httim.csv'), index_col = 0)
timedf = pd.read_csv(os.path.join('..', 'Results','tdf.csv'), index_col = 0)

fig, ax = plt.subplots(tight_layout = True)

for no,col in enumerate(moddf.columns):
        color = cm(int(0.5*no)/(len(moddf.columns)-2)*2)
        ax.plot(timedf[col],moddf[col], label = f'mf {col}', color = color)
        ax.plot(timedf[col],ttimdf[col], label = f'ttim {col}', linestyle = '--', color = color)
        ax.set_ylabel (r'$s_d = \frac{K_2 4\pi s}{Q}$')
        ax.set_xlabel(r'$t_d= \frac{K_2t}{S_2r^2}$')

fig.savefig(os.path.join('..', 'Images', 'Validation.pdf'))