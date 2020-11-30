import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


fpath = '/Users/smd4193/Documents/hx_ratefit_gabe/HHH_rd4_0518.pdb_hx_rate_1500iter.csv'

df = pd.read_csv(fpath, skiprows=5)

plt.plot(df['num'].values, df['fitted_rates'].values, ls='-', marker='o', color='black', mfc='red')
plt.xlabel('Residues (from most to least protected)')
plt.ylabel('Rate')
plt.grid()
plt.savefig(str(fpath).split('.csv')[0]+'.png')
plt.close()

print('heho')