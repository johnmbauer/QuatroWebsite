#!/usr/bin/env python

"""
Fit Hill Equation
-----------------

Script to fit gene expression data to the Hill equation for an 
inducible promoter controlled by a repressor protein.

The script outputs a CSV file containing the fit parameters (also printed
to the terminal) and a PNG file of the data and fit line.

Usage: python fit_hill_equation_final.py DATA.CSV Kd_BOUND n_BOUND 
       ymin_i ymax_i OUTPUT_FILENAME

    DATA.csv        = CSV file containing the data (see standard CSV file)
    Kd_BOUND        = Upper bound for Kd, the apparent dissociation constant. 
                      Write "inf" for no upper bound
    n_BOUND         = Upper bound for n, the Hill coefficient. 
                      Write "inf" for no upper bound
    ymin_i          = Initial guess for ymin. Write "None" for no initial guess
    ymax_i          = Initail guess for ymax. Write "None" for no initial guess
    OUTPUT_FILENAME = Base prefix for output filenames

Notes - 
    - The data are first transformed to a log10 scale before fitting to the Hill
      equation using the linear least squares method
    - The Hill equation used is - 
            y = ymin + (ymax - ymin) * (x^n) / (Kd^n + x^n),
      where y = output, ymin = minimal (basal) expression, 
      ymax = maximum expression, x = inducer concentration, 
      Kd = apparent dissociation constant, and n = Hill coefficient
"""

__author__ = "Stephanie Call"
__version__ = '1.1'

# ////////////////////////////////////////
# Python script to fit the activator Hill equation to promoter characterization data
# with given constraints using the SciPy optimization functions.
# This is done to convert the inducer concentrations for the 
# promoter characterization data into RPU units.
# This script also explores the differences between the curve_fit() and 
# least_squares() functions in SciPy

# Inputs -
#       1 - Name of the CSV file containing the data to be fitted. Should have a header line with a column of 
#           'Inducer_conc' and 'RPU_avg' for proper parsing of the input data.
#       2 - Upper bound for K_d, the apparent dissociation constant. If there are no bounds, put 'inf'
#       3 - Upper bound for n, the Hill coefficient. If there are no bounds, put 'inf'
#       4 - Initial guess for ymin. If not needed, put 'None'
#       5 - Initial guess for ymax. If not needed, put 'None'
#       6 - Name of the output CSV file where the values of the fit coefficients will be saved (no ext.)

# Outputs - 
#       1 - A graph of the data and resulting fit line for inspection. You can save it in the popup 
#           window, if desired. Note that it is saved as a PNG as well.
#       2 - The values of the fit coefficients are printed to the terminal and to an output CSV file.

# Notes - 
#       1 - The Kd is the apparent dissociation constant and, therefore, must be positive. 
#           The n is the Hill coefficient and must be positive. >1: positive cooperative binding of ligands, 
#           <1: negative cooperative binding of ligands, =1: no cooperative binding of ligands.
#       2 - All parameters should be above 0, so lower bounds are all set to 0.
#       3 - Currently, this script only supports fitting the Hill equation for repressor conditions. 
#           For example, Ptet and Ptac. Support for fitting activator inducible promoters may be added later
#       4 - This version of the fit_hill_equation.py script contains only the linear fitting method on 
#           a log scale since this was determined to be the best fit for the data and for a more concise output.

# ////////////////////////////////////////
import seaborn as sns
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import itertools
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
import csv
import matplotlib.pyplot as plt
import scipy as scipy
from scipy.optimize import curve_fit
from scipy.optimize import least_squares

#%%Plotting parameters
plt.rcParams["errorbar.capsize"] = 1.0
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['svg.fonttype'] = 'none'
sns.set(font_scale=120)
sns.set_style("ticks")
plt.rcParams['axes.xmargin'] = 0.05  # use 0.05 for the large graph
plt.rcParams['axes.ymargin'] = 1
plt.rcParams["mathtext.default"] = 'regular'
plt.rcParams['axes.linewidth'] = 0.5

plt.tick_params(
    axis='both',          # changes apply to the y-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelright=False,
    direction='out',
    length=2,
    width=2,
    colors='k',
    labelcolor='k') # labels along the bottom edge are off

#For controlling the appearance of ticklabels
params = {'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7}
plt.rcParams.update(params)


#%%
if sys.argv[1] not in ['help', '-h', '-help', '--help']:
    data_file = sys.argv[1]
else:
    print(__doc__)
    sys.exit(0)
Kd_ubound = int(sys.argv[2]) if sys.argv[2].lower() != 'inf' else np.inf
n_ubound = int(sys.argv[3]) if sys.argv[3].lower() != 'inf' else np.inf
ymin_guess = int(sys.argv[4]) if sys.argv[4].lower() != 'none' else 1
ymax_guess = int(sys.argv[5]) if sys.argv[5].lower() != 'none' else 1
outfile = sys.argv[6] + '.csv'

 
# Define the activator Hill equation in a function (** is for exponentials, not ^ (the exclusive OR operator))
# for the SciPy curve_fit function (takes the independent variable first, then parameters)
def hill(x, Kd, n, ymin, ymax):
    return ymin + (ymax - ymin) * (x**n) / (Kd**n + x**n) # Linear scale

def hill_log(x, Kd, n, ymin, ymax):
    return np.log10(ymin + (ymax - ymin) * (x**n) / (Kd**n + x**n)) # Log scale

# Define an equation to calculate the residuals between Hill equation and measured data 
# for the SciPy least_squares equation (takes the parameters first in a vector, then anything else)
def hill_res(params, ind_conc, exp_val):
    """
    Calculate the Hill Equation residuals on a linear scale.
    """
    Kd = params[0]
    n = params[1]
    ymin = params[2]
    ymax = params[3]
    rpu = ymin + (ymax - ymin) * (ind_conc**n) / (Kd**n + ind_conc**n)
    return rpu - exp_val

def hill_res_log(params, ind_conc, exp_val):
    """
    Calculate the Hill Equation residuals on a log scale.
    """
    Kd = params[0]
    n = params[1]
    ymin = params[2]
    ymax = params[3]
    log_rpu = np.log10(ymin + (ymax - ymin) * (ind_conc**n) / (Kd**n + ind_conc**n))
    log_exp = np.log10(exp_val)
    return log_rpu - log_exp

# Read the averaged experimental x and y data from the input data file 
data = pd.read_csv(data_file)
inducer_conc = data['Inducer_conc'] # x data
try:
    RPU = data['RPU_avg']               # y data (RPU mean values)
    sd_data = data['RPU_sd']            # Standard deviation values
except KeyError:
    RPU = data['REU_avg']
    sd_data = data['REU_sd']

# Use the least_squares functions to fit the data to the hill equation (see original script for curve_fit() method)

# The least_squares function takes in the function to calculate the residuals, an array of initial guesses for 
# the parameters, and bounds, if desired. The callable function must be in the form of the vector/list
# of parameters to adjust first, then the other values that may be necessary for the function, 
# including x values and the experimental data.
params0 = [0.1, 1, 1, 1]
res_lsq_lin_log10 = least_squares(hill_res_log, params0, args = (inducer_conc, RPU), bounds = (0, [Kd_ubound, n_ubound, np.inf, np.inf]))
print(f'For the linear least_squares function using a log scale for residuals,\nKd = {res_lsq_lin_log10.x[0]}, n = {res_lsq_lin_log10.x[1]}, ymin = {res_lsq_lin_log10.x[2]}, ymax = {res_lsq_lin_log10.x[3]}')



# Plot the averages of the raw data and the fit curve 
x_minor = mpl.ticker.LogLocator(base=10.0, subs="auto", numticks=999)

fig, ax = plt.subplots(1, figsize=(8, 6))  #fig.7 1.8, 1.6
ax.tick_params(which='both', width=0.5)
ax.tick_params(which='major', length=5) 
ax.tick_params(which='minor', length=3)
ax.tick_params(which='both', direction='out')
ax.tick_params(axis='both', which='major', pad=5)
#axis bounds
ax.set_ylim(10**-4, 10**2)
ax.set_xlim(10**-3, 10**1)
ax.set_yscale('log')
ax.set_xscale('log')
ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))
ax.yaxis.set_minor_locator(x_minor)
ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))
ax.xaxis.set_minor_locator(x_minor)
ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
# plt.setp(ax.get_xticklabels()[::2], visible=False)

#fig.suptitle(sys.argv[6])
for ax in fig.axes:
    ax.loglog(inducer_conc, RPU, 'ro')
    ax.errorbar(inducer_conc, RPU, yerr=sd_data, fmt='none', color='k', linestyle="None", elinewidth=0.5, capthick=0.5, alpha=0.8, zorder=3)
    ax.errorbar(inducer_conc, RPU, yerr=None, fmt='o', mfc = 'None', mec='k', alpha=0.8, markersize=3, zorder=2, mew=0.5)
    ax.tick_params(labelsize=10)
    ax.set_ylabel('Output [RPU]', fontsize=20)
    ax.set_xlabel('aTc [ng/mL]', fontsize=20)
#x_vals = np.logspace(-8, 4, num = 1000) # x values to create model line. calcs values evenly spaced on log scale
#ax.errorbar(x_vals, hill(x_vals, *res_lsq_lin_log10.x), lw=1, alpha=0.8, zorder=1, color = '#6F3996')
#ax.set_title('least_squares\nLinear, log10')



#ax.set_ylim(10**-4, 10**2)
#ax.set_xlim(10**-4, 10**1)
#ax.tick_params(labelsize=10)
#ax.set_ylabel('Output [RPU]', fontsize=20)
#ax.set_xlabel('aTc [ng/mL]', fontsize=20)

plt.tight_layout()
fig.savefig('Pqac'+'.svg')
#plt.savefig('Pqac.svg',transparent=True)
plt.show()




# Write the results of each fit to an output CSV file for inspection
with open(outfile, 'w') as ofile:
    ofile.write('Fit,K_d,n,ymin,ymax\n')
    ofile.write(f'least_squares - linear+log10,{res_lsq_lin_log10.x[0]},{res_lsq_lin_log10.x[1]},{res_lsq_lin_log10.x[2]},{res_lsq_lin_log10.x[3]}\n')
