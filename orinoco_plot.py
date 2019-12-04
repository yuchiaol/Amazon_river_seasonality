# ================================================================
# Yu-Chiao @ WHOI Dec 20, 2018
# Surface air temperature plots
# ================================================================

# ================================================================
# Import functions
# ================================================================
import argparse
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab
from math import isnan, radians
from IPython import get_ipython
import sys, os, ast
import xlrd
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from pylab import setp, genfromtxt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import griddata
from datetime import datetime
from scipy import stats
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.io.img_tiles as cimgt
from cartopy.io.img_tiles import StamenTerrain
import matplotlib.path as mpath
from scipy.signal import detrend

sys.path.append('/home/yliang/lib/python_tools/python_functions/whoi_projects/')
import whoi_data_process_f
sys.path.append('/home/yliang/lib/python_tools/python_functions/data_process/')
import data_process_f
import ERA_interim_data_process_f, ORAS4_data_process_f
sys.path.append('/home/yliang/lib/python_tools/python_functions/statistics')
import statistical_f, MCA_f, EOF_f
sys.path.append('/home/yliang/lib/python_tools/python_functions')
import plot_functions

# ================================================================
# Define functions
# ================================================================
def read_river_mask():

    print('read river mask')
    river_mask_tmp = np.loadtxt('/home/yliang/whoi_projects/Amazon_river/mask_new/obidosglobalnew.txt')
    river_mask_tmp = np.flipud(river_mask_tmp)

    river_mask = np.zeros((river_mask_tmp.shape))
    river_mask[:,0:361] = river_mask_tmp[:,360:].copy()
    river_mask[:,361:] = river_mask_tmp[:,:360].copy()
    river_mask[river_mask<=0.] = np.nan

    return river_mask

def read_hybam_river(year_st,year_ed):

    print('read hybam river')
    dirname = '/home/yliang/whoi_projects/Amazon_river/observations/streamflow/'
    filename = 'Amazon_river_obidos_1969_2015.xls'

    yearN = year_ed - year_st + 1
    year_river = np.linspace(year_st,year_ed,yearN)
    year_ref = 1969
    t0 = (year_st-year_ref)*12
    t1 = (year_ed-year_ref)*12+12

    wb = xlrd.open_workbook(dirname+filename)
    sheet = wb.sheet_by_index(0)

    river_discharge = np.zeros(((2015-1969+1)*12))

    for II in range((2015-1969+1)*12):
        river_discharge[II] = sheet.cell_value(II+2,3)

    return river_discharge[t0:t1]

def read_hybam_Orinoco_river(year_ed):

    print('read hybam Orinoco river')
    dirname = '/home/yliang/whoi_projects/Amazon_river/observations/streamflow/'
    filename = 'Orinoco_hybam.xls'

    year_ref = 2003
    yearN = year_ed - year_ref + 1
    year_river = np.linspace(year_st,year_ed,yearN)

    wb = xlrd.open_workbook(dirname+filename)
    sheet = wb.sheet_by_index(0)

    river_discharge = np.zeros((yearN*12))

    for II in range(yearN*12):
        river_discharge[II] = sheet.cell_value(II+2,3)

    return river_discharge

def read_Dai_Trenberth_river_data(year_st, year_ed, river_index):
    dirname = '/stormtrack/data4/yliang/observations/streamflow/Dai_Trenberth/'
    filename = 'coastal-stns-Vol-monthly.updated-Aug2014.nc'

    yearN = year_ed - year_st + 1
    year_river_dai = np.linspace(year_st,year_ed,yearN)
    year_ref = 1900
    t0 = (year_st-year_ref)*12
    t1 = (year_ed-year_ref)*12+12

    f = Dataset(dirname + filename, mode='r')
    flow = f.variables['FLOW'][t0:t1,river_index].data
    flow[flow<0.] = np.nan
    river_name = f.variables['riv_name'][river_index].data
    print(river_name)
#    time = f.variables['time'][t0:t1].data
    f.close()

    return flow

def max_min_select(ts_monthly,year_N):

    print('max min select')
    ts_max = np.zeros((year_N))*np.nan
    ts_min = np.zeros((year_N))*np.nan
    for II in range(year_N-2):
        ts_max[II+1] = ts_monthly[II+1,:].max()
        ts_min[II+1] = ts_monthly[II+1,:].min()

    return ts_max, ts_min

def max_min_month_select(ts_monthly,year_N,t0_max,t1_max,t0_min,t1_min):

    print('max min select')
    ts_max = np.zeros((year_N))*np.nan
    ts_min = np.zeros((year_N))*np.nan
    for II in range(year_N-2):
        ts_max[II+1] = ts_monthly[II+1,t0_max:t1_max+1].mean()
        ts_min[II+1] = ts_monthly[II+1,t0_min:t1_min+1].mean()

    return ts_max, ts_min


def trend_calculate_here(year, ts_temp):

    coeff = np.polyfit(np.linspace(1,len(year[6:-6]),len(year[6:-6])), ts_temp, 1)
    [t, prob] = stats.ttest_ind(np.linspace(1,len(year[6:-6]),len(year[6:-6]))*0, ts_temp)
    sig_label=''
    if prob<0.05: sig_label = '*'

    return coeff, sig_label


def trend_calculate_here_new(year, ts_temp):

    coeff = np.polyfit(np.linspace(1,len(year),len(year)), ts_temp, 1)
    [t, prob] = stats.ttest_ind(np.linspace(1,len(year),len(year))*0, ts_temp)
    sig_label=''
    if prob<0.05: sig_label = '*'

    return coeff, sig_label

# ================================================================
# main starts
# ================================================================
def main(inargs,year_st,year_ed,flag_detrend):
    """Run the program."""

    print('stop here!!!')

# ================================================================
# Describe this script and control flags
# ================================================================
if __name__ == '__main__':
    description='code template'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("run_flag", type=str, help="False: do not run main program, just plotting figure; True: run the main program")

    args = parser.parse_args()

    year_st = 2003
    year_ed = 2014

    flag_detrend = 0

    sig_level = 0.05 

    plt.close('all')

    year_N = year_ed - year_st + 1
    year = np.linspace(year_st,year_ed,year_N)

    month_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    if ast.literal_eval(args.run_flag) == True:


       ts_river_dai = read_Dai_Trenberth_river_data(1979, year_ed-15, 2)/100000.

# Read Hybam Orinoco river data
       ts_river3 = read_hybam_Orinoco_river(year_ed)/100000.
       print(ts_river3)

# Combine Orinoco river
#       ts_river3_dai[-144:] = ts_river3.copy()
#       print(ts_river3_dai)

# Perform 3-month running average
       N = 3
       ts_river_rmean = (ts_river3).copy()
       ts_river_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve((ts_river3), np.ones((N,))/N, mode='valid')
       ts_river3_rmean = ts_river_dai.copy()
       ts_river3_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_river_dai, np.ones((N,))/N, mode='valid')

# Arrange to monthly data
       river_monthly = ts_river_rmean.reshape((year_N,12))
       river3_monthly = ts_river3_rmean.reshape((1999-1979+1,12))

# Arrange max and min values
       [river_max, river_min] = max_min_select(river_monthly,year_N)
       print(river_max.shape) 
       print(year[1:-1])
       print(river_max[1:-1])  
       ts_temp = river_max.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('hybam max', str(reg*10), sig)
       ts_temp = river_min.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('hybam min', str(reg*10), sig)
       ts_temp = (river_max-river_min).copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('hybam max-min', str(reg*10), sig)

       year3 = np.linspace(1979,1999,1999-1979+1)
       year_N = len(year3)

       [river3_max, river3_min] = max_min_select(river3_monthly,year_N)
       print(river3_max.shape)
       print(year3[1:-1])
       print(river3_max[1:-1])
       ts_temp = river3_max.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year3[1:-1], ts_temp[1:-1], sig_level)
       print('Dai max', str(reg*10), sig)
       ts_temp = river3_min.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year3[1:-1], ts_temp[1:-1], sig_level)
       print('Dai min', str(reg*10), sig)
       ts_temp = (river3_max-river3_min).copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year3[1:-1], ts_temp[1:-1], sig_level)
       print('Dai max-min', str(reg*10), sig)

# ================================================================
# Plot figures
# ================================================================
       fig = plt.figure()
       fig.set_size_inches(10, 10, forward=True)
       window = 21
       sig_level = 0.05

       plt.axes([0.13, 0.6, 0.35, 0.2])
       plt.plot(year,river_max,'b-', linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],river_max[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'b-', linewidth=1)
       plt.plot(year,river_min,'b--', linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],river_min[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'b--', linewidth=1)
       plt.ylabel('(m$^3$/s x10$^5$)', fontsize=10)
#       plt.axis([1979, 2014, 0.5, 3])
       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010],['1980','1985','1990','1995','2000','2005','2010'], fontsize=10)
       plt.title('(a) HYBAM Orinoco River Discharge Max. and Min.', fontsize=10)

       plt.axes([0.593, 0.6, 0.35, 0.2])
       plt.plot(year,river_max-river_min, 'b-', linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],river_max[1:-1]-river_min[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'b-', linewidth=1, label='HYBAM')
       plt.legend(loc='upper left', ncol=1)
       plt.ylabel('(m$^3$/s x10$^5$)', fontsize=10)
#       plt.axis([1979, 2014, 0.7, 2])
       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010],['1980','1985','1990','1995','2000','2005','2010'], fontsize=10)
       plt.title('(b) HYBAM Orinoco River Discharge Seasonality (Max. minus Min.)', fontsize=10)


       plt.axes([0.13, 0.3, 0.35, 0.2])
       plt.plot(year3,river3_max,'b-', linewidth=0.5)
       corrcoef = np.polyfit(year3[1:-1],river3_max[1:-1],1)
       plt.plot(year3[1:-1], corrcoef[0]*year3[1:-1]+corrcoef[1], 'b-', linewidth=1)
       plt.plot(year3,river3_min,'b--', linewidth=0.5)
       corrcoef = np.polyfit(year3[1:-1],river3_min[1:-1],1)
       plt.plot(year3[1:-1], corrcoef[0]*year3[1:-1]+corrcoef[1], 'b--', linewidth=1)
       plt.ylabel('(m$^3$/s x10$^5$)', fontsize=10)
#       plt.axis([1979, 2014, 0.5, 3])
       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010],['1980','1985','1990','1995','2000','2005','2010'], fontsize=10)
       plt.title('(c) Dai-Trenberth Orinoco River Discharge Max. and Min.', fontsize=10)

       plt.axes([0.593, 0.3, 0.35, 0.2])
       plt.plot(year3,river3_max-river3_min, 'b-', linewidth=0.5)
       corrcoef = np.polyfit(year3[1:-1],river3_max[1:-1]-river3_min[1:-1],1)
       plt.plot(year3[1:-1], corrcoef[0]*year3[1:-1]+corrcoef[1], 'b-', linewidth=1, label='HYBAM')
       plt.legend(loc='upper left', ncol=1)
       plt.ylabel('(m$^3$/s x10$^5$)', fontsize=10)
#       plt.axis([1979, 2014, 0.7, 2])
       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010],['1980','1985','1990','1995','2000','2005','2010'], fontsize=10)
       plt.title('(d) Dai-Trenberth Orinoco River Discharge Seasonality', fontsize=10)



       plt.savefig('Orinoco_time_series_all_tmp_plot.jpg', format='jpeg', dpi=200)

       plt.show()

    else:

       print('hahahah')


