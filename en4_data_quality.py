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

    year_st = 1979
    year_ed = 2014

    flag_detrend = 0

    sig_level = 0.05 

    plt.close('all')

    year_N = year_ed - year_st + 1
    year = np.linspace(year_st,year_ed,year_N)

    month_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    if ast.literal_eval(args.run_flag) == True:

       year_ref = 1979
       t0 = (year_st - year_ref)*12
       t1 = (year_ed - year_ref)*12 + 12

       ttt = np.linspace(year_st,year_ed+1-1./12.,t1)

       filename = '/stormtrack/data4/yliang/observations/EN4/EN4_2_0_197901_201612.nc'
       f = Dataset(filename, 'r')
       lon_en4 = f.variables['lon'][:].data
       lat_en4 = f.variables['lat'][:].data
       sss_uncertainty = f.variables['salinity_uncertainty'][t0:t1,0,:,:].data
       sss_obs_weight = f.variables['salinity_observation_weights'][t0:t1,0,:,:].data
       sss_uncertainty[sss_uncertainty<-10000] = np.nan
       sss_obs_weight[sss_obs_weight<-10000] = np.nan
       f.close()

# Interpolate to precipitation grid
#       lon_amz1 = 295
#       lon_amz2 = 311
#       lat_amz1 = -3
#       lat_amz2 = 15

       lon_amz1 = 0
       lon_amz2 = 360
       lat_amz1 = -65
       lat_amz2 = 65

       mask_sss = sss_uncertainty[111,:,:]/sss_uncertainty[111,:,:]

       area_so = data_process_f.area_calculate_nonuniform(lon_en4, lat_en4)
       [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(lon_amz1,lon_amz2,lat_amz1,lat_amz2,lon_en4,lat_en4)

       ts_uncertainty = np.zeros((sss_uncertainty.shape[0]))
       ts_obs_weight = np.zeros((sss_obs_weight.shape[0]))
       for NT in range(len(ts_uncertainty)):
           ts_uncertainty[NT] = np.nansum(sss_uncertainty[NT,y1:y2+1,x1:x2+1]*area_so[y1:y2+1,x1:x2+1]*mask_sss[y1:y2+1,x1:x2+1])/np.nansum(area_so[y1:y2+1,x1:x2+1]*mask_sss[y1:y2+1,x1:x2+1])
           ts_obs_weight[NT] = np.nansum(sss_obs_weight[NT,y1:y2+1,x1:x2+1]*area_so[y1:y2+1,x1:x2+1]*mask_sss[y1:y2+1,x1:x2+1])/np.nansum(area_so[y1:y2+1,x1:x2+1]*mask_sss[y1:y2+1,x1:x2+1])

# ================================================================
# Plot figures
# ================================================================
       fig = plt.figure()
       fig.set_size_inches(10, 10, forward=True)

       plt.axes([0.15, 0.55, 0.65, 0.35])
       plt.plot(ttt, ts_obs_weight, 'r-', linewidth=0.7)
       plt.axis([1979, 2014, 0, 1])
       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010],['1980','1985','1990','1995','2000','2005','2010'], fontsize=10)
       plt.title('(a) 5-m Salinity Observation Weight')

       plt.axes([0.15, 0.1, 0.65, 0.35])
       plt.plot(ttt, ts_uncertainty, 'r-', linewidth=0.7)
#       plt.axis([1979, 2014, 0.58, 0.67])
       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010],['1980','1985','1990','1995','2000','2005','2010'], fontsize=10)
       plt.title('(b) 5-m Salinity Uncertainty')

       plt.savefig('en4_data_quality_tmp_plot.jpg', format='jpeg', dpi=200)

       plt.show()

    else:

       print('hahahah')


