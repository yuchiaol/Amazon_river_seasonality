flag_run = 1
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
import shapefile
from scipy.io import loadmat

sys.path.append('/home/yliang/lib/python_tools/python_functions/whoi_projects/')
import whoi_data_process_f
sys.path.append('/home/yliang/lib/python_tools/python_functions/data_process/')
import data_process_f
import ERA_interim_data_process_f
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

    print(river_mask.shape)

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

def plot_box(ax,lon1,lon2,lat1,lat2, color_txt):

    ax.plot(np.linspace(lon1,lon1,100), np.linspace(lat1, lat2, 100), transform=ccrs.PlateCarree(), color=color_txt, linewidth=1)
    ax.plot(np.linspace(lon2,lon2,100), np.linspace(lat1, lat2, 100), transform=ccrs.PlateCarree(), color=color_txt, linewidth=1)
    ax.plot(np.linspace(lon1,lon2,100), np.linspace(lat1, lat1, 100), transform=ccrs.PlateCarree(), color=color_txt, linewidth=1)
    ax.plot(np.linspace(lon1,lon2,100), np.linspace(lat2, lat2, 100), transform=ccrs.PlateCarree(), color=color_txt, linewidth=1)

def plot_figure_here(ax1,lon,lat,map_diff,map_var_diff,map_ttest,clevel,clevel_line,cmap):

    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    factor = 1.

    ax1.set_extent([-180, 180, 30, 90], ccrs.PlateCarree())
    map_2D = map_diff.copy()
    map_2D[0:43,:] = np.nan
    [map_ext, lon_ext] = data_process_f.extend_longitude(map_2D,lat,lon)
    im2      = ax1.contourf(lon_ext, lat, map_ext*factor, levels=clevel, extend='both',cmap=cmap, transform=ccrs.PlateCarree())
    map_2D = map_var_diff.copy()
    map_2D[0:43,:] = np.nan
    [map_ext, lon_ext] = data_process_f.extend_longitude(map_2D,lat,lon)
    im      = ax1.contour(lon_ext, lat, map_ext/100, levels=clevel_line, extend='both', transform=ccrs.PlateCarree(), colors='k', linewidths=0.8)

    var_sig = map_ttest.copy()
    var_sig[0:43,:] = np.nan
    for II in range(len(lat[::2])):
        ax1.plot(lon, lat[II*2]*var_sig[II*2,:], 'ko', markersize=0.5, markerfacecoloralt='k', transform=ccrs.PlateCarree())
#    plot_box(ax1,50,120,50,70, 'c')
    ax1.coastlines('110m')
    ax1.set_boundary(circle, transform=ax1.transAxes)
    ax1.set_aspect('auto')

    return im2

def output_fields_here(filename,river_max,river_min,river3_max,river3_min,prcp_max,prcp_min,sss_max,sss_min,year_N):

    files = os.listdir(os.getcwd())
    for file in files:
        if file == filename:
           print('Delete ' + filename)
           os.system("rm -rf " + file)

    f = Dataset(filename, 'w',format='NETCDF4')
    f.description = 'temporary output fields'
    f.createDimension('dt', year_N)

    river_max_out = f.createVariable('river_max','f4',('dt'))
    river_min_out = f.createVariable('river_min','f4',('dt'))   
    river3_max_out = f.createVariable('river3_max','f4',('dt'))
    river3_min_out = f.createVariable('river3_min','f4',('dt'))
    prcp_max_out = f.createVariable('prcp_max','f4',('dt'))
    prcp_min_out = f.createVariable('prcp_min','f4',('dt'))   
    sss_max_out = f.createVariable('sss_max','f4',('dt'))
    sss_min_out = f.createVariable('sss_min','f4',('dt'))   

    river_max_out[:] = river_max[:]
    river_min_out[:] = river_min[:]
    river3_max_out[:] = river3_max[:]
    river3_min_out[:] = river3_min[:]
    prcp_max_out[:] = prcp_max[:]
    prcp_min_out[:] = prcp_min[:]
    sss_max_out[:] = sss_max[:]
    sss_min_out[:] = sss_min[:]
    f.close()

def read_fields_here(filename):

    f = Dataset(filename, 'r')
    river_max = f.variables['river_max'][:].data
    river_min = f.variables['river_min'][:].data
    river3_max = f.variables['river3_max'][:].data
    river3_min = f.variables['river3_min'][:].data
    prcp_max = f.variables['prcp_max'][:].data
    prcp_min = f.variables['prcp_min'][:].data
    sss_max = f.variables['sss_max'][:].data
    sss_min = f.variables['sss_min'][:].data
    f.close()

    return river_max,river_min,river3_max,river3_min,prcp_max,prcp_min,sss_max,sss_min

def bar_plot_here(ax, tt, ts_input, color_txt):

    for NM in range(12):
        std_value = ts_input[:,NM].std()
        mean_value = np.nanmean(ts_input[:,NM])
        ax.plot([tt[NM], tt[NM]], [-std_value, std_value]+mean_value, color_txt + '-', linewidth=0.5)
        ax.plot([tt[NM]-0.1, tt[NM]+0.1], [-std_value, -std_value]+mean_value, color_txt + '-', linewidth=0.5)
        ax.plot([tt[NM]-0.1, tt[NM]+0.1], [std_value, std_value]+mean_value, color_txt + '-', linewidth=0.5)

def read_hybam_river_here(dirname,filename,year1,year2,year_ref,tag_num):

    N = 3
    year_tmp = np.linspace(year1,year2,year2-year1+1)
    year_N = len(year_tmp)
    t0 = (year_st-year_ref)*12
    t1 = (year_ed-year_ref)*12+12
    wb = xlrd.open_workbook(dirname+filename)
    sheet = wb.sheet_by_index(0)
    river_discharge_tmp = np.zeros(((year_N)*12))
    for II in range(len(river_discharge_tmp)):
        river_discharge_tmp[II] = sheet.cell_value(II+tag_num,3)
    river_discharge_tmp[river_discharge_tmp<-900] = np.nan
    ts_rmean = river_discharge_tmp.copy()
    ts_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(river_discharge_tmp, np.ones((N,))/N, mode='valid')
    river_monthly_tmp = ts_rmean.reshape((year_N,12))             
 
    return year_tmp, river_monthly_tmp


# ================================================================
# main starts
# ================================================================
if flag_run == 1:

    plt.close('all')
    year_st = 1979
    year_ed = 2014


    year_N = year_ed - year_st + 1
    year = np.linspace(year_st,year_ed,year_N)

    month_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    N = 3

    if True:

# Read river mask
#       river_mask = read_river_mask()

#       x = loadmat('river_basin_mask.mat')
#       river_tmp = x['A'].T
#       river_mask = np.zeros((river_tmp.shape))
#       river_mask[:,0:180] = river_tmp[:,180:]
#       river_mask[:,180:] = river_tmp[:,0:180]
#       del river_tmp
#       river_mask[river_mask>1] = np.nan 


       sf = shapefile.Reader('lineaire_1km.shp')
       sf2 = shapefile.Reader('amazlm_1608.shp')

       fig = plt.figure()
       fig.set_size_inches(10, 7, forward=True)

       ax1 = plt.axes([0.28, 0.14, 0.42, 0.65], projection=ccrs.PlateCarree())
       ax1.set_extent([-80, -45, -21, 6])

       for shape in sf.shapeRecords():
           x = [i[0] for i in shape.shape.points[:]]
           y = [i[1] for i in shape.shape.points[:]]
           ax1.plot(x,y,'b-', transform=ccrs.PlateCarree(), linewidth=0.7)
       for shape in sf2.shapeRecords():
           x = [i[0] for i in shape.shape.points[:-10]]
           y = [i[1] for i in shape.shape.points[:-10]]
           ax1.plot(x,y,'k-', transform=ccrs.PlateCarree(), linewidth=0.8)

       ax1.stock_img()
       ax1.coastlines('110m')
#       mask_plot = river_mask.copy()
#       mask_plot[np.isnan(river_mask)==True] = 0.
#       ax1.contour(np.linspace(0,360,360), np.linspace(-90,90,180), mask_plot, [0, 1], colors='k', linewidths=1, transform=ccrs.PlateCarree())

       ax1.set_xticks([-80, -70, -60, -50], crs=ccrs.PlateCarree())
       ax1.set_yticks([-20, -15, -10, -5, 0, 5], crs=ccrs.PlateCarree())
       lon_formatter = LongitudeFormatter(zero_direction_label=True)
       lat_formatter = LatitudeFormatter()
       ax1.xaxis.set_major_formatter(lon_formatter)
       ax1.yaxis.set_major_formatter(lat_formatter)

       ax1.plot([-73.8230, -73.8230], [-10.6637, -10.6637], 'ko', markersize=6, transform=ccrs.PlateCarree())
       ax1.text(-76,-10,'Atalaya Aval', transform=ccrs.PlateCarree(), color='k')
       ax1.plot([-77.5500, -77.5500], [-4.4333, -4.4333], 'mo', markersize=6, transform=ccrs.PlateCarree())
       ax1.text(-78,-4,'Borja', transform=ccrs.PlateCarree(), color='m')
       ax1.plot([-61.1242, -61.1242], [1.8136, 1.8136], 'mo', markersize=6, transform=ccrs.PlateCarree())
       ax1.text(-62,2.3,'Caracarai', transform=ccrs.PlateCarree(), color='m')
       ax1.plot([-59.6225, -59.6225], [-4.3778, -4.3778], 'ko', markersize=6, transform=ccrs.PlateCarree())
       ax1.text(-61,-5.6,'Fazenda Vista Alegre', transform=ccrs.PlateCarree(), color='k')
       ax1.plot([-76.9808, -76.9808], [-0.4753, -0.4753], 'mo', markersize=6, transform=ccrs.PlateCarree())
       ax1.text(-78,0,'Francisco de Orellana', transform=ccrs.PlateCarree(), color='m')
       ax1.plot([-55.9958, -55.9958], [-4.2875, -4.2875], 'mo', markersize=6, transform=ccrs.PlateCarree())
       ax1.text(-55,-4.1875,'Itaituba', transform=ccrs.PlateCarree(), color='m')
       ax1.plot([-64.8111, -64.8111], [-7.2522, -7.2522], 'mo', markersize=6, transform=ccrs.PlateCarree())
       ax1.text(-64.5,-6.8,'Labrea', transform=ccrs.PlateCarree(), color='m')
       ax1.plot([-60.6303, -60.6303], [-3.3122, -3.3122], 'mo', markersize=6, transform=ccrs.PlateCarree())
       ax1.text(-61,-2.8,'Manacapuru', transform=ccrs.PlateCarree(), color='m')
       ax1.plot([-70.0358, -70.0358], [-4.1208, -4.1208], 'ko', markersize=6, transform=ccrs.PlateCarree())
       ax1.text(-71,-5.3,'Nazareth', transform=ccrs.PlateCarree(), color='k')
       ax1.plot([-63.9460, -63.9460], [-8.7997, -8.7997], 'mo', markersize=6, transform=ccrs.PlateCarree())
       ax1.text(-63,-9,'Porto Velho', transform=ccrs.PlateCarree(), color='m')
       ax1.plot([-67.5351, -67.5351], [-14.4410, -14.4410], 'mo', markersize=6, transform=ccrs.PlateCarree())
       ax1.text(-68,-14,'Rurrenabaque', transform=ccrs.PlateCarree(), color='m')
       ax1.plot([-64.8289, -64.8289], [-0.4850, -0.4850], 'mo', markersize=6, transform=ccrs.PlateCarree())
       ax1.text(-64.5,-.45,'Serrinha', transform=ccrs.PlateCarree(), color='m')
       ax1.plot([-69.9622, -69.9622], [-4.2181, -4.2181], 'mo', markersize=6, transform=ccrs.PlateCarree())
       ax1.text(-70,-3.9,'Tabatinga', transform=ccrs.PlateCarree(), color='m')



       ax1.set_aspect('auto')


       plt.savefig('suppl_geographic_map_plot.jpg', format='jpeg', dpi=200)

       plt.show()

       sys.exit()

# Read Borja river data
       print('read hybam Borja river')
       dirname = '/home/yliang/whoi_projects/Amazon_river/observations/streamflow/'
       filename = 'Borja.xls'
       [year_borja, river_monthly_borja] = read_hybam_river_here(dirname,filename,1987,2014,1987,13)
       river_max_borja = np.zeros((len(year_borja)))
       river_min_borja = np.zeros((len(year_borja)))
       for II in range(len(year_borja)):
           river_max_borja[II] = river_monthly_borja[II,:].max()
           river_min_borja[II] = river_monthly_borja[II,:].min()
       index_array = np.isnan(river_min_borja)==False
       corrcoef_borja = np.polyfit(year_borja[index_array],river_max_borja[index_array]-river_min_borja[index_array],1)

       print('read hybam Atalaya river')
       filename = 'Atalaya.xls'
       [year_atalaya, river_monthly_atalaya] = read_hybam_river_here(dirname,filename,2010,2014,2010,12)
       river_max_atalaya = np.zeros((len(year_atalaya)))
       river_min_atalaya = np.zeros((len(year_atalaya)))
       for II in range(len(year_atalaya)):
           river_max_atalaya[II] = river_monthly_atalaya[II,:].max()
           river_min_atalaya[II] = river_monthly_atalaya[II,:].min()
       index_array = np.isnan(river_min_atalaya)==False
       corrcoef_atalaya = np.polyfit(year_atalaya[index_array],river_max_atalaya[index_array]-river_min_atalaya[index_array],1)

       print('read hybam Caracari river')
       filename = 'Caracarai.xls'
       [year_caracarai, river_monthly_caracarai] = read_hybam_river_here(dirname,filename,1979,2014,1979,146)
       river_max_caracarai = np.zeros((len(year_caracarai)))
       river_min_caracarai = np.zeros((len(year_caracarai)))
       for II in range(len(year_caracarai)):
           river_max_caracarai[II] = river_monthly_caracarai[II,:].max()
           river_min_caracarai[II] = river_monthly_caracarai[II,:].min()
       index_array = np.isnan(river_min_caracarai)==False
       corrcoef_caracarai = np.polyfit(year_caracarai[index_array],river_max_caracarai[index_array]-river_min_caracarai[index_array],1)

       print('read hybam Fazenda river')
       filename = 'Fazenda.xls'
       [year_fazenda, river_monthly_fazenda] = read_hybam_river_here(dirname,filename,1979,2014,1979,142)
       river_max_fazenda = np.zeros((len(year_fazenda)))
       river_min_fazenda = np.zeros((len(year_fazenda)))
       for II in range(len(year_fazenda)):
           river_max_fazenda[II] = river_monthly_fazenda[II,:].max()
           river_min_fazenda[II] = river_monthly_fazenda[II,:].min()
       index_array = np.isnan(river_min_fazenda)==False
       corrcoef_fazenda = np.polyfit(year_fazenda[index_array],river_max_fazenda[index_array]-river_min_fazenda[index_array],1)

       print('read hybam Francisco river')
       filename = 'Francisco.xls'
       [year_francisco, river_monthly_francisco] = read_hybam_river_here(dirname,filename,2000,2014,2000,4)
       river_max_francisco = np.zeros((len(year_francisco)))
       river_min_francisco = np.zeros((len(year_francisco)))
       for II in range(len(year_francisco)):
           river_max_francisco[II] = river_monthly_francisco[II,:].max()
           river_min_francisco[II] = river_monthly_francisco[II,:].min()
       index_array = np.isnan(river_min_francisco)==False
       corrcoef_francisco = np.polyfit(year_francisco[index_array],river_max_francisco[index_array]-river_min_francisco[index_array],1)

       print('read hybam Itaituba river')
       filename = 'Itaituba.xls'
       [year_itaituba, river_monthly_itaituba] = read_hybam_river_here(dirname,filename,1979,2014,1979,133)
       river_max_itaituba = np.zeros((len(year_itaituba)))
       river_min_itaituba = np.zeros((len(year_itaituba)))
       for II in range(len(year_itaituba)):
           river_max_itaituba[II] = river_monthly_itaituba[II,:].max()
           river_min_itaituba[II] = river_monthly_itaituba[II,:].min()
       index_array = np.isnan(river_min_itaituba)==False
       corrcoef_itaituba = np.polyfit(year_itaituba[index_array],river_max_itaituba[index_array]-river_min_itaituba[index_array],1)

       print('read hybam Labrea river')
       filename = 'Labrea.xls'
       [year_labrea, river_monthly_labrea] = read_hybam_river_here(dirname,filename,1979,2014,1979,140)
       river_max_labrea = np.zeros((len(year_labrea)))
       river_min_labrea = np.zeros((len(year_labrea)))
       for II in range(len(year_labrea)):
           river_max_labrea[II] = river_monthly_labrea[II,:].max()
           river_min_labrea[II] = river_monthly_labrea[II,:].min()
       index_array = np.isnan(river_min_labrea)==False
       corrcoef_labrea = np.polyfit(year_labrea[index_array],river_max_labrea[index_array]-river_min_labrea[index_array],1)

       print('read hybam Manacapuru river')
       filename = 'Manacapuru.xls'
       [year_manacapuru, river_monthly_manacapuru] = read_hybam_river_here(dirname,filename,1979,2014,1979,81)
       river_max_manacapuru = np.zeros((len(year_manacapuru)))
       river_min_manacapuru = np.zeros((len(year_manacapuru)))
       for II in range(len(year_manacapuru)):
           river_max_manacapuru[II] = river_monthly_manacapuru[II,:].max()
           river_min_manacapuru[II] = river_monthly_manacapuru[II,:].min()
       index_array = np.isnan(river_min_manacapuru)==False
       corrcoef_manacapuru = np.polyfit(year_manacapuru[index_array],river_max_manacapuru[index_array]-river_min_manacapuru[index_array],1)

       print('read hybam Nazareth river')
       filename = 'Nazareth.xls'
       [year_nazareth, river_monthly_nazareth] = read_hybam_river_here(dirname,filename,1990,2005,1990,2)
       river_max_nazareth = np.zeros((len(year_nazareth)))
       river_min_nazareth = np.zeros((len(year_nazareth)))
       for II in range(len(year_nazareth)):
           river_max_nazareth[II] = river_monthly_nazareth[II,:].max()
           river_min_nazareth[II] = river_monthly_nazareth[II,:].min()
       index_array = np.isnan(river_min_nazareth)==False
       corrcoef_nazareth = np.polyfit(year_nazareth[index_array],river_max_nazareth[index_array]-river_min_nazareth[index_array],1)

       print('read hybam Porto Velho river')
       filename = 'Porto_Velho.xls'
       [year_porto, river_monthly_porto] = read_hybam_river_here(dirname,filename,1979,2014,1979,143)
       river_max_porto = np.zeros((len(year_porto)))
       river_min_porto = np.zeros((len(year_porto)))
       for II in range(len(year_porto)):
           river_max_porto[II] = river_monthly_porto[II,:].max()
           river_min_porto[II] = river_monthly_porto[II,:].min()
       index_array = np.isnan(river_min_porto)==False
       corrcoef_porto = np.polyfit(year_porto[index_array],river_max_porto[index_array]-river_min_porto[index_array],1)

       print('read hybam Rurrenabaque river')
       filename = 'Rurrenabaque.xls'
       [year_rurrenabaque, river_monthly_rurrenabaque] = read_hybam_river_here(dirname,filename,1979,2014,1979,139)
       river_max_rurrenabaque = np.zeros((len(year_rurrenabaque)))
       river_min_rurrenabaque = np.zeros((len(year_rurrenabaque)))
       for II in range(len(year_rurrenabaque)):
           river_max_rurrenabaque[II] = river_monthly_rurrenabaque[II,:].max()
           river_min_rurrenabaque[II] = river_monthly_rurrenabaque[II,:].min()
       index_array = np.isnan(river_min_rurrenabaque)==False
       corrcoef_rurrenabaque = np.polyfit(year_rurrenabaque[index_array],river_max_rurrenabaque[index_array]-river_min_rurrenabaque[index_array],1)

       print('read hybam Serrinha river')
       filename = 'Serrinha.xls'
       [year_serrinha, river_monthly_serrinha] = read_hybam_river_here(dirname,filename,1979,2014,1979,136)
       river_max_serrinha = np.zeros((len(year_serrinha)))
       river_min_serrinha = np.zeros((len(year_serrinha)))
       for II in range(len(year_serrinha)):
           river_max_serrinha[II] = river_monthly_serrinha[II,:].max()
           river_min_serrinha[II] = river_monthly_serrinha[II,:].min()
       index_array = np.isnan(river_min_serrinha)==False
       corrcoef_serrinha = np.polyfit(year_serrinha[index_array],river_max_serrinha[index_array]-river_min_serrinha[index_array],1)

       print('read hybam Tabatinga river')
       filename = 'Tabatinga.xls'
       [year_tabatinga, river_monthly_tabatinga] = read_hybam_river_here(dirname,filename,1983,2014,1983,8)
       river_max_tabatinga = np.zeros((len(year_tabatinga)))
       river_min_tabatinga = np.zeros((len(year_tabatinga)))
       for II in range(len(year_tabatinga)):
           river_max_tabatinga[II] = river_monthly_tabatinga[II,:].max()
           river_min_tabatinga[II] = river_monthly_tabatinga[II,:].min()
       index_array = np.isnan(river_min_tabatinga)==False
       corrcoef_tabatinga = np.polyfit(year_tabatinga[index_array],river_max_tabatinga[index_array]-river_min_tabatinga[index_array],1)

# Read Obidos river data
       ts_river = read_hybam_river(year_st, year_ed)/100000.
#       print(ts_river)

# Read Hybam Orinoco river data
#       ts_river3 = read_hybam_Orinoco_river(year_ed)
#       print(ts_river3)

# Arrange to monthly data
       river_monthly = ts_river.reshape((year_N,12))

# ================================================================
# Plot figures
# ================================================================
if True:
       font_size = 10

       plt.close('all')

       fig = plt.figure()
       fig.set_size_inches(10, 10, forward=True)

       ax1 = plt.axes([0.06, 0.85, 0.25, 0.12])
       ax1.plot(year_atalaya, river_max_atalaya-river_min_atalaya, 'b-', label='Atalaya', markersize=3)
       ax1.plot(year_atalaya, corrcoef_atalaya[0]*year_atalaya+corrcoef_atalaya[1], 'r-', markersize=3)
       plt.title('(a) Atalaya Aval, ' + str(round(corrcoef_atalaya[0]*10,2)) + ' m$^3$/s/decade', fontsize=10)

       ax1 = plt.axes([0.06, 0.65, 0.25, 0.12])
       ax1.plot(year_borja, river_max_borja-river_min_borja, 'b-', label='Borja', markersize=3)
       ax1.plot(year_borja, corrcoef_borja[0]*year_borja+corrcoef_borja[1], 'r-', markersize=3)
       plt.title('(b) Borja, ' + str(round(corrcoef_borja[0]*10,2)) + ' m$^3$/s/decade', fontsize=10)

       ax1 = plt.axes([0.06, 0.45, 0.25, 0.12])
       ax1.plot(year_caracarai, river_max_caracarai-river_min_caracarai, 'b-', label='Caracarai', markersize=3)
       ax1.plot(year_caracarai, corrcoef_caracarai[0]*year_caracarai+corrcoef_caracarai[1], 'r-', markersize=3)
       plt.title('(c) Caracarai, ' + str(round(corrcoef_caracarai[0]*10,2)) + ' m$^3$/s/decade', fontsize=10)

       ax1 = plt.axes([0.06, 0.25, 0.25, 0.12])
       ax1.plot(year_fazenda, river_max_fazenda-river_min_fazenda, 'b-', label='Fazenda', markersize=3)
       ax1.plot(year_fazenda, corrcoef_fazenda[0]*year_fazenda+corrcoef_fazenda[1], 'r-', markersize=3)
       plt.title('(d) Fazenda Vista Alegre, ' + str(round(corrcoef_fazenda[0]*10,2)) + ' m$^3$/s/decade', fontsize=10)

       ax1 = plt.axes([0.06, 0.05, 0.25, 0.12])
       ax1.plot(year_francisco, river_max_francisco-river_min_francisco, 'b-', label='Francisco', markersize=3)
       ax1.plot(year_francisco, corrcoef_francisco[0]*year_francisco+corrcoef_francisco[1], 'r-', markersize=3)
       plt.title('(e) Francisco de Orellana, ' + str(round(corrcoef_francisco[0]*10,2)) + ' m$^3$/s/decade', fontsize=10)

       ax1 = plt.axes([0.38, 0.85, 0.25, 0.12])
       ax1.plot(year_itaituba, river_max_itaituba-river_min_itaituba, 'b-', label='Itaituba', markersize=3)
       ax1.plot(year_itaituba, corrcoef_itaituba[0]*year_itaituba+corrcoef_itaituba[1], 'r-', markersize=3)
       plt.title('(f) Itaituba, ' + str(round(corrcoef_itaituba[0]*10,2)) + ' m$^3$/s/decade', fontsize=10)

       ax1 = plt.axes([0.38, 0.65, 0.25, 0.12])
       ax1.plot(year_labrea, river_max_labrea-river_min_labrea, 'b-', label='Labrea', markersize=3)
       ax1.plot(year_labrea, corrcoef_labrea[0]*year_labrea+corrcoef_labrea[1], 'r-', markersize=3)
       plt.title('(g) Labrea, ' + str(round(corrcoef_labrea[0]*10,2)) + ' m$^3$/s/decade', fontsize=10)

       ax1 = plt.axes([0.38, 0.45, 0.25, 0.12])
       ax1.plot(year_manacapuru, river_max_manacapuru-river_min_manacapuru, 'b-', label='Manacapuru', markersize=3)
       ax1.plot(year_manacapuru, corrcoef_manacapuru[0]*year_manacapuru+corrcoef_manacapuru[1], 'r-', markersize=3)
       plt.title('(h) Manacapuru, ' + str(round(corrcoef_manacapuru[0]*10,2)) + ' m$^3$/s/decade', fontsize=10)

       ax1 = plt.axes([0.38, 0.25, 0.25, 0.12])
       ax1.plot(year_nazareth, river_max_nazareth-river_min_nazareth, 'b-', label='Nazareth', markersize=3)
       ax1.plot(year_nazareth, corrcoef_nazareth[0]*year_nazareth+corrcoef_nazareth[1], 'r-', markersize=3)
       plt.title('(i) Nazareth, ' + str(round(corrcoef_nazareth[0]*10,2)) + ' m$^3$/s/decade', fontsize=10)

       ax1 = plt.axes([0.38, 0.05, 0.25, 0.12])
       ax1.plot(year_porto, river_max_porto-river_min_porto, 'b-', label='Porto', markersize=3)
       ax1.plot(year_porto, corrcoef_porto[0]*year_porto+corrcoef_porto[1], 'r-', markersize=3)
       plt.title('(j) Porto Velho, ' + str(round(corrcoef_porto[0]*10,2)) + ' m$^3$/s/decade', fontsize=10)

       ax1 = plt.axes([0.70, 0.85, 0.25, 0.12])
       ax1.plot(year_rurrenabaque, river_max_rurrenabaque-river_min_rurrenabaque, 'b-', label='Rurrenabaque', markersize=3)
       ax1.plot(year_rurrenabaque, corrcoef_rurrenabaque[0]*year_rurrenabaque+corrcoef_rurrenabaque[1], 'r-', markersize=3)
       plt.title('(k) Rurrenabaque, ' + str(round(corrcoef_rurrenabaque[0]*10,2)) + ' m$^3$/s/decade', fontsize=10)

       ax1 = plt.axes([0.70, 0.65, 0.25, 0.12])
       ax1.plot(year_serrinha, river_max_serrinha-river_min_serrinha, 'b-', label='Serrinha', markersize=3)
       ax1.plot(year_serrinha, corrcoef_serrinha[0]*year_serrinha+corrcoef_serrinha[1], 'r-', markersize=3)
       plt.title('(l) Serrinha, ' + str(round(corrcoef_serrinha[0]*10,2)) + ' m$^3$/s/decade', fontsize=10)

       ax1 = plt.axes([0.70, 0.45, 0.25, 0.12])
       ax1.plot(year_tabatinga, river_max_tabatinga-river_min_tabatinga, 'b-', label='Tabatinga', markersize=3)
       ax1.plot(year_tabatinga, corrcoef_tabatinga[0]*year_tabatinga+corrcoef_tabatinga[1], 'r-', markersize=3)
       plt.title('(m) Tabatinga, ' + str(round(corrcoef_tabatinga[0]*10,2)) + ' m$^3$/s/decade', fontsize=10)

       plt.savefig('suppl_other_river.jpg', format='jpeg', dpi=200)

       plt.show()

       sys.exit()

# Figure 1(a)
       lon_obidos = -55.6753
       lat_obidos = -1.9225
       lon_bolivar = -63.5361
       lat_bolivar = 8.1536
       ax1 = plt.axes([0.08, 0.14, 0.42, 0.76], projection=ccrs.PlateCarree())
       ax1.set_extent([-80, -45, -20, 20])
       im = ax1.contourf(lon_salt-360,lat_salt,salt_rg, levels=np.linspace(34.5,37,15), projection=ccrs.PlateCarree(), cmap='viridis', extend='both')
       ax1.contour(lon_salt-360,lat_salt,salt_rg, [35.399, 35.4], colors='m', linewidths=0.6, projection=ccrs.PlateCarree())
       Q = ax1.quiver(lon_u[::scale_in]-360, lat_u[::scale_in], u_rg[::scale_in,::scale_in], v_rg[::scale_in,::scale_in], transform=ccrs.PlateCarree(), scale=v_scale, headwidth=7, headlength=7, width=0.004)
       ax1.quiverkey(Q, 0.46, 0.49, 0.5, str(0.5) + ' m/s', labelpos='N', coordinates='figure', transform=ccrs.PlateCarree())
       ax1.stock_img()
       ax1.coastlines('110m')
       ax1.add_feature(cfeature.RIVERS)
       ax1.plot(lon_obidos, lat_obidos, 'm*', transform=ccrs.PlateCarree())
       ax1.text(-69, -2, 'Obidos Gauge Station', fontsize=font_size, transform=ccrs.PlateCarree(),color='m')
       ax1.plot(lon_bolivar, lat_bolivar, 'b*', transform=ccrs.PlateCarree())
       ax1.text(-76, 6, 'Ciudad Bolivar Gauge Station', fontsize=font_size, transform=ccrs.PlateCarree(), color='b')
       mask_plot = river_mask.copy()
       mask_plot[np.isnan(river_mask)==True] = 0.
       ax1.contour(np.linspace(0,360,721), np.linspace(-90,84,348), mask_plot, [0, 1], colors='k', linewidths=1, transform=ccrs.PlateCarree())
#       ax1.contour(lon_sss, lat_sss-1, so_mask, [0, 1], colors='k', linewidths=1, transform=ccrs.PlateCarree())
       ax1.plot([-49, -49], [-0, 15], 'r-', linewidth=1, transform=ccrs.PlateCarree())
       ax1.plot([-65, -49], [15, 15], 'r-', linewidth=1, transform=ccrs.PlateCarree())
       ax1.plot([-65, -65], [10, 15], 'r-', linewidth=1, transform=ccrs.PlateCarree())
       ax1.set_xticks([-80, -70, -60, -50], crs=ccrs.PlateCarree())
       ax1.set_yticks([-25, -20, -15, -10, -5, 0, 5, 10, 15, 20], crs=ccrs.PlateCarree())
       lon_formatter = LongitudeFormatter(zero_direction_label=True)
       lat_formatter = LatitudeFormatter()
       ax1.xaxis.set_major_formatter(lon_formatter)
       ax1.yaxis.set_major_formatter(lat_formatter)
       ax1.set_aspect('auto')
       plt.title('(a) Geography of Amazon River Basin and Amazon Plume Region', fontsize=font_size)

       cbaxes = fig.add_axes([0.08, 0.08, 0.42, 0.01])
       cbar = plt.colorbar(im, cax=cbaxes, orientation='horizontal', ticks=[34.5, 35, 35.5, 36, 36.5, 37])
       cbar.set_clim(34.5,37)
       cbar.set_label('(PSU)',fontsize=font_size-1)
       setp(cbar.ax.yaxis.get_ticklabels(), fontsize=font_size-1)


# Figure 1(b)
       ax1 = plt.axes([0.58, 0.7, 0.35, 0.2])
       ax1.plot(tt, np.nanmean(prcp_gpcp_monthly, axis=0), 'ro-', label='GPCP', markersize=3)
       bar_plot_here(ax1, tt, prcp_gpcp_monthly, 'r')
       ax1.plot(tt, np.nanmean(prcp_precl_monthly, axis=0), 'mo-', label='PREC/L', markersize=3)
       bar_plot_here(ax1, tt, prcp_precl_monthly, 'm')
       ax1.plot(tt, np.nanmean(prcp_gpcc_monthly, axis=0), 'yo-', label='GPCC', markersize=3)
       bar_plot_here(ax1, tt, prcp_gpcc_monthly, 'y') 
       ax1.plot(tt, np.nanmean(trmm_monthly, axis=0), 'go-', label='TRMM', markersize=3)
       bar_plot_here(ax1, tt, trmm_monthly, 'g')
       ax1.plot(tt, np.nanmean(cmap_monthly, axis=0), 'co-', label='CMAP', markersize=3)
       bar_plot_here(ax1, tt, cmap_monthly, 'c')
       ax1.grid('on')
       plt.ylabel('(mm/day)', fontsize=10)
       plt.xticks(tt, month_txt, fontsize=9)
       plt.legend(loc='upper right', ncol=3, fontsize=7)
       plt.title('(b) Precipitation Climatoloty', fontsize=font_size)

# Figure 1(c)
       ax1 = plt.axes([0.58, 0.42, 0.35, 0.2])
       ax1.plot(tt, np.nanmean(river3_monthly, axis=0), 'co-', label='Dai-Trenberth', markersize=3)
       ax1.plot(tt, np.nanmean(river_monthly, axis=0), 'bo-', label='HYBAM', markersize=3)
       bar_plot_here(ax1, tt, river_monthly, 'b')
       bar_plot_here(ax1, tt, river3_monthly, 'c')
       ax1.grid('on')
       plt.ylabel('(m$^3$/s x10$^5$)', fontsize=10)
       plt.xticks(tt, month_txt, fontsize=9)
       plt.legend(loc='upper right', ncol=1, fontsize=7)
       plt.title('(c) Amazon River Discharge Climatoloty', fontsize=font_size)       

# Figure 1(d)
       ax1 = plt.axes([0.58, 0.14, 0.35, 0.2])
#       ax1.plot(tt, np.nanmean(sss_aquarius_monthly, axis=0), 'ko-', label='Aquarius', markersize=3)
#       bar_plot_here(ax1, tt, sss_aquarius_monthly, 'k')
#       ax1.plot(tt, np.nanmean(sss_smos_monthly, axis=0), 'bo-', label='SMOS', markersize=3)
#       bar_plot_here(ax1, tt, sss_smos_monthly, 'b')
       ax1.plot(tt, np.nanmean(sss_ecco_monthly, axis=0), 'go-', label='ECCO', markersize=3)
       bar_plot_here(ax1, tt, sss_ecco_monthly, 'g')
       ax1.plot(tt, np.nanmean(sss_monthly, axis=0), 'bo-', label='GECCO2', markersize=3)
       bar_plot_here(ax1, tt, sss_monthly, 'b')
       ax1.plot(tt, np.nanmean(sss_oras5_monthly, axis=0), 'yo-', label='ORAS5', markersize=3)
       bar_plot_here(ax1, tt, sss_oras5_monthly, 'y')
       ax1.plot(tt, np.nanmean(sss_soda_monthly, axis=0), 'mo-', label='SODA3.1.1', markersize=3)
       bar_plot_here(ax1, tt, sss_soda_monthly, 'm')
       ax1.plot(tt, np.nanmean(sss_en4_monthly, axis=0), 'co-', label='EN4', markersize=3)
       bar_plot_here(ax1, tt, sss_en4_monthly, 'c')
       ax1.grid('on')
       plt.ylabel('(PSU)', fontsize=10)
       plt.xticks(tt, month_txt, fontsize=9)
       plt.legend(loc='lower right', ncol=1, fontsize=6)
       plt.title('(d) 5-meter Ocean Salinity Climatoloty', fontsize=font_size)
       plt.savefig('geographic_tmp_plot.jpg', format='jpeg', dpi=200)

       fig = plt.figure()
       ax1 = plt.axes([0.3, 0.3, 0.5, 0.5]) 
       ax1.plot(tt, np.nanmean(sss_aquarius_monthly, axis=0), 'ko-', label='Aquarius', markersize=3)
       bar_plot_here(ax1, tt, sss_aquarius_monthly, 'k')
       ax1.plot(tt, np.nanmean(sss_smos_monthly, axis=0), 'bo-', label='SMOS', markersize=3)
       bar_plot_here(ax1, tt, sss_smos_monthly, 'b')
#       ax1.plot(tt, np.nanmean(sss_ecco_monthly, axis=0), 'co-', label='ECCO', markersize=3)
#       bar_plot_here(ax1, tt, sss_ecco_monthly, 'c')
#       ax1.plot(tt, np.nanmean(sss_monthly, axis=0), 'go-', label='GECCO2', markersize=3)
#       bar_plot_here(ax1, tt, sss_monthly, 'g')
       ax1.grid('on')
       plt.ylabel('(PSU)', fontsize=10)
       plt.xticks(tt, month_txt, fontsize=9)
       plt.legend(loc='lower right', ncol=1, fontsize=6)
       plt.title('5-meter Ocean Salinity Climatoloty', fontsize=font_size)
       plt.savefig('suppl_geographic_tmp_plot.jpg', format='jpeg', dpi=200)

       plt.show()

#    else:

#       print('hahahah')

# ================================================================
# Describe this script and control flags
# ================================================================


