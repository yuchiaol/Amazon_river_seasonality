run_flag = 1

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import ast, os, sys
from pylab import setp
from netCDF4 import Dataset
from scipy.io import loadmat
from scipy.interpolate import griddata
from scipy import stats
from numpy import unravel_index

sys.path.append('/home/yliang/lib/python_tools/python_functions/data_process/')
import data_process_f
sys.path.append('/home/yliang/lib/python_tools/python_functions/statistics')
import statistical_f

def read_here1(dirname, filename):

    print(filename)
    f = Dataset(dirname + filename, mode='r')
    ts_max = f.variables['ts_max'][:].data
    ts_min = f.variables['ts_min'][:].data
    f.close()

    return ts_max, ts_min

def max_min_select(ts_monthly,year_N):

    print('max min select')
    ts_max = np.zeros((year_N))*np.nan
    ts_min = np.zeros((year_N))*np.nan
    for II in range(year_N-2):
        ts_max[II+1] = ts_monthly[II+1,:].max()
        ts_min[II+1] = ts_monthly[II+1,:].min()

    return ts_max, ts_min

def read_river_mask():

    print('read river mask')
    river_mask_tmp = np.loadtxt('/home/yliang/whoi_projects/Amazon_river/mask_new/obidosglobalnew.txt')
    river_mask_tmp = np.flipud(river_mask_tmp)

    river_mask = np.zeros((river_mask_tmp.shape))
    river_mask[:,0:361] = river_mask_tmp[:,360:].copy()
    river_mask[:,361:] = river_mask_tmp[:,:360].copy()
    river_mask[river_mask<=0.] = np.nan

    return river_mask

def sss_calculate_here(TLONG,TLAT,sss_test):

# Regrid
    lon = np.linspace(0,360, 360)
    lat = np.linspace(-79,89,180)
    [gridx_out, gridy_out] = np.meshgrid(lon, lat)
    ny = 180
    nx = 360
    [gridx_in, gridy_in] = [TLONG, TLAT]
    nt = sss_test.shape[0]

    sss_test_interp = np.zeros((nt,ny,nx))
    for NT in range(nt):
        print(NT)
        sss_test_interp[NT,:,:] = griddata((gridy_in.ravel(), gridx_in.ravel()), np.squeeze(sss_test[NT,:,:]).ravel(), (gridy_out, gridx_out), method='nearest')

# Select region
    lon_amz1 = 295
    lon_amz2 = 311
    lat_amz1 = -3
    lat_amz2 = 15

    mask_sss = (sss_test_interp[11,:,:]/sss_test_interp[11,:,:]).copy()

    [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(lon_amz1,lon_amz2,lat_amz1,lat_amz2,lon,lat)
    area_so = data_process_f.area_calculate_nonuniform(lon, lat)
    ts_sss_test = np.zeros((nt))
    for NT in range(nt):
        ts_sss_test[NT] = np.nansum(sss_test_interp[NT,y1:y2+1,x1:x2+1]*area_so[y1:y2+1,x1:x2+1])/np.nansum(area_so[y1:y2+1,x1:x2+1]*mask_sss[y1:y2+1,x1:x2+1])

    ts_rmean = ts_sss_test.copy()
    ts_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_sss_test.copy(), np.ones((N,))/N, mode='valid')
    sss_monthly = ts_rmean.reshape((int(nt/12),12))
    [sss_max_test, sss_min_test] = max_min_select(sss_monthly,int(nt/12))

    [sss_test_trend, sss_test_sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(np.linspace(1,nt/12,nt/12)[0:-2], sss_max_test[1:-1]-sss_min_test[1:-1], sig_level)

    return sss_test_trend, sss_test_sig, pvalue 

def rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,year_ref):

    ts_rain_test = np.zeros((year_N*12))
    ts_runoff_test = np.zeros((year_N*12))
    NN = 0
    for NY in range(year_N):
        for NM in range(12):
            filename = 'P_R_' + str(year_ref+NY+31) + '-' + str(NM+1).zfill(2) + '.nc'
            print(filename)
            f = Dataset(dirname + filename, mode='r')
            P_in = f.variables['RAIN'][:,:,:].data.squeeze() + f.variables['SNOW'][:,:,:].data.squeeze()
            P_in[P_in>1.e11] = np.nan
            R_in = f.variables['QCHANR'][:,178,256].data.squeeze()
#            R_in[R_in>1.e11] = np.nan
            f.close()
            ts_rain_test[NN] = np.nansum(P_in[:,:]*mask_rg*area_in)/np.nansum(mask_rg*area_in)
#            ts_runoff_test[NN] = np.nansum(R_in[:,:]*mask_rg*area_in)#/np.nansum(mask_rg*area_in)
            ts_runoff_test[NN] = R_in
            NN = NN + 1

    ts_rmean = ts_rain_test.copy()
    ts_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_rain_test.copy(), np.ones((N,))/N, mode='valid')
    rain_monthly = ts_rmean.reshape((year_N,12))
    [rain_max_test, rain_min_test] = max_min_select(rain_monthly,year_N)

    [rain_test_trend, rain_test_sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(np.linspace(1,year_N,year_N)[0:-2], rain_max_test[1:-1]-rain_min_test[1:-1], sig_level)

    ts_rmean = ts_runoff_test.copy()
    ts_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_runoff_test.copy(), np.ones((N,))/N, mode='valid')
    runoff_monthly = ts_rmean.reshape((year_N,12))
    [runoff_max_test, runoff_min_test] = max_min_select(runoff_monthly,year_N)

    [runoff_test_trend, runoff_test_sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(np.linspace(1,year_N,year_N)[0:-2], runoff_max_test[1:-1]-runoff_min_test[1:-1], sig_level)

    return rain_test_trend, rain_test_sig, runoff_test_trend, runoff_test_sig, rain_max_test-rain_min_test

if run_flag == 1:
# ================================================================
# Basic setting
# ================================================================
   sig_level = 0.05

   dirname = '/home/yliang/whoi_projects/Amazon_river/model/CORE-II/'
   argo_case_txt = ['argo_ciso_tmp.nc','argo_iprc_tmp.nc','argo_isas_tmp.nc','argo_jamstec_tmp.nc','argo_scripps_tmp.nc']
   cesm2_case_txt = ['cesm2_control_tmp.nc']
#   core2_case_txt = ['core2_awi_tmp.nc','core2_bergen_tmp.nc','core2_cerfacs_tmp.nc','core2_cmcc_tmp.nc','core2_cnrm_tmp.nc','core2_fsu_tmp.nc','core2_gfdlg_tmp.nc','core2_gfdlm_tmp.nc','core2_ictp_tmp.nc','core2_inmom_tmp.nc','core2_mria_tmp.nc','core2_mrif_tmp.nc','core2_ncar_tmp.nc','core2_nocs_tmp.nc']
   core2_case_txt = ['core2_bergen_tmp.nc','core2_cerfacs_tmp.nc','core2_cmcc_tmp.nc','core2_cnrm_tmp.nc','core2_fsu_tmp.nc','core2_gfdlg_tmp.nc','core2_gfdlm_tmp.nc','core2_ictp_tmp.nc','core2_inmom_tmp.nc','core2_mria_tmp.nc','core2_mrif_tmp.nc','core2_ncar_tmp.nc']


   sss_argo_trends = np.zeros((len(argo_case_txt)+1))
   sss_argo_sig = np.zeros((len(argo_case_txt)+1))
   sss_cesm2_trends = np.zeros((len(cesm2_case_txt)+1))
   sss_core2_trends = np.zeros((len(core2_case_txt)+1))
   sss_core2_sig = np.zeros((len(core2_case_txt)+1))

   ts_sss_core2 = np.zeros((len(core2_case_txt)+1,60))

   run_trend_core2 = np.zeros((len(core2_case_txt), 58)) 
   run_sig_core2 = np.zeros((len(core2_case_txt), 58))

# ================================================================
# Read Argo data
# ================================================================
   for II in range(len(argo_case_txt)):
       [ts_max, ts_min] = read_here1(dirname, argo_case_txt[II])
       nt = len(ts_max)
#       sss_argo_trends[II] = np.polyfit(np.linspace(1,nt,nt)[1:-1],ts_max[1:-1]-ts_min[1:-1],1)[0]
#       sss_argo_sig[II] = stats.ttest_ind(np.linspace(1,nt,nt)[1:-1]*0.,ts_max[1:-1]-ts_min[1:-1])[1]
       [sss_argo_trends[II], sss_argo_sig[II], pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(np.linspace(1,nt,nt)[1:-1], ts_max[1:-1]-ts_min[1:-1], sig_level)
       if II == 0: ts_argo_ciso = ts_max - ts_min
       if II == 1: ts_argo_iprc = ts_max - ts_min
       if II == 2: ts_argo_isas = ts_max - ts_min
       if II == 3: ts_argo_jamstec = ts_max - ts_min
       if II == 4: ts_argo_scripps = ts_max - ts_min
   sss_argo_trends[-1] = sss_argo_trends[:-1].mean()
   sss_argo_trends[:] =sss_argo_trends[:]*10

   for II in range(len(core2_case_txt)):
       [ts_max, ts_min] = read_here1(dirname, core2_case_txt[II])
       nt = len(ts_max[-31:])
#       nt = len(ts_max[:])

       [sss_core2_trends[II], sss_core2_sig[II], pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(np.linspace(1,nt,nt)[:-1], ts_max[-31:-1], sig_level)
       print('sss core2 max', str(sss_core2_trends[II]*10), sss_core2_sig[II])

       [sss_core2_trends[II], sss_core2_sig[II], pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(np.linspace(1,nt,nt)[:-1], ts_min[-31:-1], sig_level)
       print('sss core2 min', str(sss_core2_trends[II]*10), sss_core2_sig[II])

       [sss_core2_trends[II], sss_core2_sig[II], pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(np.linspace(1,nt,nt)[:-1], ts_max[-31:-1]-ts_min[-31:-1], sig_level)
       print('sss core2 max-min', str(sss_core2_trends[II]*10), sss_core2_sig[II])

       ts_sss_core2[II,:] = ts_max - ts_min

#       [y_slope, y_sig] = statistical_f.running_linear_trends(np.linspace(1,nt,nt)[1:-1], ts_max[1:-1]-ts_min[1:-1],21,sig_level)
#       run_trend_core2[II,:] = y_slope.copy()
#       run_sig_core2[II,:] = y_sig.copy()
                
   sss_core2_trends[-1] = np.nanmean(sss_core2_trends[:-1])
   sss_core2_trends[:] = sss_core2_trends[:]*10
   ts_sss_core2[-1,:] = np.nanmean(ts_sss_core2[0:-1,:], axis=0)

# ================================================================
# Plot figures
# ================================================================
   plt.close('all')

   fig = plt.figure()

   fig.set_size_inches(10, 10, forward=True)

   plt.axes([0.2, 0.6, 0.45, 0.3])
   tt = np.linspace(1,len(core2_case_txt)+1,len(core2_case_txt)+1)
   plt.bar(tt, sss_core2_trends, width=0.5, color='m', alpha=0.5)
   plt.plot([tt[-1], tt[-1]], [sss_core2_trends[-1], sss_core2_trends[-1]+np.nanstd(sss_core2_trends[:-1])], 'k-')
   plt.plot([tt[-1], tt[-1]], [sss_core2_trends[-1], sss_core2_trends[-1]-np.nanstd(sss_core2_trends[:-1])], 'k-')
#   for II in range(len(sss_core2_sig)):
#       if sss_core2_sig[II] == 1:
#          plt.plot([tt[II], tt[II]], [sss_core2_trends[II], sss_core2_trends[II]], 'k*')

   plt.ylabel('(PSU per decade)')
   plt.xticks(tt[:], ['BERGEN','CERFAC','CMCC','CNRM','FSU','GFDL-GOLD','GFDL-MOM','ICTP','INMOM','MRI-A','MRI-F','NCAR','Mean'], rotation=60)
   plt.title('(a) CORE-II Models')

   plt.axes([0.2, 0.13, 0.45, 0.3])
   tt = np.linspace(1,len(argo_case_txt)+1,len(argo_case_txt)+1)
   plt.bar(tt, sss_argo_trends, width=0.5, color='m', alpha=0.5)
   plt.plot([tt[-1], tt[-1]], [sss_argo_trends[-1], sss_argo_trends[-1]+np.nanstd(sss_argo_trends[:-1])], 'k-')
   plt.plot([tt[-1], tt[-1]], [sss_argo_trends[-1], sss_argo_trends[-1]-np.nanstd(sss_argo_trends[:-1])], 'k-')
#   for II in range(len(sss_argo_sig)):
#       if sss_argo_sig[II] == 1:
#          plt.plot([tt[II], tt[II]], [sss_argo_trends[II], sss_argo_trends[II]], 'k*')

   plt.ylabel('(PSU per decade)')
   plt.xticks(tt[:], ['CISO','IPRC','ISAS','JAMSTEC','SCRIPPS','Mean'], rotation=60)
   plt.title('(b) Interpolated ARGO')


   plt.savefig('core2_tmp_plot.jpg', format='jpeg', dpi=200)

   plt.show()

