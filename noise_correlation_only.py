
''' This script does only the cross correlation aprt and the stacking. It pulls form an existing raw data file.

Can be run from the command line, with the station pair and the channel as overriding arguments. An example here:

 nohup python -u ./processing_noise_correlations.py SE14 SW15 HHZ > SE14_SW15_HHZ.log &


Yoram Terleth - JUNE 2024 

requires the seisgo python package. 
This script is heavily based on the dv/v example provided with the seisgo documentation.
'''


#%% import functions
import sys,obspy,os,glob,time

# make sure the parent directory is in the Python path
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

import numpy as np
import pandas as pd
from seisgo.utils import split_datetimestr
from seisgo import downloaders,noise_v2,utils,monitoring
from obspy import clients
# from seisgo_v2 import noise as noise_new

#%% define parameters 
# max_tries = 3                                                  #maximum number of tries when downloading, in case the server returns errors.
# flag      = False                                               # print progress when running the script; recommend to use it at the begining
# samp_freq = 5                                                  # targeted sampling rate at X samples per seconds
# rmresp   = True
# rmresp_out = 'DISP'
# pressure_chan = [None]  #Added by Xiaotao Yang. This is needed when downloading some special channels, e.g., pressure data. VEL output for these channels.
# freqmin   = 0.001                                              # pre filtering frequency bandwidth
# freqmax   =   0.5*samp_freq
# note this cannot exceed Nquist freq

# targeted region/station information: only needed when use_down_list is False
lamin,lamax,lomin,lomax= 59.9605,60.1538,-139.9209,-139.4938                # regional box: min lat, min lon, max lat, max lon (-114.0)
chan_list = ["HHZ"]
net_list  =  ["YG"] #                                             # network list
sta_list  = ["SW2","SE7"]                                       # station (using a station list is way either compared to specifying stations one by one)                               # end date of download
#inc_hours  = 24*30                                                 # length of data for each request (in hour)
maxseischan = 2                                                  # the maximum number of seismic channels, excluding pressure channels for OBS stations.
ncomp      = maxseischan #len(chan_list)


# paths and filenames
skip_if_H5_files_exist = False # SET TO TRUE ONLY IF THERE IS NO CHANGE IN PARAMETERS!

# dir to get the mseed data from 
basic_data_dir = '/data/stor/basic_data/seismic_data/day_vols/TURNER/'

# dir to get the response information from 
resp_file = '/data/stor/basic_data/seismic_data/day_vols/TURNER/resp/xml_files/Turner_station_20241018.xml'
#response_dir = '/data/stor/basic_data/seismic_data/day_vols/TURNER/resp/'

# save H5 format dir
datapath = "/data/stor/proj/Turner/NOISE/data_SK_2025/"
rootpath = "/data/stor/proj/Turner/NOISE/data_SK_2025/"                                 # roothpath for the project


# for xcorr
DATADIR  = os.path.join(datapath,'Raw',sta_list[0]+'_'+sta_list[1],chan_list[0])    # directory to save H5 files to                 
down_list  = os.path.join(DATADIR,'station.txt')                                    # CSV file for station location info 
CCFDIR    = os.path.join(rootpath,'cross_correlations',sta_list[0]+'_'+sta_list[1],chan_list[0])   # dir to store CORRELATION data
MERGEDIR  =  os.path.join(rootpath,'merged_stations',sta_list[0]+'_'+sta_list[1],chan_list[0])  # dor to store merged stations data   

# for autocorr
# DATADIR  = os.path.join(datapath,'Raw',sta_list[0]+'_'+sta_list[1],chan_list[0])    # directory to save H5 files to                 
# down_list  = os.path.join(DATADIR,'station.txt')  
# CCFDIR    = os.path.join(rootpath,'auto_correlations',sta_list[0]+'_'+sta_list[0],chan_list[0])   # dir to store CORRELATION data
# MERGEDIR  =  os.path.join(rootpath,'autoCorr_merged_stations',sta_list[0]+'_'+sta_list[0],chan_list[0])  # dor to store merged stations data   

print(DATADIR)
print(CCFDIR)
print(MERGEDIR)



#################################################################################################################
# NOW the cross correlation part of the script
# some control parameters
freq_norm   = 'rma'                                                         # 'no' for no whitening, or 'rma' for running-mean average, 'phase' for sign-bit normalization in freq domain
time_norm   = 'no'                                                          # 'no' for no normalization, or 'rma', 'one_bit' for normalization in time domain
cc_method   = 'xcorr'                                                       # 'xcorr' for pure cross correlation, 'deconv' for deconvolution; FOR "COHERENCY" PLEASE set freq_norm to "rma" and time_norm to "no"
acorr_only  = False                                                         # only perform auto-correlation
xcorr_only  = True                                                         # only perform cross-correlation or not

# pre-processing parameters
cc_len    = 3600*2                                                          # basic unit of data length for fft (sec)
step      = 3600*1                                                           # overlapping between each cc_len (sec)

# cross-correlation parameters
maxlag         = 50                                                        # lags of cross-correlation to save (sec)
substack       = True                                                      # sub-stack daily cross-correlation or not
substack_len   = cc_len*6                                            # how long to stack over (for monitoring purpose): need to be multiples of cc_len

freqmin   = 0.01
freqmax   = 5



##############################################################################################################
# make a dictionary to store all variables: also for later cc
fc_para={'cc_len':cc_len,'step':step,'freqmin':freqmin,'freqmax':freqmax,
        'freq_norm':freq_norm,'time_norm':time_norm,'cc_method':cc_method,
        'substack':substack,'substack_len':substack_len,'maxlag':maxlag}
# save fft metadata for future reference
fc_metadata  = os.path.join(CCFDIR,'fft_cc_data.txt')
if not os.path.isdir(CCFDIR):os.makedirs(CCFDIR)
# save metadata
fout = open(fc_metadata,'w')
fout.write(str(fc_para));fout.close()

# set variables to broadcast
tdir = sorted(glob.glob(os.path.join(DATADIR,'*.h5')))
#%%
#loop through all data files.
for i in range(len(tdir)):
    sfile=tdir[i]
    #call the correlation wrapper.
    noise_v2.do_correlation(sfile,cc_len,step,maxlag,cc_method=cc_method,
                         acorr_only=acorr_only,xcorr_only=xcorr_only,substack=substack,
                         substack_len=substack_len,freqmin=freqmin,freqmax=freqmax,
                         time_norm=time_norm,freq_norm=freq_norm,outdir=CCFDIR)
    


##############################################################################################################
#%% Additionally, we merge the station pirs into a new directory
'''
Stacking script of SeisGo to:
    1) load cross-correlation data for each station pair
    2) merge all time chuncks
    3) save outputs in ASDF;
'''



if not os.path.isdir(MERGEDIR):os.makedirs(MERGEDIR)
# cross-correlation files
ccfiles   = sorted(glob.glob(os.path.join(CCFDIR,'*.h5')))
pairs_all,netsta_all=noise_v2.get_stationpairs(ccfiles,False)

#
for s in netsta_all:
    tmp = os.path.join(MERGEDIR,s)
    if not os.path.isdir(tmp):os.mkdir(tmp)

# MPI loop: loop through each user-defined time chunck
for ipair in range(len(pairs_all)):
    pair=pairs_all[ipair]
    print('station-pair %s'%(pair))
    noise_v2.merge_pairs(ccfiles,outdir=MERGEDIR,verbose=True)




# %%
# Note: I think the error is in the corr function. 
