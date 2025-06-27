
''' This script takes local data and transforms it to the correct format, then performs the noise correlations and makes station pairs.
The results of this script can be investigated with the associated jupyter notebook.

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

#%% define parameters 
max_tries = 3                                                  #maximum number of tries when downloading, in case the server returns errors.
flag      = False                                               # print progress when running the script; recommend to use it at the begining
samp_freq = 10                                                  # targeted sampling rate at X samples per seconds
rmresp   = True
rmresp_out = 'DISP'
pressure_chan = [None]  #Added by Xiaotao Yang. This is needed when downloading some special channels, e.g., pressure data. VEL output for these channels.
freqmin   = 0.001                                              # pre filtering frequency bandwidth
freqmax   =   0.5*samp_freq
# note this cannot exceed Nquist freq

# targeted region/station information: only needed when use_down_list is False
lamin,lamax,lomin,lomax= 59.0,61.0,-141.9,-139.0                # regional box: min lat, min lon, max lat, max lon (-114.0)
chan_list = ["BHZ","BHZ"]
net_list  = ["AT","AK"] #                                             # network list
sta_list  = ["YKU2","PNL"]                                       # station (using a station list is way either compared to specifying stations one by one)
start_date = "2020_09_15_0_0_0"                                 # start date of download
end_date   = "2024_09_01_0_0_0"                                 # end date of download
inc_hours  = 24*10                                                 # length of data for each request (in hour)
maxseischan = 10                                                  # the maximum number of seismic channels, excluding pressure channels for OBS stations.
ncomp      = maxseischan #len(chan_list)

back_azimuth = None # set to None to actually look at the set channel!! 


# paths and filenames
skip_if_H5_files_exist = False # SET TO TRUE ONLY IF THERE IS NO CHANGE IN PARAMETERS!

# dir to get the mseed data from 
basic_data_dir = '/data/stor/basic_data/seismic_data/day_vols/TURNER/'

# dir to get the response information from 
# resp_file = '/data/stor/basic_data/seismic_data/day_vols/TURNER/resp/xml_files/Turner_station_20241018.xml'
#response_dir = '/data/stor/basic_data/seismic_data/day_vols/TURNER/resp/'

# save H5 format dir
rootpath = "/data/stor/proj/Turner/NOISE/data_SK_2025_IRIS/"                                                         # roothpath for the project
DATADIR  = os.path.join(rootpath,'Raw',sta_list[0]+'_'+sta_list[1],chan_list[0])    # directory to save H5 files to                 
down_list  = os.path.join(DATADIR,'station.txt')                                    # CSV file for station location info 

CCFDIR    = os.path.join(rootpath,'cross_correlations',sta_list[0]+'_'+sta_list[1],chan_list[0])   # dir to store CORRELATION data
MERGEDIR  =  os.path.join(rootpath,'merged_stations',sta_list[0]+'_'+sta_list[1],chan_list[0])  # dor to store merged stations data   

print(DATADIR)
print(CCFDIR)
print(MERGEDIR)

# define user credentials fro accessing embargoed data 
source = 'IRIS'
user = 'username'
password = 'password'
credentials = [user,password]

tt0=time.time()

##################################################
# we expect no parameters need to be changed below
# assemble parameters used for pre-processing waveforms in downloading
prepro_para = {'rmresp':rmresp,'rmresp_out':rmresp_out,'freqmin':freqmin,'freqmax':freqmax,\
                'samp_freq':samp_freq}

downlist_kwargs = {"source":source, 'net_list':net_list, "sta_list":sta_list, "chan_list":chan_list, \
                    "starttime":start_date, "endtime":end_date, "maxseischan":maxseischan, "lamin":lamin, "lamax":lamax, \
                    "lomin":lomin, "lomax":lomax, "pressure_chan":pressure_chan, "fname":down_list}



if not os.path.isdir(DATADIR):os.makedirs(DATADIR)
stalist=downloaders.get_sta_list(**downlist_kwargs) # saves station list to "down_list" file

print(sta_list)                                        # here, file name is "station.txt"
# save parameters for future reference
metadata = os.path.join(DATADIR,'download_info.txt')
fout = open(metadata,'w')
fout.write(str({**prepro_para,**downlist_kwargs,'inc_hours':inc_hours,'ncomp':ncomp}));fout.close()

all_chunk = split_datetimestr(start_date,end_date,inc_hours)
if len(all_chunk)<1:
    raise ValueError('Abort! no data chunk between %s and %s' % (start_date,end_date))

########################################################
#################DOWNLOAD SECTION#######################
########################################################
# loop through each time chunk
for ick in range(len(all_chunk)-1):
    s1= all_chunk[ick]
    s2=all_chunk[ick+1]
    
    print('time segment:'+s1+' to '+s2)
    
    download_kwargs = {"source":source,"rawdatadir": DATADIR, "starttime": s1, "endtime": s2, \
              "stationinfo": down_list,"credentials":credentials,"verbose":True, **prepro_para}

    # Download for ick
    downloaders.download(**download_kwargs)

tt1=time.time()
print('downloading step takes %6.2f s' %(tt1-tt0))
