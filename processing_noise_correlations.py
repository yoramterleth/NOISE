
''' This script takes local data and transforms it to the correct format.

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
lamin,lamax,lomin,lomax= 59.9605,60.1538,-139.9209,-139.4938                # regional box: min lat, min lon, max lat, max lon (-114.0)
chan_list = ["HHN"]
net_list  = ["YG"] #                                             # network list
sta_list  = ["SW2","SE2R"]                                       # station (using a station list is way either compared to specifying stations one by one)
start_date = "2020_09_15_0_0_0"                                 # start date of download
end_date   = "2024_09_01_0_0_0"                                 # end date of download
inc_hours  = 24*10                                                 # length of data for each request (in hour)
maxseischan = 10                                                  # the maximum number of seismic channels, excluding pressure channels for OBS stations.
ncomp      = maxseischan #len(chan_list)

back_azimuth = None # set to None to actually look at the set channel!! 

# back azimuth from google earth, from true north ! this is the angle normal to the transverse one we consider. 
# SW8_SE9 : 193 
# Sw14_SE15: 205
# SW2_SE4: 173 : 
# SW14 SE14: 222
# SW7_SE7: 203
# SW8_SE7: 236
# SW2 SE2R: 209

# if len(sys.argv) == 3:
#     sta_list = [str(sys.argv[1]),str(sys.argv[2])]
#     print(["Station list overidden by command line input. sta_list =" + str(sta_list)])
# elif len(sys.argv)==4:
#     sta_list = [str(sys.argv[1]),str(sys.argv[2])]
#     chan_list = [str(sys.argv[3])]
#     print(["Station list and channel overidden by command line input. sta_list =" + str(sta_list) + ' and channel =' + str(chan_list)])
# else:
#     print('No changes to input parameters.')

# print(['sys argv length=' + str(len(sys.argv))])
# print(["Station list and channel  =" + str(sta_list) + ' and channel =' + str(chan_list)])

# storage details 
# paths and filenames
skip_if_H5_files_exist = False # SET TO TRUE ONLY IF THERE IS NO CHANGE IN PARAMETERS!

# dir to get the mseed data from 
basic_data_dir = '/data/stor/basic_data/seismic_data/day_vols/TURNER/'

# dir to get the response information from 
resp_file = '/data/stor/basic_data/seismic_data/day_vols/TURNER/resp/xml_files/Turner_station_20241018.xml'
#response_dir = '/data/stor/basic_data/seismic_data/day_vols/TURNER/resp/'

# save H5 format dir
rootpath = "/data/stor/proj/Turner/NOISE/data_SK_2025/"                                                         # roothpath for the project
DATADIR  = os.path.join(rootpath,'Raw',sta_list[0]+'_'+sta_list[1],chan_list[0])    # directory to save H5 files to                 
down_list  = os.path.join(DATADIR,'station.txt')                                    # CSV file for station location info 

CCFDIR    = os.path.join(rootpath,'cross_correlations',sta_list[0]+'_'+sta_list[1],chan_list[0])   # dir to store CORRELATION data
MERGEDIR  =  os.path.join(rootpath,'merged_stations',sta_list[0]+'_'+sta_list[1],chan_list[0])  # dor to store merged stations data   

print(DATADIR)
print(CCFDIR)
print(MERGEDIR)


#%% actual downloading 
tt0=time.time()

if not os.path.isdir(DATADIR):os.makedirs(DATADIR)
#stalist=downloaders.get_sta_list(**downlist_kwargs) # saves station list to "down_list" file
                                          # here, file name is "station.txt"
# save parameters for future reference
metadata = os.path.join(DATADIR,'download_info.txt')
fout = open(metadata,'w')
#fout.write(str({**prepro_para,**downlist_kwargs,'inc_hours':inc_hours,'ncomp':ncomp}));fout.close()

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
    sdatetime = obspy.UTCDateTime(s1)
    edatetime = obspy.UTCDateTime(s2)
  
    fname_test = os.path.join(DATADIR,str(sdatetime).replace(':', '-') + 'T' + str(edatetime).replace(':', '-') + '.h5')
    print(fname_test)
    if os.path.isfile(fname_test) and skip_if_H5_files_exist:
        print(fname_test + ' already exists in directory. Skipping.')
        continue
    else:
        # this is the actual conversion step
        try:
                downloaders.download(s1, s2, stationinfo=None, network=net_list, station=sta_list,channel=chan_list,
                source='IRIS',rawdatadir=DATADIR,sacheader=False, getstainv=True, max_tries=3,
                savetofile=False,pressure_chan=None,samp_freq=samp_freq,freqmin=freqmin,freqmax=freqmax,
                rmresp=True, rmresp_out='DISP',respdir=resp_file,qc=True,event=None,verbose=flag,credentials=None,download_local=True,backazimuth=back_azimuth,mseed_datadir=basic_data_dir)
        except:
            print('Conversion step failed. Moving on.')
            continue
        

    tt1=time.time()
    print('downloading step takes %6.2f s' %(tt1-tt0))


    # Up to here, the traces are resampled to "samp_freq" response is removed, data is demeaned and detrended, 
    # and bandpassed between freq_min and freqmax (the ones given above).

#################################################################################################################
