''' This script performs the actual dv/v computation.'''
#%%

import sys,obspy,os,glob,time
import numpy as np
import pandas as pd
from seisgo.utils import split_datetimestr
from seisgo import downloaders,noise,utils,monitoring
from obspy import clients
import pickle 
import pyasdf
import matplotlib.pyplot as plt

#%% parameters 
chan_list = ["HHZ"]
net_list  = ["YG"] #                                             # network list
sta_list  = ["SW2","SE2R"]      
cc_comp = 'ZZ'       

rootpath = "/data/stor/proj/Turner/NOISE/data_SK_2025/"     
pickle_path = rootpath+'dvv_pickles_final/'

save = True                                # Save to pickle?
cc_step = False                             # option to turn off cc_step if unnecessary 
shift_corrdata = True                       # shift array based on pahse arrival 

# lwin_test = 60*30 #  np.arange(15, 60, 5)
#var = np.zeros(len(lwin_test))

#for i in range(len(lwin_test)):
#    print(i)

# ''' Here are the parmeters for the dvv function: 
    
#     freq: minimum and maximum frequencies for dvv measurements.
freq = [0.7,1.5]                              #targeted frequency band for monitoring

#     win: window length for dvv measurements.
lwin = 20 # 5/freq[0]                              #window length in seconds for the coda waves

#     ref: reference trace/stack. default is None (will get by stacking all in the data)
#     dvmax: maximum dvv searching range. default: 0.05 (5%).\
dv_max = 10/100                              #limit for the dv/v

#     vmin: minimum velocity for the main phase, only considered in cross-correlations between
#         two stations. default 1.0 km/s.
vmin = 1.0                              #minimum velocity of the direct waves -> 

#     offset: offset from 0.0 seconds (for autocorr) and from the maximum arrival time of the main
#         phase (for xcorr)
offset=6                                  # seconds from zero

#     resolution: in seconds, specifying the temporal resolution (resampling/substacking) before measuring
#         dvv.
win_len = 3600*24*3                          #Window length i.e. Resolution

#     stack_method: stacking method to get the reference trace (if not given) and the short-window substacks. default is 'linear'.
#     normalize: Ture or False for data normalization in measuring dvv.
#     whiten='no',whiten_smooth=20, whiten_pad=100: parameters for whitening the trace before measuring dv/v.
#                 whiten: default is 'no', could be 'phase_only' or 'rma'.
whiten = 'no'

#     method: dvv measuring method.
method = 'ts'
#     subfreq: keep all frequencies in the dvv result. default is True. Otherwise, only get one dvv result.'''


#% derived params - might need modification!
#sta_list.sort()
source= net_list[0]+'.'+ sta_list[0]                        #source station
receiver = net_list[0]+'.'+ sta_list[1]                      #reciever station
                                                #Which components to use (ZZ, EZ, EN, EE, NZ, NN)

#% dirs 
CCFDIR    = os.path.join(rootpath,'cross_correlations',sta_list[0]+'_'+sta_list[1],chan_list[0])   # dir to store CORRELATION data
MERGEDIR  =  os.path.join(rootpath,'merged_stations',sta_list[0]+'_'+sta_list[1],chan_list[0])  # dor to store merged stations data   
print(CCFDIR)
print(MERGEDIR)

if cc_step:
    #% merging the corrfiles
    if not os.path.isdir(MERGEDIR):os.makedirs(MERGEDIR)
    # cross-correlation files
    ccfiles   = sorted(glob.glob(os.path.join(CCFDIR,'*.h5')))
    pairs_all,netsta_all=noise.get_stationpairs(ccfiles,False)

    #
    for s in netsta_all:
        tmp = os.path.join(MERGEDIR,s)
        if not os.path.isdir(tmp):os.mkdir(tmp)

    # MPI loop: loop through each user-defined time chunck
    for ipair in range(len(pairs_all)):
        pair=pairs_all[ipair]
        print('station-pair %s'%(pair))
        try:
            noise.merge_pairs(ccfiles,pairlist=pair,outdir=MERGEDIR,verbose=True,split=False)
            flag=False
        except:
            print('Header was too long. Processing a subset of corrfiles. Will need to manually adjust and do part 2!')
            noise.merge_pairs(ccfiles[0:100],pairlist=pair,outdir=MERGEDIR,verbose=True,split=False)
            flag=True


# % file handling#
key=source+'_'+receiver
ccfile=sorted(glob.glob(os.path.join(MERGEDIR,source,'*'+receiver+'*.h5')))[0]

###### Loading Waveform Data #######
corrdata=noise.extract_corrdata(ccfile,pair=source+'_'+receiver,comp=cc_comp)
cdatashifted=corrdata[key][cc_comp].copy()
cdata=corrdata[key][cc_comp].copy()


#%% corrdata correction 
if shift_corrdata:
    pad_with = 100 # padd around so as to allow row shift
    twin = [-6,6] # window for phase arrival, to compute shift in corr - should be different than coda!
    # make a padded array of zeros 
    cdata_padded = np.pad(cdata.data,pad_width=((0,0),(pad_with,pad_with)),mode='constant',constant_values=0)

    # find the time window to run corss correlation over 
    tvec = np.arange(-cdata.lag,cdata.lag+0.5*cdata.dt,cdata.dt)
    
    # Find indices where tvec is within twin
    indices = np.where((tvec >= twin[0]) & (tvec <= twin[1]))[0] + pad_with

    # define reference function megastack
    ref = np.median(cdata_padded[:,indices],axis=0)

    #% find lags at which correlation between each row and the ref function is highest
    from scipy.signal import correlate, correlation_lags
    lags = correlation_lags(len(ref),len(cdata_padded[0,indices]),mode='full')
    lag = np.zeros(cdata_padded.shape[0])
    for i in range(cdata_padded.shape[0]):
        corr = correlate(ref,cdata_padded[i,indices])
        lag[i] = lags[np.argmax(corr)]

    #% shift each row according to the respective lag
    shifted_data = np.zeros_like(cdata_padded)
    for j in range(cdata_padded.shape[0]):
        shift = int(lag[j])
        row = cdata_padded[j,:]
        if shift == 0:
            shifted_data[j,:] = row
            continue
        if shift>0: # shift to the right
            shifted_data[j,shift:] = row[:-shift]
        else:
            shifted_data[j, :shift] = row[-shift:]

    # implement into cdata object
    cdatashifted.data = shifted_data[:,100:-100]

    #% plot results 
    plt.figure()
    plt.imshow(cdata_padded,vmin=-.01,vmax=.01)
    plt.axvline(indices[0],color='r')
    plt.axvline(indices[-1],color='r')

    plt.figure()
    plt.imshow(shifted_data,vmin=-.01,vmax=.01)
    plt.axvline(indices[0],color='r')
    plt.axvline(indices[-1],color='r')

#%%
t1=time.time()
nproc=8 #number of processors to use.
vdata=monitoring.get_dvv(cdatashifted,freq,lwin,resolution=win_len,vmin=vmin,dvmax=dv_max,offset=offset,\
                        method=method,plot=True,nproc=nproc)
print("with %d processors, took %f seconds"%(nproc,time.time()-t1))

#%
vdata.plot(ylim=(-2,2))

#%%
if save:
    if not os.path.isdir(pickle_path):os.makedirs(pickle_path)
    if shift_corrdata:
        with open(pickle_path+'dvv_output_'+sta_list[0]+'_'+sta_list[1]+'_s_'+chan_list[0]+'.pkl', 'wb') as f:
            pickle.dump(vdata, f)
    else:
        with open(pickle_path+'dvv_output_'+sta_list[0]+'_'+sta_list[1]+'_'+chan_list[0]+'.pkl', 'wb') as f:
            pickle.dump(vdata, f)


# %%

# corrdata['YG.SW7_YG.SE9']['EE'].plot_for_paper(lag=35,freqmin=0.01,freqmax=5,coda_times=vdata.window,vlim=0.6)
# cdatashifted.plot_for_paper(lag=35,freqmin=0.01,freqmax=5,coda_times=vdata.window,vlim=0.6)



# %%
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from seisgo import stacking

# find the time window to run corss correlation over 
tvec = np.arange(-cdata.lag,cdata.lag+0.5*cdata.dt,cdata.dt)
times = cdatashifted.time
timestamp =[]
tmarks =[]
for t in times:
    #timestamp.append = obspy.UTCDateTime(t)
    tmarks.append(obspy.UTCDateTime(t).strftime('%Y-%m-%d'))

# find where there are holes in the times.
date_series = pd.Series(pd.to_datetime(tmarks))
gaps = date_series.diff().dt.days
np.mean(gaps)
gap_indices = []
for i in range(1, len(gaps)):
    if gaps[i] > 5:
        gap_indices.append(i)



xlim = 35 # seconds
vlim = 0.05 
coda_times = vdata.window

# Create the figure and GridSpec
fig = plt.figure(figsize=(12, 10))  # Adjust figure size as needed
gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[6,1], width_ratios=[4,1]) # 2 rows, 2 cols

# Subplots sharing x-axis (left column)
ax0 = fig.add_subplot(gs[0, 0])  # Top left
ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)  # Bottom left, sharing x with ax0

# Subplots sharing y-axis (top row)
#ax2 = fig.add_subplot(gs[0, 1], sharey=ax0)  # Top right, sharing y with ax0

import matplotlib.patheffects as pe

mesh = ax0.pcolormesh(tvec,tmarks,cdatashifted.data,cmap='RdBu_r',vmin = -vlim, vmax=vlim)
ax0.yaxis.set_major_locator(ticker.IndexLocator(150, offset=0))  # Show every 50th index

ax0.axvline(0,color='k',linestyle=':')
ax0.axvspan(coda_times[0], coda_times[1], color='g', alpha=0.1) 
ax0.axvspan(-coda_times[1], -coda_times[0], color='g', alpha=0.1,label='dv/v windows')
ax0.set_xlim(-xlim,xlim)
for gi in gap_indices:
    #ax0.axhline(tmarks[gi],xmin=-xlim-2,xmax=xlim+2,color='white',linewidth=2,zorder=3,clip_on=False)
    ax0.plot(
    [ax0.get_xlim()[0] - 2, ax0.get_xlim()[1] + 2], [tmarks[gi], tmarks[gi]], color='white', linewidth=2,
    zorder=3,clip_on=False)
    ax0.scatter(xlim,tmarks[gi-1],marker='_',s=200,color='black',zorder=4,clip_on=False)
    ax0.scatter(xlim,tmarks[gi+10],marker='_',s=200,color='black',zorder=4,clip_on=False)
    ax0.scatter(-xlim,tmarks[gi-1],marker='_',s=200,color='black',zorder=4,clip_on=False)
    ax0.scatter(-xlim,tmarks[gi+10],marker='_',s=200,color='black',zorder=4,clip_on=False)
    # ax0.axhline(tmarks[gi],xmin=-xlim-2,xmax=xlim+5,color='white',linewidth=3,zorder=6)
    ax0.text(xlim-10,tmarks[gi-30],str(tmarks[gi-1]),fontsize=10)
    ax0.text(xlim-10,tmarks[gi+60],str(tmarks[gi+1]),fontsize=10)
ax0.invert_yaxis()
dstack = stacking.seisstack(cdatashifted.data,method='linear',par=None)
ax1.plot(tvec,dstack,color='k',linewidth=1)
ax1.axvline(0,color='k',linestyle=':')
ax1.axvspan(coda_times[0], coda_times[1], color='g', alpha=0.1) 
ax1.axvspan(-coda_times[1], -coda_times[0], color='g', alpha=0.1,label='dv/v windows')
ax1.set_xlim(-xlim,xlim)
ax1.set_xlabel('time (s)')
# Remove ticks
#ax2.tick_params(axis='y', which='both', left=False, right=False, labelleft=False, labelright=False)
ax0.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False, labeltop=False)
ax1.tick_params(axis='y', which='both', left=False, right=False, labelleft=False, labelright=False)

ax0.set_title(sta_list[0]+' - '+sta_list[1] +' '+ chan_list[0] + ' cross correlation functions')
# cbar = plt.colorbar(mesh, ax=ax0)
# cbar.set_label("Colorbar Label")

plt.subplots_adjust(hspace=0.01, wspace=0.1)
#fig.tight_layout()  # Adjust these values as needed
plt.show()

#%%

plt.figure()
plt.scatter(vdata.time,vdata.maxcc2)
# %%
plt.figure(figsize=(50,10))
plt.plot(tmarks,lag)
ax = plt.gca()
ax.xaxis.set_major_locator(ticker.IndexLocator(20, offset=0))
# gap_indices = []
# for i in range(1, len(gaps)):
#     if gaps[i] > expected_delta_days:
#         gap_indices.append(i)


# %%
