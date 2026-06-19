
"""this file contains the various helper functions used to make the main plots, including snow correction and three component averaging of the signals"""

import numpy as np 
import pickle 
from datetime import datetime 
import pandas as pd
import matplotlib.dates as mdates 


def load_dvv(direc_pickle,sp,channel,window_size,cc_threshold, get_neg_values=False):
    
    file = direc_pickle + 'dvv_output_' + sp + '_' + channel + '.pkl'
    with open(file, 'rb') as f:
        vdata = pickle.load(f)


    # correct HHZ direction error for SE15
    if sp == 'SW14_SE15_s' and channel =='HHZ':
        get_neg_values=True

    if get_neg_values:
        data = vdata.data1
        print('getting negative side values')
    else:
        data = vdata.data2
    cc = vdata.maxcc2
    date_range = [datetime.utcfromtimestamp(ts) for ts in vdata.time]

    # Step 1: Filter by correlation coefficient
    filtered_data = [d for d, c in zip(data, cc) if c > cc_threshold and d != 0]
    filtered_date_range = [date for date, d, c in zip(date_range, data, cc) if c > cc_threshold and d != 0]
            
    moving_mean = np.convolve(filtered_data, np.ones(window_size) / window_size, mode='valid')
    moving_mean_masked = np.ma.masked_invalid(moving_mean)
    date_range_moving_mean = filtered_date_range[window_size - 1:]

    return filtered_data,filtered_date_range,moving_mean_masked,date_range_moving_mean,data,date_range,cc



def load_velocity(vel_read_file):
    #vel_read_file = velocity_file_km15
    velocity_data = pd.read_csv(vel_read_file, parse_dates=['mid_date'], dayfirst=True)
    # Assuming velocity_data is a pandas DataFrame with columns 't' (time) and 's' (speed)
    velocity_data['t'] = pd.to_datetime(velocity_data['mid_date'])  # Ensure 't' is a datetime column

    # Sort the data by time
    velocity_data.sort_values(by='t', inplace=True)
    # take mean if there are multiple vlaues on same mid_date
    velocity_data = velocity_data.groupby('t', as_index=False)['s'].mean()
    # Set time as the index
    velocity_data.set_index('t', inplace=True)
    # 10-day moving mean
    velocity_data['moving_mean'] = velocity_data['s'].rolling('10D').mean()
    # Reset index if needed
    velocity_data.reset_index(inplace=True)

    if not vel_read_file[-8:-6] == '15':
        velocity_data['moving_mean'] = velocity_data['moving_mean']/365

    return velocity_data

def seasonal_dvv_correction(dates,data,snowmasskm=7,output_snow=False):
 
    # slope and intercept from snowmass analysis based on east side stations 
    slope = 0.32456273096865684 # east stations only
    intercept = -0.2847229798205404 # east

    if snowmasskm < 5:
        slope =  0.1910130789135306 # east side stations only 
        intercept = -0.23607644341872175 # east side stations only 

    # load the appropriate snow 
    if snowmasskm > 10:
        snowmass = pd.read_csv('./timeseries_data/snowmass_2025_timeseries_km15.csv',parse_dates=['time'], dayfirst=True)
    elif snowmasskm < 5: 
        snowmass = pd.read_csv('./timeseries_data/snowmass_2025_timeseries_km5.csv',parse_dates=['time'], dayfirst=True)
    else:
        snowmass = pd.read_csv('./timeseries_data/snowmass_2025_timeseries_km7.csv',parse_dates=['time'], dayfirst=True)
    snow_mwe = snowmass['runoff']
    snow_time = snowmass['time']
    snow_time_num = mdates.date2num(snow_time)
    # synthetic dvv based on linear fit
    dvv_synth = snow_mwe * slope + intercept
    raw_dates_num = mdates.date2num(dates) # adjust dates so we can interpolate
    dvv_c = np.interp(raw_dates_num, snow_time_num, dvv_synth) # interpolate
    raw_data = data - dvv_c # subtract synthetic portion from signal to hopefully be left with the rest of the signal
    if output_snow:
        return raw_data, dvv_c, snow_time,snow_mwe
    else:
        return raw_data,dvv_c
    



# function to average the three component signal
def average_components(sp, comp_list, direc_pickle,window_size=1,cc_threshold=0.5,get_neg=False):
    
  # Initialize lists to hold data
    all_m_data = []
    all_m_dates = []
    all_cc = []
    
    for channel in comp_list:
        data, dates, m_data, m_dates, raw_data, raw_dates, cc = load_dvv(direc_pickle, sp, channel, window_size, cc_threshold, get_neg_values=get_neg)

        m_data = np.array(raw_data)
        m_dates = np.array([d.toordinal() for d in raw_dates])  # Convert to ordinal dates
        cc = np.array(cc)

        # remove vals with very weird cc values
        m_data[cc<0.2] = np.nan
        m_data[np.abs(m_data)>1.5] = np.nan
        
        all_m_data.append(m_data)
        all_m_dates.append(m_dates)
        all_cc.append(cc)
    
    unique_dates = np.unique(np.concatenate(all_m_dates))
    
    # interplotat the values
    interpolated_data = []
    interpolated_cc = []
    for i in range(len(all_m_data)):
        interp_values = np.interp(unique_dates, all_m_dates[i], all_m_data[i], left=np.nan, right=np.nan)
        interp_cc = np.interp(unique_dates, all_m_dates[i], all_cc[i], left=np.nan, right=np.nan)
        interpolated_data.append(interp_values)
        interpolated_cc.append(interp_cc)
    
    threec_data = np.array(interpolated_data)
    threec_cc = np.array(interpolated_cc)
    
    #finally, average vertically
    averaged_data = np.nanmean(threec_data, axis=0)
    averaged_cc = np.nanmean(threec_cc, axis=0)
    
    # convert to datetime
    averaged_dates = [datetime.fromordinal(int(d)) for d in unique_dates]
    
    return averaged_data, averaged_dates, averaged_cc
