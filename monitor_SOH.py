"""This is a short function to pull and plot the SOH data from a seismic instrument. Allows to verify that the periods of clock 
drift apparent in the cross correlation stacks are indeed associated with poor data quality.

Yoram Terleth, March 2024"""

import obspy 

soh = obspy.read('/data/stor/basic_data/seismic_data/RAW/Turner/sv2108/SW14/2021/01/SW14.YG..SOH.2021.030')
print(soh)
gst = soh.select(channel='GST')

print(gst[0].data)
gst.plot()
# %%
import obspy
import glob
import os

# Define the root directory containing the files
root_directory = '/data/stor/basic_data/seismic_data/RAW/Turner/sv2108/SW14/'

# Find all files in the root directory and subdirectories that contain "SOH" in their filename
soh_files = glob.glob(os.path.join(root_directory, '**', '*SOH*'), recursive=True)

# Initialize an empty Stream object to hold all GST data
gst_stream = obspy.Stream()

# Loop through each matching file and read in the GST channels
for file in soh_files:
    try:
        soh = obspy.read(file)
        gst = soh.select(channel='GPL')
        gst_stream += gst
    except Exception as e:
        print(f"Error reading {file}: {e}")

# Merge the GST streams if they contain overlapping traces
gst_stream.merge(method=1)

# Print the final merged stream
print(gst_stream)
gst_stream.plot()

# %%
from obspy import UTCDateTime
start_time = UTCDateTime("2020-12-08T00:00:00")
end_time = UTCDateTime("2020-12-15T23:59:59")
gst_stream2 = gst_stream.copy()
# Trim the stream to the desired time range
gst_stream2.trim(start_time, end_time)

# Plot the trimmed stream
gst_stream2.plot()

# %%
