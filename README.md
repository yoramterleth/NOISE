# NOISE
Seismic noise correlation processing and analysis. Built for Terleth et al (in prep.). These scripts rely heavily on the python package seisgo authored by Yang et al. (2022) and available here: https://github.com/xtyangpsp/SeisGo. Some of the seisgo functions weere modified slightly to work within this workflow. These are the noise_v2.py and downloaders.py included here. They are mostly copy pasted from seisgo, with some modified/added functions. To run these scripts, the easiest solution is to download seisgo and extract it in the working directory, then call the modified packages when needed. 

The workflow used here is also heavily inspired by one of the working example notebooks provided with seisgo. 

## Processing steps for each station pair:
- Run processing_noise_correlations.py: this loads the data from mseed dayvols, does some preprocessing, and saves it as station pairs in ASDF. Best run from the command line.
- Run noise_correlation_only.py: performs the cross correlation for the station pairs. Best run from the command line.
- Run dv_v_script.py Performs the dvv computations, plots figures, and saves csvs for more plotting. Best run as a notebook.
- Script_for_main_figure.py can be used to re-create figure 3 in the paper. This calls on the functions in helpers that implement the averaging and snow corrections.
