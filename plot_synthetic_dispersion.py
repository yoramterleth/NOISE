
""" This is the script used to produce the synthetic dispersion used to select a frequency range targeting the base of the glacier.
It is heavily inspiried by the examples that goe with the package disba: https://github.com/keurfonluu/disba.

Experiments heavily inspired by Zhan et al.: Zhan, Z., 2019. Seismic noise interferometry reveals transverse drainage configuration 
beneath the surging Bering Glacier. Geophysical Research Letters, 46(9), pp.4747-4756. """
#%%
import numpy as np
from disba import PhaseDispersion, GroupDispersion
import matplotlib.pyplot as plt


#%% plot the experiments 


velocity_model_SK = np.array([ # from Zhan et al (almost similar)
   [1.0, 3400.00, 1500.00, 850.00], # firn layer 
   [50.0, 3700.00, 1820.00, 917.00], # ice layer top
   [270.0, 3700.00, 1820.00, 917.00], # ice layer main
   [30.0, 3700.00, 1820.00, 917.00], # ice layer bottom 
   [10.0, 1700.00, 250, 1800.00], # till layer
   [5000.0, 3650.00, 1980.00, 2500.00], # sediment layer
   [1000.0,4500.00,2500.00,3000.00] # bedrock layer
])

velocity_model_experiments = [
    np.array([ # from Zhan et al (almost similar)
   [10.0, 3400.00, 1500.00, 600.00], # firn layer increased to 10 m
   [50.0, 3700.00, 1820.00, 917.00], # ice layer top
   [270.0, 3700.00, 1820.00, 917.00], # ice layer
   [30.0, 3700.00, 1820.00, 917.00], # ice layer bottom 
   [10.0, 1700.00, 250, 1800.00], # till 
   [5000.0, 3650.00, 1980.00, 2500.00], # sediment layer
   [1000.0,4500.00,2500.00,3000.00] # bedrock layer
]),
 np.array([ # from Zhan et al (almost similar)
   [1.0, 3400.00, 1500.00, 850.00], # firn layer 
   [50.0, 3700.00, 1820.00, 917.00], # ice layer top
   [320.0, 3700.00, 1820.00, 917.00], # ice layer thickened by 50m
   [30.0, 3700.00, 1820.00, 917.00], # ice layer bottom 
   [10.0, 1700.00, 250, 1800.00], # till layer
   [5000.0, 3650.00, 1980.00, 2500.00], # sediment layer
   [1000.0,4500.00,2500.00,3000.00] # bedrock layer
]),
 np.array([ # from Zhan et al (almost similar)
   [1.0, 3400.00, 1500.00, 850.00], # firn layer 
   [50.0, 3700.00*.97, 1820.00*.97, 917.00], # ice layer top slowed down by 3%
   [270.0, 3700.00, 1820.00, 917.00], # ice layer
   [30.0, 3700.00, 1820.00, 917.00], # ice layer bottom 
   [10.0, 1700.00, 250, 1800.00], # till layer
   [5000.0, 3650.00, 1980.00, 2500.00], # sediment layer
   [1000.0,4500.00,2500.00,3000.00] # bedrock layer
]),
 np.array([ # from Zhan et al (almost similar)
   [1.0, 3400.00, 1500.00, 850.00], # firn layer 
   [50.0, 3700.00, 1820.00, 917.00], # ice layer top
   [270.0, 3700.00, 1820.00, 917.00], # ice layer
   [30.0, 3700.00*.7, 1820.00*.97, 917.00], # ice layer bottom slowed down by 3%
   [10.0, 1700.00, 250, 1800.00], # till layer 
   [5000.0, 3650.00, 1980.00, 2500.00], # sediment layer
   [1000.0,4500.00,2500.00,3000.00] # bedrock layer
]),
 np.array([ # from Zhan et al (almost similar)
   [1.0, 3400.00, 1500.00, 850.00], # firn layer 
   [50.0, 3700.00, 1820.00, 917.00], # ice layer top
   [270.0, 3700.00, 1820.00, 917.00], # ice layer
   [30.0, 3700.00, 1820.00, 917.00], # ice layer bottom 
   [10.0, 1700.00, 250*.75, 1800.00], # till layer slowed down 25%
   [5000.0, 3650.00, 1980.00, 2500.00], # sediment layer
   [1000.0,4500.00,2500.00,3000.00] # bedrock layer
]),
],
experiment_names = ['Firn/snow layer increased from 1m to 10m ', 'Ice thickness increased by 50m','Upper 50 m of ice slowed by 3%','Lower 30m of ice slowed down by 3%','10m Till layer - Vs slowed by 25%']



fig,ax = plt.subplots(5,1,figsize=(5,10),sharex=True)

for i in range(5):

    f = np.arange(.01,5,0.001)
    t = np.sort(1/f)
    # model
    velocity_model = velocity_model_SK/1000
    pd = PhaseDispersion(*velocity_model.T)
    cpr = pd(t, mode=0, wave="rayleigh")
    cpl = pd(t, mode=0, wave="love")

    velocity_model_experiment = velocity_model_experiments[0][i]/1000
    pd_experiment = PhaseDispersion(*velocity_model_experiment.T)
    cpr_experiment = pd_experiment(t, mode=0, wave="rayleigh")
    cpl_experiment = pd_experiment(t, mode=0, wave="love")

    # calculate differences in percentages
    rayleigh_exp = ((cpr[1]-cpr_experiment[1])/cpr[1])*100
    love_exp = ((cpl[1]-cpl_experiment[1])/cpl[1])*100
    
    # plot 
    ax[i].plot(1/cpr[0],rayleigh_exp,label='Rayleigh-waves')
    ax[i].plot(1/cpl[0],love_exp,label='Love-waves')
    ax[i].axvspan(.7,1.5,color='grey',alpha=0.2)
    ax[i].set_ylabel('velocity reduction (%)')
    ax[i].set_ylim(-.5,2.5)
    if i ==4:
        ax[i].set_xlabel('Frequency (Hz)')
    if i ==0:
        ax[i].legend()

    ax[i].grid()
    ax[i].set_title(experiment_names[i],loc='right',verticalalignment='top')
    plt.tight_layout()

 # Extract depth, Vp, and Vs
depth = velocity_model_SK[:, 0]
Vp = velocity_model_SK[:, 1]
Vs = velocity_model_SK[:, 2]

# Create step-like profiles
vs_profile = np.repeat(Vs, 2)[1:]  # Duplicate each Vs value and skip the first
vp_profile = np.repeat(Vp, 2)[1:]  # Duplicate each Vp value and skip the first
depth_profile = np.repeat(np.cumsum(depth), 2)[:-1]  # Cumulative depth, repeat, and trim the last

# Plot the seismic wave model
plt.figure(figsize=(5, 10))

['firn','upper glacier','glacier','bottom glacier','till','sediment','bedrock']

# Plot step-like Vs and Vp profiles
plt.plot(vs_profile/1000, depth_profile, label='$V_s$', color='red', linestyle='--')
plt.plot(vp_profile/1000, depth_profile, label='$V_p$', color='blue', linestyle='--')


# Format the plot
plt.gca().invert_yaxis()  # Depth increases downward
plt.yscale('log')  # Use a linear scale for depth (log is optional)
plt.xlabel('Wave Speed (km/s)', fontsize=12)
plt.ylabel('Depth (m)', fontsize=12)
plt.title('Seismic Wave Speed Reference Model', fontsize=14)
plt.grid(True, which='both', linestyle='--', alpha=0.7)
plt.legend(fontsize=10)

# Show the plot
plt.show()
