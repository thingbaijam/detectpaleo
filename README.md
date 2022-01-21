# A Simple Python Implemention of the UCERF3 workflow


See: Weldon RJ, Biasi GP. Appendix I: Probability of detection of ground rupture at paleoseismic sites. 
 US Geol. Surv. Open‐File Rept. 2013‐1165‐I, and California Geol. Surv. Special Rept. 228‐I. 2013.

#### Usage

from detectpaleo import dpaleo

#### Examples and analysis

See the Jupyter notebooks   
- workflow_general.ipynb : A general setup for detectability of ground rupture at paleoseismic sites 
- workflow_informalUCERF3.ipynb : An informal implementation of UCERf3 models 
- workout_gev_surfslipdist.ipynb : An analysis for Generalized Extreme Value Distribution Model for Surface Slip Distribution   


#### Key Functions

<i> Note: The units used here are: slip in metres,and magntude Mw </i>

#####  Probability of Surface Rupture 
- prob_detectsurfrup(magnitude=-1, model="UCERf3", doplot = False)
  returns probability 
- plot_prob_detectsurfrup(model="UCERF3") 
   - Model can be also a name of a file (along with it's relative or absolute location). 
   - File should contain list of comma separated values: magnitude, probability

##### Average Surface Slip 
- mag2avg_surfslip(magnitude=-1, model="UCERF3")
   - Alternately, model = NZNSHM2022 (derived from simplied mag-area scaling, dominant strike-slip faulting) or TMG2017 (pending implementation) 
   - returns average surface slip
- plot_mag2avg_surfslip(model="UCERF3")

##### Slip Profile
- slip_profile(avg_surfslip, x_by_RL=np.linspace(0.05, 0.5,50), model='sinesqrt')
    - returns slip values at locations x_by_R, based on a sinesqrt functionL. The location is that normalized by rupture length
- sinesqrt(x)
    - returns slip values at locations x_by_RL, based on a sinesqrt function. The location is that normalized by rupture length
- plot_sinesqrt()

##### Detectability of Paleo-Slip 
- prob_detectpaleoslip(sampledslip, prob_sampledslip = 1, model="wrightwood2013", slipfactor=1)
    - slipfactor (typically, <=1.0) is to enable modulations on detectability of slip  
    - returns probability
- plot_prob_detect_paleoslip(model="wrightwood2013")

##### Slip PDF at a Point of the Fault 
- prob_slip_profpoint (slip_x, xi=-1, model = "UCERF3", normalized = False)
    - Alternately, model = "GEV" 
    - return sampledslip, slipprob

- surfslipdist_gev(fRL)
    - returns the GEV parameters (shape, location, scale) for a given location fRL on rupture profile
