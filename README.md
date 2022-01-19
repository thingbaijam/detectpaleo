# A simple independent implemention of the UCERF3 workflow


See: Weldon RJ, Biasi GP. Appendix I: Probability of detection of ground rupture at paleoseismic sites. 
 US Geol. Surv. Open‐File Rept. 2013‐1165‐I, and California Geol. Surv. Special Rept. 228‐I. 2013.


Note: The units used here are: slip in metres,and magntude Mw 

#### Key functions

##### prob_detectsurfrup(magnitude=-1, model="UCERf3", doplot = False)

Also, <b> plot_prob_detectsurfrup(model="UCERF3") </b>

Model can be also a name of a file (along with it's relative or absolute location). 
File should contain list of comma separated values: magnitude, probabilty

mag2avg_surfslip(mag, model='UCERF3', doplot = False)

Alternately, model = NZNSHM2022 (derived from simplied mag-area scaling, dominant strike-slip faulting) 
or TMG2017 (pending implementation) 

normslip_prof(avg_surfslip, x_by_RL= np.linspace(0.05, 0.5,50), model="sinesqrt")

##### prob_slip_profpoint (avg_surfslip, model = "UCERF3")
   
##### prob_detectpaleoslip(sampledslip, prob_sampledslip = 1, model="wrightwood")

#### Key inputs

mag_surfslipprob.csv
