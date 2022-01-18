# A simple independent implemention of the UCERF3 workflow.
# ------------------------------------------------------------------------------------
# See: Weldon RJ, Biasi GP. Appendix I: Probability of detection of ground rupture at paleoseismic sites. 
# US Geol. Surv. Open‐File Rept. 2013‐1165‐I, and California Geol. Surv. Special Rept. 228‐I. 2013.
# ------------------------------------------------------------------------------------
# Units used here are: slip metres, magntude Mw 
#
# Key functions:
#
# prob_detectsurfrup(magnitude=-1, model="UCERf3", doplot = False)
# Also, plot_prob_detectsurfrup(model="UCERF3")
# Model can be also a name of a file (along with it's relative or absolute location). 
# File contains list of comma separated values: magnitude, probabilty

# 
# (2) mag2avg_surfslip(mag, model='UCERF3', doplot = False)
#      note: model = NZNSHM2022 or TMG2017  ... 
# This could in be two step => mag2area, 
# area2avgslip assuming that averagelsip and # average surface slip concides reasonably, 
# espcailly for strikeslip events 
#     
# (3) normslip_prof(avg_surfslip, x_by_RL= np.linspace(0.05, 0.5,50), model='sinesqrt')
#   default x_by_RL = np.linspace(0.05, 0.5,50)
#       sinesqrt(x)* avgslip 
#  
# (4) prob_slip_profpoint (avg_surfslip, model = 'UCERF3')
#    = lognormal 
#    == I prefer extrme value distribution
#   
# (5) prob_detectpaleoslip(sampledslip, prob_sampledslip = 1, model='wrightwood2013')

# What do you need?
# (1) mag_surfslipprob.csv


import numpy as np
import csv
import matplotlib.pyplot as plt
#from typing import List
import statistics as stat
import scipy.stats as statdist
import numpy as np

#------------------------------------------------------------------------------------------------------
def prob_detectsurfrup(magnitude=-1, model="UCERf3", doplot = False):
    # from Appendix I UCERF3 model - WC1994 based model
    mags = []
    probs = []
    
    if model == "UCERF3":
       # from Appendix I UCERF3 model - WC1994 based model
       myfile = "files/UCERF3_mag_surfslipprob.csv"
    else:
      myfile = model
    
    with open(myfile, "r") as file:
            # assuming no header
            csvr = csv.reader(file)
            for row in csvr:
                mags.append(float(row[0]))
                probs.append(float(row[1]))
    if not doplot:
       if magnitude<=min(mags):
          p = 0.0
       elif magnitude>=max(mags):
          p = 1.0
       else:
          p = np.interp(magnitude, mags, probs)
           
    if doplot is True:
        plt.plot(mags, probs, 'o')
        plt.xlabel("Magnitude (Mw)");
        plt.ylabel("Probability of surface rupture, Psr");
        p = None
   
    return p


def plot_prob_detectsurfrup(model="UCERF3"):
   prob_detectsurfrup(model=model, doplot = True)

#-------------------------------------------------------------------------------------------

def mag2avg_surfslip(magnitude=-1, model="UCERF3"):
    if model=="UCERF3":
       # For the record: I think this scaling relation is awkward. 
       # from the appendix I - Figure I8
       avg_surfslip = 2.32*(magnitude-6.12)
       if avg_surfslip<=0.0:
          avg_surfslip = 0.007
    elif model=="TMG2017":
       # a link to eqsrcpy package .. perhaps, not now
       pass
    elif model=="NZNSHM22":
       # 
       pass
    return avg_surfslip


def plot_mag2avg_surfslip(model="UCERF3"):
    mags = np.linspace(5,8.5, 50)

    avg_surfslip = []
    for m in mags:
       avg_surfslip.append(mag2avg_surfslip(m))

    plt.plot(mags, avg_surfslip, 'o')
    plt.xlabel('Magnitude (Mw)');
    plt.ylabel('Average surface slip (m)');
    plt.rcParams.update({'font.size': 14})
    # fig = plt.gcf()
    # fig.set_size_inches(7, 4)

#-------------------------------------------------------------------------------------------

# (3) normslip_prof(avg_surfslip, x_by_RL= np.linspace(0.05, 0.5,50), model='sinesqrt')
#   default x_by_RL = np.linspace(0.05, 0.5,50)
#       sinesqrt(x)* avgslip 

def slip_profile(avg_surfslip, x_by_RL=np.linspace(0.05, 0.5,50), model='sinesqrt'):
    # As of now- we do not have any model other than sinesqrt
    
    # slip at the point on the fault (x_by_RL) is distrbuted according to sinesqrt function
    slip_x = sinesqrt(x_by_RL)*avg_surfslip
    return slip_x

def sinesqrt(x):
    # sinesqrt function 
    # note that the angle(s) is in radians
    normslip_x = 1.311*np.sqrt(np.sin(np.pi*x))
    return normslip_x

def plot_sinesqrt():
    x = np.linspace(0.0,1,100)
    y = []
    for val in x:
        y.append(sinesqrt(val))
    plt.plot(x, y, 'o-')
    plt.xlabel('normalized rupture length (x/L)');
    plt.ylabel('Normalized average surface slip');

#-------------------------------------------------------------------------------------------

def prob_detectpaleoslip(sampledslip, prob_sampledslip = 1, model="wrightwood2013"):
    #  As of now, we have no other model except - The Wrightwood model
    #  Table from UCERF Appendix I
    
    slips = [0, 0.10, 0.20, 0.30, 1.00, 2.00, 4.00]
    prob_detect = [0, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99]

    # This fits an exponential model
    # prob = 1-np.exp(-1.6*slip)
    # But UCERF3 uses 1-D interpolation.
 
    prob = np.interp(sampledslip, slips, prob_detect)
    
    if (prob<0.0) | (prob>1.0):
        print("*** Probability of detecting paleo-slip is not correct! ")
    else:
        prob = prob*prob_sampledslip
    return prob

def plot_prob_detect_paleoslip(model="wrightwood2013"):
    slips = [0, 0.10, 0.20, 0.30, 1.00, 2.00, 4.00]
    y = []
    for s in slips:
        y.append(prob_detectpaleoslip(s))
    plt.plot(slips, y, 'o-')
    plt.xlabel('Paleo-slip (m)');
    plt.ylabel('Detection probability');

#-------------------------------------------------------------------------------------------

def prob_slip_profpoint (slip_x, normx=None, model = "UCERF3"):
    # normx is location (normalized such that RL = 1, where RL is rupture length) 
    # along the rupture profile. 
    
    if model=="UCERF3":
        # UCERF3 model is logrnormal 
        # An  approximation for normalized slip is applied here with as \mu = 0.9 and sigma = 0.5
        # \mu is given by log(m) where m is scale parameter 
        # This approximation is fixed - no use for normx
        # UCERF3 appendix I provides no info, but figures to caliberate. not doing that - okay
        
        mu = 0.9  # mean is given by log(m) where m is scale parameter
        sigma = 0.5      
        nsampling = 100  # idealy 10,000
        dist=statdist.lognorm(s= sigma,scale=mu);
        maxslp = 4   # again a rough estimate
        sampledslip = np.linspace(0.001,maxslp,nsampling) 
        sprob = dist.pdf(sampledslip)
        #dx = (maxslp-0.001)/nsampling .. not doing this
        slipprob = sprob/sum(sprob) # normalized so that summation of probabilies = 1
     
    sampledslip = sampledslip*slip_x
    return sampledslip, slipprob









