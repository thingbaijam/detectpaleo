import numpy as np
import csv
import matplotlib.pyplot as plt
from typing import List

def sinesqrt(x):
    # sinesqrt function - note the angle is in radians
    slip_x = 1.311*np.sqrt(np.sin(np.pi*x))
    return slip_x

def plot_sinesqrt():
    x = np.linspace(0.0,1,100)
    y = []
    for val in x:
        y.append(sinesqrt(val))
    plt.plot(x, y, 'o-')
    
def prob_detect_paleoslip(slip):
    #  The Wrightwood model
    #  Table from UCERF Appendix I
    slips = [0, 0.10, 0.20, 0.30, 1.00, 2.00, 4.00]
    prob_detect = [0, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99]

    #  I like  This fits an exponential model
    # prob = 1-np.exp(-1.6*slip)

    # But UCERF3 uses  1-D interpolation .
    # check prob is within bounds

    prob = np.interp(slip, slips, prob_detect)


    if (prob<0.0) | (prob>1.0):
        print("*** Probability of detecting paleo-slip is not correct! ")
    return prob

def plot_prob_detect_paleoslip(slip = [-1]):
    if slip[0]==-1:
        slip = [0, 0.10, 0.20, 0.30, 1.00, 2.00, 4.00]
    y = []
    for s in slip:
        y.append(prob_detect_paleoslip(s))
    plt.plot(slip, y, 'o-')

def avg_surf_slip(mag):
    # get surface-average slip from mag
    # For the record: I think this UCERF approach is total wrong!!
    # from the appendix I - Figure I8
    avgslip = 2.32*(mag-6.12)
    if avgslip<=0.0:
        avgslip = 0.007
    return avgslip

def prob_detectsurfrupture(mag, doplot=False, model='UCERF3'):
    # from Appendix I UCERF3 model - WC1994 based model
    mags = []
    probs = []
    if model == "UCERF3":
        with open('mag_surfslipprob.csv', 'r') as file:
            csvr = csv.reader(file)
            for row in csvr:
                mags.append(float(row[0]))
                probs.append(float(row[1]))

    if mag<=5.5:
        p = 0.0
    elif mag>=8.0:
        p = 1.0
    else:
        p = np.interp(mag, mags, probs)
    if doplot is True:
        plt.plot(mags, probs, 'o')
    return p

def plot_prob_detect_surfslip(model='UCERF3'):
    # from Appendix I UCERF3 model - WC1994 based model
    M = []
    prob = []
    if model=='UCERF3':
        with open('mag_surfslipprob.csv', 'r') as file:
           csvr = csv.reader(file)
           for row in csvr:
               M.append(float(row[0]))
               prob.append(float(row[1]))
    plt.plot(M, prob, 'o')


"""
function [slips] =  slip_distribution(mean_slip)
% This is surely dependend on location of falult
% and also should not be lognormal but extreme valu distirbution
% but for this test - I am taking it lognormal with
%

sigma= 0.6;
mu=log(mean_slip);
slips = lognrnd(mu,sigma,1000,1);

% To get back the mean slip ----------
% modeled_mean_slip = exp(mean(log(slips)));
"""
