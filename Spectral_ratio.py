#!/home/meichen/anaconda3/bin/python

import numpy as np
import pickle
import os
import obspy
from Event_Trace import *
from mtspec import mtspec
from scipy.signal import tukey
from scipy.interpolate import interp1d

def spectrum(data,delta):
    data = data * tukey(len(data),alpha=0.1)
    nfft = 2**np.ceil(np.log2(10 * len(data)))
    spec, freq = mtspec(data=data,delta=delta,time_bandwidth=3.5,number_of_tapers=5,nfft=nfft)
    spec = np.sqrt(spec)
    return freq, spec

def commonfreq(masterfreq, egffreq):
    return np.max(masterfreq[0],egffreq[0]), np.min(masterfreq[-1],egffreq[-1])

def resample(freq, spec, pointnum, freqmin, freqmax):
    freqlog = np.log10(freq)
    speclog = np.log10(spec)
    f = interp1d(freqlog, speclog)
    freqlog_output = np.linspace(np.log10(freqmin),np.log10(freqmax), num=pointnum)
    speclog_output = f(speclog_output)
    return 10**freqlog_output, 10**speclog_output
    
def fit(freq,spec):


def fit_spectrum(master, egf, stn):
    mastertr = obspy.read("{}/{}/{}.{}.{}.{}.*.SAC".format(dataset_path,master.id,stn.netw,stn.stn,stn.loc,stn.chn))[0]
    egftr = obspy.read("{}/{}/{}.{}.{}.{}.*.SAC".format(dataset_path,egf.id,stn.netw,stn.stn,stn.loc,stn.chn))[0]

    winnum = 5
    pointnum = 800
    masterarrivaltime = mastertr.stats.starttime + master.origintime + tt
    masterwindowlen = (master.timeafter + master.timebefore) / winnum
    egfarrivaltime = egftr.stats.starttime + egf.origintime + tt
    egfwindowlen = (egf.timeafter + egf.timebefore) / winnum
    spec = np.zeros(pointnum)
    for i in np.arange(winnum):
        masterdata = mastertr.slice(masterarrivaltime - master.timebefore + masterwindowlen / 2 * i, masterarrivaltime - master.timebefore + masterwindowlen / 2 * (i + 2))
        masterfreq, masterspec = spectrum(masterdata,mastertr.delta)
        egfdata = egftr.slice(egfarrivaltime - egf.timebefore + egfwindowlen / 2 * i, egfarrivaltime - egf.timebefore + egfwindowlen / 2 * (i + 2))
        egffreq, egfspec = spectrum(egfdata,egftr.delta)
        commonfreqmin, commonfreqmax = commonfreq(masterfreq, egffreq)
        egffreq, egfspec = resample(egffreq,egfspec,pointnum, commonfreqmin, commonfreqmax)
        masterfreq, masterspec = resample(masterfreq,masterspec,pointnum,commonfreqmin,commonfreqmax)
        spec = spec + masterspec / egfspec
    spec = spec / winnum

    fcmaster, fcegf, momentratio = fit(spec)


def main(filename):
    global dataset_path
    dataset_path = "/home/meichen/work1/RidgeCrest"
    current_path = "/home/meichen/Research/RidgeCrest"

    with open(filename,"rb") as f:
        eventpairs_with_stations = pickle.load(f)
    f.close()

    for master, egf_stations in eventpairs_with_stations.items():
        for egf, staions in egf_stations.items():
            for stn in stations:
                fcmaster, fcegf = fit_spectrum(master, egf, stn)

main("eventpairs_with_stations_S.pkl")
