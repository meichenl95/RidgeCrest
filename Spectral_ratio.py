#!/home/meichen/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import obspy
from Event_Trace import *
from mtspec import mtspec
from scipy.signal import tukey
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

def spectrum(data,delta):
    data = data * tukey(len(data), alpha=0.1)
    nfft = int(2**np.ceil(np.log2(10 * len(data))))
    spec, freq = mtspec(data=data,delta=delta,time_bandwidth=3.5,number_of_tapers=5,nfft=nfft)
    spec = np.sqrt(spec)
    return freq, spec

def commonfreq(master, egf, masterfreq, egffreq):
    return np.max([masterfreq[1],egffreq[1]]), np.min([master.freqmax,egf.freqmax])

def resample(freq, spec, pointnum, freqmin, freqmax):
    freq = freq[1::]
    spec = spec[1::]
    freqlog = np.log10(freq)
    speclog = np.log10(spec)
    f = interp1d(freqlog, speclog)
    freqlog_output = np.linspace(np.log10(freqmin),np.log10(freqmax), num=pointnum)
    speclog_output = f(freqlog_output)
    return 10**freqlog_output, 10**speclog_output
    
def func(x,a,b,c):
    return np.log10(a) + np.log10(1 + x**2 / b**2) - np.log10(1 + x**2 / c**2)

def bootstrapfit(freq,spec,iternum):
    spec = np.log10(spec)
    poptlist = []
    try:
        popt, pcov = curve_fit(func, freq, spec, bounds=([0,0.0001,0.0001],[1e5,1000,1000]),method='trf',loss='huber',f_scale=0.1)
        res = spec - func(freq, *popt)
        for i in np.arange(iternum):
            random_index = np.random.randint(0,len(spec),size=len(spec))
            new_spec = func(freq, *popt) + [res[j] for j in random_index]
            try:
                new_popt, new_pcov = curve_fit(func, freq, new_spec, bounds=([0,0.0001,0.0001],[1e5,1000,1000]), method='trf', loss='huber',f_scale=0.1)
                poptlist.append(new_popt)
            except:
                print("curve_fit failed")
        poptlist = np.array(poptlist)
        meanmomentratio = np.mean(poptlist[:,0])
        meanfcegf = 10**np.mean(np.log10(poptlist[:,1]))
        meanfcmaster = 10**np.mean(np.log10(poptlist[:,2]))
        stdfcmaster = np.std(np.log10(poptlist[:,2]),ddof=1)
        
        return meanmomentratio, meanfcegf, meanfcmaster, stdfcmaster
    except:
        print("curve_fit failed")
        return False,False,False,False

def fit_spectrum(master, egf, stn):
    mastertr = obspy.read("{}/{}/{}.{}.{}.{}.*.SAC".format(dataset_path,master.id,stn.netw,stn.stn,stn.loc,stn.chn))[0]
    mastertr.detrend('demean')
    egftr = obspy.read("{}/{}/{}.{}.{}.{}.*.SAC".format(dataset_path,egf.id,stn.netw,stn.stn,stn.loc,stn.chn))[0]
    egftr.detrend('demean')

    winnum = 5
    pointnum = 800
    masterarrivaltime = mastertr.stats.starttime + master.origintime + stn.masterarrivaltime
    masterwindowlen = (master.timeafter + master.timebefore) / winnum
    egfarrivaltime = egftr.stats.starttime + egf.origintime + stn.egfarrivaltime
    egfwindowlen = (egf.timeafter + egf.timebefore) / winnum
    spec = np.zeros(pointnum)
    for i in np.arange(winnum):
        masterdata = mastertr.slice(masterarrivaltime - master.timebefore + masterwindowlen / 2 * i, masterarrivaltime - master.timebefore + masterwindowlen / 2 * (i + 2)).data
        masterfreq, masterspec = spectrum(masterdata,mastertr.stats.delta)
        egfdata = egftr.slice(egfarrivaltime - egf.timebefore + egfwindowlen / 2 * i, egfarrivaltime - egf.timebefore + egfwindowlen / 2 * (i + 2)).data
        egffreq, egfspec = spectrum(egfdata,egftr.stats.delta)
        commonfreqmin, commonfreqmax = commonfreq(master, egf, masterfreq,egffreq)
        egffreq, egfspec = resample(egffreq,egfspec,pointnum, commonfreqmin, commonfreqmax)
        masterfreq, masterspec = resample(masterfreq,masterspec,pointnum,commonfreqmin,commonfreqmax)
        spec = spec + masterspec / egfspec
    spec = spec / winnum

    meanmomentratio, meanfcegf, meanfcmaster, stdfcmaster = bootstrapfit(masterfreq,spec,100)
#    print(master.id,master.mag, egf.id,egf.mag,meanmomentratio,meanfcegf,meanfcmaster,stdfcmaster)

    return meanmomentratio, meanfcegf, meanfcmaster, stdfcmaster


def pairs_spectralratio(filename):
    global dataset_path
    dataset_path = "/home/meichen/work1/RidgeCrest"
    current_path = "/home/meichen/Research/RidgeCrest"

    with open(filename,"rb") as f:
        eventpairs_with_stations = pickle.load(f)
    f.close()

    pairs_with_spectralratio = {}
    for master, egf_stations in eventpairs_with_stations.items():
        masterdic = {}
        for egf, stations in egf_stations.items():
            egfdic = {}
            for stn in stations:
                momentratio, fcegf, fcmaster, stdmaster = fit_spectrum(master, egf, stn)
                if (fcmaster):
                    egfdic[stn] = {"momentratio": momentratio, "fcegf": fcegf, "fcmaster": fcmaster, "stdmaster": stdmaster}
            masterdic[egf] = egfdic
        pairs_with_spectralratio[master] = masterdic

    return pairs_with_spectralratio

def main():
    pairs_with_spectralratio_S = pairs_spectralratio("eventpairs_with_stations_S.pkl")
    with open("eventpairs_with_spectralratio_S.pkl","wb") as f:
        pickle.dump(pairs_with_spectralratio_S, f, pickle.HIGHEST_PROTOCOL)
    f.close()
    pairs_with_spectralratio_P = pairs_spectralratio("eventpairs_with_stations_P.pkl")
    with open("eventpairs_with_spectralratio_P.pkl","wb") as f:
        pickle.dump(pairs_with_spectralratio_P, f, pickle.HIGHEST_PROTOCOL)
    f.close()

main()
