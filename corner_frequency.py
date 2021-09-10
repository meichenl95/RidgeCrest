#!/home/meichen/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
from Event_Trace import *
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import pickle
import pandas as pd

def smooth(specl,n):
    if (n%2 != 1):
        print("Wrong smooth factor")
        return False
    countsmallindex = 0
    while(specl[countsmallindex] == 0):
        countsmallindex += 1
    countlargeindex = specl.shape[0] - 1
    while(specl[countlargeindex] == 0):
        countlargeindex -= 1
    speclog = list(specl[countsmallindex:countlargeindex])
    speclog_shift = speclog.copy()
    speclog_output = np.array(speclog.copy())
    halfwin = int(n/2)
    for i in np.arange(halfwin):
        speclog_shift = speclog_shift[1:] + [0]
        speclog_output += np.array(speclog_shift)
    speclog_shift = speclog.copy()
    for i in np.arange(halfwin):
        speclog_shift = [0] + speclog_shift[:-1]
        speclog_output += np.array(speclog_shift)
    for i in np.arange(halfwin):
        speclog_output[i] = speclog_output[i] / (halfwin + i + 1)
        speclog_output[-i-1] =  speclog_output[-i-1] / (halfwin + i + 1)
    speclog_output[halfwin:len(speclog)-halfwin] = speclog_output[halfwin:len(speclog)-halfwin] / n

    return np.array([0]*countsmallindex + list(speclog_output) + [0]*(specl.shape[0]-countlargeindex))

def resample(commonfreqmin, commonfreqmax, spec, freqlog_to_align):
    freqlog = np.linspace(np.log10(commonfreqmin),np.log10(commonfreqmax),num=spec.shape[0])
    speclog = np.log10(spec)
    countsmallindex = np.sum(freqlog[0] > freqlog_to_align)
    countlargeindex = freqlog_to_align.shape[0] - np.sum(freqlog[-1] < freqlog_to_align)
    f = interp1d(freqlog, speclog)
    speclog_output = np.zeros(freqlog_to_align.shape)
    speclog_output[countsmallindex:countlargeindex] = f(freqlog_to_align[countsmallindex:countlargeindex])
    count_output = np.zeros(freqlog_to_align.shape)
    count_output[countsmallindex:countlargeindex-1] = 1.0

    return speclog_output, count_output

def func(x,a,b,c):
    return np.log10(a) + np.log10(1 + x**2 / b**2) - np.log10(1 + x**2 / c**2)

def bootstrapfit(freq,spec,iternum):
    spec = np.log10(spec)
    poptlist = []
    popt, pcov = curve_fit(func, freq, spec, bounds=([0,0.0001,0.0001],[    1e5,1000,1000]),method='trf',loss='huber',f_scale=0.1)
    try:
        popt, pcov = curve_fit(func, freq, spec, bounds=([0,0.0001,0.0001],[    1e5,1000,1000]),method='trf',loss='huber',f_scale=0.1)
        res = spec - func(freq, *popt)
        for i in np.arange(iternum):
            random_index = np.random.randint(0,len(spec),size=len(spec))
            new_spec = func(freq, *popt) + [res[j] for j in random_index]
            try:
                new_popt, new_pcov = curve_fit(func, freq, new_spec, bounds=    ([0,0.0001,0.0001],[1e5,1000,1000]), method='trf', loss='huber',f_scale=0.1)
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

def combine_corner_frequency(filename):
    with open(filename,"rb") as f:
        eventpairs_with_spectralratio = pickle.load(f)
    f.close()
    
    minmastersnr = 2
    minegfsnr = 2
    maxstdmaster = 0.1
    minstngcarc = 0.1
    minstnnum = 5
    iternum = 50
    smoothfactor = 5
    
    # save to csv
    dataframedic = {}
    masterid = []
    mastermag = []
    egfid = []
    egfmag = []
    stnnumber = []
    cornerfrequency_master = []
    list_uniq_masterid = []
    list_uniq_mastermag = []
    list_uniq_fc = []
    list_uniq_egfnum = []

    masteriddic = {}
    freqlog = np.linspace(-2,2,num=1000)
    for master, egfs in eventpairs_with_spectralratio.items():
        egfiddic = {}
        uniq_fc = 0
        uniq_egfnum = 0
        for egf, stations in egfs.items():
            stnnum = 0
            fcmaster = 0
            speclog = np.zeros(freqlog.shape)
            countsum = np.zeros(freqlog.shape)
            for stn, stninfo in stations.items():
                if (stn.mastersnr1 < minmastersnr or stn.mastersnr2 < minmastersnr or stn.egfsnr1 < minegfsnr or stn.egfsnr2 < minegfsnr):
                    continue
                elif (stn.mastergcarc < minstngcarc or stn.egfgcarc < minstngcarc):
                    continue
                else:
                    spec = stninfo['spectral_ratio']
                    specl, count = resample(stninfo['common_freq_min'],stninfo['common_freq_max'],spec,freqlog)
                    specl = smooth(specl,smoothfactor)
                    speclog += specl
                    countsum += count
                    stnnum += 1
            if (stnnum > 0):
                speclog = np.array([speclog[i]/countsum[i] if countsum[i]>0 else 0 for i in np.arange(speclog.shape[0])])
                freq = 10**freqlog[countsum>0]
                spec = 10**speclog[countsum>0]
                momentratio, fcegf, fcmaster, stdfcmaster = bootstrapfit(freq, spec, iternum)
                if (stdfcmaster < maxstdmaster and fcegf > fcmaster):
                    masterid.append(master.id)
                    mastermag.append(master.mag)
                    egfid.append(egf.id)
                    egfmag.append(egf.mag)
                    stnnumber.append(stnnum)
                    cornerfrequency_master.append(fcmaster)
                    egfiddic[egf.id] = {'station_number': stnnum, 'fcmaster': fcmaster}
                    if (stnnum >= minstnnum):
                        uniq_fc += fcmaster
                        uniq_egfnum += 1
        masteriddic[master.id] = egfiddic
        if (uniq_egfnum > 0):
            list_uniq_masterid.append(master.id)
            list_uniq_mastermag.append(master.mag)
            list_uniq_fc.append(uniq_fc / uniq_egfnum)
            list_uniq_egfnum.append(uniq_egfnum)

    dataframedic['masterid'] = masterid
    dataframedic['mastermag'] = mastermag
    dataframedic['egfid'] = egfid
    dataframedic['egfmag'] = egfmag
    dataframedic['station_number'] = stnnumber
    dataframedic['fcmaster'] = cornerfrequency_master
    df = pd.DataFrame.from_dict(dataframedic)
    print(df)
    df.to_csv("{}_separateegf.csv".format(filename.replace('.','_').split('_')[-2]))

    uniqdataframedic = {}
    uniqdataframedic['masterid'] = list_uniq_masterid
    uniqdataframedic['master_magnitude'] = list_uniq_mastermag
    uniqdataframedic['fc'] = list_uniq_fc
    uniqdataframedic['egf_number'] = list_uniq_egfnum
    df = pd.DataFrame.from_dict(uniqdataframedic)
    print(df)
    df.to_csv("{}_uniqmaster.csv".format(filename.replace('.','_').split('_')[-2]))

    plt.scatter(list_uniq_mastermag,list_uniq_fc)
    plt.yscale('log')
    plt.xlabel("Magnitude")
    plt.ylabel("Corner frequency")
    plt.title(filename.replace('.','_').split('_')[-2])
    plt.savefig("{}_uniqmaster.pdf".format(filename.replace('.','_').split('_')[-2]))
    plt.close()

def main():
    combine_corner_frequency("eventpairs_with_spectralratio_S.pkl")
    combine_corner_frequency("eventpairs_with_spectralratio_P.pkl")
    
main()
