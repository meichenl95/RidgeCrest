#!/home/meichen/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
from Event_Trace import *
import pickle
import pandas as pd

def combine_corner_frequency(filename):
    with open(filename,"rb") as f:
        eventpairs_with_spectralratio = pickle.load(f)
    f.close()
    
    minmastersnr = 1
    minegfsnr = 1
    maxstdmaster = 0.01
    minstnnum = 5
    
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
    for master, egfs in eventpairs_with_spectralratio.items():
        egfiddic = {}
        uniq_fc = 0
        uniq_egfnum = 0
        for egf, stations in egfs.items():
            stnnum = 0
            fcmaster = 0
            for stn, stninfo in stations.items():
                if (stn.mastersnr < minmastersnr or stn.egfsnr < minegfsnr):
                    continue
                elif (stninfo['fcegf'] <= stninfo['fcmaster']):
                    continue
                elif (stninfo['stdmaster'] > maxstdmaster):
                    continue
                else:
                    fcmaster += stninfo['fcmaster']
                    stnnum += 1
            if (stnnum > 0):
                fcmaster = fcmaster / stnnum
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

def main():
    combine_corner_frequency("eventpairs_with_spectralratio_S.pkl")
    combine_corner_frequency("eventpairs_with_spectralratio_P.pkl")
    
main()
