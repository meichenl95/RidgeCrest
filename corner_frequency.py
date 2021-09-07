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
    
    # save to csv
    dataframedic = {}
    masterid = []
    mastermag = []
    egfid = []
    egfmag = []
    stnnumber = []
    cornerfrequency_master = []

    masteriddic = {}
    for master, egfs in eventpairs_with_spectralratio.items():
        egfiddic = {}
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
        masteriddic[master.id] = egfiddic

    dataframedic['masterid'] = masterid
    dataframedic['mastermag'] = mastermag
    dataframedic['egfid'] = egfid
    dataframedic['egfmag'] = egfmag
    dataframedic['station_number'] = stnnumber
    dataframedic['fcmaster'] = cornerfrequency_master
    df = pd.DataFrame.from_dict(dataframedic)
    print(df)
    df.to_csv("{}_separateegf.csv".format(filename.replace('.','_').split('_')[-2]))

def main():
    combine_corner_frequency("eventpairs_with_spectralratio_S.pkl")
    combine_corner_frequency("eventpairs_with_spectralratio_P.pkl")
    
main()
