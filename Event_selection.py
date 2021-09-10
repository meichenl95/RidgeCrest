#!/home/meichen/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import os
import subprocess
import datetime
import obspy
import glob
from Event_Trace import *

def dateTime2datetime(date,time):
    date = date.split('-')
    time = time.replace('.',':').split(':')
    timelist = date + time
    timelist = [int(time) for time in timelist]
    dt = datetime.datetime(timelist[0],timelist[1],timelist[2],timelist[3],timelist[4],timelist[5],timelist[6])
    return dt

def select_eventpairs(eventlist):
    eventpairs = {}
    for i in np.arange(len(eventlist)):
        egflist = []
        masterevent = eventlist[i]
        if (masterevent.mag > 3.5 or masterevent.mag < 3):
            continue
        for j in np.arange(i + 1,len(eventlist)):
            egfevent = eventlist[j]
            if (masterevent.time_to(egfevent) >= masterevent.mintimediff and masterevent.mag_to(egfevent) >= masterevent.minmagdiff and masterevent.distance_to(egfevent) <= masterevent.maxdistdiff and egfevent.mag < 3 and egfevent.mag > 2):
                egflist.append(egfevent)
        eventpairs[eventlist[i]] = egflist
    return eventpairs

def signal2noise_and_gcarc(event,netw,stn,chn,loc,tt):
    os.chdir("{}/{}".format(dataset_path,event.id))

    filename = glob.glob("{}.{}.{}.{}.*.SAC".format(netw,stn,loc,chn))
    if (filename == []):
        return False,False,False
    tr = obspy.read(filename[0])[0]

    # frequency range 1
    tr_snr1 = tr.copy()
    print(event.id,netw,stn,loc,chn)
    tr_snr1.filter('bandpass',freqmin=event.freqsnr1low,freqmax=event.freqsnr1high,zerophase=True,corners=2)
    arrivaltime = tr_snr1.stats.starttime + event.origintime + tt
    signal = tr_snr1.slice(arrivaltime-event.timebefore, arrivaltime+event.timeafter).data
    noise = tr_snr1.slice(arrivaltime - tt - 2*(event.timebefore + event.timeafter), arrivaltime - tt - (event.timebefore + event.timeafter)).data
    if (signal.size == 0  or noise.size == 0):
        return False, False, False
    maxsignal1 = np.max(np.abs(signal))
    maxnoise1 = np.max(np.abs(noise))

#    fig,ax = plt.subplots(2,1,figsize=[5,8])
#    ax[0].plot(signal)
#    ax[0].plot(noise)
#    ax[0].set_title(chn)
#    ax[1].loglog(np.fft.rfftfreq(len(signal),d=tr.stats.delta),np.abs(np.fft.rfft(signal)))
#    ax[1].loglog(np.fft.rfftfreq(len(noise),d=tr.stats.delta),np.abs(np.fft.rfft(noise)))
#    plt.show()
#    plt.close()

    # frequency range 2
    tr_snr2 = tr.copy()
    tr_snr2.filter('bandpass',freqmin=event.freqsnr2low,freqmax=event.freqsnr2high,zerophase=True,corners=2)
    arrivaltime = tr_snr2.stats.starttime + event.origintime + tt
    signal = tr_snr2.slice(arrivaltime-event.timebefore, arrivaltime+event.timeafter).data
    noise = tr_snr2.slice(arrivaltime - tt - 2*(event.timebefore+event.timeafter), arrivaltime - tt - (event.timebefore+event.timeafter)).data
    if (signal.size == 0 or noise.size == 0):
        return False,False,False
    maxsignal2 = np.max(np.abs(signal))
    maxnoise2 = np.max(np.abs(noise))

#    fig,ax = plt.subplots(2,1,figsize=[5,8])
#    ax[0].plot(signal)
#    ax[0].plot(noise)
#    ax[0].set_title(chn)
#    ax[1].loglog(np.fft.rfftfreq(len(signal),d=tr.stats.delta),np.abs(np.fft.rfft(signal)))
#    ax[1].loglog(np.fft.rfftfreq(len(noise),d=tr.stats.delta),np.abs(np.fft.rfft(noise)))
#    plt.show()
#    plt.close()
    
    return maxsignal1/maxnoise1, maxsignal2/maxnoise2, tr.stats.sac['gcarc']

def common_stations(master,egf,stationdic):
    masterlist = stationdic.get(master)
    egflist = stationdic.get(egf)
    
    commonlist = []
    for station in masterlist:
        if station in egflist:
            tracesnr = Tracesnr(station.netw,station.stn,station.loc,station.chn)
            index = egflist.index(station)
            tracesnr.set_masteregfsnr(station.snr1,station.snr2,egflist[index].snr1,egflist[index].snr2)
            tracesnr.set_masteregfarrivaltime(station.arrivaltime,egflist[index].arrivaltime)
            tracesnr.set_masteregfgcarc(station.gcarc,egflist[index].gcarc)
            commonlist.append(tracesnr)
    return commonlist

def select_stations(event):
    os.chdir("{}".format(dataset_path))
    phaseinfo = pd.read_csv("{}.phase".format(event.id),skipinitialspace=True,sep="\s+",names=['network','stations','channel','location','phase','signal_onset_quality','pick_quality','traveltime'],usecols=[0,1,2,3,7,9,10,12],skiprows=1,dtype={'network':str,'stations':str,'channel':str,'location':str,'phase':str,'signal_onset_quality':str,'pick_quality':np.float64,'traveltime':np.float64})

    stationlist_S = []
    stationlist_P = []
    for i in np.arange(phaseinfo.shape[0]):
        if ((phaseinfo['channel'][i][1::] == 'HN' or phaseinfo['channel'][i][1::] == 'HE') and phaseinfo['phase'][i] == 'S'):
            snr1, snr2, gcarc = signal2noise_and_gcarc(event,phaseinfo['network'][i],phaseinfo['stations'][i],phaseinfo['channel'][i],phaseinfo['location'][i].replace("-",""),phaseinfo['traveltime'][i])
            tr = Trace(phaseinfo['network'][i],phaseinfo['stations'][i],phaseinfo['location'][i].replace("-",""),phaseinfo['channel'][i])
            if (snr1 != False and snr2 != False):
                tr.set_snr(snr1,snr2)
                tr.set_arrivaltime(phaseinfo['traveltime'][i])
                tr.set_gcarc(gcarc)
                stationlist_S.append(tr)
        elif (phaseinfo['channel'][i][1::] == 'HZ' and phaseinfo['phase'][i] == 'P'):
            snr1, snr2, gcarc = signal2noise_and_gcarc(event,phaseinfo['network'][i],phaseinfo['stations'][i],phaseinfo['channel'][i],phaseinfo['location'][i].replace("-",""),phaseinfo['traveltime'][i])
            tr = Trace(phaseinfo['network'][i],phaseinfo['stations'][i],phaseinfo['location'][i].replace("-",""),phaseinfo['channel'][i])
            if (snr1 != False and snr2 != False):
                tr.set_snr(snr1,snr2)
                tr.set_arrivaltime(phaseinfo['traveltime'][i])
                tr.set_gcarc(gcarc)
                stationlist_P.append(tr)

    return stationlist_S,stationlist_P            
    

def main():
    catalog = pd.read_csv("catalog_Trugman.txt",skipinitialspace=True,sep=" ",names=['eventid','date','time','magnitude','longitude','latitude','depth'],usecols=[0,1,2,3,4,5,6],skiprows=28,dtype={'eventid':str,'date':str,'time':str,'magnitude':np.float64,'longitude':np.float64,'latitude':np.float64,'depth':np.float64})
    
    # make a list of class event
    eventlist = []
    for i in np.arange(catalog.shape[0]):
        dt = dateTime2datetime(catalog['date'][i],catalog['time'][i])
        event = Event(catalog['eventid'][i],catalog['latitude'][i],catalog['longitude'][i],catalog['depth'][i],catalog['magnitude'][i],dt)
        event.mag_to_freqtimewin()
        if (event.mag >= 2 and event.mag < 3.5):
            eventlist.append(event)

    # sort the list based on magnitude in descending order
    eventlist.sort(key=lambda x: x.mag,reverse=True)

    # find pairs of events
    eventpairs = select_eventpairs(eventlist)
    
    with open("eventpairs.pkl","wb") as f:
        pickle.dump(eventpairs, f, pickle.HIGHEST_PROTOCOL)
    print("complete saving eventpairs")
    f.close()
#    with open("eventpairs.pkl","rb") as f:
#        temp = pickle.load(f)

    # compute snr of each station for events
    global dataset_path
    dataset_path = "/home/meichen/work1/RidgeCrest"
    current_path = "/home/meichen/Research/RidgeCrest"
    stationdic_S = {}
    stationdic_P = {}
    for event in eventlist:
        stationlist_S, stationlist_P = select_stations(event)
        stationdic_S[event] = stationlist_S
        stationdic_P[event] = stationlist_P

    os.chdir(current_path)
    with open("stationdic_S.pkl","wb") as f:
        pickle.dump(stationdic_S, f, pickle.HIGHEST_PROTOCOL)
    f.close()
    with open("stationdic_P.pkl","wb") as f:
        pickle.dump(stationdic_P, f, pickle.HIGHEST_PROTOCOL)
    f.close()
    print("complete saving stationdic")

    # find same stations for each pair of master/egf
    eventpairs_with_stations_S = {}
    eventpairs_with_stations_P = {}
    for master, egfs in eventpairs.items():
        egf_stations_S = {}
        egf_stations_P = {}
        for egf in egfs:
            stations_S = common_stations(master, egf, stationdic_S)
            stations_P = common_stations(master, egf, stationdic_P)
            egf_stations_S[egf] = stations_S
            egf_stations_P[egf] = stations_P
        eventpairs_with_stations_S[master] = egf_stations_S
        eventpairs_with_stations_P[master] = egf_stations_P
    with open("eventpairs_with_stations_S.pkl","wb") as f:
        pickle.dump(eventpairs_with_stations_S, f, pickle.HIGHEST_PROTOCOL)
    f.close()
    with open("eventpairs_with_stations_P.pkl","wb") as f:
        pickle.dump(eventpairs_with_stations_P, f, pickle.HIGHEST_PROTOCOL)
    f.close()
    print("complete saving eventpairs with stations")
main()
