#!/home/meichen/anaconda3/bin/python3

import numpy as np
import pandas as pd
import os
import subprocess
import datetime
import obspy
import glob

class Trace:
    def __init__(self,init_netw,init_stn,init_loc,init_snr):
        self.netw = init_netw
        self.stn = init_stn
        self.loc = init_loc
        self.snr = init_snr

class Event:
    freqmin = 0
    freqmax = 50
    timebefore = 0
    timeafter = 0
    def __init__(self,init_id,init_lat,init_lon,init_dep,init_mag,init_datetime):
        self.id = init_id
        self.lat = init_lat
        self.lon = init_lon
        self.dep = init_dep
        self.mag = init_mag
        self.datetime = init_datetime # in datetime format

    # Distance to event2 in km
    def distance_to(self,event2):
        myradius = 6371 - self.dep
        radius2 = 6371 - event2.dep
        arg = np.cos(np.deg2rad(self.lat)) * np.cos(np.deg2rad(event2.lat)) * np.cos(np.deg2rad(self.lon) - np.deg2rad(event2.lon)) + np.sin(np.deg2rad(self.lat)) * np.sin(np.deg2rad(event2.lat))
        dist = np.sqrt(myradius**2 + radius2**2 - 2 * myradius * radius2 * arg)
        return dist
        
    # Difference of origin time in seconds
    def time_to(self,event2):
        return np.abs((self.datetime - event2.datetime).seconds)

    # Difference of moment magnitude
    def mag_to(self,event2):
        return np.abs(self.mag - event2.mag)

    # compute the frequency range and time window
    def mag_to_freqtimewin(self):
        if (self.mag >= 6):
            self.freqmin = 0.05
            self.freqmax = 5
            self.timebefore = 20
            self.timeafter = 100
            self.origintime = 90
        elif (self.mag >= 5 and self.mag < 6):
            self.freqmin = 0.1
            self.freqmax = 5
            self.timebefore = 10
            self.timeafter = 50
            self.origintime = 75
        elif (self.mag >= 4 and self.mag < 5):
            self.freqmin = 1
            self.freqmax = 10
            self.timebefore = 5
            self.timeafter = 25
            self.origintime = 60
        elif (self.mag >= 3 and self.mag < 4):
            self.freqmin = 2
            self.freqmax = 20
            self.timebefore = 3
            self.timeafter = 12
            self.origintime = 45
        elif (self.mag >= 2 and self.mag < 3):
            self.freqmin = 2
            self.freqmax = 50
            self.timebefore = 2
            self.timeafter = 4
            self.origintime = 30
        else:
            self.freqmin = 5
            self.freqmax = 50
            self.timebefore = 0.5
            self.timeafter = 2.5
            self.origintime = 15

    # string representation
    def __str__(self):
        return "id: {}, mag: {}".format(self.id, self.mag)


def dateTime2datetime(date,time):
    date = date.split('-')
    time = time.replace('.',':').split(':')
    timelist = date + time
    timelist = [int(time) for time in timelist]
    dt = datetime.datetime(timelist[0],timelist[1],timelist[2],timelist[3],timelist[4],timelist[5],timelist[6])
    return dt

def select_eventpairs(eventlist,minmagdiff,mintimediff,maxdistdiff):
    eventpairs = {}
    for i in np.arange(len(eventlist)):
        egflist = []
        masterevent = eventlist[i]
        for j in np.arange(i + 1,len(eventlist)):
            egfevent = eventlist[j]
            masterevent.time_to(egfevent)
            if (eventlist[i].time_to(eventlist[j]) >= mintimediff and eventlist[i].mag_to(eventlist[j]) >= minmagdiff and eventlist[i].distance_to(eventlist[j]) <= maxdistdiff):
                egflist.append(eventlist[j])
        eventpairs[eventlist[i]] = egflist
    return eventpairs

def signal2noise(event,netw,stn,chn,loc,tt):
    os.chdir("{}/{}".format(dataset_path,event.id))

    filename = glob.glob("{}.{}.{}.{}.*.SAC".format(netw,stn,loc,chn))
    if (filename == []):
        return False
    tr = obspy.read(filename)[0]
    tr.filter('bandpass',freqmin=event.freqmin,freqmax=event.freqmax,zerophase=True,corners=2)
    arrivaltime = tr.stats.starttime + event.origintime + tt
    signal = tr.slice(arrivaltime-event.timebefore, arrivaltime+event.timeafter).data
    maxsignal = np.max(np.abs(signal))
    noise = tr.slice(event.origintime/3,event.origintime*2/3).data
    maxnoise = np.max(np.abs(noise))
    
    return maxsignal/maxnoise

def select_stations(event):
    os.chdir("{}".format(dataset_path))
    phaseinfo = pd.read_csv("{}.phase".format(event.id),skipinitialspace=True,sep="\s+",names=['network','stations','channel','location','phase','signal_onset_quality','pick_quality','traveltime'],usecols=[0,1,2,3,7,9,10,12],skiprows=1,dtype={'network':str,'stations':str,'channel':str,'location':str,'phase':str,'signal_onset_quality':str,'pick_quality':np.float64,'traveltime':np.float64})

    stationlist_S = []
    stationlist_P = []
    for i in np.arange(phaseinfo.shape[0]):
        if (phaseinfo['channel'][i][1::] == 'HN' and phaseinfo['phase'][i] == 'S'):
            snr = signal2noise(event,phaseinfo['network'][i],phaseinfo['stations'][i],phaseinfo['channel'][i],phaseinfo['location'][i].replace("-",""),phaseinfo['traveltime'][i])
            if (snr != False):
                stationlist_S.append(Trace(phaseinfo['network'][i],phaseinfo['stations'][i],phaseinfo['location'][i],snr))
        elif (phaseinfo['channel'][i][1::] == 'HZ' and phaseinfo['phase'][i] == 'P'):
            snr = signal2noise(event,phaseinfo['network'][i],phaseinfo['stations'][i],phaseinfo['channel'][i],phaseinfo['location'][i].replace("-",""),phaseinfo['traveltime'][i])
            if (snr != False):
                stationlist_P.append(Trace(phaseinfo['network'][i],phaseinfo['stations'][i],phaseinfo['location'][i],snr))

    return stationlist_S,stationlist_P            
    

def main():
    catalog = pd.read_csv("catalog_test.txt",skipinitialspace=True,sep=" ",names=['eventid','date','time','magnitude','longitude','latitude','depth'],usecols=[0,1,2,3,4,5,6],skiprows=28,dtype={'eventid':str,'date':str,'time':str,'magnitude':np.float64,'longitude':np.float64,'latitude':np.float64,'depth':np.float64})
    
    # make a list of class event
    eventlist = []
    for i in np.arange(catalog.shape[0]):
        dt = dateTime2datetime(catalog['date'][i],catalog['time'][i])
        event = Event(catalog['eventid'][i],catalog['latitude'][i],catalog['longitude'][i],catalog['depth'][i],catalog['magnitude'][i],dt)
        event.mag_to_freqtimewin()
        eventlist.append(event)

    # sort the list based on magnitude in descending order
    eventlist.sort(key=lambda x: x.mag,reverse=True)

    # find pairs of events
    minmagdiff = 1.0
    mintimediff = 10 # seconds
    maxdistdiff = 100 # km
    eventpairs = select_eventpairs(eventlist,minmagdiff,mintimediff,maxdistdiff)
    
    # compute snr of each station for events
    global dataset_path
    dataset_path = "/home/meichen/work1/RidgeCrest"
    stationdic_S = {}
    stationdic_P = {}
    for event in eventlist:
        print(event)
        stationlist_S, stationlist_P = select_stations(event)
        stationdic_S[event] = stationlist_S
        stationdic_P[event] = stationlist_P

main()
