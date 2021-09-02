#!/home/meichen/anaconda3/bin/python3

import numpy as np
import pandas as pd
import os
import subprocess
import datetime
import obspy

class Event:
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

def select_stations(event):
    dataset_path = "/home/meichen/work1/RidgeCrest"
    os.chdir("{}".format(dataset_path))
    phaseinfo = pd.read_csv("{}.phase".format(event.id),skipinitialspace=True,sep=" ",names=['network','stations','channel','location','phase','signal_onset_quality','pick_quality','traveltime'],usecols[0,1,2,3,7,9,10,12],skiprows=1,dtype={'network':str,'stations':str,'channel':str,'location':str,'phase':str,'signal_onset_quality':str,'pick_quality':np.float64,'traveltime':np.float64})
    for i in np.arange(phaseinfo.shape[0]):
        if
    

def main():
    catalog = pd.read_csv("catalog_test.txt",skipinitialspace=True,sep=" ",names=['eventid','date','time','magnitude','longitude','latitude','depth'],usecols=[0,1,2,3,4,5,6],skiprows=28,dtype={'eventid':str,'date':str,'time':str,'magnitude':np.float64,'longitude':np.float64,'latitude':np.float64,'depth':np.float64})
    
    # make a list of class event
    eventlist = []
    for i in np.arange(catalog.shape[0]):
        dt = dateTime2datetime(catalog['date'][i],catalog['time'][i])
        event = Event(catalog['eventid'][i],catalog['latitude'][i],catalog['longitude'][i],catalog['depth'][i],catalog['magnitude'][i],dt)
        eventlist.append(event)

    # sort the list based on magnitude in descending order
    eventlist.sort(key=lambda x: x.mag,reverse=True)
    print(eventlist)

    # find pairs of events
    minmagdiff = 1.0
    mintimediff = 10 # seconds
    maxdistdiff = 100 # km
    eventpairs = select_eventpairs(eventlist,minmagdiff,mintimediff,maxdistdiff)
    
    # select stations for events
    for event in eventlist:
        stationlist = select_stations(event)

main()
