#!/home/meichen/anaconda3/bin/python3

import numpy as np

class Tracename:
    def __init__(self,init_netw,init_stn,init_loc):
        self.netw = init_netw
        self.stn = init_stn
        self.loc = init_loc
    def __str__(self):
        return "netw: {}, stn: {}, loc: {}".format(self.netw,self.stn,self.loc)

class Trace(Tracename):
    snr = 0;
    def set_snr(self,init_snr):
        self.snr = init_snr
    def __eq__(self, trace2):
        return ((self.netw == trace2.netw) and (self.stn == trace2.stn) and (self.loc == trace2.loc))

class Tracesnr(Tracename):
    mastersnr = 0
    egfsnr = 0
    def set_masteregfsnr(self,master,egf):
        self.mastersnr = master
        self.egfsnr = egf

class Event:
    freqmin = 0
    freqmax = 45
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
            self.freqmax = 45
            self.timebefore = 2
            self.timeafter = 4
            self.origintime = 30
        else:
            self.freqmin = 5
            self.freqmax = 45
            self.timebefore = 0.5
            self.timeafter = 2.5
            self.origintime = 15

    # string representation
    def __str__(self):
        return "id: {}, mag: {}".format(self.id, self.mag)
