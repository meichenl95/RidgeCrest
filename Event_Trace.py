#!/home/meichen/anaconda3/bin/python3

import numpy as np

class Tracename:
    def __init__(self,init_netw,init_stn,init_loc,init_chn):
        self.netw = init_netw
        self.stn = init_stn
        self.loc = init_loc
        self.chn = init_chn
    def __str__(self):
        return "netw: {}, stn: {}, loc: {}, chn: {}".format(self.netw,self.stn,self.loc,self.chn)

class Trace(Tracename):
    snr1 = 0;
    snr2 = 0;
    arrivaltime = 0
    gcarc = 0
    def set_snr(self,init_snr1,init_snr2):
        self.snr1 = init_snr1
        self.snr2 = init_snr2
    def set_arrivaltime(self,init_arrivaltime):
        self.arrivaltime = init_arrivaltime
    def set_gcarc(self, init_gcarc):
        self.gcarc = init_gcarc
    def __eq__(self, trace2):
        return ((self.netw == trace2.netw) and (self.stn == trace2.stn) and (self.loc == trace2.loc) and (self.chn == self.chn))

class Tracesnr(Tracename):
    mastersnr1 = 0
    mastersnr2 = 0
    egfsnr1 = 0
    egfsnr2 = 0
    masterarrivaltime = 0
    egfarrivaltime = 0
    mastergcarc = 0
    egfgcarc = 0
    def set_masteregfsnr(self,mastersnr1,mastersnr2,egfsnr1,egfsnr2):
        self.mastersnr1 = mastersnr1
        self.mastersnr2 = mastersnr2
        self.egfsnr1 = egfsnr1
        self.egfsnr2 = egfsnr2
    def set_masteregfarrivaltime(self, mastertime,egftime):
        self.masterarrivaltime = mastertime
        self.egfarrivaltime = egftime
    def set_masteregfgcarc(self,mastergcarc,egfgcarc):
        self.mastergcarc = mastergcarc
        self.egfgcarc = egfgcarc

class Event:
    freqmin = 0
    freqmax = 45
    timebefore = 0
    timeafter = 0
    minmagdiff = 0
    mintimediff = 0
    maxdistdiff = 0
    freqsnr1low = 1
    freqsnr1high = 10
    freqsnr2low = 1
    freqsnr2high = 10
    
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
            self.minmagdiff = 1
            self.mintimediff = 120
            self.maxdistdiff = 50
        elif (self.mag >= 5 and self.mag < 6):
            self.freqmin = 0.1
            self.freqmax = 5
            self.timebefore = 10
            self.timeafter = 50
            self.origintime = 75
            self.minmagdiff = 1
            self.mintimediff = 60
            self.maxdistdiff = 20
        elif (self.mag >= 4 and self.mag < 5):
            self.freqmin = 1
            self.freqmax = 15
            self.timebefore = 5
            self.timeafter = 25
            self.origintime = 60
            self.minmagdiff = 1
            self.mintimediff = 30
            self.maxdistdiff = 10
        elif (self.mag >= 3.5 and self.mag < 4):
            self.freqmin = 1
            self.freqmax = 20
            self.timebefore = 0.5
            self.timeafter = 2.5
            self.origintime = 45
            self.minmagdiff = 0.5
            self.mintimediff = 5
            self.maxdistdiff = 1
        elif (self.mag >= 3 and self.mag < 3.5):
            self.freqmin = 1
            self.freqmax = 35
            self.timebefore = 0
            self.timeafter = 1.0
            self.origintime = 45
            self.minmagdiff = 0.5
            self.mintimediff = 5
            self.maxdistdiff = 0.5
            self.freqsnr1low = 1
            self.freqsnr1high = 10
            self.freqsnr2low = 10
            self.freqsnr2high = 20
        elif (self.mag >= 2 and self.mag < 3):
            self.freqmin = 1
            self.freqmax = 35
            self.timebefore = 0
            self.timeafter = 0.5
            self.origintime = 30
            self.minmagdiff = 0.5
            self.mintimediff = 3
            self.maxdistdiff = 0.2
            self.freqsnr1low = 2
            self.freqsnr1high = 20
            self.freqsnr2low = 20
            self.freqsnr2high = 35
        else:
            self.freqmin = 5
            self.freqmax = 45
            self.timebefore = 0.5
            self.timeafter = 2.5
            self.origintime = 15
            self.minmagdiff = 1
            self.mintimediff = 3
            self.maxdistdiff = 1

    # string representation
    def __str__(self):
        return "id: {}, mag: {}".format(self.id, self.mag)
