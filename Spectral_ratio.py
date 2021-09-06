#!/home/meichen/anaconda3/bin/python

import numpy as np
import pickle
import os
import obspy
from Event_Trace import Trace
from Event_Trace import Event

def main():
    with open("eventpairs.pkl","rb") as f:
        eventpairs = pickle.load(f)
    f.close()
    with open("stationdic_S.pkl","rb") as f:
        stationdic_S = pickle.load(f)
    f.close()
    with open("stationdic_P.pkl","rb") as f:
        stationdic_P = pickle.load(f)
    f.close()

    for master, egfs in eventpairs.items():
        for egf in egfs:
            print(master,egf)

    for 

main()
