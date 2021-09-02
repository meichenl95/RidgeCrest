#!/home/meichen/anaconda3/bin/python3

import numpy as np
import obspy
import os
import subprocess

dataset_path = "/home/meichen/work1/RidgeCrest"

def mseed2sac(eventid,evla,evlo,evdp):
    isExist = os.path.exists('{}/{}'.format(dataset_path,eventid))
    if not isExist:
        print("Folder {}/{} not found".format(dataset_path,eventid))
        exit(-1);

    os.chdir("{}/{}".format(dataset_path,eventid))
    if os.path.isfile('{}.station.txt'.format(eventid)):
        subprocess.call(['rm *.SAC'],shell=True)
    subprocess.call(['echo "" >> {}.iris.txt'.format(eventid)],shell=True)
    subprocess.call(['echo "" >> {}.scedc.txt'.format(eventid)],shell=True)
    subprocess.call(['cat {}.iris.txt {}.scedc.txt {}.ncedc.txt > {}.station.txt'.format(eventid,eventid,eventid,eventid)],shell=True)
    subprocess.call(['mseed2sac {}/{}/*.?H[NEZ].*ms -m {}.station.txt -E "2019,,/{}/{}/{}"'.format(dataset_path,eventid,eventid,evla,evlo,evdp)],shell=True)

def main():
    catalog = np.genfromtxt("catalog_test.txt",skip_header=28,usecols=(0,4,5,6),)
    eventid = catalog[:,0]
    evlo = catalog[:,1]
    evla = catalog[:,2]
    evdp = catalog[:,3]
    
    for i in np.arange(len(eventid)):
        mseed2sac(int(eventid[i]),evla[i],evlo[i],evdp[i])

main()
