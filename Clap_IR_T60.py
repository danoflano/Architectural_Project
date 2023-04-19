#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 13:45:28 2023

@author: dantegarcia
"""
from scipy.io import wavfile
import numpy as np
from scipy import fftpack as ft
import matplotlib.pyplot as plt
global fs

def ImpulseToT60(time_sig,t1,t2,slope,adj_Lp,title,b_wn,xmax,ymin):
    # adj_Lp - Amplitude offset [dB]
    # b_wn - boolean - 0=IIR, 1=INM
    global fs
    i1,i2 = round(t1*fs),round(t2*fs)
    time_sig = time_sig[i1:i2]

    dt = 1/fs
    N = len(time_sig)
    counts = np.arange(N) # array from 0 to N-1
    time = counts/fs # time array
    freq = np.linspace(0.0, fs/2, N//2) #use // to return int
    freq_sig = ft.fft(time_sig)
    # plot time-series and frequency response of unfiltered signal
    if 1:
        fig, ax1 = plt.subplots(2)
        ax1[0].plot(time,time_sig/max(time_sig))
        ax1[0].set(title=title+' - IR vs time(s)')#,xlim=[0,0.002])
        ax1[1].semilogx(freq, 2.0/N * np.abs(freq_sig[:N//2])/max(2.0/N * np.abs(freq_sig[:N//2])))
        ax1[1].set(title=title+' - IR vs frequency (Hz)',xlim=[0,1000])
        plt.tight_layout()
        plt.show()
        fig.savefig('_'.join(title.split())+'_IR.png', format='png', dpi=500)
    xfit = np.arange(0,2,0.01)
    yfit = xfit*slope
    Lp = np.cumsum(np.flip(time_sig **2))*dt
    Lp = np.flip(10*np.log10(Lp)) - max(np.flip(10*np.log10(Lp))) + adj_Lp
    if b_wn:
        Lp = 20*np.log10(time_sig) - max(np.flip(20*np.log10(time_sig))) + adj_Lp
    fig, ax1 = plt.subplots()
    ax1.plot(time,Lp,label='Data')
    ax1.plot(xfit,yfit,label='Linear Fit')
    ax1.hlines(y=-60, xmin=0, xmax=2, color='r',label='-60 dB')
    ax1.grid(visible=None, which='major', axis='both')
    ax1.legend(fontsize=8,loc=1)

    T60 = -60/slope
    print(T60)  # 1.54 s
    ax1.set(xlim=[0,xmax],ylim=[ymin,0],title=title+', T60 = '+str(round(T60,3))+' s',ylabel='Magintude [dB]',xlabel='time [s]')
    #xlim=[0,1],ylim=[-61,0]
    fig.savefig('_'.join(title.split())+'_T60_fit.png', format='png', dpi=500)
    return T60


#%% Untreated
fs,ht = wavfile.read('Clap_Untreated_Middle_251.wav') 
time = np.arange(len(ht))/fs
plt.close('all')
# fig,ax1 = plt.subplots()
# ax1.plot(time,ht)
ImpulseToT60(ht,0,0.11,-91,-4.5,'Studio without Treatment, WA-251',0,1,-61) # 0.214 s start = 0.0178

#%% Treated
fs,ht = wavfile.read('Clap_Treated_Middle_251.wav') 
time = np.arange(len(ht))/fs
# plt.close('all')
# fig,ax1 = plt.subplots()
# ax1.plot(time,ht)
ImpulseToT60(ht,0,0.08,-250,-11.5,'Studio with Treatment, WA-251',0,0.5,-61) # 0.214 s start = 0.0178
# 0 to 0.08


#%% WA-84 FRF
if 0:
    import csv
    
    # opening the CSV file
    with open('WA84_FRF.csv', mode ='r')as file:
       
        # reading the CSV file
        csvFile = csv.reader(file)
         
        freq,mag = [],[]
        for lines in csvFile:
            print(lines)
            freq.append(float(lines[0]))
            mag.append(float(lines[1]))
    
    fig, ax1 = plt.subplots()
    ax1.semilogx(freq,mag)
    ax1.set(xlim=[20,20000],ylim=[-25,25])
    
    






