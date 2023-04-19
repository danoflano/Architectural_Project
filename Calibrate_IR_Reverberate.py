#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 15:40:45 2023

@author: dantegarcia
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 11:10:15 2023

@author: dantegarcia
"""

from scipy.io import wavfile
import numpy as np
from scipy import fftpack as ft
import matplotlib.pyplot as plt
global fs

def AverageTimeSignals(s_wav_file,chirp_length):
    global fs
    
    fs,time_sig = wavfile.read(s_wav_file)
    time_sig = np.array(time_sig)
    N = len(time_sig)
    
    datablock = chirp_length*fs
    n = 10 #round(N/datablock)
    time_sigs = np.empty((n,datablock))
    time_sigs[:] = np.nan
    for i in range(0,n):
        time_sigs[i,:] = time_sig[i*datablock:i*datablock+datablock]
    avg_t_sig = np.average(time_sigs, axis=0)
    avg_t_sig = avg_t_sig/max(avg_t_sig)
    N = len(avg_t_sig)
    counts = np.arange(N) # array from 0 to N-1
    time = counts/fs # time array
    freq = np.linspace(0.0, fs/2, N//2) #use // to return int
    freq_sig = ft.fft(avg_t_sig)/max(np.abs(ft.fft(avg_t_sig)))

    return time,freq,avg_t_sig,freq_sig

def Convolve(A,B):
    # convolve 2 time signals
    if len(A) > len(B):
        B = np.append(B,np.zeros(len(A)-len(B)))
    elif len(B) > len(A):
        A = np.append(A,np.zeros(len(B)-len(A)))
    C = ft.ifft(ft.fft(A)*ft.fft(B))
    C = np.int16(np.real(C)/max(np.real(C))*32767)
    return C

def Deconvolve(A,B):
    # deconvolve 2 time signals
    # Order matters here!
    if len(A) > len(B):
        B = np.append(B,np.zeros(len(A)-len(B)))
    elif len(B) > len(A):
        A = np.append(A,np.zeros(len(B)-len(A)))
    C = ft.ifft(ft.fft(A)/ft.fft(B))
    C = np.int16(np.real(C)/max(np.real(C))*32767)
    return C

if __name__ == "__main__":
    
    #%% Find IR of Microphones
    # 3 second chirps - Find IR's
    time,freq,t_sig_chirp1,f_sig_chirp1 =  AverageTimeSignals('exponential_chirp.wav',3)
    time,freq,t_sig23,f_sig23 = AverageTimeSignals('Chamber/Earthworks M23 chirp1.wav',3)
    time,freq,t_sig57,f_sig57 = AverageTimeSignals('Chamber/SM57 chirp1.wav',3)
    time,freq,t_sig251,f_sig251 = AverageTimeSignals('Chamber/251 chirp1.wav',3)
    time,freq,t_sig84_dr,f_sig84_dr = AverageTimeSignals('Chamber/84 drum5 chirp1.wav',3)
    time,freq,t_sig84_de,f_sig84_de = AverageTimeSignals('Chamber/84 desk6 chirp1.wav',3)
    
    speaker_response = f_sig23/f_sig_chirp1
    IR57 = f_sig57/speaker_response/f_sig_chirp1
    IR251 = f_sig251/speaker_response/f_sig_chirp1
    IR84_dr = f_sig84_dr/speaker_response/f_sig_chirp1
    IR84_de = f_sig84_de/speaker_response/f_sig_chirp1
    
    
    
    #%% Obtain IR of BEDROOM using a given mic (251 in this case)
    time,freq,t_sig251_BR,f_sig251_BR = AverageTimeSignals('Bedroom Chirps/WA-251_Bedroom.wav',3)
    # Plot averaged time signal
    if 0 :
        fig,ax1 = plt.subplots()
        ax1.plot(time,t_sig251_BR)
    
    # H, h are impulse responses of bedroom seeing using 251
    H = f_sig251_BR/speaker_response/IR251/f_sig_chirp1 # rec_sig/speaker_IR/mic_IR/inp_chp
    h = ft.ifft(H)
    
    # Plot IR vs time and freq
    time = np.arange(len(h))/fs
    freq = np.linspace(0.0, fs/2, len(h)//2) #use // to return int
    fig,ax1 = plt.subplots(2)
    ax1[0].plot(time,h)
    ax1[1].semilogx(freq,20*np.log10(np.abs(H[:len(h)//2])/max(np.abs(H[:len(h)//2]))))
    ax1[0].set(title='Room Impulse Response vs Time (s)')
    ax1[1].set(title='Room Impulse Response vs Frequency (Hz)',xlim=[1e2,2e4])
    
    #%% Reverberate (Convolve) Audio Track recorded in chamber
    
    fs,moonrocks = wavfile.read('Chamber/talkingheads_moonrocks.wav')
    reverberated_sig = Convolve(h,moonrocks)
    time = np.arange(len(reverberated_sig))/fs
    
    # plot reverberated time signal and save wav file
    if 0 :
        fig, ax1 = plt.subplots()
        ax1.plot(time,reverberated_sig) #h
        ax1.set(title='Reverbed Moon')#,xlim=[0,0.002])
    # wavfile.write('reverbed_raw_moonrocks_BR.wav', fs, reverberated_sig)


    #%% Obtain original signal (Deconvolve) Audio Track
    # moonrocksfft = ft.fft(moonrocks)
    # time = np.arange(len(moonrocks))/fs
    # freq = np.linspace(0.0, fs/2, len(moonrocks)//2) #use // to return int
    
    # deconv_moonrocks = Deconvolve(moonrocksfft,IR84_dr)
    # deconv_moonrocks = Deconvolve(deconv_moonrocks,speaker_response)
    # reverberated_deconv_sig = Convolve(h,deconv_moonrocks)
    # wavfile.write('deconv_moonrocks_BR.wav', fs, deconv_moonrocks)
    # wavfile.write('reverbed_deconv_moonrocks_BR.wav', fs, reverberated_deconv_sig)
    
    
    
    
    
    
    
    
    
    