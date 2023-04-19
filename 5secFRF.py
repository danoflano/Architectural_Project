#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 16:52:40 2023

@author: dantegarcia
"""
from scipy.io import wavfile
import numpy as np
from scipy import fftpack as ft
import matplotlib.pyplot as plt
import Project_Final as pf
global fs

#%% Obtaining Mic responses from 5 second chirp and PLOTTING

# Obtain Response for each mic
time,freq,t_sig_chirp,f_sig_chirp = pf.AverageTimeSignals('exponential_chirp2.wav',5)
time,freq,t_sig23,f_sig23 = pf.AverageTimeSignals('Chamber/Earthworks M23 chirp2#01.wav',5)
time,freq,t_sig57,f_sig57 = pf.AverageTimeSignals('Chamber/SM57 chirp 2#01.wav',5)
time,freq,t_sig251,f_sig251 = pf.AverageTimeSignals('Chamber/251 chirp2#03.wav',5)
time,freq,t_sig84_dr,f_sig84_dr = pf.AverageTimeSignals('Chamber/84 drum5 chirp2#01.wav',5)
time,freq,t_sig84_de,f_sig84_de = pf.AverageTimeSignals('Chamber/84 desk5 chirp2#01.wav',5)

speaker_response = f_sig23/f_sig_chirp
IR57 = f_sig57/speaker_response/f_sig_chirp
IR251 = f_sig251/speaker_response/f_sig_chirp
IR84_dr = f_sig84_dr/speaker_response/f_sig_chirp
IR84_de = f_sig84_de/speaker_response/f_sig_chirp


if 0:
    plt.close('all')
    N = len(time)
    fig, ax1 = plt.subplots()
    ax1.semilogx(freq, 20*np.log10(IR57[:N//2]),label='SM57')
    ax1.semilogx(freq, 20*np.log10(IR251[:N//2]),label='251')
    ax1.semilogx(freq, 20*np.log10(IR84_dr[:N//2]),label='84-1')
    ax1.semilogx(freq, 20*np.log10(IR84_de[:N//2]),label='84-2')
    ax1.set(title='Frequency Response',xlim=[100,2e4])
    ax1.legend(fontsize=8,loc=4)
    ax1.grid(visible=None, which='major', axis='both')


# Slice where f > 100 and plot (3 second chirps)

idx1 = np.argmax(freq>=100)
idx2 = np.argmax(freq>=2e4)
idx3 = np.argmax(freq>=1e3)
freq_tr = freq[idx1:idx2]
IR57_tr = IR57[idx1:idx2]
IR251_tr = IR251[idx1:idx2]
IR84_dr_tr = IR84_dr[idx1:idx2]
IR84_de_tr = IR84_de[idx1:idx2]


plt.close('all')
fig, ax1 = plt.subplots()
ax1.semilogx(freq_tr, 20*np.log10(IR57_tr/np.mean(np.abs(IR57_tr[:idx3]))),label='SM57')
ax1.semilogx(freq_tr, 20*np.log10(IR251_tr/np.mean(np.abs(IR251_tr[:idx3]))),label='251')
ax1.semilogx(freq_tr, 20*np.log10(IR84_dr_tr/np.mean(np.abs(IR84_dr_tr[:idx3]))),label='84-1')
ax1.semilogx(freq_tr, 20*np.log10(IR84_de_tr/np.mean(np.abs(IR84_de_tr[:idx3]))),label='84-2')
ax1.set(title='Frequency Response',xlim=[100,2e4])
ax1.legend(fontsize=8,loc=4)
ax1.grid(visible=None, which='major', axis='both')