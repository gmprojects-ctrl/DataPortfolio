'''
SiPM single event modelling
(a) Make a pulse, always 50 ns after the start of the event.
(b) Randomly decide with the dark count rate how many
    additional pulses should be created.
(c) for more than zero, add additional pulses with random amplitude
    in discrete units of scale and random position in time
    according to exponential distribution and dark count rate.
(d) Analyse event: first filter with matched filter and find peaks.
(e) Use results on peak positions and number to fit all peaks,
    especially if there are more than one.
(f) Draw data and fit.
'''

import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit


def matchedFilter(data, template):
    return np.correlate(data, template, mode='full')


def pulseSequence(t, *pos):
    total = np.zeros_like(t)
    for idx in range(0, len(pos) - 3, 4):  # superposition of pulses
        total += pulse(t, pos[idx], pos[idx + 1], pos[idx + 2], pos[idx + 3])
    return total


def pulse(time, amplitude, start, rt, dt):
    singlepulse = np.exp(-(time - start) / rt) - np.exp(-(time - start) / dt)
    singlepulse[np.where(
        time < start)] = 0.0  # not defined before onset time, set 0
    return -1 * amplitude * singlepulse


def makeTemplate(rt, dt):
    timevalues = np.linspace(0, 50, 101)  # 0.05 mus, 0.5 unit step size
    scale = 1.0  # some scale factor giving reasonable values
    onset = 0.0  # pulse start [ns]
    dummy = pulse(timevalues, scale, onset, rt, dt)
    template = dummy / np.trapz(dummy, x=timevalues)  # normalized
    return timevalues, template


def dataProduction(time, cfg):
    amp = cfg[0]  # some scale factor giving reasonable values
    start = cfg[1]  # pulse start [ns]
    rtime = cfg[2]  # realistic rise time
    dtime = cfg[3]  # realistic decay time
    noiselevel = cfg[4]  # noise level scale
    dcr = cfg[5]  # [1/ns] => 2 MHz

    framestop = time[-1]  # final time bin

    Npulses = np.random.poisson(framestop * dcr)
    print('n pulses: ', Npulses)

    pp = pulse(time, amp, start, rtime, dtime)
    noisy = np.random.normal(pp, scale=noiselevel)
    frame = noisy  # first triggered pulse at onset
    for _ in range(Npulses):  # additional pulses
        npe = np.random.poisson(1.0)  # n photo electrons given DCR
        print('npe: ', npe)
        pretrigger = start
        triggertime = random.expovariate(dcr)  # rate parameter
        start = pretrigger + triggertime
        if start > framestop - 300:
            break
        if npe > 0:
            print('next onset: ', start)
            pp = pulse(time, npe * amp, start, rtime, dtime)
            frame += pp
    return frame


def analysis(tvalues, data, cfg):
    scale = cfg[0]  # some scale factor giving reasonable values
    rtime = cfg[2]  # realistic rise time
    dtime = cfg[3]  # realistic decay time

    # prepare the analysis with the matched filter - get a template pulse
    time, tplt = makeTemplate(rtime, dtime)
    time -= time[-1]
    filtered = matchedFilter(result, tplt)  # filter
    responsetime = np.concatenate((time[:-1], tvalues), axis=None)

    # search the filtered result for peaks in response
    peakfilt, _ = find_peaks(filtered, height=6.0, distance=6.0)
    print('in filtered: ', responsetime[peakfilt])

    # fit the pulse, try similar initial values to construction
    # and the identified peak number and positions
    if responsetime[peakfilt].size == 0:
        return None  # failed peak finding
    init = []

    # construct the initial parameter array
    for val in peakfilt:
        init.append([scale, val, rtime, dtime])
    try:
        fitParams, _ = curve_fit(pulseSequence,
                                 tvalues,
                                 result,
                                 p0=np.array(init))
        print(fitParams)
    except (RuntimeError, ValueError):
        fitParams = None  # failed fit

    return fitParams


# Start main script
# make a pulse, consider times in nano seconds [ns]
timevalues = np.linspace(0, 1500, 3001)  # 1.5 mus, 0.5 unit step size
scale = 10.0  # some scale factor giving reasonable values
risetime = 2.0  # realistic rise time
decaytime = 150.0  # realistic decay time
onset = 50.0  # pulse start [ns]
nlevel = 0.2  # noise level scale
darkCountRate = 0.002  # [1/ns] => 2 MHz
config = [scale, onset, risetime, decaytime, nlevel, darkCountRate]

# Data production first
result = dataProduction(timevalues, config)

# then analyse the event
bestfit = analysis(timevalues, result, config)

# finally plotting
plt.plot(timevalues, result, 'r-')
if bestfit is not None:
    plt.plot(timevalues, pulseSequence(timevalues, bestfit), 'b-')
plt.title('SiPM pulse', size=12)
plt.xlabel('Time [ns]', size=12)
plt.ylabel('Amplitude [mV]', size=12)

plt.show()
