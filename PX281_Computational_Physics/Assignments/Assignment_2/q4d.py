'''
Use moving average filter for a simple but noisy pulse

This example contains programming faults!
'''

import numpy as np
import matplotlib.pyplot as plt

def moving_average(inputdata, filterorder):
    ''' implement a moving average filter of odd order
        inputdata: to be filtered, NumPy array
        filterorder: integer number as filter order,
                     must be made odd if even, i.e. even+1.
    '''
    if filterorder % 2 == 0:
        filterorder+=1
    response = 1 / filterorder  # a float normalization factor
    output = []

    # padding the input for averaging data borders
    leftextension = np.flip(inputdata[:filterorder // 2], 0)
    rightextension = np.flip(inputdata[-filterorder // 2 + 1:], 0)
    padded = np.concatenate([leftextension, inputdata, rightextension])

    output.append(np.sum(padded[:filterorder]) * response) # first average
    n = filterorder // 2 + 1
    for idx in range(n, len(inputdata) + n - 2):
        term = output[-1]  # previous
        term += padded[idx + filterorder // 2] * response
        term -= padded[idx - filterorder // 2 - 1] * response
        output.append(term)
    return np.array(output)

def simple_pulse(length, amplitude, risetime, decaytime):
    '''Simple Pulse function'''
    time = np.linspace(0, length, length + 1)
    onset = 0.3 * length # start 30% into the length
    pulse = np.exp(-(time - onset)/risetime) - np.exp(-(time - onset)/decaytime)
    pulse[np.where(time < onset)] = 0.0 # not defined before onset time, set 0
    return -amplitude * pulse

# make the data
amp = 10.0
risetime_ = 3.0
decaytime_ = 300.0
data = simple_pulse(2000, amp, risetime_, decaytime_)
noisy = data + 0.4 * amp * np.random.normal(size = len(data))

# filter and plot
smooth = moving_average(noisy, 12)
plt.plot(noisy)
plt.plot(smooth, 'r-')
plt.title('pulse with noise and smoothing')
plt.ylabel('amplitude', fontsize = 14)
plt.xlabel('samples', fontsize = 14)
plt.show()
