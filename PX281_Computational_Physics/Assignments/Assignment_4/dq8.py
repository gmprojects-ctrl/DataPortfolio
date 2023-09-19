# YOUR CODE HERE
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.signal import peak_widths
from scipy.stats import norm
import matplotlib.pyplot as plt


def simple_pulse(time, onset, amplitude, risetime, decaytime):
    pulse = np.exp(-(time - onset) / risetime) - np.exp(
        -(time - onset) / decaytime)
    pulse[np.where(time < onset)] = 0.0
    return -amplitude * pulse


def oscillation(length, amplitude, freq, decaytime):
    time = np.arange(length)
    return amplitude * np.sin(2 * np.pi * freq * time) * np.exp(
        -time / decaytime)


num_pulse = 2000
len_pulse = 1000
time = np.linspace(0, len_pulse, len_pulse)
lower_half_A_ch = (11, 16, 31, 82)
lower_half_A = np.array([
    np.random.choice(lower_half_A_ch) for _ in range(0, int(0.45 * num_pulse))
])
upper_half_A = np.array(
    [np.random.uniform(0, 100) for _ in range(0, int(0.45 * num_pulse))])
bottom_A = np.array(
    [np.random.uniform(1, 20) for _ in range(0, int(0.10 * num_pulse))])
rise_time = 6
onset_time = 250
decay_time = 200


def gen_pul(i):
    noise = np.sqrt(i) * np.random.normal(size=len_pulse)
    return simple_pulse(time, onset_time, i, rise_time, decay_time) + noise


first_45 = np.array(
    [gen_pul(upper_half_A[i]) for i in range(0, int(0.45 * num_pulse))])
second_45 = np.array(
    [gen_pul(lower_half_A[i]) for i in range(0, int(0.45 * num_pulse))])
third_10 = np.array([
    oscillation(len_pulse, bottom_A[i], 1 / 80, 500)
    for i in range(0, int(0.10 * num_pulse))
])
data = np.vstack((first_45, second_45, third_10))
ampli = np.concatenate((upper_half_A, lower_half_A, bottom_A))
fitted_ampli = np.array([])
relative_error = np.empty((0, 4))
for i in range(0, data.shape[0]):
    try:
        guess = curve_fit(simple_pulse,
                          time,
                          data[i],
                          p0=[onset_time, 50, rise_time, decay_time],
                          bounds=((0, 0, 0, 0),
                                  (onset_time * 2, 100 * 2, rise_time * 2,
                                   decay_time * 2)))
        fitted_ampli = np.append(fitted_ampli, guess[0][1])
        actual = np.array([onset_time, ampli[i], rise_time, decay_time])
        relative_error = np.vstack((relative_error, guess[0] - actual))
    except ValueError or RuntimeError:
        pass
