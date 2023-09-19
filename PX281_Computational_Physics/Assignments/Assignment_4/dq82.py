fit2 = fitted_ampli[fitted_ampli < 100]
fit2 = fit2[np.where(np.abs(relative_error[:, 0]) < 1)]
amplitudes = np.array(np.unique(np.sort(np.round(fit2)), return_counts=True))
amplitude_r = np.linspace(0, 100, 1000)
peak_index = find_peaks(amplitudes[1], height=25)[0]
width_size = peak_widths(amplitudes[1], peak_index)[0]
plt.figure(figsize=(16, 9))
plt.title("Scatter Plot of Estimated Ampltitude against Onset Relative Error",
          size=18)
plt.grid(True)
plt.scatter(relative_error[:, 0], fitted_ampli)
plt.xlabel("Onset Relative Error (ns)", size=12)
plt.ylabel("Estimated Amplitude", size=12)
plt.figure(figsize=(16, 9))
plt.title("Histogram of Estimated Amplitudes", size=18)
plt.grid(True)
plt.xlabel("Estimated Amplitude", size=12)
plt.ylabel("Frequency Density", size=12)
plt.hist(fit2, bins=100, density=True)
for i in range(0, len(peak_index)):
    plt.plot(amplitude_r,
             norm.pdf(amplitude_r, amplitudes[0][peak_index[i]],
                      width_size[i]),
             label=f"Peak at {amplitudes[0][peak_index[i]]}")

plt.legend(fontsize=12)
plt.figure(figsize=(16, 9))
plt.title("Fitted Peak Width against Fitted Peaks", size=18)
plt.grid(True)
plt.xlabel("Amplitude of Peak", size=12)
plt.ylabel("Peak Width", size=12)
plt.scatter(amplitudes[0, peak_index], width_size, label="")
plt.show()
