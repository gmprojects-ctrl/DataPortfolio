""" Our code"""

# INSERT YOUR CODE HERE

import numpy as np
import matplotlib.pyplot as plt


class Model:
    """Our Model"""

    def __init__(self, pathtodata, binamount=20):
        # Number of bins
        self.bins = binamount
        # Storing the data in the class obj
        data_obj = np.genfromtxt(pathtodata)
        # Extracting time
        self.time = data_obj[:, 0]
        # Extracting flux
        self.flux = data_obj[:, 1]
        # Extracting flux_error
        self.flux_err = data_obj[:, 2]
        # Creating the title
        self.title = (pathtodata.rsplit("/")[-1]).rsplit("_")[0]
        # Creating a tuple to hold our Period results:
        self.presults = None
        # Presults = (time,s,bestperiod)
        # Used classes, as it creates a unqiue object for each experiment
        # Note that self.VARIABLENAME means that the variable is owned
        # by that classs

    def PERIODIZE(self, period):
        """Periodize the time array"""
        # Periodize time  (t mod p)/ p
        return np.mod(self.time, period) / period

    def BINNED_ARRAYS_OP(self, bin_edges, point_count, p_time):
        """Operating on binned arrays to get A & E"""
        # Get the mean flux
        meanflux = np.mean(self.flux)
        # Set E_STAT to zero
        E_STAT = 0
        # Set A_STAT to zero
        A_STAT = 0
        for i in range(0, self.bins):
            # Get the fluxes in each mean
            # In particular, the code will find the find the corresponding
            # fluxes of the times contained in the bin periodised by our
            # time p
            forw = bin_edges[i + 1]
            back = bin_edges[i]
            curr_arr = self.flux[(back <= p_time) & (p_time < forw)]
            # Some times, we get an empty array
            if curr_arr.size != 0:
                # Find the mean of each bin array
                mean_flux_curr = np.mean(curr_arr)
                # Find the variance of each bin and added to E_STAT
                # E  = \sum_{1}^{M} \sum_{1}^{nj} (f_j - mean_flux_bin )**2
                # Note tha the calculation is vectorised
                E_STAT += np.sum((curr_arr - mean_flux_curr) ** 2)
                # Find the A_STAT for each bin and added
                # A  = \sum_{1}^{M} n_j * (mean_flux_bin - mean_fux)**2
                mean_bin_diff = (mean_flux_curr - meanflux) ** 2
                A_STAT += point_count[i] * mean_bin_diff

            # Return a tuple with (E_Stat,A_Stat)
        return (E_STAT, A_STAT)

    def S_STAT(self, period):
        """Getting S_stat"""
        # Get the len of flux
        Nf = len(self.flux)
        # Periodize the time according to the period
        p_time = self.PERIODIZE(period)
        # Find the points and edges of each bin
        point_count, bin_edges = np.histogram(p_time, bins=self.bins)
        # print(bin_edges) # TESTING PURPOSES
        # Caclulate the A and E statistics
        E_stat, A_stat = self.BINNED_ARRAYS_OP(bin_edges, point_count, p_time)
        # Return S = (N-M)*A / ( (M-1)*E )
        return ((Nf - self.bins) * A_stat) / ((self.bins - 1) * E_stat)

    def FIND_PERIOD(self, search):
        """Finding our Period"""
        # Creates our time search,
        # From 0.1 of a day to 13 days as specified
        time_search = np.linspace(search[0], search[1], search[2])
        # Using Vectorise to apply S_STAT to our time array
        S_array = np.vectorize(self.S_STAT)(time_search)
        # print(S_array)
        # finding the best period
        best_period = time_search[np.argmax(S_array)]
        # Assign it to our p results variable
        self.presults = (time_search, S_array, best_period)

    def PERIODOGRAM(self):
        """Creating a Periodogram"""
        # If no self.presults, it will return with an Value Error
        if self.presults is None:
            raise ValueError(
                """
                Cannot generate periodogram without data
                Please run object.FIND_PERIOD before running PERIODGRAM"""
            )
        # print(self.presults)
        plt.figure(figsize=(16, 9))
        plt.grid(True)
        plt.xlabel("Time (Days)", fontsize=18)
        plt.ylabel("S(t)", fontsize=18)
        plt.title(f"{self.title} Periodogram", fontsize=18)
        plt.plot(self.presults[0], self.presults[1])
        plt.show()

    def PHASE_FOLD(self, bperiod):
        """Defining our Phase Fold"""
        # If no self.presults, it will return with an Value Error
        if self.presults is None:
            raise ValueError(
                """
                Cannot generate periodogram without data
                Please run object.FIND_PERIOD before running PHASE_FOLD"""
            )
        # Time periodised by the best period
        bp_time = self.PERIODIZE(bperiod)
        plt.figure(figsize=(16, 9))
        plt.grid(True)
        plt.xlabel("Time periodized by P", fontsize=18)
        plt.ylabel("Flux", fontsize=18)
        plt.title(f"{self.title} Phasefold at P={bperiod:.4f}", fontsize=18)
        plt.scatter(bp_time, self.flux)
        plt.show()

    def TESS_CURVE(self):
        """Creating Tesss Curve"""
        # Gets the minimum time
        min_t = np.min(self.time)
        # Minuses the minimum time from our time
        time = self.time - min_t
        plt.figure(figsize=(16, 9))
        plt.grid(True)
        # Creates the error plot
        plt.errorbar(
            time,
            self.flux,
            self.flux_err,
            fmt="o",
            linestyle="none",
            markersize=5,
            elinewidth=2.5,
            ecolor="red",
            label="Data +- Error",
        )
        plt.xlabel(f"Time (+{min_t: 0.3E}s)", fontsize=18)
        plt.ylabel("Flux", fontsize=18)
        plt.title(f"{self.title} TESS curve", fontsize=18)
        plt.legend(fontsize=18)
        plt.show()
        # Plots the TESS Curve

    def MAIN(self, true_value=None, ov_p=None, search=None):
        """Main Function"""
        # Concatentates all the previously defined functions into a single
        # main code
        # Defaulting Search
        if search is None:
            search = [0.1, 13, 3000]
        self.FIND_PERIOD(search)
        self.TESS_CURVE()
        self.PERIODOGRAM()
        # If you want to set a custom period
        # you can be setting
        if ov_p is None:
            bp = self.presults[2]
        else:
            if ov_p == -1:
                bp = self.presults[2] * 2
            elif ov_p > 0:
                bp = ov_p
            else:
                raise ValueError("Period can't be negative")
        self.PHASE_FOLD(bp)
        print(f"Best Period : {bp}")
        if true_value is not None:
            diff = np.abs(bp - true_value) / true_value
            print(f"Relative Error : {100*diff}%")


# Please HASH OUT any Testing code here
X = Model("./data/WASP-7_TOI-2197_FLUX.dat", 30)
# print(X.S_STAT(1.6284246))
# print(X.FIND_PERIOD()[1])
# X.FIND_PERIOD()
# X.PERIODOGRAM()
# X.PHASE_FOLD()
# plt.plot(graph[0],graph[1][0])
# plt.xlim(3,4)
# X.TESSCurve()
# print(__name__)
X.MAIN()
