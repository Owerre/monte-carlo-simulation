####################################
# Author: S. A. Owerre
# Date modified: 29/11/2021
# Class: Persistent Random Walk 1D
####################################
import numpy as np
from numpy.random import rand


class PersistentRandomWalk1D:
    """monte carlo simulation of one-dimensional persistent random walk.

    inputs:
        (integer) ntrials: number of times to repeat the walk
        (integer) nsteps: number of steps to take in one trial
        (float) p: probability to step in the same direction as the previous step
        1-p: probability to step in the opposite direction to the previous step

    outputs:
        (1d array) x_arr: displacement
        (dictionary) visited_sites: visisted sites and their count
        (1d array) x_avg: average displacement
        (1d array) sigma2: displacement variance
        (integer) count: count the number of distinct sites visited
        (1d array) mean_count: mean count the number of distinct sites visited
    """

    def __init__(self, nsteps, ntrials, p):
        """define parameters of the model."""
        self.nsteps = nsteps
        self.ntrials = ntrials
        self.p = p

    def monte_carlo(self):
        """monte carlo simulation."""
        x_arr = np.zeros((self.ntrials, self.nsteps + 1))
        x2_arr = np.zeros((self.ntrials, self.nsteps + 1))
        visited_sites = {}  # track visited sites and their count
        for i in range(self.ntrials):
            dir = np.zeros(self.nsteps + 1)   # direction of the last step
            x = 0   # initial position
            dir[0] = 1   # initial (previous) direction
            for j in range(self.nsteps):
                # step in the same direction as the previous step
                if rand() <= self.p:
                    if dir[j] == 1:
                        x += 1   # step right
                        dir[j + 1] = 1   # update direction of step
                    else:
                        x -= 1   # step left
                        dir[j + 1] = -1   # update direction of step

                else:  # step in the opposite direction to the prevoius step
                    if dir[j] == 1:
                        x -= 1
                        dir[j + 1] = -1
                    else:
                        x += 1
                        dir[j + 1] = 1

                # update the position array at each n steps
                x_arr[i][j + 1] = x
                x2_arr[i][j + 1] = x * x

            # update visited sites after n steps
            visited_sites[x] = visited_sites.get(x, 0) + 1

        # average over ntrials
        x_avg = np.mean(x_arr, axis=0)
        x2_avg = np.mean(x2_arr, axis=0)
        sigma2 = x2_avg - x_avg * x_avg
        return x_arr, visited_sites, x_avg, sigma2

    def sites_visited(self):
        """count the number of distinct sites visited
        during the course of n steps."""
        x = 0   # initial position
        dir = np.zeros(self.nsteps + 1)   #  direction of the last step
        count = np.zeros(self.nsteps + 1)   #  no. of distinct visited sites
        visited_sites = {}   #  visited sites
        dir[0] = 1   # initial (previous) direction
        count[0] = 1   # initial position counted as 1
        visited_sites[0] = 1   # initial position already visited once
        for i in range(self.nsteps):
            # step in the same direction as the previous step
            if rand() <= self.p:
                if dir[i] == 1:
                    x += 1   # step right
                    dir[i + 1] = 1   # update direction of step
                else:
                    x -= 1   # step left
                    dir[i + 1] = -1   # update direction of step

            else:  # step in the opposite direction to the prevoius step
                if dir[i] == 1:
                    x -= 1
                    dir[i + 1] = -1
                else:
                    x += 1
                    dir[i + 1] = 1

            if x in visited_sites:
                count[i + 1] = count[i]   # the same as the previous count
            else:
                count[i + 1] = 1 + count[i]   # update count by 1

            # update visited sites
            visited_sites[x] = visited_sites.get(x, 0) + 1
        return count

    def average_sites_visited(self):
        """compute the average number of distinct sites visited during
        the course of n steps over n trials or walkers."""
        arr = np.zeros((self.ntrials, self.nsteps + 1))
        for i in range(self.ntrials):
            count = self.sites_visited()
            arr[i, :] = count
        mean_count = np.mean(arr, axis=0)
        return mean_count
