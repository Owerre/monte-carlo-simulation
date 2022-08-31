#############################
# Author: S. A. Owerre
# Date modified: 02/01/2021
# Class: 1D Random Walk
############################
import math
import numpy as np
from numpy.random import rand


class RandomWalk1D:
    """monte carlo simulation of one-dimensional random walk.

    inputs:
        (integer) ntrials: number of times to repeat the walk
        (integer) nsteps: number of steps to take in one trial
        (float) p: probability to step to the right and 1-p to step left

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
        visited_sites = {}  # map visited sites to count after n steps
        for i in range(self.ntrials):
            x = 0   # initial position
            for j in range(self.nsteps):
                if rand() <= self.p:
                    x += 1
                else:
                    x -= 1

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
        during the course of n steps.
        """
        count = np.zeros(self.nsteps + 1)   # number of distinct visited sites
        visited_sites = {}   # # track visited sites and their count
        x = 0   # initial position
        count[0] = 1   # initial position counted as 1
        visited_sites[0] = 1   # initial position already visited once
        for i in range(self.nsteps):
            if rand() <= self.p:
                x += 1
            else:
                x -= 1

            if x in visited_sites:
                count[i + 1] = count[i]   # the same as the previous count
            else:
                count[i + 1] = 1 + count[i]   # update count by 1

            # update visited sites
            visited_sites[x] = visited_sites.get(x, 0) + 1
        return count

    def average_sites_visited(self):
        """compute the average number of distinct sites visited during
        the course of n steps over n trials or walkers.
        """
        arr = np.zeros((self.ntrials, self.nsteps + 1))
        for i in range(self.ntrials):
            count = self.sites_visited()
            arr[i, :] = count
        mean_count = np.mean(arr, axis=0)
        return mean_count

    def exact_dist(self):
        """the exact probability P(x,n) that the displacement
        of the walker is x after n steps, given by the
        binomial distribution.
        """
        _, hashmap, _, _ = self.monte_carlo()
        x_pos = list(hashmap.keys())
        prob = [0] * len(x_pos)

        for i in range(len(x_pos)):
            k = (self.nsteps + x_pos[i]) / 2
            dem = math.factorial(self.nsteps - k) * math.factorial(k)
            coeff = math.factorial(self.nsteps) / dem
            prob[i] = coeff * (self.p**k) * (1 - self.p) ** (self.nsteps - k)
        return x_pos, prob

    def gaussian_approx(self, const):
        """gaussian approximation of the position distribution
        for large n steps.
        """
        _, hashmap, _, _ = self.monte_carlo()
        x_pos = list(hashmap.keys())

        xbar = self.nsteps * (2 * self.p - 1)   # exact average for large n
        sigma2 = (
            4 * self.nsteps * self.p * (1 - self.p)
        )   # exact average for large n
        gau_prob = [0] * len(x_pos)

        for i in range(len(x_pos)):
            coeff = const / np.sqrt(2 * np.pi * sigma2)
            gau_prob[i] = coeff * np.exp(
                -((x_pos[i] - xbar) ** 2) / (2 * sigma2)
            )
        return x_pos, gau_prob

    def exact_enumeration(self):
        """evaluate averages by exactly enumerating all the possible walks
        for a given n steps at p = 0.5; there are 2**n possible walks.
        note: this should be done for very small n steps
        """
        all_pos = [c for c in range(-self.nsteps, self.nsteps + 1)]
        actual_pos = []   # x positions after n steps (equiprobable events)
        for x in all_pos:
            k = (self.nsteps + x) / 2   # number of steps to the right
            if 2 * k % 2 == 0:
                # the binomial coefficients count how many walks ended at x
                dem = math.factorial(self.nsteps - int(k)) * math.factorial(
                    int(k)
                )
                bino_coeff = math.factorial(self.nsteps) / dem
                for _ in np.arange(bino_coeff):
                    actual_pos.append(x)

        # compute averages
        # 0 for p = 0.5
        x_avg = sum([x for x in actual_pos]) / 2**self.nsteps
        # nsteps for p = 0.5
        x2_avg = sum([x * x for x in actual_pos]) / 2**self.nsteps
        return actual_pos, x_avg, x2_avg


# if __name__ == "__main__":
#     walk = RandomWalk1D(nsteps=2, ntrials=100, p=0.5)
#     print(walk.exact_enumeration())
