#############################
# Author: S. A. Owerre
# Date modified: 02/01/2021
# Class: 2D Random Walk
#############################
import numpy as np
import random


class RandomWalk2D:
    """monte carlo simulation of two-dimensional random walk

    inputs:
        (integer) nwalkers: number of walkers
        (integer) nsteps: number of steps to take for each walkers

    outputs:
        (1d array) x_arr: displacement in x direction
        (1d array) y_arr: displacement in y direction
        (1d array) sigmax: displacement variance in x direction
        (1d array) sigmay: displacement variance in y direction
        (1d array) r: total displacement variance
        (integer) count: count the number of distinct sites visited
        (1d array) mean_count: mean count the number of distinct sites visited
    """

    def __init__(self, nsteps, nwalkers):
        """define parameters of the model."""
        self.nsteps = nsteps
        self.nwalkers = nwalkers

    def monte_carlo(self):
        """monte carlo simulation."""
        # measured quantities
        x_arr = np.zeros((self.nwalkers, self.nsteps + 1))
        x2_arr = np.zeros((self.nwalkers, self.nsteps + 1))
        y_arr = np.zeros((self.nwalkers, self.nsteps + 1))
        y2_arr = np.zeros((self.nwalkers, self.nsteps + 1))

        for i in range(self.nwalkers):
            x = 0
            y = 0
            for j in range(self.nsteps):
                nearest_neighbors = np.array(
                    [[x + 1.0, y], [x - 1.0, y], [x, y + 1.0], [x, y - 1.0]]
                )

                # choose a random direction and move x and y there
                k = random.randint(0, len(nearest_neighbors) - 1)
                x = nearest_neighbors[k, 0]
                y = nearest_neighbors[k, 1]

                # update the position arrays
                x_arr[i][j + 1] = x
                y_arr[i][j + 1] = y
                x2_arr[i][j + 1] = x**2
                y2_arr[i][j + 1] = y**2

        # average over n walkers
        x_avg = np.mean(x_arr, axis=0)
        y_avg = np.mean(y_arr, axis=0)
        x2_avg = np.mean(x2_arr, axis=0)
        y2_avg = np.mean(y2_arr, axis=0)

        sigma2x = x2_avg - x_avg * x_avg
        sigma2y = y2_avg - y_avg * y_avg
        r2 = sigma2x + sigma2y
        return x_arr, y_arr, sigma2x, sigma2y, r2

    def sites_visited(self):
        """count the number of distinct sites visited
        during the course of n steps."""
        visited_sites = {}   # track visited sites
        count = np.zeros(
            self.nsteps + 1
        )   # track the no. of distinct visited sites
        x, y = 0, 0   # starting position
        visited_sites[(0, 0)] = 1   # starting position already visited once
        count[0] = 1   # starting position counted as 1
        for i in range(self.nsteps):
            nearest_neighbors = np.array(
                [[x + 1.0, y], [x - 1.0, y], [x, y + 1.0], [x, y - 1.0]]
            )
            # choose a random direction and move x and y there
            k = random.randint(0, len(nearest_neighbors) - 1)
            x = nearest_neighbors[k, 0]
            y = nearest_neighbors[k, 1]
            if (x, y) in visited_sites:
                count[i + 1] = count[i]  # same as the previous count
            else:
                count[i + 1] = 1 + count[i]   # update count by 1

            # update visited sites
            visited_sites[(x, y)] = visited_sites.get((x, y), 0) + 1
        return count

    def average_sites_visited(self):
        """compute the average number of distinct sites visited during
        the course of n steps over n walkers.
        """
        arr = np.zeros((self.nwalkers, self.nsteps + 1))
        for i in range(self.nwalkers):
            count = self.sites_visited()
            arr[i, :] = count
        mean_count = np.mean(arr, axis=0)
        return mean_count
