####################################
# Author: S. A. Owerre
# Date modified: 4/12/2021
# Class: True Self-Avoiding Walk 1D
####################################
import numpy as np
from numpy.random import rand

class TrueSelfAvoidingWalk1D:
    """monte carlo simulation of one-dimensional true self-avoiding walk.

    inputs:
        (integer) ntrials: number of times to repeat the walk
        (integer) nsteps: number of steps to take in one trial
        (float) g: strength with which the walk avoids itself
    
    outputs:
        (1d array) x_arr: displacement
        (dictionary) visited_sites: visisted sites and their count
        (1d array) x_avg: average displacement
        (1d array) sigma2: displacement variance
        (integer) count: count the number of distinct sites visited
        (1d array) mean_count: mean count the number of distinct sites visited
    """
    def __init__(self, nsteps, ntrials, g):
        """define parameters of the model."""
        self.nsteps = nsteps
        self.ntrials = ntrials
        self.g = g

    def monte_carlo(self):
        """monte carlo simulation."""
        x_arr = np.zeros((self.ntrials, self.nsteps+1))
        x2_arr = np.zeros((self.ntrials, self.nsteps+1))
        visited_sites = {}  # map visited sites to count after n steps
        for i in range(self.ntrials):
            nv = np.zeros(2*self.nsteps+1) # track the number of visits to x
            x = 0 # initial position
            nv[0] = 1 # initial position counted as 1
            for j in range(self.nsteps):
                dem = ( np.exp(-self.g*nv[x+1]) + np.exp(-self.g*nv[x-1]) )
                p = np.exp(-self.g*nv[x+1])/ dem # probability to jump to x+1
                if rand() <= p: 
                    x += 1 # step right
                    nv[x] += 1 # update the number of visits to x
                else:
                    x -= 1 # step left
                    nv[x] += 1 # update the number of visits to x

                # update the position array at each n steps
                x_arr[i][j+1] = x  
                x2_arr[i][j+1] = x*x

            # update visited sites after n steps
            visited_sites[x] = visited_sites.get(x,0) + 1

        # average over ntrials
        x_avg = np.mean(x_arr, axis =0) 
        x2_avg = np.mean(x2_arr, axis =0) 
        sigma2 = x2_avg - x_avg*x_avg 
        return x_arr, visited_sites, x_avg, sigma2

    def sites_visited(self):
        """count the number of distinct sites visited 
        during the course of n steps.
        """        
        nv = np.zeros(2*self.nsteps + 1) # track number of visits to x
        count = np.zeros(self.nsteps + 1) # track no. of distinct visited sites 
        visited_sites = {} # track visited sites

        x = 0 # initial position
        nv[0] = 1 # initial position already visited once
        count[0] = 1 # initial position counted as 1
        visited_sites[0] = 1 # initial position already visited once

        for i in range(self.nsteps):
            dem = ( np.exp(-self.g*nv[x+1]) + np.exp(-self.g*nv[x-1]) )
            p = np.exp(-self.g*nv[x+1])/ dem # probability to jump to x+1
            if rand() <= p: 
                x += 1 # step right
                nv[x] += 1 # update the number of visits to x
            else:
                x -= 1 # step left
                nv[x] += 1 # update the number of visits to x

            if x in visited_sites:
                count[i+1] = count[i] # the same as the previous count
            else:
                count[i+1] = 1 + count[i] # update count by 1
            
            # update visited sites
            visited_sites[x] = visited_sites.get(x,0) + 1
        return count

    def average_sites_visited(self):
        """compute the average number of distinct sites visited during 
        the course of n steps over n trials or walkers.
        """
        arr = np.zeros((self.ntrials, self.nsteps +1))
        for i in range(self.ntrials):
            count = self.sites_visited()
            arr[i,:] = count
        mean_count = np.mean(arr, axis = 0)
        return mean_count