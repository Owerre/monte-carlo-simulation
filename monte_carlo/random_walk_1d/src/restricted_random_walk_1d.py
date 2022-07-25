####################################
# Author: S. A. Owerre
# Date modified: 05/12/2021
# Class: 1D Restricted Random Walk
####################################
import numpy as np
from numpy.random import rand

class RestrictedRandomWalk1D:
    """monte carlo simulation of one-dimensional restricted random walk.

    inputs:
        (integer) ntrials: number of times to repeat the walk
        (integer) nsteps: number of steps to take in one trial
        (float) p: probability to step to the right and 1-p to step to the left
        (integer) L: lattice size
    
    outputs:
        (1d array) x_arr: displacement
        (dictionary) visited_sites: visisted sites and their count
        (1d array) x_avg: average displacement
        (1d array) sigma2: displacement variance
        (integer) count: count the number of distinct sites visited
        (1d array) mean_count: mean count the number of distinct sites visited
    """
    def __init__(self, nsteps, ntrials, p, L):
        """define parameters of the model."""
        self.nsteps = nsteps
        self.ntrials = ntrials
        self.p = p
        self.L = L

    def step_count_b4_trap(self):
        """count the number of steps before being trapped at x = 0 and x = L."""
        lat = [t for t in range(self.L+1)] # lattice sites
        traj = [] # trajectory of the walker
        x =  lat[self.L//2] # starting position
        count = 0 # count total number of steps before being trapped
        for _ in range(self.nsteps):
            # walk terminates at the trap sites
            if x == lat[0] or x == lat[self.L]: 
                break
            
            # random walk
            if rand() <= self.p:
                x += 1
            else:
                x -= 1
            
            # count steps and track trajectory
            count += 1
            traj.append(x)
        return traj, count

    def average_nsteps_trap(self):
        """mean number of steps for the walker to be trapped;
        this is called the mean first passage time.
        """
        step_count = np.zeros(self.ntrials)
        for i in range(len(step_count)):
            _, count = self.step_count_b4_trap()
            step_count[i] = count
        return np.mean(step_count)

    def step_count_b4_trap0(self):
        """count when walker gets trapped at x = 0."""
        lat = [t for t in range(self.L+1)] # lattice sites
        x =  lat[self.L//2] # starting position
        count = 0 # count the number of times the walker hits x = 0
        for _ in range(self.nsteps):
            # walk terminates at the trap sites
            if x == lat[0]: 
                count += 1 
                break
            elif x == lat[self.L]:
                break 
            
            # random walk
            if rand() <= self.p:
                x += 1
            else:
                x -= 1
        return count

    def prob_trap(self):
        """probability of the walker being trapped at x = 0."""
        count = 0
        for _ in range(self.ntrials):
            count += self.step_count_b4_trap0()
        return count/self.ntrials
    
    def reflecting_boundaries(self):
        """monte carlo simulation of reflected boundary random walk."""
        x_arr = np.zeros((self.ntrials, self.nsteps+1))
        x2_arr = np.zeros((self.ntrials, self.nsteps+1))
        visited_sites = {}  # track visited sites and their count
        lat = [t for t in range(-self.L, self.L+1)] # lattice sites
        for i in range(self.ntrials):
            x =  0 # initial position
            for j in range(self.nsteps):
                if rand() <= self.p:
                    if x == lat[2*self.L]:  # right reflection site
                        x -= 1 
                    else:
                        x += 1
                else:
                    if x == lat[0]: # left reflection site
                        x += 1
                    else:
                        x -= 1
                
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
        count = np.zeros(self.nsteps + 1) # number of distinct visited sites 
        visited_sites = {} # # track visited sites and their count
        lat = [x for x in range(-self.L, self.L+1)] # lattice sites

        x = 0 # initial position
        count[0] = 1 # initial position counted as 1
        visited_sites[0] = 1 # initial position already visited once
        
        for i in range(self.nsteps):
            if rand() <= self.p:
                    if x == lat[2*self.L]:  # right reflection site
                        x -= 1 
                    else:
                        x += 1
            else:
                if x == lat[0]: # left reflection site
                    x += 1
                else:
                    x -= 1

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