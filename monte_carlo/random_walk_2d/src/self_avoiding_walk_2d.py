####################################
# Author: S. A. Owerre
# Date modified: 02/01/2021
# Class: 2D Self-Avoiding Walk
####################################
import numpy as np
import random

class SelfAvoidingWalk2D:
    """
    monte carlo simulation of two-dimensional self-avoiding random walk;
    the constraint is that no lattice site can be visited more
    than once in each walk.

    inputs:
        (integer) nwalkers: number of walkers
        (integer) nsteps: number of steps to take for each walkers
    """
    def __init__(self, nsteps, nwalkers):
        """
        define parameters
        """
        self.nsteps = nsteps
        self.nwalkers = nwalkers

    def saw(self):
        """
        self-avoiding walk
        """
        x, y = 0,0 # starting position
        visited = np.zeros((self.nsteps+1, 2)) # track visited sites
        x_arr =  np.zeros(self.nsteps+1) # accumulate x positions
        y_arr =  np.zeros(self.nsteps+1) # accumulate y positions
        for i in range(self.nsteps):
            nearest_neighbors = np.array([
                    [x+1.0, y], 
                    [x-1.0, y],
                    [x, y+1.0],
                    [x, y-1.0]
                ])
            # extract neighbors that have not been visited
            unvisited_neighors = set([tuple(t) for t in nearest_neighbors])
            visited_sites = set([tuple(t) for t in visited])
            dest = np.array([d for d in unvisited_neighors if d not in visited_sites])
            if len(dest) != 0:
                k = random.randint(0,len(dest)-1)
                x = dest[k,0] # move x
                y = dest[k,1] # move y
            visited[i+1,:] = [x,y] # update trajectory already visited
            x_arr[i+1] = x # update the position arrays
            y_arr[i+1] = y
        return visited, x_arr, y_arr, 

    def monte_carlo(self):
        """
        monte carlo simulation for computing mean-squared displacement
        """
        # measured quantities 
        x_arr =  np.zeros((self.nwalkers, self.nsteps+1))
        x2_arr =  np.zeros((self.nwalkers, self.nsteps+1))
        y_arr =  np.zeros((self.nwalkers, self.nsteps+1))
        y2_arr =  np.zeros((self.nwalkers, self.nsteps+1))

        for i in range(self.nwalkers):
            _, x, y = self.saw()
            x_arr[i,:] = x
            y_arr[i,:] = y
            x2_arr[i,:] = x*x
            y2_arr[i,:] = y*y

        # average over n walkers
        x_avg = np.mean(x_arr, axis =0)
        y_avg = np.mean(y_arr, axis =0)
        x2_avg = np.mean(x2_arr, axis =0)
        y2_avg = np.mean(y2_arr, axis =0) 

        sigma2x = x2_avg - x_avg*x_avg 
        sigma2y = y2_avg - y_avg*y_avg 
        r2 = sigma2x + sigma2y
        return sigma2x, sigma2y, r2