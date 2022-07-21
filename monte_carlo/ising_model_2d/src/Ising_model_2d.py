####################################
# Author: S. A. Owerre
# Date modified: 01/01/2021
# Class: Ising model
####################################
import random
import numpy as np
from numpy.random import rand

class Ising:
    """
    monte carlo simulation of 2D Ising model
    """
    def __init__(self, L, nwarmup, nsteps):
        """
        define default parameters
        """
        self.L = L # lattice size
        self.N = L*L # number of spins
        self.nwarmup = nwarmup # number of warm up step
        self.nsteps = nsteps # number of mc step

    def initialize(self):
        """
        random initialization of spin on the 2D square lattice 
        """
        state = [[0]*self.L for _ in range(self.L)]
        for i in range(self.L):
            for j in range(self.L):
                if rand() < 0.5:
                    state[i][j] = -1
                else:
                    state[i][j] = 1
        return state

    def neighbor_pos(self):
        """
        set up array for nearest neigbour positions 
        with periodic boundary condition
        """
        np_ = [0]*self.L
        nm = [0]*self.L
        for i in range(self.L):
            np_[i] = i+1
            nm[i] = i-1
            if i == self.L-1:
                np_[i] = 0
            if i == 0:
                nm[i] = self.L-1
        return np_, nm

    def precom_expo(self, T):
        """
        precompute energy change when spin is flipped
        """
        res = [0]*17
        for de in range(-8,9,4):
            res[de+8] = np.exp(-de/T)
        return res

    def energy(self, spinconf):
        """
        total energy of a given spin
        configuration
        """
        ene = 0
        np_, nm = self.neighbor_pos()
        for i in range(self.L):
            for j in range(self.L):
                ene += -1*spinconf[i][j]*(spinconf[i][np_[j]] +\
                     spinconf[i][nm[j]] +spinconf[np_[i]][j] + \
                          spinconf[nm[i]][j] )
        return ene/2 # to compensate for over-counting
        
    def magnetization(self, spinconf):
        """
        total magnetization of a given spin
        configuration
        """
        mag = 0
        for i in range(self.L):
            for j in range(self.L):
                mag += spinconf[i][j]
        return mag
        
    def metropolis(self, spinconf, pw):
        """
        metropolis algorithm
        """
        ene = self.energy(spinconf) # initial energy state
        mag = self.magnetization(spinconf) # initial magnetization
        np_, nm = self.neighbor_pos() # periodic boundary condition

        for _ in range(self.N):
            ix = random.randint(0,self.L-1) # select random position
            iy = random.randint(0,self.L-1) # select random position

            # compute energy change
            de = 2*spinconf[ix][iy]*(spinconf[ix][np_[iy]] + spinconf[ix][nm[iy]] +\
                    spinconf[np_[ix]][iy] + spinconf[nm[ix]][iy] )
            if de <= 0 or rand() < pw[de + 8]:
                spinconf[ix][iy] *= -1  # flip spin and retain the configuration
                ene += de # update energy
                mag += 2*spinconf[ix][iy] # update magnetization
        return ene, mag

    def mcsweeps(self, spinconf, pw):
        """
        perform monte-carlo sweeps
        """
        avg = np.zeros(6) # initialize averages to zero
        for _ in range(self.nwarmup):   # equilibrate by warm up
            self.metropolis(spinconf, pw)

        for _ in range(self.nsteps):
            ene, mag = self.metropolis(spinconf, pw)
            avg[0] += ene 
            avg[1] += ene*ene
            avg[2] += mag
            avg[3] += mag*mag
            avg[4] += np.sqrt(mag*mag)
            avg[5] += mag*mag*mag*mag
        return avg