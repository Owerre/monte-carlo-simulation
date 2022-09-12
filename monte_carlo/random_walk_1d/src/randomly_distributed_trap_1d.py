#######################################
# Author: S. A. Owerre
# Date modified: 25/03/2022
# Class: 1D Randomly Distributed Traps
#######################################
import numpy as np
from numpy.random import rand
import random


class RandomlyDistributedTrap:
    """monte carlo simulation of one-dimensional random walk
    on lattices with randomly distributed traps.

    inputs:
        (integer) ntrials: number of times to repeat the walk
        (integer) nsteps: number of steps to take in one trial
        (float) p: probability to step to the right and 1-p to step to the left
        (float) rho: trap concentration - probability that a site is a trap site
        (integer) L: lattice size

    outputs:
        (1d array) sites: lattice sites
        (1d array) site_label: site label representing trap or non-trap sites
        (1d array) traj: walker's trajectory
        (integer) count: number of steps before being trapped
        (1d array) proba_trap_config: probability configuration
        (1d array) mean_proba: mean survival probability
    """

    def __init__(self, nsteps, ntrials, p, rho, L):
        """define parameters of the model."""
        self.nsteps = nsteps
        self.ntrials = ntrials
        self.p = p
        self.L = L
        self.rho = rho

    def lattice(self, dtype=None):
        """lattice with randomly distributed sites and site label:
        traps = 0 and nontraps = 1.
        """
        sites = np.zeros(self.L, dtype=dtype)   # lattice sites
        site_label = np.zeros(self.L, dtype=dtype)   # site labels
        for i in range(len(sites)):
            if rand() <= self.rho:
                site_label[i] = 0   # trap label
                sites[i] = i   # trap sites
            else:
                site_label[i] = 1   # nontrap label
                sites[i] = i   # nontrap sites
        return sites, site_label

    def step_count_b4_trap(self):
        """count the number of steps before being trapped
        at site_label[x] = 0.
        """
        traj = []   # trajectory of the walker
        sites, site_label = self.lattice(dtype=int)
        x = random.choice(sites)   # random starting position
        count = 0   # count total number of steps before being trapped
        for _ in range(self.nsteps):
            # walk terminates at the trap site
            if site_label[x] == 0:
                break

            # random walk
            if rand() <= self.p:
                x += 1
            else:
                x -= 1

            # apply periodic boundary condition
            if x > sites[-1]:
                x = sites[0]
            elif x < sites[0]:
                x = sites[-1]

            # count steps and track trajectory
            count += 1
            traj.append(x)
        return traj, count

    def average_nsteps_trap(self):
        """mean number of steps before the walker is trapped;
        this is called the mean survival time or the mean first passage time.
        """
        step_count = np.zeros(self.ntrials)
        for i in range(len(step_count)):
            _, count = self.step_count_b4_trap()
            step_count[i] = count
        return np.mean(step_count)

    def neighbor_pos(self, dtype=None):
        """set up array for nearest neigbour positions
        with periodic boundary condition.
        """
        rp = np.zeros(self.L, dtype=dtype)
        lp = np.zeros(self.L, dtype=dtype)
        for i in range(self.L):
            rp[i] = i + 1
            lp[i] = i - 1
            if i == self.L - 1:
                rp[i] = 0
            elif i == 0:
                lp[i] = self.L - 1
        return rp, lp

    def exact_enumeration_proba(self):
        """survival probability per trap configuration
        using exact enumeration.
        """
        rp, lp = self.neighbor_pos(dtype=int)   # periodic boundary condition
        _, site_label = self.lattice()   # random trap configuration
        proba_trap_config = np.zeros(self.nsteps + 1)   # survival probability
        proba_trap_config[0] = 1   # for initial trap configuration
        w = np.zeros(self.L)   # number of walkers
        for i in range(self.nsteps):
            for j in range(self.L):
                if site_label[j] == 0:
                    w[j] = 0   # traps do not move
                else:
                    w[j] = (site_label[rp[j]] + site_label[lp[j]]) / 2
            proba_trap_config[i + 1] = sum(w) / np.count_nonzero(w)
            site_label = w  # start a new configuration for the next step
            w = np.zeros(self.L)   # use the same number of walkers
        return proba_trap_config

    def average_survival_proba(self):
        """mean survival probability over several trap configurations."""
        arr = np.zeros((self.ntrials, self.nsteps + 1))
        for i in range(self.ntrials):
            proba_trap_config = self.exact_enumeration_proba()
            arr[i, :] = proba_trap_config
        mean_proba = np.mean(arr, axis=0)
        return mean_proba
