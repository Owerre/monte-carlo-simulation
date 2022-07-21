#################################################
# Author: S. A. Owerre
# Date modified: 21/07/2022
# Function: Test for 1D Randomly Distributed Traps
##################################################

import sys
base_path = '/Users/sowerre/Documents/python/'
sys.path.append(base_path + 'monte-carlo/monte_carlo/random_walk_1d/src/')
import randomly_distributed_trap_1d as walk

def test_randomly_distributed_trap():
    # test default values
    model = walk.RandomlyDistributedTrap(10, 100, 0.5, 0.1, 20)
    assert model.nsteps == 10
    assert model.ntrials == 100
    assert model.p == 0.5
    assert model.rho == 0.1
    assert model.L == 20

    # test lattice method
    sites, site_label = model.lattice()
    assert len(sites) == model.L
    assert len(site_label) == model.L

    # test step_count_b4_trap method
    traj, count = model.step_count_b4_trap()
    assert len(traj) != 0
    assert count != 0

    # test average_nsteps_trap method
    mean_step_count = model.average_nsteps_trap()
    assert mean_step_count != 0

    # test exact_enumeration_proba method
    proba_trap_config = model.exact_enumeration_proba()
    assert len(proba_trap_config) == model.nsteps +1

    # test average_survival_proba method
    mean_proba = model.average_survival_proba()
    assert len(mean_proba) ==  model.nsteps +1


