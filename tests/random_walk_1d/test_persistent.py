#############################################
# Author: S. A. Owerre
# Date modified: 21/07/2022
# Function: Test for Persistent Random Walk 1D
#############################################

import sys
base_path = ''
sys.path.append(base_path + 'monte-carlo/monte_carlo/random_walk_1d/src/')
import persistent_random_walk_1d as walk

def test_persistent():
    # test default values
    model = walk.PersistentRandomWalk1D(10, 100, 0.5)
    assert model.nsteps == 10
    assert model.ntrials == 100
    assert model.p == 0.5

    # test monte carlo method
    x_arr, visited_sites, x_avg, sigma2 = model.monte_carlo()
    assert x_arr.shape == (model.ntrials, model.nsteps+1)
    assert len(visited_sites) != 0
    assert len(x_avg) == model.nsteps+1
    assert len(sigma2) == model.nsteps+1

    # test sites visited method
    count = model.sites_visited()
    assert len(count) == model.nsteps+1

    # test average sites visited method
    mean_count = model.average_sites_visited()
    assert len(mean_count) == model.nsteps+1

