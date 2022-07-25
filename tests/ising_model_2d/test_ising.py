#######################################
# Author: S. A. Owerre
# Date modified: 21/07/2022
# Function: Test for 2D Ising Model
#####################################

import sys
base_path = ''
sys.path.append(base_path + 'monte-carlo/monte_carlo/ising_model_2d/src/')
import Ising_model_2d as ising

def test_ising():
    # test default values
    model = ising.Ising(5, 10, 100)
    assert model.L == 5
    assert model.N == 5*5
    assert model.nwarmup == 10
    assert model.nsteps == 100

    # test initialize method
    spinconf = model.initialize()
    assert len(spinconf) == model.L

    # test neighbor_pos method
    np_, nm = model.neighbor_pos()
    assert len(np_) == model.L
    assert len(nm) == model.L

    # test precom_expo method
    pw = model.precom_expo(0.1)
    assert len(pw) == 17

    # test energy method
    ene = model.energy(spinconf)
    assert ene != 0

    # test magnetization method
    mag = model.magnetization(spinconf)
    assert mag != 0

    # test metropolis method
    ene, mag = model.metropolis(spinconf, pw)
    assert ene != 0
    assert mag != 0

    # test mcsweeps method
    avg = model.mcsweeps(spinconf, pw)
    assert len(avg) == 6