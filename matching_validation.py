from origin.ephemeris import Ephemeris
from origin.screening import Screen
from origin.data_parser import debris_parser
import numpy as np
import time
from bayes_opt import BayesianOptimization

# Load training data
out = debris_parser(deb_id=25, is_training=True)
ephem = out["ephemeris"]
sat_id_true = out["sat_id"]
print("True satellite id: {}".format(sat_id_true))

# Setup screening parameters
screen = Screen()
R = np.diag([
    50.0**2, 0.001**2, 0.003**2, 0.1**2, 0.1**2, 0.1**2,
])
cram0 = out["cram"]
print("log10 of cram: {}".format(np.log10(cram0*1e6)))
orb_elems0 = ephem.get_keplerian()[0]

# Bayesian optimization setup
obj_func = lambda a, e, i, W, w, E, log_cram: -screen.argmin_rms(
    ephem.time[0], a, e, i, W, w, E, 10**log_cram*1e-6, R
)[0]
state_bounds = {
    "a": tuple((orb_elems0[0] + np.array([-50, 50])).tolist()),
    "e": tuple((orb_elems0[1] + np.array([-0.001, 0.001])).tolist()),
    "i": tuple((orb_elems0[2] + np.array([-0.003, 0.003])).tolist()),
    "W": tuple((orb_elems0[3] + np.array([-0.1, 0.1])).tolist()),
    "w": tuple((orb_elems0[4] + np.array([-0.1, 0.1])).tolist()),
    "E": tuple((orb_elems0[5] + np.array([-0.1, 0.1])).tolist()),
    "log_cram": tuple((np.log10(cram0*1e6) + np.array([-0.02, 0.02])).tolist()),
}
optimizer = BayesianOptimization(
    f=obj_func,
    pbounds=state_bounds,
    random_state=1,
)

# Run
time_start = time.time()
optimizer.maximize(
    init_points=10,
    n_iter=10,
)
time_stop = time.time()
print("Computation finished in {} sec".format(time_stop - time_start))

opt_result = optimizer.max
state_opt = opt_result["params"]
screen_result = screen.argmin_rms(
    ephem.time[0], 
    state_opt["a"],
    state_opt["e"],
    state_opt["i"],
    state_opt["W"], 
    state_opt["w"], 
    state_opt["E"], 
    10**state_opt["log_cram"]*1e-6,
    R
)

print("")
print("-----Summary-----")
print("True satellite id: {}".format(sat_id_true))
print("Matched satellite id: {}".format(screen_result[1]))
print("Estimated Cr(A/m): {}".format(10**state_opt["log_cram"]*1e-6))
print("True Cr(A/m): {}".format(cram0))
print("Estimated log10[Cr(A/m)]: {}".format(state_opt["log_cram"]))
print("True log10[Cr(A/m)]: {}".format(np.log10(cram0*1e6)))
print("----------")

print("done")