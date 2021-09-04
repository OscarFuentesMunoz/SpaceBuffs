from origin.ephemeris import Ephemeris
import numpy as np
import argparse
import csv
import os
from bayes_opt import BayesianOptimization
import time

from numpy.lib.utils import info
from origin.screening import Screen
from origin.data_parser import debris_parser

# Using Bayesian optimization, perform matching and output the solution into csv
def main(fname_out, fname_est, n_iter0, n_iter1):

    # Setup
    time_very_beginning = time.time()
    solution = []  # array of [debris_id, parent_id, cram]
    logdir = "data/solution"
    if not os.path.isdir(logdir):
            os.makedirs(logdir)
    
    # Read existing estimation results
    ids_estimated = []  # array of int
    info_estimated = []  # array of dict
    with open(os.path.join("data/results", fname_est)) as f:
        reader = csv.reader(f)
        for row in reader:
            ids_estimated.append(float(row[0]))
            epoch_in_days = float(row[1])
            state = np.array([
                float(row[2]), float(row[3]), float(row[4]),
                float(row[5]), float(row[6]), float(row[7]),
                float(row[8]),
            ])
            sigmas = np.array([
                float(row[9]), float(row[10]), float(row[11]),
                float(row[12]), float(row[13]), float(row[14]),
                float(row[15]),
            ])
            info_estimated.append({
                "epoch": float(row[1]),
                "state": state,
                "sigmas": sigmas,
            })

    # If the output file already has entries, read these and skip computation
    id_already_processed = []
    try:
        with open(os.path.join(logdir, fname_out)) as f:
            reader = csv.reader(f)
            for row in reader:
                id = int(float(row[0]))
                cost = float(row[1])
                sat_id = int(float(row[2]))
                cram = float(row[3])
                solution.append([id, cost, sat_id, cram])
                id_already_processed.append(id)
        print("The entries for the following debris already exist. The computation will be skipped for those.")
        print(id_already_processed)
    except:
        print("The output file does not contain any data.")
        print("Starting from the debris id 1.")

    # Setup screening parameters
    screen = Screen()
    R = np.diag([
        50.0**2, 0.001**2, 0.003**2, 0.1**2, 0.1**2, 0.1**2,
    ])

    # Loop over 100 debris
    for i in range(100):
        # Check if the entry is alrady in the output file
        id = i + 1
        if id in id_already_processed:
            continue  # => skip the entry
        
        print("Debris ID: {}".format(id))

        # Define Bayesian optimization spec depending on whether the estimates are available or not. Need (obj_func, state_bounds, n_iter, init_points)
        # Check if PF estimate exists
        if id in ids_estimated:
            init_points = 0
            n_iter = n_iter0
            ind = ids_estimated.index(id)
            epoch_s = info_estimated[ind]["epoch"] * 86400.0
            state = info_estimated[ind]["state"]  # cram is in m^2/kg
            sigmas = info_estimated[ind]["sigmas"]
            state_upper = np.zeros(7)
            state_upper[0] = state[0] + sigmas[0]/2
            state_upper[1] = state[1] + sigmas[1]/2
            state_upper[2] = state[2] + sigmas[2]*0.75
            state_upper[3] = state[3] + sigmas[3]/2
            state_upper[4] = state[4] + sigmas[4]/2
            state_upper[5] = state[5] + sigmas[5]/4
            state_upper[6] = np.log10(state[6] + sigmas[6])
            state_lower = np.zeros(7)
            state_upper[0] = state[0] - sigmas[0]/2
            state_upper[1] = state[1] - sigmas[1]/2
            state_upper[2] = state[2] - sigmas[2]*0.75
            state_upper[3] = state[3] - sigmas[3]/2
            state_upper[4] = state[4] - sigmas[4]/2
            state_upper[5] = state[5] - sigmas[5]/4
            state_lower[6] = np.log10(state[6] - sigmas[6])
            state_probe = np.zeros(7)
            state_probe[:6] = state[:6]
            state_probe[6] = np.log10(state[6])
        else:
            n_iter = n_iter1
            init_points = 0
            debris_data = debris_parser(deb_id=id, is_training=False)
            ephem0 = debris_data["ephemeris"]
            epoch_s = ephem0.time[0]
            keplerian = ephem0.get_keplerian()[0]
            sigmas_kep = np.array([
                70.0, 0.001, 0.004, 0.8, 0.1, 0.1, 
            ])
            state_upper = np.zeros(7)
            state_upper[:6] = keplerian + sigmas_kep
            state_upper[6] = 1.8
            state_lower = np.zeros(7)
            state_lower[:6] = keplerian - sigmas_kep
            state_lower[6] = -0.5
            state_probe = np.zeros(7)
            state_probe[:6] = keplerian
            state_probe[6] = np.log10(0.65)

        obj_func = lambda a, e, i, W, w, E, log_cram: \
                -screen.argmin_rms(
                epoch_s, a, e, i, W, w, E, 10**log_cram*1e-6, R
            )[0]
        state_bounds = {
                "a": (state_lower[0], state_upper[0]),
                "e": (state_lower[1], state_upper[1]),
                "i": (state_lower[2], state_upper[2]),
                "W": (state_lower[3], state_upper[3]),
                "w": (state_lower[4], state_upper[4]),
                "E": (state_lower[5], state_upper[5]),
                "log_cram": (state_lower[6], state_upper[6]),
            }


        # Run Bayesian optimization
        optimizer = BayesianOptimization(
            f=obj_func,
            pbounds=state_bounds,
            random_state=1,
        )

        optimizer.probe(
            params={
                "a": state_probe[0],
                "e": state_probe[1],
                "i": state_probe[2],
                "W": state_probe[3],
                "w": state_probe[4],
                "E": state_probe[5],
                "log_cram": state_probe[6],
            },
            lazy=True,
        )
        time_start = time.time()
        optimizer.maximize(
            init_points=init_points,
            n_iter=n_iter,
        )
        time_stop = time.time()
        print("Computation finished in {} min".format((time_stop - time_start)/60))
        print("Cumulative computation time: {} min".format((time_stop - time_very_beginning)/60))


        # Get corresponding satellite ID
        opt_result = optimizer.max
        cost = -opt_result["target"]
        state_opt = opt_result["params"]
        screen_result = screen.argmin_rms(
            epoch_s, 
            state_opt["a"],
            state_opt["e"],
            state_opt["i"],
            state_opt["W"], 
            state_opt["w"], 
            state_opt["E"], 
            10**state_opt["log_cram"]*1e-6,
            R
        )
        sat_id = screen_result[1]
        cram = 10**state_opt["log_cram"]

        print("Matched satellite id: {}".format(screen_result[1]))
        print("Estimated Cr(A/m): {}".format(10**state_opt["log_cram"]*1e-6))
        print("Estimated log10[Cr(A/m)]: {}".format(state_opt["log_cram"]))
        print("Optimal cost: {}".format(cost))


        # Save output to csv  (debris_id, parent_id, cram)
        with open(os.path.join(logdir, fname_out), "a") as f:
            writer = csv.writer(f)
            writer.writerow([id, cost, sat_id, cram])
        solution.append([id, cost, sat_id, cram])

    # Save another text file for submission
    solution_np = np.array([np.array(row) for row in solution])
    solution_np = solution_np[solution_np[:, 0].argsort()]
    np.savetxt(os.path.join(logdir, fname_out.replace(".csv", "_for_submission.txt")), solution_np[:, 2:])


    print("Total cost: {}".format(np.sum(solution_np[:, 1])))
    return solution


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Perform matching and output solution.')
    parser.add_argument('--name_out', type=str, default="solution_bayes.csv",
                        help='Name of the output txt file.')
    parser.add_argument('--name_source', type=str, default="pf_result.csv",
                        help='Name of the csv file in data/results directory that contains state estimates.')
    parser.add_argument('--n_iter0', type=int, default="12",
                        help='Number of iteration for debris with at least 2 measurements.')
    parser.add_argument('--n_iter1', type=int, default="30",
                        help='Number of iteration for debris with only one measurement.')
    
    args = parser.parse_args()
    main(args.name_out, args.name_source, args.n_iter0, args.n_iter1)
