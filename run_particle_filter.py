import csv
import time
import os
import numpy as np
from origin.particle_filter import PF_Bootstrap_OE
from origin.data_parser import debris_parser

def main():
    # Setup particle filter
    Qpf = np.zeros((7, 7))
    # Measurement noise (this is important for the performance of PF!)
    Rpf0 = np.diag([
        100.0**2, 0.01**2, 0.01**2, 0.1**2, 0.1**2, 0.5**2,
    ]) * 1.0
    Rpf1 = np.diag([50**2, 0.001**2, 0.003**2, 0.1**2, 0.1**2, 0.1**2])
    Npf = 400
    pf = PF_Bootstrap_OE(Qpf, Rpf0, Rpf1, Npf, switch_R_after=3)
    cov_oe0 = np.diag([
        50.0**2, 0.001**2, 0.003**2, 0.1**2, 0.1**2, 0.1**2,
    ])

    logdir = "data/results"
    if not os.path.isdir(logdir):
            os.makedirs(logdir)

    # For each debris with more than 2 data points, run the PF filter and save the estimates in csv.
    for i in range(100):
        id = i+1
        debris_data = debris_parser(id, is_training=False)
        ephem = debris_data["ephemeris"]
        time_obs = np.flip(ephem.time)

        if len(ephem.time) >= 2:  # then run particle filter
            oe0 = ephem.get_keplerian()[-1]
            obs_his = ephem.get_keplerian()[1:]
            time_start = time.time()
            pf_out = pf(time_obs, obs_his, oe0, cov_oe0, is_parallel=True)
            time_end = time.time()
            print("Total time running PF for debris id {0}: {1} min".format(id, (time_end - time_start)/60))

            # Output in csv
            with open(os.path.join(logdir, "pf_result.csv"), "W") as f:
                writer = csv.writer(f)
                row = []
                row.append(id)
                row.extend(pf_out[1][-1, :].tolist())
                row.extend(np.sqrt(pf_out[2][-1, :]).tolist())
                writer.writerow(row)
        else:
            print("Skipping debris id {} as there is only one data point.".format(id))



if __name__ == "__main__":
    main()