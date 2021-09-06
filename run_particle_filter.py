import csv
import time
import os
import argparse
import numpy as np
from origin.particle_filter import PF_Bootstrap_OE
from origin.data_parser import debris_parser

def main(fname):
    # Setup particle filter
    Qpf = np.zeros((7, 7))
    # Measurement noise (this is important for the performance of PF!)
    Rpf0 = np.diag([
        50.0**2, 0.001**2, 0.003**2, 0.1**2, 0.1**2, 0.1**2,
    ]) * 10
    # Rpf1 = np.diag([50**2, 0.001**2, 0.003**2, 0.1**2, 0.1**2, 0.1**2])
    Rpf1 = Rpf0
    Npf = 600
    cov_oe0 = np.diag([
        50.0**2, 0.001**2, 0.003**2, 0.1**2, 0.1**2, 0.1**2,
    ])

    logdir = "data/results"
    if not os.path.isdir(logdir):
            os.makedirs(logdir)
    id_already_processed = []
    try:
        with open(os.path.join(logdir, fname)) as f:
            reader = csv.reader(f)
            for row in reader:
                id = int(float(row[0]))
                id_already_processed.append(id)
        print("The entries for the following debris already exist. The computation will be skipped for those.")
        print(id_already_processed)
    except:
        print("The output file does not contain any data.")
        print("Starting from the debris id 1.")

    # For each debris with more than 2 data points, run the PF filter and save the estimates in csv.
    for i in range(100):
        id = i+1
        if id in id_already_processed:
            continue

        # Otherwise, do the following
        debris_data = debris_parser(id, is_training=False)
        ephem = debris_data["ephemeris"]
        time_obs = np.flip(ephem.time)
        print("Debris ID {}".format(id))

        if len(ephem.time) >= 2:  # then run particle filter
            pf = PF_Bootstrap_OE(Qpf, Rpf0, Rpf1, Npf, switch_R_after=3)
            oe0 = ephem.get_keplerian()[-1]
            obs_his = np.flip(ephem.get_keplerian(), axis=0)[1:, :]
            time_start = time.time()
            pf_out = pf(time_obs, obs_his, oe0, cov_oe0, is_parallel=True)
            time_end = time.time()
            print("Total time running PF for debris id {0}: {1} min".format(id, (time_end - time_start)/60))

            # Output in csv
            with open(os.path.join(logdir, fname), "a") as f:
                writer = csv.writer(f)
                row = []
                row.append(id)
                row.append(time_obs[-1] / 86400.0)
                row.extend(pf_out[1][-1, :].tolist())
                row.extend(np.sqrt(pf_out[2][-1, :]).tolist())
                writer.writerow(row)
        else:
            print("Skipping debris id {} as there is only one data point.".format(id))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Perform particle filtering for test debris.')
    parser.add_argument('--name', type=str, default="pf_result.csv",
                        help='Name of the output csv file. Must be name.csv.')
    args = parser.parse_args()
    main(args.name)