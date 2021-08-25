import os
import origin.ephemeris as ephemeris
import numpy as np
import kepler

# function that reads the contents of the debris data, given debris number
def debris_parser(deb_id:int, is_training:bool):
    out = dict()
    if is_training:
        base_dir = "data/space-debris-the-origin/deb_train"
        fname = "eledebtrain{:03}".format(deb_id) + ".dat"
    else:
        base_dir = "data/space-debris-the-origin/deb_test"
        fname = "eledebnewfd{:03}".format(deb_id) + ".dat"
    path = os.path.join(base_dir, fname)

    # Read the dat file
    with open(path) as f:
        data_all = f.readlines()
        data_all_np = np.zeros((len(data_all), 7))
        # data: epoch, a, e, i, M, w, W (angles are in deg)
        # what we need: epoch, a, e, i, W, w, E (angles are in rad)
        indices = [0, 1, 2, 3, 6, 5, 4]
        for i, data in enumerate(data_all):
            data_all_np[i, :] = np.array(list(map(float, data.split())))[indices]
            data_all_np[i, 0] = data_all_np[i, 0] * 86400.0  # convert from day to sec
            data_all_np[i, (3,4,5,6)] = data_all_np[i, (3,4,5,6)] * np.pi/180 # conver from deg to rad
            eccentric_anomaly, _, _ = kepler.kepler(data_all_np[i, 6], data_all_np[i, 2])
            data_all_np[i, 6] = eccentric_anomaly


    time = data_all_np[:, 0]
    keplerian = data_all_np[:, 1:]
    ephem = ephemeris.Ephemeris(time, keplerian, state_type="keplerian")

    out["deb_id"] = deb_id
    out["path"] = path
    out["time"] = time
    out["ephemeris"] = ephem

    if is_training:
        sol_path = os.path.join(base_dir, "../", "labels_train.dat")
        with open(sol_path) as f:
            sol_info = f.readlines()[deb_id-1].split()
            
        out["sat_id"] = int(sol_info[0])
        out["cram"] = float(sol_info[1]) * 1e-6
        time_event = np.array([float(sol_info[2])]) * 86400.0  # convert from day to sec
        keplerian_event_orig = np.array([list(map(float, sol_info[3:]))])
        keplerian_event = keplerian_event_orig[:, (0, 1, 2, 5, 4, 3)]
        keplerian_event[0, 2:] = keplerian_event[0, 2:] * np.pi/180
        eccentric_anomaly, _, _ = kepler.kepler(keplerian_event[0, 5], keplerian_event[0, 1])
        keplerian_event[0, 5] = eccentric_anomaly
        ephem_event = ephemeris.Ephemeris(
            time_event, keplerian_event, state_type="keplerian"
        )
        out["event_ephemeris"] = ephem_event


    return out
    

# function that reads the contents of the satellite data, given satellite number
def satellite_parse(sat_id:int):
    out = dict()
    base_dir = "data/space-debris-the-origin/sat"
    fname = "elesat{:03}".format(sat_id) + ".dat"
    path = os.path.join(base_dir, fname)
    with open(path) as f:
        data_all = f.readlines()
        data_all_np = np.zeros((len(data_all), 7))
        indices = [0, 1, 2, 3, 6, 5, 4]
        for i, data in enumerate(data_all):
            data_all_np[i, :] = np.array(list(map(float, data.split())))[indices]
            data_all_np[i, 0] = data_all_np[i, 0] * 86400.0  # convert from day to sec
            data_all_np[i, (3,4,5,6)] = data_all_np[i, (3,4,5,6)] * np.pi/180 # conver from deg to rad
            eccentric_anomaly, _, _ = kepler.kepler(data_all_np[i, 6], data_all_np[i, 2])
            data_all_np[i, 6] = eccentric_anomaly

    time = data_all_np[:, 0]
    keplerian = data_all_np[:, 1:]
    ephem = ephemeris.Ephemeris(time, keplerian, state_type="keplerian")
    out["sat_id"] = sat_id
    out["ephemeris"] = ephem

    return out


if __name__ == "__main__":
    debris_parser(1, True)
    print("done")