from origin.data_parser import debris_parser
from origin.propagation import Propagator
from origin.ephemeris import Ephemeris
from origin.util import normalize_diff
from origin.particle_filter import PF_Bootstrap_OE
import numpy as np
import time
from matplotlib import pyplot as plt

if __name__ == "__main__":
    np.set_printoptions(linewidth=150)

    out = debris_parser(deb_id=25, is_training=True)
    ephem = out["ephemeris"]
    cram = out["cram"]
    print("True Cr(A/m): {}".format(cram))

    # Construct PF object
    # Process noise
    Qpf = np.zeros((7, 7))
    # Measurement noise (this is important for the performance of PF!)
    Rpf0 = np.diag([
        50.0**2, 0.001**2, 0.003**2, 0.1**2, 0.1**2, 0.1**2,
    ]) * 10
    # Rpf1 = np.diag([50**2, 0.001**2, 0.003**2, 0.1**2, 0.1**2, 0.1**2])
    Rpf1 = Rpf0
    Npf = 600
    pf = PF_Bootstrap_OE(Qpf, Rpf0, Rpf1, Npf, switch_R_after=3)

    # Call PF
    time_obs = np.flip(ephem.time)
    oe0 = ephem.get_keplerian()[-1]
    cov_oe0 = np.diag([
        50.0**2, 0.001**2, 0.003**2, 0.1**2, 0.1**2, 0.1**2,
    ])
    obs_his = np.flip(ephem.get_keplerian(), axis=0)[1:, :]
    time_start = time.time()
    pf_out = pf(time_obs, obs_his, oe0, cov_oe0, is_parallel=True)
    time_end = time.time()
    print("Total time running PF: {} min".format((time_end - time_start)/60))

    # Compute truth
    print("Propagating the truth")
    ephemeris_event = out["event_ephemeris"]
    x0 = ephemeris_event.get_cartesian()[-1]
    t0 = ephemeris_event.time[-1]
    prop = Propagator()\
        .state(x0)\
        .time(t0)\
        .cram(out["cram"])

    t_grid = ephem.time
    ephem_truth = prop(t_grid)
    diff = np.flip(pf_out[1], axis=0) - \
        np.concatenate(
            (
                ephem_truth.get_keplerian(), 
                (cram*1e6) * np.ones(len(t_grid)).reshape((len(t_grid), 1)),
            ),
            axis=1
        )
    normalize_diff(diff[:, 2])
    normalize_diff(diff[:, 3])
    normalize_diff(diff[:, 4])
    normalize_diff(diff[:, 5])

    diff_obs = ephem.get_keplerian() - ephem_truth.get_keplerian()
    normalize_diff(diff_obs[:, 2])
    normalize_diff(diff_obs[:, 3])
    normalize_diff(diff_obs[:, 4])
    normalize_diff(diff_obs[:, 5])

    sigmas = np.sqrt(np.flip(pf_out[2], axis=0))

    # Plot state errors
    fig, axs = plt.subplots(3, 3)
    axs = axs.reshape(9,)
    for i in range(7):
        axs[i].scatter(t_grid / 86400.0, diff[:, i], label="estimation error")
        axs[i].plot(t_grid / 86400.0, sigmas[:, i], linestyle="dashed", color="gray")
        axs[i].plot(t_grid / 86400.0, -sigmas[:, i], linestyle="dashed", color="gray")
        if i in (0,1,2,3,4,5):
            axs[i].scatter(t_grid / 86400.0, diff_obs[:, i], label="obs error")
        axs[i].grid()
    axs[0].set_ylabel("a (km)")
    axs[1].set_ylabel("e")
    axs[2].set_ylabel("i (rad)")
    axs[3].set_ylabel("W (rad)")
    axs[4].set_ylabel("w (rad)")
    axs[5].set_ylabel("E (rad)")
    axs[6].set_ylabel("Cr(A/m)")
    axs[0].legend()
    axs[7].hist(pf_out[0][-1, :, 6]*1e-6); axs[7].plot((out["cram"], out["cram"]), axs[7].get_ylim())
    axs[7].set_ylabel("Count"); axs[7].set_xlabel("Cr(A/m)")
    plt.tight_layout()
    plt.show()


    print("done")

