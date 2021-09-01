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
    Rpf = np.diag([
        100.0**2, 0.01**2, 0.01**2, 0.1**2, 0.1**2, 2.0**2,
    ]) * 1.0
    Npf = 200
    pf = PF_Bootstrap_OE(Qpf, Rpf, Npf)

    # Call PF
    time_obs = np.flip(ephem.time)
    oe0 = ephem.get_keplerian()[-1]
    cov_oe0 = np.diag([
        50.0**2, 0.001**2, 0.003**2, 0.05**2, 0.05**2, 0.05**2,
    ])
    obs_his = ephem.get_keplerian()[1:]
    time_start = time.time()
    pf_out = pf(time_obs, obs_his, oe0, cov_oe0, is_parallel=True)
    time_end = time.time()
    print("Total time running PF: {} min".format((time_end - time_start)/60))

    # TODO need bug fix
    # # Compute truth
    # print("Propagating the truth")
    # ephemeris_event = out["event_ephemeris"]
    # x0 = ephemeris_event.get_cartesian()[-1]
    # t0 = ephemeris_event.time[-1]
    # prop = Propagator()\
    #     .state(x0)\
    #     .time(t0)\
    #     .cram(out["cram"])

    # t_grid = ephem.time
    # ephem_truth = prop(t_grid)
    # diff = np.flip(pf_out[1], axis=0) - \
    #     np.concatenate(
    #         (
    #             ephem_truth.get_keplerian(), 
    #             np.log10(cram*1e6) * np.ones(len(t_grid)).reshape((len(t_grid), 1)),
    #         ),
    #         axis=1
    #     )
    # normalize_diff(diff[:, 2])
    # normalize_diff(diff[:, 3])
    # normalize_diff(diff[:, 4])
    # normalize_diff(diff[:, 5])

    # diff_obs = ephem.get_keplerian()

    # # Plot state errors
    # fig, axs = plt.subplots(4, 2)
    # axs = axs.reshape(8,)
    # for i in range(7):
    #     axs[i].scatter(t_grid / 86400.0, diff[:, i], )
    #     if i in (0,1,2,3,4,5):
    #         axs[i].scatter(t_grid / 86400.0, )
    #     axs[i].grid()
    # axs[0].set_ylabel("a (km)")
    # axs[1].set_ylabel("e")
    # axs[2].set_ylabel("i (rad)")
    # axs[3].set_ylabel("W (rad)")
    # axs[4].set_ylabel("w (rad)")
    # axs[5].set_ylabel("E (rad)")
    # axs[6].set_ylabel("log10[Cr(A/m)]")
    # plt.tight_layout()
    # plt.show()


    print("done")

