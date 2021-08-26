from origin.data_parser import debris_parser
from origin.propagation import Propagator
from origin.ephemeris import Ephemeris
from origin.util import normalize_diff
import numpy as np
import time
from matplotlib import pyplot as plt

np.set_printoptions(linewidth=150)

out = debris_parser(deb_id=25, is_training=True)
ephem = out["ephemeris"]
ephem_event = out["event_ephemeris"]
time_all = np.concatenate((ephem_event.time, ephem.time))
kep_all = np.concatenate((ephem_event.get_keplerian(), ephem.get_keplerian()))
ephem_all = Ephemeris(time_all, kep_all, state_type="keplerian")

# Setup propagation job
x0 = ephem_all.get_cartesian()[0]
t0 = ephem_all.time[0]
prop = Propagator()\
    .state(x0)\
    .time(t0)\
    .cram(out["cram"])

t_grid = ephem_all.time
start_time = time.time()
ephem_out = prop(t_grid)
end_time = time.time()
print("Propagation wall-clock time: {} sec\n".format(end_time - start_time))


# Compute diff
diff_kep = ephem_all.get_keplerian() - ephem_out.get_keplerian()
normalize_diff(diff_kep[:, -1])
normalize_diff(diff_kep[:, -2])
normalize_diff(diff_kep[:, -3])
normalize_diff(diff_kep[:, -4])
diff_cart = ephem_all.get_cartesian() - ephem_out.get_cartesian()

print("Time grid (days):")
print(ephem_all.time / 86400)
print("Keplerian (Observed):")
print(ephem_all.get_keplerian())
print("Keplerian (Propagated):")
print(ephem_out.get_keplerian())
print("Diff in Keplerian:")
print(diff_kep)

print("Cartesian (Observed):")
print(ephem_all.get_cartesian())
print("Cartesian (Propagated):")
print(ephem_out.get_cartesian())
print("Diff in Cartesian:")
print(diff_cart)


# Plot difference in the state to better understand the noise characteristics
# Diff in Keplerian
fig, axs = plt.subplots(3, 2)
axs = axs.reshape(6,)
for i in range(6):
    axs[i].scatter(ephem_all.time / 86400.0, diff_kep[:, i], )
    axs[i].grid()
axs[0].set_ylabel("a (km)")
axs[1].set_ylabel("e")
axs[2].set_ylabel("i (rad)")
axs[3].set_ylabel("W (rad)")
axs[4].set_ylabel("w (rad)")
axs[5].set_ylabel("E (rad)")

# Diff in Cartesian
fig, axs = plt.subplots(3, 2)
axs = axs.reshape(6,)
for i in range(6):
    axs[i].scatter(ephem_all.time / 86400.0, diff_cart[:, i], )
    axs[i].grid()
axs[0].set_ylabel("X (km)")
axs[1].set_ylabel("Y (km)")
axs[2].set_ylabel("Z (km)")
axs[3].set_ylabel("Vx(km/s)")
axs[4].set_ylabel("Vy (km/s)")
axs[5].set_ylabel("Vz (km/s)")


plt.tight_layout()
plt.show()

print("done")

