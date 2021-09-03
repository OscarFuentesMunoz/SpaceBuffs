from numpy.lib.function_base import diff
from origin.propagation import Propagator
from origin.ephemeris import Ephemeris
from origin.data_parser import satellite_parse
from origin.util import normalize_diff
import numpy as np

class Screen():
    def __init__(self) -> None:
        # Store ephemeris of every satellites
        ephemeris = np.empty((100, 1096, 6))  # (sat_id, time_index, orbital_elem)
        for i in range(100):
            sat_id = i+1
            sat_data = satellite_parse(sat_id)
            ephem = sat_data["ephemeris"]
            ephemeris[i, :, :] = ephem.get_keplerian()
        self.sat_ephemeris = ephemeris
        self.time = ephem.time
    

    def argmin_rms(self, epoch, a, e, i, W, w, E, cram, R):
        """
        Given a debris ephemeris at a single epoch, cram, and covariance, returns the minimum Mahalanobis distance, corresponding satellite id, 
        and event epoch.
        """
        # Backward propagate the debris trajectory to match the time step of the satellites' ephemeris.
        orb_elems = [a,e,i,W,w,E]
        ephem_debris0 = Ephemeris(
            time=np.array([epoch]),
            states=np.array([orb_elems]),
            state_type="keplerian",
        )
        x0 = ephem_debris0.get_cartesian()[0]  # assuming the input only has state at an epoch
        t0 = ephem_debris0.time[0]
        prop = Propagator()\
            .state(x0)\
            .time(t0)\
            .cram(cram)
        is_valid = self.time <= t0
        t_grid = self.time[is_valid]
        ephem_debris = prop(t_grid)

        # Handle a case where the debris collide with the Earth
        if ephem_debris.get_keplerian().size == 0:  # which means current cram is super infeasible
            return 1e8, None, None
            
        elif ephem_debris.get_cartesian().shape[0] < len(t_grid):  
            # then redefine the time array and propagate again
            len_t = ephem_debris.get_cartesian().shape[0]
            t_min = t_grid[-len_t+1]  # to be safe +1
            is_valid = is_valid and t_min < self.time
            t_grid = self.time[is_valid]
            prop = Propagator()\
                .state(x0)\
                .time(t0)\
                .cram(cram)
            ephem_debris = prop(t_grid)

        
        # Compute weighted RMS
        kep_target = self.sat_ephemeris[:, is_valid, :]
        kep_debris = ephem_debris.get_keplerian()
        diff_kep = kep_debris - kep_target  # (sat_id, time_index, orbital_elem)
        normalize_diff(diff_kep[:, :, 2])
        normalize_diff(diff_kep[:, :, 3])
        normalize_diff(diff_kep[:, :, 4])
        normalize_diff(diff_kep[:, :, 5])
        loss = np.sum(np.dot(diff_kep, np.linalg.pinv(R)) * diff_kep, axis=2)

        # Find argmin of the RMS
        ind = np.unravel_index(np.argmin(loss), loss.shape)
        sat_id = ind[0] + 1
        ephem_event = Ephemeris(np.array([t_grid[ind[1]]]), np.array([kep_target[ind[0], ind[1], :]]), state_type="keplerian")
        loss_min = loss[ind]
        
        return loss_min, sat_id, ephem_event
