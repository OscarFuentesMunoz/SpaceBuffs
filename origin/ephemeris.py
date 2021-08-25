import numpy as np
import pykep as pk
from origin.const import GMe


# Define ephemeris object 
class Ephemeris:
    def __init__(self, time, states, state_type="keplerian") -> None:
        """
        time: np.array (n_time, )
        states: np.array (n_time, 6)
        state_type: either "keplerian" or "cartesian"
        """
        if time.ndim != 1:
            raise ValueError("Time array must be 1-D numpy array")
        if states.ndim != 2:
            raise ValueError("State array must be 2-D numpy array")
        self.time = np.array(time, dtype=float)
        self.size = np.shape(states)
        if state_type.lower() == "keplerian":
            self._keplerian = np.array(states, dtype=float)
        elif state_type.lower() == "cartesian":
            self._cartesian = np.array(states, dtype=float)
        else:
            raise ValueError("state type not recognized")

    
    def get_keplerian(self):
        try:
            value = self._keplerian
        except AttributeError:
            # Compute keplerian elements
            size = self.size
            keplerian = np.zeros(size)
            for i in range(size[0]):
                keplerian[i] = pk.ic2par(self._cartesian[i, 0:3], self._cartesian[i, 3:], GMe)
            self._keplerian = keplerian
            value = keplerian
        return value


    def get_cartesian(self):
        try:
            value = self._cartesian
        except AttributeError:
            # Compute cartesian elements
            size = self.size
            cartesian = np.zeros(size)
            for i in range(size[0]):
                r, v = pk.par2ic(self._keplerian[i], GMe)
                cartesian[i] = np.append(r, v)
            self._cartesian = cartesian
            value = cartesian
        return value