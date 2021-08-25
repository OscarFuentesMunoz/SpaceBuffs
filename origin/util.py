import pykep as pk



def normalize_diff(radians):
    mask  = radians < -180 * pk.DEG2RAD
    radians[mask] = radians[mask] + 360* pk.DEG2RAD
    mask  = radians > 180 * pk.DEG2RAD
    radians[mask] = radians[mask] - 360* pk.DEG2RAD