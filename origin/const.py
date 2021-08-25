import numpy as np

#constants
GMe = 3.986004407799724e+5 # [km^3/sec^2]
GMo = 1.32712440018e+11 #[km^3/sec^2]
GMm = 4.9028e+3 #[km^3/sec^2]
Re = 6378.1363 #[km]
C20 = -4.84165371736e-4
C22 = 2.43914352398e-6
S22 = -1.40016683654e-6
theta_g = (np.pi/180)*280.4606 #[rad]
nu_e = (np.pi/180)*(4.178074622024230e-3) #[rad/sec]
nu_o = (np.pi/180)*(1.1407410259335311e-5) #[rad/sec]
nu_ma = (np.pi/180)*(1.512151961904581e-4) #[rad/sec]
nu_mp = (np.pi/180)*(1.2893925235125941e-6) #[rad/sec]
nu_ms = (np.pi/180)*(6.128913003523574e-7) #[rad/sec]
alpha_o = 1.49619e+8 #[km]
epsilon = (np.pi/180)*23.4392911 #[rad]
phi_o = (np.pi/180)*357.5256 #[rad]
Omega_plus_w = (np.pi/180)*282.94 #[rad]
PSRP = 4.56e-3 #[kg/(km*sec^2)]