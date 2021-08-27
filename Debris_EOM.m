function dx = Debris_EOM(t,state,pars)
cram = pars.cram;

dx = zeros(6,1);
X = state(1);
Y = state(2);
Z = state(3);
VX = state(4);
VY = state(5);
VZ = state(6);

%constants
GMe = 3.986004407799724e+5; % [km^3/sec^2]
GMo = 1.32712440018e+11; %[km^3/sec^2]
GMm = 4.9028e+3; %[km^3/sec^2]
Re = 6378.1363; %[km]
C20 = -4.84165371736e-4;
C22 = 2.43914352398e-6;
S22 = -1.40016683654e-6;
theta_g = (pi/180)*280.4606; %[rad]
nu_e = (pi/180)*(4.178074622024230e-3); %[rad/sec]
nu_o = (pi/180)*(1.1407410259335311e-5); %[rad/sec]
nu_ma = (pi/180)*(1.512151961904581e-4); %[rad/sec]
nu_mp = (pi/180)*(1.2893925235125941e-6); %[rad/sec]
nu_ms = (pi/180)*(6.128913003523574e-7); %[rad/sec]
alpha_o = 1.49619e+8; %[km]
epsilon = (pi/180)*23.4392911; %[rad]
phi_o = (pi/180)*357.5256; %[rad]
Omega_plus_w = (pi/180)*282.94; %[rad]
PSRP = 4.56e-3; %[kg/(km*sec^2)]

%Sun's position
lo = phi_o + nu_o*t;
lambda_o = Omega_plus_w + lo + (np.pi/180)*( (6892/3600)*sin(lo) + (72/3600)*sin(2*lo) );
ro = (149.619 - 2.499*cos(lo) - 0.021*cos(2*lo))*(10^6);

Xo = ro*cos(lambda_o);
Yo = ro*sin(lambda_o)*np.cos(epsilon);
Zo = ro*sin(lambda_o)*np.sin(epsilon);

%Moon's position
phi_m = nu_o*t;
phi_ma = nu_ma*t;
phi_mp = nu_mp*t;
phi_ms = nu_ms*t;
L0 = phi_mp + phi_ma + (np.pi/180)*218.31617;
lm = phi_ma + (np.pi/180)*134.96292;
llm = phi_m + (np.pi/180)*357.5256;
Fm = phi_mp + phi_ma + phi_ms + (np.pi/180)*93.27283;
Dm = phi_mp + phi_ma - phi_m  + (np.pi/180)*297.85027;

rm = 385000 - 20905*cos(lm) - 3699*cos(2*Dm - lm) - 2956*cos(2*Dm) -...
     570*cos(2*lm) + 246*cos(2*lm - 2*Dm) - 205*cos(llm - 2*Dm) -...
     171*cos(lm + 2*Dm) - 152*cos(lm + llm - 2*Dm);
     
lambda_m = L0 + (np.pi/180)*( (22640/3600)*sin(lm) + (769/3600)*sin(2*lm) - (4856/3600)*sin(lm - 2*Dm) +...
     (2370/3600)*sin(2*Dm) - (668/3600)*sin(llm) - (412/3600)*sin(2*Fm) -...
     (212/3600)*sin(2*lm - 2*Dm) - (206/3600)*sin(lm + llm - 2*Dm) +...
     (192/3600)*sin(lm + 2*Dm) - (165/3600)*sin(llm - 2*Dm) +...
     (148/3600)*sin(lm - llm) - (125/3600)*sin(Dm) - (110/3600)*sin(lm + llm) -...
     (55/3600)*sin(2*Fm - 2*Dm) );
     
bm = (np.pi/180)*( (18520/3600)*sin(Fm + lambda_m - L0 + (np.pi/180)*((412/3600)*sin(2*Fm) + (541/3600)*sin(llm)) ) -...
     (526/3600)*sin(Fm - 2*Dm) + (44/3600)*sin(lm + Fm - 2*Dm) - (31/3600)*sin(-lm + Fm -2*Dm) -...
     (25/3600)*sin(-2*lm + Fm) - (23/3600)*sin(llm + Fm - 2*Dm) + (21/3600)*sin(-lm + Fm) +...
     (11/3600)*sin(-llm + Fm - 2*Dm) );
     
Xm =  cos(bm)*cos(lambda_m)*rm;
Ym = -np.sin(epsilon)*sin(bm)*rm + np.cos(epsilon)*cos(bm)*sin(lambda_m)*rm;
Zm =  np.cos(epsilon)*sin(bm)*rm + cos(bm)*np.sin(epsilon)*sin(lambda_m)*rm;

%Earth's Keplerian terms
magR2 = X^2 + Y^2 + Z^2;
fKepX = -GMe*X/(magR2^(3./2));
fKepY = -GMe*Y/(magR2^(3./2));
fKepZ = -GMe*Z/(magR2^(3./2));

%Earth's J2 terms
J2term1 = GMe*(Re^2)*np.sqrt(5)*C20/(2*magR2^(1./2));
J2term2 = 3/(magR2^2);
J2term3 = 15*(Z^2)/(magR2^3);
fJ2X = J2term1*X*(J2term2 - J2term3);
fJ2Y = J2term1*Y*(J2term2 - J2term3);
fJ2Z = J2term1*Z*(3*J2term2 - J2term3);

%Earth's C22 and S22 terms
x =  X*cos(theta_g + nu_e*t) + Y*sin(theta_g + nu_e*t);
y = -X*sin(theta_g + nu_e*t) + Y*cos(theta_g + nu_e*t);
z = Z;
magr2 = x^2 + y^2 + z^2;

C22term1 = 5*GMe*(Re^2)*np.sqrt(15)*C22/(2*magr2^(7./2));
C22term2 = GMe*(Re^2)*np.sqrt(15)*C22/(magr2^(5./2));
fC22x = C22term1*x*(y^2 - x^2) + C22term2*x;
fC22y = C22term1*y*(y^2 - x^2) - C22term2*y;
fC22z = C22term1*z*(y^2 - x^2);

S22term1 = 5*GMe*(Re^2)*np.sqrt(15)*S22/(magr2^(7./2));
S22term2 = GMe*(Re^2)*np.sqrt(15)*S22/(magr2^(5./2));
fS22x = -S22term1*(x^2)*y + S22term2*y;
fS22y = -S22term1*x*(y^2) + S22term2*x;
fS22z = -S22term1*x*y*z;

fC22X = fC22x*cos(theta_g + nu_e*t) - fC22y*sin(theta_g + nu_e*t);
fC22Y = fC22x*sin(theta_g + nu_e*t) + fC22y*cos(theta_g + nu_e*t);
fC22Z = fC22z;

fS22X = fS22x*cos(theta_g + nu_e*t) - fS22y*sin(theta_g + nu_e*t);
fS22Y = fS22x*sin(theta_g + nu_e*t) + fS22y*cos(theta_g + nu_e*t);
fS22Z = fS22z;

%Sun's gravity
magRo2 = Xo^2 + Yo^2 + Zo^2;
magRRo2 = (X - Xo)^2 + (Y - Yo)^2 + (Z - Zo)^2;
fSunX = -GMo*( (X - Xo)/(magRRo2^(3./2)) + Xo/(magRo2^(3./2)) );
fSunY = -GMo*( (Y - Yo)/(magRRo2^(3./2)) + Yo/(magRo2^(3./2)) );
fSunZ = -GMo*( (Z - Zo)/(magRRo2^(3./2)) + Zo/(magRo2^(3./2)) );

%Moon's gravity 
magRm2 = Xm^2 + Ym^2 + Zm^2;
magRRm2 = (X - Xm)^2 + (Y - Ym)^2 + (Z - Zm)^2;
fMoonX = -GMm*( (X - Xm)/(magRRm2^(3./2)) + Xm/(magRm2^(3./2)) );
fMoonY = -GMm*( (Y - Ym)/(magRRm2^(3./2)) + Ym/(magRm2^(3./2)) );
fMoonZ = -GMm*( (Z - Zm)/(magRRm2^(3./2)) + Zm/(magRm2^(3./2)) );

%Sun's radiation pressure
SRPterm = cram*PSRP*(alpha_o^2)/(magRRo2^(3./2));
fSRPX = SRPterm*(X - Xo);
fSRPY = SRPterm*(Y - Yo);
fSRPZ = SRPterm*(Z - Zo);

dx(1) = VX;
dx(2) = VY;
dx(3) = VZ;
dx(4) = fKepX + fJ2X + fC22X + fS22X + fSunX + fMoonX + fSRPX;
dx(5) = fKepY + fJ2Y + fC22Y + fS22Y + fSunY + fMoonY + fSRPY;
dx(6) = fKepZ + fJ2Z + fC22Z + fS22Z + fSunZ + fMoonZ + fSRPZ;

end