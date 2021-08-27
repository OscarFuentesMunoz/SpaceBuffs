function [position_vector,velocity_vector] = orbel2rv(varargin)
% Orbel2RV - Orbital Elements to position and velocity
% Method 1:
% function [[position_vector,velocity_vector] = orbel2rv (orbel)
% orbel           = Structure of orbital elements: orbel.semi, orbel.ecc,
%                   orbel.incl, orbel.peri_arg, orbel.node_arg, orbel.nu.
% position_vector = Position vector in kilometers.
% velocity_vector = Velocity vector in kilometers per second.
% Method 2:
% function [[position_vector,velocity_vector] = orbel2rv (semimajor_axis,eccentricity,inclination,argument_of_perigee,right_ascension,nu)
% semimajor_axis      = Semimajor axis in kilometers.
% eccentricity        = Eccentrcity, Unitless.
% inclination         = Inclination, radians.
% right_ascention     = RAAN, radians.
% argument_of_perigee = Argument of perigee, radians.
% nu                  = True anomaly, radians.
% position_vector     = Position vector in kilometers.
% velocity_vector     = Velocity vector in kilometers per second.

% Adopted from Vallado, rewritten by Alex Sizemore
% Written 7/13/2016
% Modified 7/20/2016
% Dependencies:

if nargin == 1 || nargin == 2
    orbel = varargin{1};
    semimajor_axis = orbel.semi;
    eccentricity = orbel.ecc;
    inclination = orbel.incl;
    right_ascension = orbel.node_arg;
    argument_of_perigee = orbel.peri_arg;
    general_nu = orbel.nu;
    if nargin == 1
        mu = 398600;
    else
        mu = varargin{2};
    end
    
elseif nargin == 6 || nargin == 7
    semimajor_axis = varargin{1};
    eccentricity = varargin{2};
    inclination = varargin{3};
    right_ascension = varargin{4};
    argument_of_perigee = varargin{5};
    general_nu = varargin{6};
    if nargin == 6
        mu = 398600;
    else
        mu = varargin{7};
    end
else
    error('Either input a structure, or orbital elements individually.')
end

semilatus_rectum = semimajor_axis*(1-eccentricity^2);


%% Determine orbit type.
if eccentricity < 1e-9

    if inclination < 1e-9 || abs(inclination-pi) < 1e-9
        % Circular equatorial
        argument_of_perigee = 0;
        right_ascension = 0;
    else % Circular inclined
        argument_of_perigee= 0;
    end
    
else % Elliptical equatorial
    
    if inclination < 1e-9 || abs(inclination-pi) < 1e-9
        right_ascension = 0.0;
    end %if ( ( incl<1e-9) || (abs(incl-pi)<1e-9) )
    
end %if eccentricity < 1e-9

%% Form pqw position and velocity vectors
position_pqw = semilatus_rectum/(1+eccentricity*cos(general_nu))*[cos(general_nu) sin(general_nu) 0]';

if abs(semilatus_rectum) < 0.0001
    semilatus_rectum = 0.0001;
end

velocity_pqw = sqrt(mu)/sqrt(semilatus_rectum)*[-sin(general_nu) (eccentricity + cos(general_nu)) 0]';

%% Perform transformation to ijk
rotation_matrix = [[cos(right_ascension)*cos(argument_of_perigee)-sin(right_ascension)*sin(argument_of_perigee)*cos(inclination), -cos(right_ascension)*sin(argument_of_perigee)-sin(right_ascension)*cos(argument_of_perigee)*cos(inclination), sin(right_ascension)*sin(inclination)];
                  [sin(right_ascension)*cos(argument_of_perigee)+cos(right_ascension)*sin(argument_of_perigee)*cos(inclination), -sin(right_ascension)*sin(argument_of_perigee)+cos(right_ascension)*cos(argument_of_perigee)*cos(inclination), -cos(right_ascension)*sin(inclination)];
                  [sin(argument_of_perigee)*sin(inclination), cos(argument_of_perigee)*sin(inclination), cos(inclination)]];

              
position_vector = rotation_matrix*position_pqw;
velocity_vector = rotation_matrix*velocity_pqw;
