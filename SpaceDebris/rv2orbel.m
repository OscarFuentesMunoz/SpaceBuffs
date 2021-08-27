% rv2orbel.m converts from cartesian coordinates to the orbital elements
%   a e i raan aop true_anomaly
%
% function orb_elem = rv2orbel(r0vect,v0vect,mu)
%
% Inputs:
%   r0vect: 3 element vector of x-y-z position coordinates
%   v0vect: 3 element vector of x-y-z velocity coordinates
%   mu: gravitational parameter of system
%
% Outputs:
%   orb_elem: [6x1] vector of orbital elements [a,e,i,raan,aop,nu] with
%       angles in radians
%
% Dependencies: none
%
% Author: David Lujan, david.lujan@colorado.edu

function orb_elem = rv2orbel(r0vect,v0vect,mu)
r0 = norm(r0vect);
v0 = norm(v0vect);
evect = 1/mu*cross(v0vect,cross(r0vect,v0vect)) - r0vect/r0;
e = norm(evect);
hvect = cross(r0vect,v0vect);
h = norm(hvect);
inc = acos(dot(hvect,[0,0,1])/h);
if inc == pi
    inc = 0;
end
ainv = -2/mu*(v0^2/2-mu/r0);
if abs(ainv) > 1e-15
    % elliptic and hyperbolic
    a = 1/ainv;
    if inc < 10^(-9) && e < 10^(-9)
        % Circular Equatorial Orbit
        % Define raan and aop to equal 0, so any angular displacement from
        % the x axis is contained within true anomaly.
        true_anom = atan2(r0vect(2),r0vect(1));
        if true_anom < 0
            true_anom = 2*pi + true_anom;
        end
        orb_elem = [a,0,0,0,0,true_anom]';
        return
    elseif inc < 10^(-9) && e > 10^(-9)
        % Elliptic Equatorial Orbit
        % Define raan to be angle between x axis and eccentricity vector, aop
        % is set equal to 0, and true anomaly is the angular displacement from
        % the eccentricity vector.
        nhat = evect./e;
        raan = atan2(nhat(2),nhat(1));
        if raan < 0
            raan = 2*pi + raan;
        end
        % angle b/t x axis and position vector
        f = atan2(r0vect(2),r0vect(1));
        if f < 0
            f = 2*pi + f;
        end
        % difference between f and raan gives true anomaly
        true_anom = f - raan;
        if true_anom < 0
            true_anom = 2*pi + true_anom;
        end
        orb_elem = [a,e,0,raan,0,true_anom]';
        return
    elseif inc > 10^(-9) && e < 10^(-9)
        % Circular Inclined Orbit
        % define eccentricity vector to point in direction of ascending node,
        % aop to be equal to 0, and true anomaly to be the angular displacement
        % b/t the ascending node and position vector.
        nhat = cross([0,0,1],hvect./h)./norm(cross([0,0,1],hvect./h));
        raan = atan2(nhat(2),nhat(1));
        if raan < 0
            raan = 2*pi + raan;
        end
        true_anom = acos(dot(nhat,r0vect./r0));
        if r0vect(3) < 0
            true_anom = 2*pi - true_anom;
        end
        orb_elem = [a,0,inc,raan,0,true_anom]';
        return
    else
        % Elliptic Inclined Orbit
        % Everything is defined normally. There are no singularities.
        nhat = cross([0,0,1],hvect./h)./norm(cross([0,0,1],hvect./h));
        raan = atan2(nhat(2),nhat(1));
        if raan < 0
            raan = 2*pi + raan;
        end
        aop = acos(dot(nhat,evect)./e);
        if evect(3) < 0
            aop = 2*pi - aop;
        end
        true_anom = acos(dot(evect./e,r0vect./r0));
        if dot(r0vect,v0vect) < 0
            true_anom = 2*pi - true_anom;
        end
        orb_elem = [a,e,inc,raan,aop,true_anom]';
        return
    end

elseif ainv < -1e-5
    %parabolic - who knows?
    x = zeros(6,1);
end

end