function [x, y, z] = ellipsoid_csml(a,b,c)

n = 20;
theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi/2;
cosphi = cos(phi); cosphi(1) = 0; cosphi(n+1) = 0;
sintheta = sin(theta); sintheta(1) = 0; sintheta(n+1) = 0;

x = a*cosphi*cos(theta);
y = b*cosphi*sintheta;
z = c*sin(phi)*ones(1,n+1);


