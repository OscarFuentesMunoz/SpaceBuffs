function Xd = planar_f2bp(t,X,mu,nu,I1,I2)

I1x = I1(1,1);
I1y = I1(2,2);
I1z = I1(3,3);
I2x = I2(1,1);
I2y = I2(2,2);
I2z = I2(3,3);

r = X(1);
theta = X(2);
phi1 = X(3);
phi2 = X(4);
rd = X(5);
thetad = X(6);
phi1d = X(7);
phi2d = X(8);

Ur = -mu/r^2 - 3*mu/(2*r^4)*(I1x+I1y+I1z+I2x+I2y+I2z) +...
    9*mu/(4*r^4)*(I1x+I1y-cos(2*phi1)*(I1y-I1x)+I2x+I2y-cos(2*phi2)*(I2y-I2x));
Uphi1 = -3*mu/(2*r^3)*(sin(2*phi1)*(I1y-I1x));
Uphi2 = -3*mu/(2*r^3)*(sin(2*phi2)*(I2y-I2x));

rdd = thetad^2*r + Ur;
thetadd = -Uphi1/r^2 - Uphi2/r^2 - 2*rd*thetad/r;
phi1dd = (1+(1-nu)*r^2/I1z)*Uphi1/r^2 + Uphi2/r^2 + 2*rd*thetad/r;
phi2dd = (1+nu*r^2/I2z)*Uphi2/r^2 + Uphi1/r^2 + 2*rd*thetad/r;

Xd = [rd;thetad;phi1d;phi2d;rdd;thetadd;phi1dd;phi2dd];
end
