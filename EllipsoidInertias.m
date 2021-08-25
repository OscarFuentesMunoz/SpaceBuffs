function I = EllipsoidInertias(m,a,b,c)

I1 = (1/5)*m*(b^2+c^2);
I2 = (1/5)*m*(a^2+c^2);
I3 = (1/5)*m*(a^2+b^2);

I = diag([I1,I2,I3]);