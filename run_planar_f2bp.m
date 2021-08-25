options = odeset('Reltol',1e-12,'AbsTol',1e-12);

mA = 3e11; % Primary mass
mB = 1.3e9; % Secondary mass
mu = 6.67430e-11*(mA+mB);
nu = mA/(mA+mB);
mDART = 700;
vDART = 5000;

Beta = 2;   % Define beta

% Initial conditions
r0 = 920;
theta0 = 0;
phi10 = 0;
phi20 = 0;
rd0 = 0;
% T = 2*pi*sqrt(r0^3/mu);
% thetad0 = 2*pi/T;
phi1d0 = 4.854616375452000e-04;
phi2d0 = 0;

% Physical parameters
a1 = 370;
c1 = 370;
a2 = 100;
c2 = 80;
I1 = EllipsoidInertias(mA,a1,a1,c1);
I2 = EllipsoidInertias(mB,a2,c2,c2);
r = r0/a1;
I1z = I1(3,3)/(mA*a1^2);
Is = I1(1,1)/(mA*a1^2);
I2z = I2(3,3)/(mB*a1^2);
I2x = I2(1,1)/(mB*a1^2);
I2y = I2(2,2)/(mB*a1^2);
n = sqrt(mu/a1^3); % Mean motion

% Get initial theta_dot of system for equilibrium circular orbit
C2 = -2*I2x+I2y+I2z;
K_eq = sqrt(nu^2/r^5*(r^6+(2*I2z/nu+(3/2)*(C2+I1z-Is))*r^4+(I2z^2/nu^2+3*I2z/nu*(C2+I1z-Is))*r^2+(3/2)*I2z^2/nu^2*(C2+I1z-Is)));

% Perturb system with beta
dK = Beta*nu*(mDART/mB)*(vDART/a1/n)*r;
K = K_eq + dK;

% System spin rate after impact
thetad0 = (K-I2z*phi2d0)/(I2z+nu*r^2);
thetad0 = thetad0*n;

X0 = [r0,theta0,phi10,phi20,rd0,thetad0,phi1d0,phi2d0];
I1 = I1/mA;
I2 = I2/mB;

tspan = 0:600:200000;

[t,X_out] = ode113(@(t,X) planar_f2bp(t,X,mu,nu,I1,I2),tspan,X0,options);

theta = X_out(:,2);
phi2 = X_out(:,4);


% %% Vectors
% R = nan(3,length(t));V = R;
% a_ana = nan(1,length(t));e_ana = a_ana;w_ana = a_ana;f_ana = a_ana;
% 
% for j = 1:length(t)
% r = X_out(j,1);
% th = X_out(j,2);
% rd = X_out(j,5);
% thetad = X_out(j,6);
% xhat = [1;0;0];
% yhat = [0;1;0];
% 
% rhat = cos(th)*xhat + sin(th)*yhat;
% theta_hat = -sin(th)*xhat + cos(th)*yhat;
% 
% R(:,j) = r*rhat;
% V(:,j) = rd*rhat + r*thetad*theta_hat;
% 
% 
% %% Orbit elements
% [a_ana(j),e_ana(j),~,~,w_ana(j),f_ana(j),~,~] = Cart2Kepler(R(:,j),V(:,j),mu);
% end

for i = 1:length(t)
    if X_out(i,2) > 2*pi
        X_out(i:end,2) = X_out(i:end,2) - 2*pi;
    end
end

% [~,loc] = findpeaks(X_out(:,2));
% period = mean(diff(t(loc)));
% a_calc = ((period^2*mu)/(4*pi^2))^(1/3);

%% Figures
figure
subplot(221)
plot(t,X_out(:,1))
title('$r$','interpreter','latex','fontsize',15)
subplot(222)
plot(t,X_out(:,2))
title('$\theta$','interpreter','latex','fontsize',15)
subplot(223)
plot(t,X_out(:,3))
title('$\phi 1$','interpreter','latex','fontsize',15)
subplot(224)
plot(t,X_out(:,4))
title('$\phi 2$','interpreter','latex','fontsize',15)

figure
subplot(221)
plot(t,X_out(:,5))
title('$\dot{r}$','interpreter','latex','fontsize',15)
subplot(222)
plot(t,X_out(:,6))
title('$\dot{\theta}$','interpreter','latex','fontsize',15)
subplot(223)
plot(t,X_out(:,7))
title('$\dot{\phi 1}$','interpreter','latex','fontsize',15)
subplot(224)
plot(t,X_out(:,8))
title('$\dot{\phi 2}$','interpreter','latex','fontsize',15)

% figure
% plot(t,a_ana)
% hold on
% plot(t,a_calc*ones(1,length(t)))

%%
% Told = 39000;
% Tnew = 43800;
% a = ((Told^2*mu)/(4*pi^2))^(1/3);
% anew = ((Tnew^2*mu)/(4*pi^2))^(1/3);
% n = 2*pi/Told;
% vold = a*2*pi/Told;
% vnew = sqrt(mu*(2/a-1/anew));
% dv = (a*n/vold)^2*(Tnew-Told)/(3*Told)*vold;
% dv = vnew-vold;
% beta = Ms*dv/(Msc*Vsc);