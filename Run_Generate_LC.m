

% Initial conditions
r0 = 860;
theta0 = 0.01;
phi10 = 0;
phi20 = 0;
rd0 = 0;
phi1d0 = 4.854616375452000e-04;
phi2d0 = 0;
tmax = 250000;

% Primary shape
a1 = 370; % fixed - do not change
c1 = 330;

% Secondary shape
a2 = 100;
a_over_c = 2.213;
c2 = a2/a_over_c;

% Impact parameter
Beta = 2.946;

profile clear
profile on

for i=1:5
% LC = Generate_LC(tmax,Beta,r0,theta0,phi10,phi20,rd0,phi1d0,phi2d0,a1,c1,a2,c2);
% LC = Generate_LC_opt(tmax,Beta,r0,theta0,phi10,phi20,rd0,phi1d0,phi2d0,a1,c1,a2,c2);
LC = Generate_LC_opt_mex(tmax,Beta,r0,theta0,phi10,phi20,rd0,phi1d0,phi2d0,a1,c1,a2,c2);
end

profile off
profile viewer

figure
plot(LC(:,1),LC(:,2),'bo','markersize',2,'markerfacecolor','b')
% hold on
% plot(LC_real(:,1),LC_real(:,2),'ro','markersize',2,'markerfacecolor','r')

