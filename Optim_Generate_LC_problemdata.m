%% Milani Challenge
% Optimizing Binary parameters given light-curve
% Trying with the problem data
clc
clear


%% Read Lightcurves
params = csvread('lightcurves/parameters.csv');

lcv = 3;
namedir  = 'lightcurves/';
namefile = [namedir 'lcvold' num2str(lcv,'%.3d') '.dat'];

LCold = load(namefile);

figure(1); clf;
scatter(LCold(:,1),LCold(:,2),10,'k','filled'); hold on; grid on;
ylabel('Mag'); xlabel('time (s)')


%% Generate Initial guesses
Beta     = params(lcv, 1);
J2       = params(lcv, 2);
a_over_c = params(lcv, 3);

tmax = LCold(end,1)*60;

% Initial conditions
r0 = 830;
theta0 = 0.05;%0.01;
phi10  = 0;
phi20  = 0;
rd0    = 0;
phi1d0 = 4.854616375452000e-04;
phi2d0 = 0;

% x_true = [Beta; J2; r0; a_over_c; theta0];
x_true = [0; J2; r0; a_over_c; theta0];

a1 = 370;
c1 = a1*(1-J2);

a2 = 130;
c2 = a2/a_over_c;
LC = Generate_LC_opt_mex(tmax,0,r0,theta0,phi10,phi20,0,0,0,a1,c1,a2,c2);

figure(1); clf;
% scatter(LCold(:,1)*60,LCold(:,2)-mean(LCold(:,2)),10,'k','filled'); hold on; grid on;
scatter(LCold(:,1)*60,LCold(:,2),10,'k','filled'); hold on; grid on;
ylabel('Mag'); xlabel('time (s)')

% scatter(LC(:,1),LC(:,2)-mean(LC(:,2)),6,'c','filled')
% scatter(LC(:,1),LC(:,2)-mean(LC(:,2)),10,'b','filled')
scatter(LC(:,1),LC(:,2),6,'c','filled')
for i=1:length(LCold(:,1))
    ids(i) = find( LC(:,1) == LCold(i,1)*60 );
end
scatter(LC(ids,1),LC(ids,2),10,'r','filled')



