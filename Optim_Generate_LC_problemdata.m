%% Milani Challenge
% Optimizing Binary parameters given light-curve
% Trying with the problem data
clc
clear
close all

%% Constants
G = 6.6743e-11;
Mp = 3e11;
Ms = 1.3e9;
mu = G*(Ms+Mp);


%% Read Lightcurves
params = csvread('lightcurves/parameters.csv');

for lcv = 3%:1:20

namedir  = 'lightcurves/';
namefile = [namedir 'lcvold' num2str(lcv,'%.3d') '.dat'];
LCold = load(namefile);
namefile = [namedir 'lcvnew' num2str(lcv,'%.3d') '.dat'];
LCnew = load(namefile);

figure(1); clf;
scatter(LCold(:,1),LCold(:,2),20,'k','filled'); hold on; grid on;
ylabel('Mag'); xlabel('time (s)')
scatter(LCnew(:,1),LCnew(:,2),20,'m','filled')

plot(LCold([1 end],1),0.95*[1 1],'r--')

legend('pre-impact','post-impact')

pause(1)
end

%% Generate Initial guesses
% 1. Exact problem solutions
Beta     = params(lcv, 1);
J2       = params(lcv, 2);
a_over_c = params(lcv, 3);
% J2 = 0;

tmax = LCold(end,1)*60;

% Initial conditions
r0 = 1356;
theta0 = -20 *pi/180;
% theta0 = 90 *pi/180;

% 2. Initial guesses from lightcurves
r0 = r0_guess_autocorr(LCold,mu);


%%
% r0 = 850;

phi10  = 0;
phi20  = 0;
rd0    = 0;
phi1d0 = 0; %4.854616375452000e-04;
phi2d0 = 0;

% x_true = [Beta; J2; r0; a_over_c; theta0];
x_true = [0; J2; r0; a_over_c; theta0];

a1 = 370;
c1 = a1*(1-J2);

% a2 = 70;
mA = 3e11; % Primary mass
mB = 1.3e9; % Secondary mass
a2 = a1*( mB/mA *(1-J2)*(a_over_c) )^(1/3);
% a2 = a2*0.4;
c2 = a2/a_over_c;

%-----------------------------------------
% Generate LC of pre-impactor
% LC = Generate_LC_opt_mex(tmax,0,r0,theta0,phi10,phi20,0,0,0,a1,c1,a2,c2);
LC = Generate_LC_opt(tmax,0,r0,theta0,phi10,phi20,0,0,0,a1,c1,a2,c2);

figure(2); clf;
% scatter(LCold(:,1)*60,LCold(:,2)-mean(LCold(:,2)),10,'k','filled'); hold on; grid on;
scatter(LCold(:,1)*60,LCold(:,2),20,'k','filled'); hold on; grid on;
ylabel('Mag'); xlabel('time (s)')

% scatter(LC(:,1),LC(:,2)-mean(LC(:,2)),6,'c','filled')
% scatter(LC(:,1),LC(:,2)-mean(LC(:,2)),10,'b','filled')
scatter(LC(:,1),LC(:,2),6,'c','filled')
for i=1:length(LCold(:,1))
    ids(i) = find( LC(:,1) == LCold(i,1)*60 );
end
scatter(LC(ids,1),LC(ids,2),10,'r','filled')
title('Pre-impactor Lightcurves')

mold_off = mean(LC(:,2))-mean(LCold(:,2));
sold = std(LCold(:,2));
LC2 = LC;
LC2(:,2) = (LC(:,2)-mold_off)*sold + mean(LCold(:,2));

% scatter(LC2(:,1),LC2(:,2),6,'c','filled')
% scatter(LC2(ids,1),LC2(ids,2),10,'r','filled')
% title('Pre-impactor Lightcurves')

%-----------------------------------------
%% Generate LC of post-impactor
phi20 = 0.1;

% LC = Generate_LC_opt_mex(tmax,0,r0,theta0,phi10,phi20,0,0,0,a1,c1,a2,c2);
LC = Generate_LC_opt(tmax,Beta,r0,theta0,phi10,phi20,0.0,0.0,0,a1,c1,a2,c2);

figure(3); clf;
% scatter(LCold(:,1)*60,LCold(:,2)-mean(LCold(:,2)),10,'k','filled'); hold on; grid on;
scatter(LCnew(:,1)*60,LCnew(:,2),20,'k','filled'); hold on; grid on;
ylabel('Mag'); xlabel('time (s)')

% scatter(LC(:,1),LC(:,2)-mean(LC(:,2)),6,'c','filled')
% scatter(LC(:,1),LC(:,2)-mean(LC(:,2)),10,'b','filled')
scatter(LC(:,1),LC(:,2),6,'c','filled')
for i=1:length(LCnew(:,1))
    ids(i) = find( LC(:,1) == LCnew(i,1)*60 );
end
scatter(LC(ids,1),LC(ids,2),10,'r','filled')
title('Post-impactor Lightcurves')

%%
mA = 3e11; % Primary mass
mB = 1.3e9; % Secondary mass

% Constant density assumption
a2_min = a1*( mB/mA *(1-0.15)*(1.3) )^(1/3)
a2_max = a1*( mB/mA *(1-0.00)*(2.5) )^(1/3)


%% Functions
function r0 = r0_guess_autocorr(LC,mu)
% Code taken from read_lightcurves.m

LC(:,1) = LC(:,1)*60; % convert to sec

fs = 1/min(diff(LC(:,1)));
LCnorm = LC(:,2)-mean(LC(:,2));
tvec = LC(1,1):1/fs:LC(end,1);
LCfilled = nan(1,length(tvec));
for i = 1:length(tvec)
    idx = find(LC(:,1)==tvec(i));
    if ~isempty(idx)
        LCfilled(i) = LCnorm(idx);
    else
        LCfilled(i) = 0;
    end
end

[autocorr,~]=xcorr(LCfilled,ceil(max(LC(:,1))/2*fs),'coeff');

[~,lclg] = findpeaks(autocorr,'minpeakprominence',0.1,'minpeakheight',0.05);
long_old = mean(diff(lclg))/fs;
Told = long_old*2;

r0 = ((Told^2*mu)/(4*pi^2))^(1/3);

end