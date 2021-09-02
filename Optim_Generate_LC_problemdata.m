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

mA = 3e11; % Primary mass
mB = 1.3e9; % Secondary mass

a1 = 370;

% Constant density assumption --> Bounds on a2, c2
a2_min = a1*( mB/mA *(1-0.15)*(1.3) )^(1/3) ;
a2_max = a1*( mB/mA *(1-0.00)*(2.5) )^(1/3) ;


%% Read Lightcurves
params = csvread('lightcurves/parameters.csv');

for lcv = 5%:1:20

namedir  = 'lightcurves/';
namefile = [namedir 'lcvold' num2str(lcv,'%.3d') '.dat'];
LCold = load(namefile);
namefile = [namedir 'lcvnew' num2str(lcv,'%.3d') '.dat'];
LCnew = load(namefile);

figure(1); clf;
scatter(LCold(:,1),LCold(:,2),20,'k','filled'); hold on; grid on;
ylabel('Mag'); xlabel('time (s)')
scatter(LCnew(:,1),LCnew(:,2),20,'m','filled')

% plot(LCold([1 end],1),0.95*[1 1],'r--')

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

% 2. Initial guesses from lightcurves
[r0, Told] = r0_guess_autocorr(LCold,mu);
% r0 = 850;
theta0 = -20 *pi/180;
% theta0 = 90 *pi/180;

%-- Rest of ICs
phi10  = 0;
phi20  = 0;
rd0    = 0;
phi1d0 = 0; %4.854616375452000e-04;
phi2d0 = 0;

% x_true = [Beta; J2; r0; a_over_c; theta0];
x_true = [0; J2; r0; a_over_c; theta0];
c1 = a1*(1-J2);

% a2 = 70;
mA = 3e11; % Primary mass
mB = 1.3e9; % Secondary mass
a2 = a1*( mB/mA *(1-J2)*(a_over_c) )^(1/3);
% a2 = a2*0.4;
c2 = a2/a_over_c;

% -----------------------------------------
%% Generate LC of pre-impactor
c1 = a1*(1-J2);
LCtheta = Generate_LC_opt(tmax*1.5,0,r0,theta0,phi10,phi20,0.0,0.0,0,a1,c1,a2,c2);
LC      = LCtheta(1:(find(LCtheta(:,1)==tmax)),:);
ids = [];
for i=1:length(LCold(:,1))
    ids(i) = find( LC(:,1) == LCold(i,1)*60 );
end
theta0 = theta_minmse(LCtheta, LCold, LC, ids);
theta0 = theta0 + pi;

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

% Testing the modification of the LC to resemble the dataset
% mold_off = mean(LC(:,2))-mean(LCold(:,2));
% sold = std(LCold(:,2));
% LC2 = LC;
% LC2(:,2) = (LC(:,2)-mold_off)*sold + mean(LCold(:,2));

% scatter(LC2(:,1),LC2(:,2),6,'c','filled')
% scatter(LC2(ids,1),LC2(ids,2),10,'r','filled')
% title('Pre-impactor Lightcurves')


% -----------------------------------------
%% Generate LC of post-impactor
[r1, Tnew] = r0_guess_autocorr(LCnew,mu)
LCtheta    = Generate_LC_opt(tmax*1.5,0,r1,theta0,phi10,phi20,0.0,0.0,0,a1,c1,a2,c2);
ids = [];
for i=1:length(LCnew(:,1))
    ids(i) = find( LC(:,1) == LCnew(i,1)*60 );
end
theta0_new = theta_minmse(LCtheta, LCnew, LC, ids);

phi20  = 0.;

Beta     = params(lcv, 1);
J2       = params(lcv, 2);
a_over_c = params(lcv, 3);
% Choose to manipulate or not
% J2 = 0;
% a_over_c = 1.3;
c2 = a2/a_over_c;
c1 = a1*(1-J2);

% LC = Generate_LC_opt_mex(tmax,0,r0,theta0,phi10,phi20,0,0,0,a1,c1,a2,c2);
% LC = Generate_LC_opt(tmax,Beta,r1_beta,theta0,phi10,phi20,0.0,0.0,0,a1,c1,a2,c2);
LC = Generate_LC_opt(tmax,0,r1,theta0_new,phi10,phi20,0.0,0.0,0,a1,c1,a2,c2);

figure(3); clf;
% scatter(LCold(:,1)*60,LCold(:,2)-mean(LCold(:,2)),10,'k','filled'); hold on; grid on;
scatter(LCnew(:,1)*60,LCnew(:,2),20,'k','filled'); hold on; grid on;
ylabel('Mag'); xlabel('time (s)')

% scatter(LC(:,1),LC(:,2)-mean(LC(:,2)),6,'c','filled')
% scatter(LC(:,1),LC(:,2)-mean(LC(:,2)),10,'b','filled')
scatter(LC(:,1),LC(:,2),6,'c','filled')
scatter(LC(ids,1),LC(ids,2),10,'r','filled')
title('Post-impactor Lightcurves')

    
%% Analytical Beta
beta_est = beta_analytical( Told, Tnew )
[r1_beta,Tnew_beta] = Tnew_analytical( Told, Beta )


%% =====================================================
%% TESTING THE OPTIMIZATION
% TRIAL 1 - Use to find a/c, assume r0, theta0

% Assumptions
NL = 40; % Noise Level
spa= 0;  % Sun phase angle

lcv = 3%:1:20

namedir  = 'lightcurves/';
namefile = [namedir 'lcvold' num2str(lcv,'%.3d') '.dat'];
LCold = load(namefile);


% 1. INFORMED INITIAL GUESSES
%  a. Radius (period) given autocorrelation
[r0, Told] = r0_guess_autocorr(LCold,mu);
pfac   = 0.5;
theta0 = 0.5;

%  b. Find theta0 close to optimal
LCtheta = Generate_LC_fit_mex(tmax*1.5,0,r0,theta0,phi10,phi20,0.0,0.0,0,a1,c1,a2,c2,pfac,NL,spa);
LC      = LCtheta(1:(find(LCtheta(:,1)==tmax)),:);

ids = []; % Find Indexes of the LC where there is a measurement
for i=1:length(LCold(:,1))
    ids(i) = find( LC(:,1) == LCold(i,1)*60 );
end
[theta0_new,idmin] = theta_minmse(LCtheta, LCold, LC, ids);
theta0_new = wrapTo2Pi(theta0_new);

%  c. Rest of initial guesses
J2_guess     = 0.07;
r0_guess     = r0;
ac_guess     = 1.3;
theta0_guess = wrapTo2Pi(theta0_new);
pfac         = 0.5;

x_guess = [0;0;0;0;0];
x_guess(1) = J2_guess;
x_guess(2) = r0_guess;
x_guess(3) = ac_guess;
x_guess(4) = theta0_guess;
x_guess(5) = pfac;

% 2. OPTIMIZATION OF PARAMETERS
%  2.a. True Values
Beta     = params(lcv, 1);
J2       = params(lcv, 2);
a_over_c = params(lcv, 3);

% x_true = [Beta; J2; r0; a_over_c; theta0];
x_true = [J2; r0; a_over_c; theta0_new; pfac];

%  2.b. Optimization setup
x0 = x_guess;

% Lower bounds and Upper bounds - OLD SETUP
% %    Beta, J2,   r0,   a/c, theta0
% lb = [1;   0;    500;  1.3];% 0     ];
% ub = [3;   0.15; 1500; 2.5];% 2*pi  ];

%     J2,   r0,  a/c,  theta0, pfac
lb = [0;    r0;  1.3;  theta0_new;   0.4 ];
ub = [0.15; r0;  2.5;  theta0_new;   0.5 ];

% Upper bounds
options = optimoptions('fmincon',...
    'Display','iter',...
    'Algorithm','active-set',...
    'TypicalX',[0.1; 800; 2; pi; 0.5],...; pi],...
    'FunctionTolerance', 1e-6,...
    'HessianApproximation','lbfgs');

[x_min,loss_min] = fmincon(@(x) loss_function(LCold(:,2),x,tmax, ids,NL,spa), x0, [],[],[],[], lb,ub,[],options);

% Plot optimized final LC
[~, LCmin] = loss_function( LCold(:,2), x_min, tmax, ids,NL,spa );
Params = {'J2';'r_0';'a/c';'theta0';'pfac'};
Results = table(Params,x_true,x_guess,x_min)

% end

F = figure(23); clf;

scatter(LCold(:,1)*60,LCold(:,2),20,'k','filled'); hold on; grid on;
ylabel('Mag'); xlabel('time (s)')

% LCguess = LCtheta(idmin:(idmin+length(LC(:,1))-1),2);
[err_guess, LCguess] = loss_function( LCold(:,2), x_guess, tmax, ids,NL,spa );

scatter(LCguess(:,1),LCguess(:,2),6,'c','filled')
scatter(LCguess(ids,1),LCguess(ids,2),10,'r','filled')

title('Pre-impactor Lightcurves')
scatter(LCmin(:,1),LCmin(:,2),6,'g','filled'); hold on; grid on;
legend('Real',...
    ['Guess - L=' num2str(err_guess, '%.2e')],...
    ['Optim - L=' num2str(loss_min,  '%.2e')])


%% Check with current knowledge...
% Beta found using read_lightcurves
acv = rand(100,1)*(1.2) + 1.3;
score = score_fxn( params, beta(1:100), acv )

figure(6)
plot( 1:100, beta(1:100)-params(1:100,1), 'r.' ); hold on
% plot( 1:100, params(1:100,1), 'g.')

%% FUNCTIONS DEFINITION

%---------------------------------------
% Optimization Loss Function
function [err, LC] = loss_function( LCgen, x, tmax, ids, NL, spa )

    % Constants of the problem
    a1 = 370;
    mA = 3e11; % Primary mass
    mB = 1.3e9; % Secondary mass
    
    % Variables of the current guess
    J2     = x(1);
    r0     = x(2);
    a_ov_c = x(3);
    theta0 = x(4);
    pfac   = x(5);
    % Beta   = x(1);

    % Variables set as constants
    phi10 = 0;
    phi20 = 0;
    rd0 = 0;
    phi1d0 = 0;%4.854616375452000e-04;
    phi2d0 = 0;
    
    % Mapping to inputs of LG gen
    a2 = a1*( mB/mA *(1-J2)*(a_ov_c) )^(1/3);
    c2 = a2/a_ov_c;
    c1 = a1*(1-J2);
    
    LC = Generate_LC_fit_mex(tmax,0,r0,theta0,phi10,phi20,rd0,phi1d0,phi2d0,a1,c1,a2,c2,pfac,NL,spa);

    % err = immse( LCgen, LC(ids,2) );
    err = norm( LCgen - LC(ids,2) )^2;
    
end


%---------------------------------------
% Generate radius and period estimates given autocorrelation
function [r0, Told] = r0_guess_autocorr(LC,mu)
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

%---------------------------------------
% Finding Beta analytically
function beta = beta_analytical( Told, Tnew )
% Code taken from read_lightcurves.m

G = 6.6743e-11;
Mp = 3e11;
Ms = 1.3e9;
mu = G*(Ms+Mp);

Msc = 700;
Vsc = 5000;

a = ((Told^2*mu)/(4*pi^2))^(1/3);
anew = ((Tnew^2*mu)/(4*pi^2))^(1/3);
n = 2*pi/Told;
vold = a*2*pi/Told;
vnew = sqrt(mu*(2/a-1/anew));
dv = (a*n/vold)^2*(Tnew-Told)/(3*Told)*vold;
% dv = vnew-vold;
beta = Ms*dv/(Msc*Vsc);
% if beta(r) > 3 || beta(r) < 1
%     beta(r) = 2*rand+1;
% end

end


%---------------------------------------
% Learning from test: period fxn of Beta
function [r1,Tnew] = Tnew_analytical( Told, Beta )
G = 6.6743e-11;
Mp = 3e11;
Ms = 1.3e9;
mu = G*(Ms+Mp);

Msc = 700;
Vsc = 5000;
dV = Beta*Vsc*Msc/Ms ;

a = ((Told^2*mu)/(4*pi^2))^(1/3);
n = 2*pi/Told;
vold = a*n;

Tnew = Told*(1 + 3*dV/vold);
r1 = ((Tnew^2*mu)/(4*pi^2))^(1/3);

end


%---------------------------------------
% Sweep over theta to find initial guess
function [theta0_est, idmin] = theta_minmse(LCtheta, LCnew, LC, ids)

    % Array of theta's from integration
    LCcut = LC;
    ntv = size(LC,1);
    nex = size(LCtheta,1) - size(LC,1); % Number of additional periods
    tv  = LCtheta(1:nex,3);
    errori = nan(nex,1);
    for t=1:nex
        ordr = t:(ntv+t-1);
        LCcut(:,2) = LCtheta(ordr,2);
        errori(t) = immse( LCcut(ids,2), LCnew(:,2) );
    end
    
    [~,idmin] = min(errori);
    theta0_est = tv(idmin);
    
%     figure(5);
%     plot( tv, errori )
%     grid on
%     xlabel('\theta (rad)')
%     ylabel('MSE')

end

%---------------------------------------
% Compute score
function score = score_fxn( params, Betav, acv )

    N = length(Betav);
    score = .5*sqrt( sum((params(1:N,1)-Betav).^2) /N ) + ...
          + .5*sqrt( sum((params(1:N,3)-acv  ).^2) /N );

end


