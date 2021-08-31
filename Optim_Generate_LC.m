%% Milani Challenge
% Optimizing Binary parameters given light-curve
% Toy problem
clc
clear

%%
% Initial conditions
r0 = 860;
theta0 = 0.5;%0.01;
phi10 = 0;
phi20 = 0;
rd0 = 0;
phi1d0 = 4.854616375452000e-04;
phi2d0 = 0;
tmax = 250000;

% Primary shape
a1 = 370; % fixed - do not change
c1 = 330;
J2 = 1 - c1/a1;

% Secondary shape
a2 = 100; % assumed, for now
a_over_c = 2.213;
c2 = a2/a_over_c;

% Impact parameter
Beta = 2.946;

% profile clear
% profile on

for i=1%:5
% LC = Generate_LC(tmax,Beta,r0,theta0,phi10,phi20,rd0,phi1d0,phi2d0,a1,c1,a2,c2);
LC = Generate_LC_opt(tmax,Beta,r0,theta0,phi10,phi20,rd0,phi1d0,phi2d0,a1,c1,a2,c2);
LC = Generate_LC_opt_mex(tmax,Beta,r0,theta0,phi10,phi20,rd0,phi1d0,phi2d0,a1,c1,a2,c2);
end

% profile off
% profile viewer

% figure
% plot(LC(:,1),LC(:,2),'bo','markersize',2,'markerfacecolor','b')
% hold on
% plot(LC_real(:,1),LC_real(:,2),'ro','markersize',2,'markerfacecolor','r')


%% Implement Loss function
figure(10); clf;

% Truth LC
x_true = [Beta; J2; r0; a_over_c; theta0];
% [err, LC] = loss_function( LCgen, x, tmax )

LCreal = Generate_LC_opt_mex(tmax,Beta,r0,theta0,phi10,phi20,rd0,phi1d0,phi2d0,a1,c1,a2,c2);
plot(LCreal(:,1),LCreal(:,2),'k--'); hold on; grid on;
ylabel('Mag'); xlabel('time (s)')

% Generate noisy measurements
NoiseLevel = 100;
LCmeas = Generate_LC_noisy(tmax,Beta,r0,theta0,phi10,phi20,rd0,phi1d0,phi2d0,a1,c1,a2,c2,NoiseLevel);
scatter(LCmeas(:,1),LCmeas(:,2),10,'b','filled')

% Minimum Loss function - nominal vs. noisy
err_nom = immse( LCmeas(:,2),LCreal(:,2) ) %MSE of the measured light-curve and the real

% Perturb initial conditions and plot
% x0 = [Beta,r0,theta0,phi10,phi20,rd0,phi1d0,phi2d0,a1,c1,a2,c2];
Beta_guess   = Beta*0.98;
J2_guess     = 0.05;
r0_guess     = r0*1.02;
ac_guess     = 2.15;
theta0_guess = 10*pi/180;

x_guess = x_true;
x_guess(1) = Beta_guess;
x_guess(2) = J2_guess;
x_guess(3) = r0_guess;
x_guess(4) = ac_guess;
x_guess(5) = theta0_guess;

[err_guess, LCguess] = loss_function( LCmeas(:,2), x_guess, tmax );
scatter(LCguess(:,1),LCguess(:,2),6,'r','filled'); hold on; grid on;
% title(['Error of guess = ' num2str(err_guess,'%.3e')])


%% Optimizer
x0 = x_guess;

% Lower bounds and Upper bounds
%    Beta, J2,   r0,   a/c, theta0
lb = [1;   0;    500;  1.3; 0     ];
ub = [3;   0.15; 1500; 2.5; 2*pi  ];

% Upper bounds
options = optimoptions('fmincon',...
    'Display','iter',...
    'Algorithm','interior-point',...
    'TypicalX',[2; 0.1; 80; 2; pi],...
    'FunctionTolerance', 1e-6,...
    'HessianApproximation','lbfgs');

[x_min,loss_min] = fmincon(@(x) loss_function(LCmeas(:,2),x,tmax), x0, [],[],[],[], lb,ub,[],options);
error_init = (x_guess - x_true);
error_min  = (x_min  - x_true);

% Plot optimized final LC
[~, LCmin] = loss_function( LCmeas(:,2), x_min, tmax );
scatter(LCmin(:,1),LCmin(:,2),6,'g','filled'); hold on; grid on;
% title(['Error of guess = ' num2str(loss_min,'%.3e')])

legend('Real',...
    ['Noisy - L=' num2str(err_nom,   '%.2e')],...
    ['Pertu - L=' num2str(err_guess, '%.2e')],...
    ['Optim - L=' num2str(loss_min,  '%.2e')])

Params = {'Beta';'J2';'r_0';'a/c';'theta0'};
Results = table(Params,x_true,x_guess,x_min);


%% Functions
function [err, LC] = loss_function( LCgen, x, tmax )

    % Constants of the problem
    a1 = 370;
    a2 = 100; % Assumed, for now
    
    % Variables of the current guess
    Beta   = x(1);
    J2     = x(2);
    r0     = x(3);
    a_ov_c = x(4);
    theta0 = x(5);
    
    % Variables set as constants
    phi10 = 0;
    phi20 = 0;
    rd0 = 0;
    phi1d0 = 4.854616375452000e-04;
    phi2d0 = 0;
    
    % Mapping to inputs of LG gen
    c1 = a1*(1-J2);
    c2 = a2/a_ov_c;
    
    LC = Generate_LC_opt_mex(tmax,Beta,r0,theta0,phi10,phi20,rd0,phi1d0,phi2d0,a1,c1,a2,c2);

    err = immse( LCgen, LC(:,2) );
end
