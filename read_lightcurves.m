testruns = strings(1,100);
for i = 1:200
    if i < 10
        num = ['00',num2str(i)];
    elseif i < 100
        num = ['0',num2str(i)];
    else
        num = i;
    end
    testruns(i) = num;
end

close all
beta = nan(200,1);

params = csvread('lightcurves/parameters.csv');

%% Loop

G = 6.6743e-11;
Mp = 3e11;
Ms = 1.3e9;
mu = G*(Ms+Mp);

Msc = 700;
Vsc = 5000;

Rp = 370;

trials = 2%1:200;
for r = trials
    
    testrun = char(testruns(r));
    
%% Old

LC = readtable(['lightcurves/lcvold',testrun,'.dat']);
LC = table2array(LC);
LC(:,1) = LC(:,1)*60; % convert to sec

% Column 1: time in minuts
% Column 2: normalized luminosity 

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

[autocorr,lags]=xcorr(LCfilled,ceil(max(LC(:,1))/2*fs),'coeff');

[pksh,lcsh] = findpeaks(autocorr,'minpeakdistance',3);
short_old = mean(diff(lcsh))/fs;

[pklg,lclg] = findpeaks(autocorr,'minpeakprominence',0.1,'minpeakheight',0.05);
long_old = mean(diff(lclg))/fs;

if length(trials) == 1
subplot(221)
plot(LC(:,1),LC(:,2),'bo','markerfacecolor','b','markersize',2)
subplot(223)
plot(lags/fs,autocorr)
hold on
pks = plot(lags(lcsh)/fs,pksh,'or',lags(lclg)/fs,pklg+0.05,'vk');
end

%% New

LC = readtable(['lightcurves/lcvnew',num2str(testrun),'.dat']);
LC = table2array(LC);
LC(:,1) = LC(:,1)*60; % convert to sec

% Column 1: time in minuts
% Column 2: normalized luminosity 

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

[autocorr,lags]=xcorr(LCfilled,ceil(max(LC(:,1))/2*fs),'coeff');


[pksh,lcsh] = findpeaks(autocorr,'minpeakdistance',3);
short_new = mean(diff(lcsh))/fs;

[pklg,lclg] = findpeaks(autocorr,'minpeakdistance',floor(long_old*fs*0.9));
long_new = mean(diff(lclg))/fs;
if long_new > 1.5*long_old
    [pklg,lclg] = findpeaks(-LCfilled,'minpeakheight',0.02);
    long_new = mean(diff(lclg))/fs;
end

if length(trials) == 1
subplot(222)
plot(LC(:,1),LC(:,2),'bo','markerfacecolor','b','markersize',2)
subplot(224)
plot(lags/fs,autocorr)
hold on
pks = plot(lags(lcsh)/fs,pksh,'or',lags(lclg)/fs,pklg+0.05,'vk');
end

%% get beta
Told = long_old;
Tnew = long_new;

a = ((Told^2*mu)/(4*pi^2))^(1/3);
anew = ((Tnew^2*mu)/(4*pi^2))^(1/3);
n = 2*pi/Told;
vold = a*2*pi/Told;
vnew = sqrt(mu*(2/a-1/anew));
dv = (a*n/vold)^2*(Tnew-Told)/(3*Told)*vold;
% dv = vnew-vold;
beta(r) = Ms*dv/(Msc*Vsc);
% if beta(r) > 3 || beta(r) < 1
%     beta(r) = 2*rand+1;
% end
end

%% Score

score = sqrt(1/200*sum((beta-params(:,1)).^2));