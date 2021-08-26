function [cram_nom,sat] = debris_mother(cram_l,cram_u,debris)
S = load('all_sats.mat'); % sat_time, sat_pos, sat_vel
GMe = 3.986004407799724e+5;
days2sec = 3600*24;

[r,v] = orbel2rv(debris(2:end),GMe);

ncram = 100;
cram = linspace(cram_l,cram_u,ncram);
min_d = zeros(1,ncram);
sat_num = zeros(1,ncram);

IC = [r;v];
opts = odeset('abstol',1e-16,'reltol',3e-14);
T = debris(1)*days2sec:10*days2sec:-7305*days2sec;

t_id = find(S.sat_time==debris(1));
sat_pos = S.sat_pos(1:t_id,:,:);

for k = 1:ncram
    pars.cram = cram(k);
    [~,y] = ode113(@(t,x) Debris_EOM(t,x,pars),T,IC,opts);
    re = sat_pos - repmat(y(:,1:3),1,1,100);
    de = sqrt(sum(re.^2,2));
    md = min(de,[],1);
    min_d(k) = min(md);
    sat_num(k) = find(md==min(md));
end
sat = sat_num(min_d==min(min_d));
cram_nom = cram(min_d==min(min_d));
end