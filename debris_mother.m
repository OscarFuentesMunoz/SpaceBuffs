function [cram_nom,sat,out] = debris_mother(cram_l,cram_u,debris)
S = load('all_sats.mat'); % sat_time, sat_pos, sat_vel, sat_oe
GMe = 3.986004407799724e+5;
days2sec = 3600*24;

[r,v] = orbel2rv(debris(2),debris(3),debris(4)*pi/180,debris(6)*pi/180,debris(7)*pi/180,M2f(debris(5)*pi/180,debris(3)),GMe);

ncram = 10;
cram = linspace(cram_l,cram_u,ncram);
min_dv = zeros(1,ncram);
min_moid = zeros(1,ncram);
sat_num = zeros(1,ncram);
t_min_moid = zeros(1,ncram);

IC = [r;v];
opts = odeset('abstol',1e-16,'reltol',3e-14);
T = debris(1)*days2sec:-10*days2sec:-7305*days2sec;

t_id = find(S.sat_time==debris(1));
sat_vel = S.sat_vel(1:t_id,:,:);
sat_oe = S.sat_oe(1:t_id,:,:);
MOID = zeros(length(T),100);
t_moid_id = zeros(100,1);
% tic
for k = 1:ncram
    [~,y] = ode113(@(t,x) Debris_EOM(t,x,cram(k)),T,IC,opts);
    for ii = 1:100
        for jj = 1:length(T)
            oe_deb = rv2orbel(y(jj,1:3),y(jj,4:6),GMe);
            MOID(jj,ii) = MOID_SDG_win(oe_deb([1 2 4 3 5]), sat_oe(jj,[1 2 4 3 5],ii));
        end
        t_moid_id(ii) = find(MOID(:,ii)==min(MOID(:,ii)));
    end
    md = min(MOID,[],1);
    dv = nan(1,100);
    for ii = 1:100
        if md(ii) < 5
            dv_vec = y(t_moid_id(ii),4:6)-sat_vel(t_moid_id(ii),:,ii);
            dv(ii) = sqrt(sum(dv_vec.^2));
        end
    end
    min_dv(k) = min(dv);
    min_dv_id = find(dv==min(dv));
    min_moid(k) = md(min_dv_id);
    sat_num(k) = min_dv_id;
    t_min_moid(k) = T(min_dv_id)/days2sec;
end
% toc
sat = sat_num(min_moid==min(min_moid));
cram_nom = cram(min_moid==min(min_moid));
out.min_dv = min_dv;
out.min_moid = min_moid;
out.sat_num = sat_num;
out.t_min_moid = t_min_moid;
out.cram = cram;
end