% inputs:
%   cram_l: lower bound on cram
%   cram_u: upper bound on cram
%   debris: [1x7] array of debris measurement
%   type: string either 'moid' or 'oe'
function [cram_nom,sat,out] = debris_mother(cram_l,cram_u,debris,type)
S = load('all_sats.mat'); % sat_time, sat_pos, sat_vel, sat_oe
GMe = 3.986004407799724e+5;
days2sec = 3600*24;

% compute cartesian states for debris
[r,v] = orbel2rv(debris(2),debris(3),debris(4)*pi/180,debris(6)*pi/180,debris(7)*pi/180,M2f(debris(5)*pi/180,debris(3)),GMe);

ncram = 10; % number of cram values to try
cram = linspace(cram_l,cram_u,ncram);

% variables to record data for each cram value
min_dv = zeros(1,ncram);
min_moid = zeros(1,ncram);
sat_num = zeros(1,ncram);
t_min_moid = zeros(1,ncram);

% set up integration stuff
IC = [r;v];
opts = odeset('abstol',1e-16,'reltol',3e-14);
T = debris(1)*days2sec:-10*days2sec:-7305*days2sec;

% get only the necessary satellite info
t_id = find(S.sat_time==debris(1));
sat_pos = S.sat_pos(1:t_id,:,:);
sat_vel = S.sat_vel(1:t_id,:,:);
sat_oe = S.sat_oe(1:t_id,:,:);

% variables to track across time for each satellite (100 is the number of
% sats)
MOID = zeros(length(T),100);
oe_dist = zeros(length(T),100);
t_moid_id = zeros(100,1);
% tic

for k = 1:ncram
    [~,y] = ode113(@(t,x) Debris_EOM(t,x,cram(k)),T,IC,opts);
    for ii = 1:100 % look at each satellite
        for jj = 1:length(T) % look at each time
            oe_deb = rv2orbel(y(jj,1:3),y(jj,4:6),GMe);
            % store moid and oe distances
            MOID(jj,ii) = MOID_SDG_win(oe_deb([1 2 4 3 5]), sat_oe(jj,[1 2 4 3 5],ii));
            oe_dist(jj,ii) = norm(oe_deb(1:5)-sat_oe(jj,1:5,ii)');
        end
        % find time that minimum occurs for each satellite
        if type == 'moid'
            t_moid_id(ii) = find(MOID(:,ii)==min(MOID(:,ii)));
        elseif type == 'oe'
            t_moid_id(ii) = find(oe_dist(:,ii)==min(oe_dist(:,ii)));
        end
    end
    if type == 'moid'
        md = min(MOID,[],1); % minimum moid for each satellite
    elseif type == 'oe'
        md = min(oe_dist,[],1); % min oe distance for each sat
    end
    
    % look at the velocity difference between debris and satellites at the
    % instance in time that the minimum distance occurs
    dv = nan(1,100);
    for ii = 1:100
        if md(ii) < 3 % threshold to investigate satellites further
            dv_vec = y(t_moid_id(ii),4:6)-sat_vel(t_moid_id(ii),:,ii);
            dv(ii) = sqrt(sum(dv_vec.^2));
        end
    end
    % pick satellite that has the minimum dv
    min_dv_id = find(dv==min(dv));
    min_dv(k) = min(dv);
    min_moid(k) = md(min_dv_id); % could be moid or oe distance
    sat_num(k) = min_dv_id;
    t_min_moid(k) = T(min_dv_id)/days2sec;
end
% toc

% assign values to function outputs
sat = sat_num(min_moid==min(min_moid));
cram_nom = cram(min_moid==min(min_moid));
out.min_dv = min_dv;
out.min_moid = min_moid;
out.sat_num = sat_num;
out.t_min_moid = t_min_moid;
out.cram = cram;
end