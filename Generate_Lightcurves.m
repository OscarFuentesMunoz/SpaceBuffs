% First run run_planar_f2bp.m to get the time history of binary dynamics

e_vec = [0,1,0];    % Vector to Earth
s_vec = [0,1,0];    % Vector to sun

[ax,ay,az] = ellipsoid(0,0,0,a1,a1,c1);
[bx,by,bz] = ellipsoid(0,0,0,a2,c2,c2);

% bxc = reshape(bx,21^2,1);
% byc = reshape(by,21^2,1);
% bzc = reshape(bz,21^2,1);
% b_verts = [bxc,byc,bzc]';

S = surf(bx,by,bz);
axis equal
sec = surf2patch(S);
close all

%% Primary brightness
% P = surf(ax,ay,az);
% prim = surf2patch(P);
% 
% plight = 0;
% for f = 1:length(prim.faces)
%     % Get surface normal
%     v1 = prim.vertices(prim.faces(f,1),:)';
%     v2 = prim.vertices(prim.faces(f,2),:)';
%     v3 = prim.vertices(prim.faces(f,3),:)';
%     v4 = prim.vertices(prim.faces(f,4),:)';
%     vec1 = v2-v1;
%     vec2 = v3-v1;
%     if norm(vec1) == 0
%         vec1 = v4-v1;
%     elseif norm(vec2) == 0
%         vec2 = v4-v1;
%     end
%     surf_norm = cross(vec1,vec2)/norm(cross(vec1,vec2));
%     
%     % Make sure this is an outer facing normal
%     if dot(surf_norm,v1) < 0
%         surf_norm = -surf_norm;
%     end
%     
%     % Check if lit
%     if dot(surf_norm,s_vec) > 0
%         % Check if visible from Earth
%         if dot(surf_norm,e_vec) > 0
%             plight = plight + dot(surf_norm,e_vec)*norm(vec1)*norm(vec2);
%         end
%     end
% end
plight = 2*8.477442968647369e+05;

%% quasi-Lightcurves
mag = zeros(1,length(t));
for i = 1:length(t)
    
    M = rotation(theta(i)+phi2(i),3);
    NB = M';
    
    vis = 0;
    for f = 1:length(sec.faces)
        % Get surface normal
        v1 = R(:,i) + NB*sec.vertices(sec.faces(f,1),:)';
        v2 = R(:,i) + NB*sec.vertices(sec.faces(f,2),:)';
        v3 = R(:,i) + NB*sec.vertices(sec.faces(f,3),:)';
        v4 = R(:,i) + NB*sec.vertices(sec.faces(f,4),:)';
        vc = [(v1(1)+v2(1)+v3(1)+v4(1))/4;
              (v1(2)+v2(2)+v3(2)+v4(2))/4;
              (v1(3)+v2(3)+v3(3)+v4(3))/4];
        vec1 = v2-v1;
        vec2 = v3-v1;
        if norm(vec1) == 0
            vec1 = v4-v1;
        elseif norm(vec2) == 0
            vec2 = v4-v1;
        end
        surf_norm = cross(vec1,vec2)/norm(cross(vec1,vec2));
        
        % Make sure this is an outer facing normal
        if dot(surf_norm,v1) < 0
            surf_norm = -surf_norm;
        end
        
        % Check if lit
        if dot(surf_norm,s_vec) > 0
            % Check if visible from Earth
            if dot(surf_norm,e_vec) > 0
                if abs(vc(1))>370 || vc(2) > 0
                    vis = vis + 1;
                    mag(i) = mag(i) + dot(surf_norm,e_vec)*norm(vec1)*norm(vec2);
                end
            end
        end   
        mag(i) = mag(i)';% + 500*randn;
    end
end
figure
lc = (mag+plight)/max(mag+plight);
tlc = t;
%% Remove some data

% num_gaps = ceil(20*rand);
% size_gaps = ceil(40*rand(1,num_gaps));
% start_gaps = floor(length(lc)*rand(1,num_gaps));
% 
% for i = 1:num_gaps
%     if start_gaps(i)+size_gaps(i) < length(lc)
%         lc(start_gaps(i):start_gaps(i)+size_gaps(i)) = nan;
%     else
%         lc(start_gaps(i):end) = nan;
%     end
% end
% lc = lc(~isnan(lc));
% tlc = t(~isnan(lc));
% 
plot(tlc,lc,'bo','markerfacecolor','b','markersize',2)


%% Autocorrelate lightcurve
% clear LC
% LC(:,1) = tlc;
% LC(:,2) = lc;
% fs = 1/min(diff(LC(:,1)));
% LCnorm = LC(:,2)-mean(LC(:,2));
% tvec = LC(1,1):1/fs:LC(end,1);
% LCfilled = nan(1,length(tvec));
% for i = 1:length(tvec)
%     idx = find(LC(:,1)==tvec(i));
%     if ~isempty(idx)
%         LCfilled(i) = LCnorm(idx);
%     else
%         LCfilled(i) = 0;
%     end
% end
% 
% [autocorr,lags]=xcorr(LCfilled,ceil(max(LC(:,1))/2*fs),'coeff');
% 
% [pksh,lcsh] = findpeaks(autocorr,'minpeakdistance',3);
% short = mean(diff(lcsh))/fs;
% 
% [pklg,lclg] = findpeaks(autocorr,'minpeakprominence',0.1,'minpeakheight',0.05);
% long = mean(diff(lclg))/fs;
% figure
% plot(lags/fs,autocorr)
% hold on
% pks = plot(lags(lcsh)/fs,pksh,'or',lags(lclg)/fs,pklg+0.05,'vk');