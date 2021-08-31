function LC = Generate_LC(tmax,Beta,r0,theta0,phi10,phi20,rd0,phi1d0,phi2d0,a1,c1,a2,c2)

% options = odeset('Reltol',1e-12,'AbsTol',1e-12);

%% Fixed system characteristics
mA = 3e11; % Primary mass
mB = 1.3e9; % Secondary mass
mu = 6.67430e-11*(mA+mB);
nu = mA/(mA+mB);
mDART = 700;
vDART = 5000;

%% Variable Physical parameters
I1 = EllipsoidInertias(mA,a1,a1,c1);
I2 = EllipsoidInertias(mB,a2,c2,c2);
r = r0/a1;
I1z = I1(3,3)/(mA*a1^2);
Is = I1(1,1)/(mA*a1^2);
I2z = I2(3,3)/(mB*a1^2);
I2x = I2(1,1)/(mB*a1^2);
I2y = I2(2,2)/(mB*a1^2);
n = sqrt(mu/a1^3);

% Get initial theta_dot of system for equilibrium circular orbit
C2 = -2*I2x+I2y+I2z;
K_eq = sqrt(nu^2/r^5*(r^6+(2*I2z/nu+(3/2)*(C2+I1z-Is))*r^4+(I2z^2/nu^2+3*I2z/nu*(C2+I1z-Is))*r^2+(3/2)*I2z^2/nu^2*(C2+I1z-Is)));

% Perturb system with beta
dK = Beta*nu*(mDART/mB)*(vDART/a1/n)*r;
K = K_eq + dK;

% System spin rate after impact
thetad0 = (K-I2z*phi2d0)/(I2z+nu*r^2);
thetad0 = thetad0*n;

% Initial condisions for integrator
X0 = [r0,theta0,phi10,phi20,rd0,thetad0,phi1d0,phi2d0];
I1 = I1/mA;
I2 = I2/mB;
tspan = 0:600:tmax;

%% Integrate

nt = length(tspan);
t = zeros(nt,1);
X_out = zeros(nt,8);

[t,X_out] = pf2bp_ode78(mu,nu,I1,I2,tspan([1 end]),X0,1e-12);
% [t,X_out] = ode78(@(t,X) planar_f2bp(t,X,mu,nu,I1,I2),tspan([1 end]),X0,1e-12);
% [t,X_out] = ode113(@(t,X) planar_f2bp(t,X,mu,nu,I1,I2),tspan,X0,options);

theta = X_out(:,2);
phi2 = X_out(:,4);


%% Vectors
R = nan(3,nt);
V = nan(3,nt);

for j = 1:length(t)
    r = X_out(j,1);
    th = X_out(j,2);
    rd = X_out(j,5);
    thetad = X_out(j,6);
    xhat = [1;0;0];
    yhat = [0;1;0];
    
    rhat = cos(th)*xhat + sin(th)*yhat;
    theta_hat = -sin(th)*xhat + cos(th)*yhat;
    
    R(:,j) = r*rhat;
    V(:,j) = rd*rhat + r*thetad*theta_hat;
end

%% Redefine theta to be between 0, 2*pi
for i = 1:length(t)
    if X_out(i,2) > 2*pi
        X_out(i:end,2) = X_out(i:end,2) - 2*pi;
    end
end

%% Set up system for lightcurve generation
e_vec = [0,1,0];    % Vector to Earth
s_vec = [0,1,0];    % Vector to sun

% [bx,by,bz] = ellipsoid(0,0,0,a2,c2,c2);
[bx,by,bz] = ellipsoid_csml(a2,c2,c2);

m=21;n=21;
v = [bx(:) by(:) bz(:)];
q = (1:m*n-m-1)';
q(m:m:end) = [];
f = [q q+m q+m+1 q+1];

sec.vertices = v;
sec.faces = f;

% S = surf(bx,by,bz);
% axis equal
% sec = surf2patch(S);
% close all

plight = 1.5*8.477442968647369e+05;   % Primary brightness


%% Pre-compute vertices
v1m = sec.vertices(sec.faces(:,1),:)';
v2m = sec.vertices(sec.faces(:,2),:)';
v3m = sec.vertices(sec.faces(:,3),:)';
v4m = sec.vertices(sec.faces(:,4),:)';

%% quasi-Lightcurves
mag = zeros(1,length(t));
for i = 1:length(t)
    
    % System rotation
    M = rotation(theta(i)+phi2(i),3);
    NB = M';
    
    vis = 0;
    for f = 1:length(sec.faces)
        % Get surface normal
        % v1 = R(:,i) + NB*sec.vertices(sec.faces(f,1),:)';
        % v2 = R(:,i) + NB*sec.vertices(sec.faces(f,2),:)';
        % v3 = R(:,i) + NB*sec.vertices(sec.faces(f,3),:)';
        % v4 = R(:,i) + NB*sec.vertices(sec.faces(f,4),:)';
        v1 = R(:,i) + NB*v1m(:,f);
        v2 = R(:,i) + NB*v2m(:,f);
        v3 = R(:,i) + NB*v3m(:,f);
        v4 = R(:,i) + NB*v4m(:,f);
        
        
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
%         surf_norm = crpd(vec1,vec2);
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
        mag(i) = mag(i)';
    end
end
%% Create lightcurve
lc = (mag+plight)/max(mag+plight);
LC = [t,lc'];
