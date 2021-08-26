function M = rotation(theta,axis)
% This function generates the rotation matrix by angle theta about a
% specified axis. This matrix will rotate the coordinate system
%   
%   theta: angle of rotation specified in radians
%   axis: the axis of rotation.
%       1: x-axis
%       2: y-axis
%       3: z-axis

if axis == 1
    M = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
elseif axis == 2
    M = [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
elseif axis == 3
    M = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
else
    disp('Enter a valid rotation axis (1:x, 2:y, 3:z)')
    M = 0;
end