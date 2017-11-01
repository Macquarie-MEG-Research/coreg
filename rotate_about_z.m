%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rotate_about_z - make a rotation matix for arbitrary rotation in degrees
% around z axis
%
% Written by Paul Sowman Oct 2017 (http://web.iitd.ac.in/~hegde/cad/lecture/L6_3dtrans.pdf - page 4)
%
% INPUTS:
% - deg        = degrees of rotation required
%
% OUTPUTS:
% - rmatx      = a 4*4 rotation matrix for deg degrees about z
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rmatx=rotate_about_z(deg)
    deg   = deg2rad(deg);
    rmatx = [cos(deg) sin(deg) 0 0;-sin(deg) cos(deg) 0 0;0 0 1 0;0 0 0 1];
end