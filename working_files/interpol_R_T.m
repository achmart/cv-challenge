function [Rint,Tint] = interpol_R_T(R,T,p)
%INTERPOLATEROTATION Summary of this function goes here
%   se E and percentage p to interpolate the E for free-viewpoint
% Code has been changed such that Quaternion.m is used instead

% Interpolate the Translation
Tint = p .* T;

% Interpolate the rotation
%rot_quat = Quaternion(R);
[s, v] = Rotation2Quaternion(R);

% interpolate rotation between 0 (no rotation) and R
%quat_interpol = rot_quat.scale(p);
[s0,v0] = Rotation2Quaternion(eye(3));
slerp = QuatSlerp([s0;v0], [s;v], p);

% Convert back to rotation matrix
%R_interpol = quat_interpol.R();
Rint = Quaternion2Rotation(slerp(1), slerp(2:4))
end

