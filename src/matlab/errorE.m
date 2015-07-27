function [ err ] = errorE( E, F )
% computes the translational and rotational error of two matries E,F in
% SE(3)
    Erot = E(1:3,1:3);
    Frot = F(1:3,1:3);
    Etrans = E(1:3,4);
    Ftrans = F(1:3,4);
err = [norm(Etrans-Ftrans), rotation_error(Erot,Frot)];
end

function [ err ] = rotation_error( R0,R1 )
%ROTATION_ERROR computes the angular distance between two rotation matrices
    err = angle_R(R0*R1');
end

function [ theta] = angle_R( R )
%ANGLE_R Compute the rotation angles theta of a rotation matrix R
    theta = real(acos((trace(R)-1.)/2.))*180./pi;
end