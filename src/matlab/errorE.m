function [ err ] = errorE( E, F )
% computes the translational and rotational error of two matries E,F in
% SE(3)
    Erot = E(1:3,1:3);
    Frot = F(1:3,1:3);
    Etrans = E(1:3,4);
    Ftrans = F(1:3,4);
    
    diff = Erot*Frot';
    angle = real(acos((trace(diff)-1.)/2.))*180./pi;
    err = [norm(Etrans-Ftrans), angle];
end


