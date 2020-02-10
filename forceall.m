%% Type of function: Forceall function and Forcepair function
%*************************************************************************
%*(c) 2020, Gert-Jan Slokker (Mechanical Engineering, TUE)               *
%* Course: 4LM30 - Multiscale Modelling for Polymer Mechanics            *
%* Exercise 2: Single polymer chain in 3D                                *
%* Sub-function goal: calculating the resulting force on each particle,  *
%* which has done by the use of two functions:                           *
%* --> Forcepair function: function which will calculate the force for   * 
%* pairs of particles.                                                   *
% --> Forceall function: function which will calculate the final         *
% resulting force on each separated particle.                            *
%*************************************************************************

function [Fnew] = forceall(pos,bond,k)

    Fnew = zeros(size(pos,1),size(pos,2));
    
    for i = 1:size(bond,1)
        [Fpair] = forcepair(pos(bond(i,1),:),pos(bond(i,2),:),bond(i,3),k);
        Fnew(bond(i,2),:) = Fnew(bond(i,2),:) + Fpair;
        Fnew(bond(i,1),:) = Fnew(bond(i,1),:) - Fpair;
    end
end

function [Fpair] = forcepair(pos1,pos2,l0,k);
    r = norm(pos2-pos1);    % Current distance between the particles
    dir = (pos2-pos1)/r;    % Normalized direction
    Fabs = k*(l0 - r);      % Absolute force
    Fpair = Fabs * dir;     % Force with direction for particle 2
end