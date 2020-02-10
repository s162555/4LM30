%% Type of function: Potential Energy function
%*************************************************************************
%*(c) 2020, Gert-Jan Slokker (Mechanical Engineering, TUE)               *
%* Course: 4LM30 - Multiscale Modelling for Polymer Mechanics            *
%* Exercise 2: Single polymer chain in 3D                                *
%* Sub-function goal: calculating the total kinetic energy for N         *
%* particles by using the provided equation: Ekin = 0.5*m*v^2            *
%*************************************************************************

function [Ekin] = KineticEnergy(vold,m);

Ekin = 0; % Setting the initial kinetic energy to 0

for i = 1:size(vold,1)
    Ekin = Ekin + 0.5*m*vold(i,:)*vold(i,:)'; % Resulting potential energy
end
end