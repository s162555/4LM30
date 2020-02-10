%% Type of function: Velocity-Verlet scheme velocity function
%*************************************************************************
%*(c) 2020, Gert-Jan Slokker (Mechanical Engineering, TUE)               *
%* Course: 4LM30 - Multiscale Modelling for Polymer Mechanics            *
%* Exercise 2: Single polymer chain in 3D                                *
%* Sub-function goal: new velocity v(n+1)based on the previous conditions*
%* for the velocity v(n), force f(n) and the new conditions for the force*
%* f(n+1) according to an Velocity-Verlet integration scheme.            *
%*************************************************************************

function  [v] = Timestepvelocity_Euler(vold,Fold,m,dt);
% Velocity calculation
v = vold+(Fold/m)*dt;
end