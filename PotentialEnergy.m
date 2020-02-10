%% Type of function: Potential Energy function
%*************************************************************************
%*(c) 2020, Gert-Jan Slokker (Mechanical Engineering, TUE)               *
%* Course: 4LM30 - Multiscale Modelling for Polymer Mechanics            *
%* Exercise 2: Single polymer chain in 3D                                *
%* Sub-function goal: calculating the total potential energy for N       *
%* particles by using the provided equation: Epot = k/2*(r-L0)^2         *
%*************************************************************************
function [Epot] = PotentialEnergy(r,bond,k);

Epot = 0; %setting the initial potential energy to 0

  for i = 1:size(bond,1)     
      L0 = bond(i,3);                                    % Initial bond distance
      rdistance = norm(r(bond(i,1),:)-r(bond(i,2),:));   % Current bond distance
      Epot = Epot + k/2*(rdistance-L0)^2;                % Resulting potential energy
  end
end


