%% Type of function: Velocity-Verlet scheme velocity function
% Description: this function describe calculated the new velocity v(n+1) 
% based on the previous conditions for the velocity v(n), force f(n) and 
% the new conditions for the force f(n+1) according to an Velocity-Verlet
% integration scheme.

function  [v] = Timestepvelocity(vold,Fnew,Fold,m,dt);
% Velocity calculation
v = vold+((Fnew+Fold)/(2*m))*dt;
end