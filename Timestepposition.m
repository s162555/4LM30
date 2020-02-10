%% Type of function: Velocity-Verlet scheme position function
% Description: this function describe calculated the new position r(n+1) 
% and based on the previous conditions for the position r(n), velocity 
% v(n) and force f(n) according to an Velocity-Verlet integration scheme.

function  [r] = Timestepposition(rold,vold,Fold,m,dt);

% Location calculation
r  = rold + dt * vold + (Fold/(2*m))*dt^2; 

end