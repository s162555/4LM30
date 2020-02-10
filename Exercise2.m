close all; clear all; clc; warning off; tic
%*************************************************************************
%*(c) 2020, Gert-Jan Slokker (Mechanical Engineering, TUE)               *
%* Course: 4LM30 - Multiscale Modelling for Polymer Mechanics            *
%* Exercise 2: Single polymer chain in 3D                                *
%* Problem description: writing an MD-code that simulates the motion of  *
%* a linear polymer chain consisting of N particles (="monomers").       *
%*************************************************************************
%% Parameters
N = 10;         % Number of particles
m = 1;          % Mass
k = 1;          % Stiffness
L0 = 1;         % Initial length
dim = 3;        % Problem dimensions

%simulation Settings
Simtime  = 10;                     % Set simulation time [s]
dt = 0.01;                         % Set timesteps [s]
    
%% Initialization
%initialization of the position
r = zeros(N,dim,Simtime/dt);
for i = 1:N
   r(i,1,1) = 0+L0*i;            % Straight line with distance L0 between neighbours
end

%initialization of the velocity
v = zeros(N,dim,Simtime/dt);
vrandom = randn(N,dim)*0.3;      % Random initial velocities [root-mean square = 0.3; average = 0]
v(:,:,1) = vrandom; 

%initialization of the force
Fnew = zeros(N,dim);             % Set initial force to zero

%bonding calculations
bond = zeros(N-1,3);
for i = 1:N-1
    bond(i,:) = [i,i+1,L0];      % General bonding information
end

%initialization of the Energy
Epot = zeros(Simtime/dt,1); 
Epot(1) = PotentialEnergy(r(:,:,1),bond,k);    % Initial potential energy

Ekin = zeros(Simtime/dt,1);
Ekin(1) = KineticEnergy(v(:,:,1),m);           % Initial kinetic energy

Etot = zeros(Simtime/dt,1);
Etot(1) = Epot(1)+Ekin(1);                     % Initial total energy


%% Simulation
% simulation part has been bluild up by using specific created functions,
% such as the 'forceall', 'Timestepposition','Timestepvelocity','Potential-
% Energy', 'KineticEnergy' which can be looked up by typing "open ...".
  
% Applying timesteps
for n = 1:(Simtime/dt)               % Starting for-loop
    
    %Smart storage
    Fold = Fnew;                     % Transporting f(n+1) to f(n) for the new timestep
    
    %Timestep - Position calculation
    %r(:,:,n+1) = r(:,:,n) + v(:,:,n)*dt+(Fold/(2*m))*(dt)^2;
    r(:,:,n+1) = Timestepposition(r(:,:,n),v(:,:,n),Fold,m,dt);
    
    %Timestep - Force calculation
    Fnew = forceall(r(:,:,n+1),bond,k);
    
    %Timestep - Velocity calculation
    %v(:,:,n+1) = v(:,:,n) + ((Fnew+Fold)/(2*m))*dt;
    %v(:,:,n+1) = Timestepvelocity_Euler(v(:,:,n),Fold,m,dt);
    v(:,:,n+1) = Timestepvelocity(v(:,:,n),Fnew,Fold,m,dt);
    
    %Timestep - Energy calculation
    Epot(n+1) = PotentialEnergy(r(:,:,n+1),bond,k);
    Ekin(n+1)= KineticEnergy(v(:,:,n+1),m);
    Etot(n+1) = Ekin(n+1)+Epot(n+1);
    
    %Plotting trajectory
    figure(1)
    %plot3(r(:,1,n),r(:,2,n),r(:,3,n),'-bo')
     plot3(r(:,1,n),r(:,2,n),r(:,3,n),'-ro',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.49 1 .63],...
    'MarkerSize',10)
    title('MD-simulation [Single polymer chain in 3D](Current position)')
    %xlabel('X-location'); ylabel('Y-location'); zlabel('Z-location'); 
    xlim([-5.5 12.5]); ylim([-5.5 5.5]); zlim([-5.5 5.5]);
    grid on
end

figure(2)
    plot3(r(:,1,1),r(:,2,1),r(:,3,1),'-ro',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.49 1 .63],...
    'MarkerSize',10)
    title('MD-simulation [Single polymer chain in 3D](Initial position)')
    xlabel('X-location'); ylabel('Y-location'); zlabel('Z-location'); 
    xlim([-5.5 12.5]); ylim([-5.5 5.5]); zlim([-5.5 5.5]);
grid on

%% End-result plotting

t = 0:dt:Simtime;     % Time vector for plotting

%resulting energy plot
figure(3)
plot(t,Ekin, 'r')
hold on
plot(t,Epot, 'b')
hold on
plot(t,Etot,'k')
xlabel('Time')
ylabel('Energy')
legend('Kinetic energy','Potential energy', 'Total energy')
grid on

toc;
