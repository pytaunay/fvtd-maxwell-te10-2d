%%% Author Pierre-Yves Taunay
%%% Nov. 2015
%%%
%%% Goal
%%% Simulate a straight section of WR112 waveguide
%%% using the FVTD method for Maxwell's equations in vac.
%%%
%%% Details
%%% The Maxwells equations are put in a conservative form, and the usual
%%% finite volume treatment is performed. 
%%% The numerical flux at the cell boundaries is calculated using a Steger
%%% and Warming flux splitting method. Time stepping is done with RK4 
%%% and explicit Euler

clear all;

%%% Define the geometry , initial power, constants
a = 28.4988e-3;
b = 12.6238e-3;
L = 4e-2;
f = 10e9; % Hz
P10 = 1; % Watts 
tmax = 4/f; % N time steps 

mu = 4*pi*1e-7;
eps = 8.851e-12;
c = 1/sqrt(mu*eps);

% Number of discretization points
Nx = 20;
Nz = 40;
Nt = 200; % Discretization of a single period 
Ns = 4; % Number of sides

dx = a/Nx;
dz = L/Nz;
dt = 1/(f*Nt); 

dS = dx*dz;

% Derived values
om = 2*pi*f; 
beta = beta_te10(om,mu,eps,a);
Y = sqrt(eps/mu);

% Get coefficient for E and H
A10 = sqrt(4*pi^2*P10/(om*mu*a^3*b)*1/real(beta));

%%% Construct initial values for E and H
Ue = cell(Nx,Nz);
Uh = cell(Nx,Nz);
[Ue{:}] = deal(zeros(1,1));
[Uh{:}] = deal(zeros(2,1));

Uen = cell(Nx,Nz);
Uhn = cell(Nx,Nz);
[Uen{:}] = deal(zeros(1,1));
[Uhn{:}] = deal(zeros(2,1));

Ueall = cell(Nx,Nz,ceil(tmax/dt));
Uhall = cell(Nx,Nz,ceil(tmax/dt));


%%% Create some matrix holders
% Matrix of surface normals
Nmat = [-1,1,0,0;
        0,0,0,0;
        0,0,-1,1];

% Assuming that mu and eps are constant everywhere, get Ap and Am for each vector direction
Apmat = cell(3,1);
Ammat = cell(3,1);
for l=1:Ns
    Apmat{l} = flux_matrix(Nmat(:,l),mu,eps);
    Ammat{l} = -flux_matrix(-Nmat(:,l),mu,eps);
end

% Alpha inverse matrix
alpha_mat = eye(3);
alpha_mat(1,1) = 1/eps;

alpha_mat(2,2) = 1/mu;
alpha_mat(3,3) = 1/mu;

am1 = alpha_mat^-1;

% Matrix for PEC BC
Tpec = zeros(3,3);
Tpec(1,1) = 2;

%%% Time stepping
t = 0;
nt = 1;
tic()

%%% Original values
for i=1:Nx
    for k=1:Nz
        Ueall{i,k,nt} = Ue{i,k};
    end
end

EHvec = zeros(3,Nx);
while t < tmax
	%t = t+dt/2;
	t = t+dt;
    half_step;

    Ue = Uen;
    Uh = Uhn;

%	t = t+dt/2;
%    half_step;

    Ue = Uen;
    Uh = Uhn;

    nt = nt+1;
    for i=1:Nx
        for k=1:Nz
            Ueall{i,k,nt} = Ue{i,k};
            Uhall{i,k,nt} = Uh{i,k};
        end
    end
end
    
toc()


