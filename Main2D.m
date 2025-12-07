
% Main File to run a element finite element problem for the AdCom project
% task b)Implement the formulation of a) in the given finite element code.
%Q4-2D-2DOF-1GP
addpath("elements", "materials", "program") 
clear
close all

%% Define System 

% define geometry 
System.nodes = [0,0; 
         1,0;
         1,1;
         0,1;
         1,2;
         0,2;
         1,3;
         0,3]; 

System.elements = [1,2,3,4;
                   4,3,5,6;
                   6,5,7,8];

Prob_pars.num_nodes = size(System.nodes,1);        % number of nodes 
Prob_pars.num_ele = size(System.elements,1);       % number of elements 
Prob_pars.ndm = size(System.nodes,2);              % number of dimensions 
Prob_pars.nnodel = size(System.elements,2);        % number of nodes per element

% define problem parameters 
Prob_pars.ndf = 2;                                 % number of degrees of freedom
Prob_pars.ngp = 4;                                 % number of Gauss Points 
Prob_pars.nh = 1;                                  % number of history variables 
Prob_pars.ntih = 0;                                % number of time-independent history variables 
Prob_pars.element_type = 'Q4_SB';                     % element type identifier 
Prob_pars.StatDyn = 'Static';                      % Switch Static analysis 'Static' and your solver

% define material  
Mat_props.type = 'St_Venant_Kirchhoff2D';            % material type identifier
Mat_props.lambda = 400;                            % Lamé Constant
Mat_props.mu = 400;                                % Lamé Constant
Mat_props.rho = 3;                                 % density

% define postprocessing --> Flag: 0 = No; 1 = Yes
PostProc.NodeNum = 1;                              % Show NodeNumbers and Location
PostProc.Support = 1;                              % Show Location of supports at the nodes
PostProc.Deform = 1;                               % Show Deformed Nodes 

%% Define boundary conditions 

% initialize matrices 
dirichlet = zeros(Prob_pars.num_nodes,Prob_pars.ndf);         % define nodes and degrees of freedom where boundary dirichlet conditions are placed
displacements = zeros(Prob_pars.num_nodes,Prob_pars.ndf);      % define values of the displacement 
forces = zeros(Prob_pars.num_nodes,Prob_pars.ndf);            % external force matrix 

% fill with specific values 
dirichlet([1,2],1) = 1; 
dirichlet([1,2],2) = 1;
dirichlet([7,8],1) = 1;


%forces(3:4,2) = 100;  
displacements([7,8],1) = 1; 

% transform into the format used in the solver 
System.etu = logical(~reshape(dirichlet.',[],1));         % equations to use - boolean array: false at the DOFS with Dirichlet boundary conditions
uBCmax = reshape(displacements.',[],1);                    % boundary displacements for displacement control
Pmax = reshape(forces.',[],1);                            % Force vector for force control 

% initialize history 
if Prob_pars.nh == 0 
    % define a pseudo history variable 
    Prob_pars.nh = 1;
end 
    history = zeros(Prob_pars.nh,Prob_pars.ngp,Prob_pars.num_ele);        % history 
    ti_history = zeros(Prob_pars.ntih,Prob_pars.ngp, Prob_pars.num_ele);   % time-independent history variables 



%% time stepping 
time = 0;
dt = 1;                          % time increment 
loadfactors = 0:0.1:1;       % loadfactors to define load stepping (has to start with zero)

% Initialize the quantities for computation
U = zeros(Prob_pars.num_nodes*Prob_pars.ndf,1);
System.U = zeros(size(U));
System.V = zeros(size(U));
System.A = zeros(size(U));

for i=2:length(loadfactors)

    time = time +dt; 

    disp('------ New Time Step -----------------------');
    disp(['time: ', num2str(time)]);
    
    % dUBC= (loadfactors(i)-loadfactors(i-1))*uBCmax; 
    U(~System.etu) = loadfactors(i)*uBCmax(~System.etu); 
    P = loadfactors(i)*Pmax;
    
    % solve timestep 
    [U,R,historyNew,ti_historyNew] = solve_timestep(System, U, P,dt,history, ti_history,Prob_pars, Mat_props, 16) ;

    history = historyNew; 
    ti_history = ti_historyNew;
end



PostProcessor2D(System,PostProc,dirichlet,displacements,U)