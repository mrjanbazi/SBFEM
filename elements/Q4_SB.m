function [k_elem,r_elem, m_elem, history1,ti_history1] = Q4_SB(x, u_elem, dt, history0, ti_history0, Prob_pars, Mat_props ) 
% INPUT:
% x: nodal coordinates (nnodel x ndf) 
% u_elem: nodal displacements (nnodel X ndf) 
% dt: time increment 
% history0: Gauss Point history array (nnodel x nh) 
% ti_history0: time-independent history array (nnodel x nh) 
% Prob_pars: Problem parameters see main.m 
% Mat_props: Material properties see main.m 

% OUTPUT:
% k_elem: element stiffness matrix (ndf*nnodel x ndf*nnodel) 
% r_elem: element residual vector (ndf*nnodel x 1) 
% history1: new trial history array 
% ti_history1: updated time-independent history array 



% Here comes your Finite Element

    % initialize the element
    k_elem = zeros(Prob_pars.ndf*Prob_pars.nnodel); 
    k_sec = zeros(Prob_pars.ndf*3);
    r_elem = zeros(Prob_pars.ndf*Prob_pars.nnodel,1); 
    r_sec = zeros(Prob_pars.ndf*3,1);
    m_elem = zeros(Prob_pars.ndf*Prob_pars.nnodel);
    m_sec = zeros(Prob_pars.ndf*3);

    % Initialize Full Stiffness Matrix (Boundary Nodes + Center Node)
    % Size: (n_nodes + 1) * ndf = 10x10
    n_dof_total = Prob_pars.ndf * (Prob_pars.nnodel + 1); 
    k_total = zeros(n_dof_total, n_dof_total);
    m_total = zeros(n_dof_total, n_dof_total);

    % Calculate Scaling Center (Geometric Center)
    X0 = mean(x, 1);
    U0 = mean(u_elem, 1);


    % Section Assembly (SBFEM Logic)
        
    % Section Coordinates: [Node1x, Node1y; Node2x, Node2y; Centerx, Centery]

    x_sec1 = [x(1,:)'; x(2,:)'; X0']; u_sec1 = [u_elem(1,:)'; u_elem(2,:)'; X0'];
    x_sec2 = [x(2,:)'; x(3,:)'; X0']; u_sec2 = [u_elem(1,:)'; u_elem(2,:)'; X0'];
    x_sec3 = [x(3,:)'; x(4,:)'; X0']; u_sec3 = [u_elem(1,:)'; u_elem(2,:)'; X0'];
    x_sec4 = [x(4,:)'; x(1,:)'; X0']; u_sec4 = [u_elem(1,:)'; u_elem(2,:)'; X0'];

    x_sec = {x_sec1, x_sec2, x_sec3, x_sec4};
    u_sec = {u_sec1, u_sec2, u_sec3, u_sec4};
        
    % Assembly Indices for K_tot [Center DOFs, Node1 DOFs, Node2 DOFs]
   
    indices1 = [1,2,3,4,9,10];
    indices2 = [3,4,5,6,9,10];
    indices3 = [5,6,7,8,9,10];
    indices4 = [7,8,1,2,9,10];
    
    indice = {indices1, indices2, indices3, indices4};


    
    % Sort
    x = [x(1,:)';x(2,:)';x(3,:)';x(4,:)'];
    u_elem = [u_elem(1,:)';u_elem(2,:)';u_elem(3,:)';u_elem(4,:)'];



    
    % Section Loop

    for i = 1:Prob_pars.nnodel

        x_se = x_sec{i};
        u_se = u_sec{i};
        indices = indice{i};
        
        % GP Loop

        for n = 1:Prob_pars.ngp


            % GP Coor
            [coor,weights] = gp01(n,Prob_pars.ngp);
        
            % Mapping [-1,1] -> [0,1] for Radial Direction
            xi_gp = coor(1);
            eta_gauss = coor(2);
            
            % Mapping [-1,1] -> [0,1] for Radial Direction
            eta_gp = 0.5 * (1 + eta_gauss); 
            Jac_map = 0.5;
            weights = weights * Jac_map;


            % Shape Functions
            [shp] = shp_SBFEM(xi_gp,eta_gp);

            % Jacobian
            [detJ,J] = jaco_SBFEM(shp,x_se);
    
            % B Operator
            [Bmat] = bmat_SBFEM(shp,J);
    
            % Strain Tensor
            [E,F,Fmat] = epsilon01(Bmat,u_se);
    
            % History Variables at GP
            h0_gp = history0(:,n);
        
            % Stress and Materialtangent
            [S,dSdE,h1_gp] = material_select(E,dt,h0_gp,Mat_props);
    
            % right hand side
            r_sec = r_sec + Bmat'*Fmat'*S*detJ*weights(1)*weights(2);
    
            % stiffness matrix
            k_sec = k_sec + Bmat'*Fmat'*dSdE*Fmat*Bmat*detJ*weights(1)*weights(2);
    
            % Geometrical stiffness
            [Smat] = smat01(S);
            k_sec = k_sec + Bmat'*Smat*Bmat*detJ*weights(1)*weights(2);
    
            % Mass Matrix
            Hmat = hmat(shp);
            m_sec = m_sec + Hmat'*Hmat*Mat_props.rho*detJ*weights(1)*weights(2); 

            % New History Variables
            history1(:,n) = h1_gp;

        end

        ti_history1 = ti_history0;
        k_total(indices, indices) = k_total(indices, indices) + k_sec;
        %r_elem = r_elem + r_sec;
        m_total(indices, indices) = m_total(indices, indices) + m_sec;
        
    end

% Static Condensation (Eliminate Center Node)
num_b_dofs = Prob_pars.ndf * Prob_pars.nnodel; 
   
K_bb = k_total(1:num_b_dofs, 1:num_b_dofs);         % Boundary-Boundary
K_b0 = k_total(1:num_b_dofs, num_b_dofs+1:end);     % Boundary-Center
K_0b = k_total(num_b_dofs+1:end, 1:num_b_dofs);     % Center-Boundary
K_00 = k_total(num_b_dofs+1:end, num_b_dofs+1:end); % Center-Center

m_bb = m_total(1:num_b_dofs, 1:num_b_dofs);         % Boundary-Boundary
m_b0 = m_total(1:num_b_dofs, num_b_dofs+1:end);     % Boundary-Center
m_0b = m_total(num_b_dofs+1:end, 1:num_b_dofs);     % Center-Boundary
m_00 = m_total(num_b_dofs+1:end, num_b_dofs+1:end); % Center-Center
    
% Schur Complement
k_elem = K_bb - K_b0 * (K_00 \ K_0b);
m_elem = m_bb - m_b0 * (m_00 \ m_0b);
    
% Residual (Linear)
r_elem = k_elem * u_elem;

end




%% Subroutines
% Gaussian points
function[coor,weights] = gp01(n,ngauss)
switch ngauss
    case 1
        gpcoor = [0;0;0];
        gpweights = [2;2;2];
     case 4
        gpcoor = 1/sqrt(3)*[-1 1 1 -1; -1 -1 1 1];
        gpweights = [1 1 1 1;1 1 1 1];
    
    case 8
        gpcoor = 1/sqrt(3)*[1 1 -1 -1 1 1 -1 -1; 1 -1 1 -1 1 -1 1 -1;-1 -1 -1 -1 1 1 1 1];
        gpweights = [1 1 1 1 1 1 1 1;1 1 1 1 1 1 1 1; 1 1 1 1 1 1 1 1];
    case 9
        gpcoor = sqrt(3/5)*[-1 -1 -1 0 0 0 1 1 1;-1 0 1 -1 0 1 -1 0 1];
        gpweights = [5/9 5/9 5/9 8/9 8/9 8/9 5/9 5/9 5/9; 5/9 8/9 5/9 5/9 8/9 5/9 5/9 8/9 5/9];
    case 16
        gpca = sqrt(3/7+2/7*sqrt(6/5));
        gpcb = sqrt(3/7-2/7*sqrt(6/5));
        gpcoor = [-gpca -gpca -gpca -gpca -gpcb -gpcb -gpcb -gpcb gpcb gpcb gpcb gpcb gpca gpca gpca gpca;
            -gpca -gpcb gpcb gpca -gpca -gpcb gpcb gpca -gpca -gpcb gpcb gpca -gpca -gpcb gpcb gpca];
        gpwa = (18-sqrt(30))/36;
        gpwb = (18+sqrt(30))/36;
        gpweights = [gpwa gpwa gpwa gpwa gpwb gpwb gpwb gpwb gpwb gpwb gpwb gpwb gpwa gpwa gpwa gpwa;
            gpwa gpwb gpwb gpwa gpwa gpwb gpwb gpwa gpwa gpwb gpwb gpwa gpwa gpwb gpwb gpwa];
end
coor = gpcoor(:,n);
weights = gpweights(:,n);
end

% Shape functions
function[shp] = shp_SBFEM(xi_gp, eta_gp)

%    shp  | 
%  matrix |  shp      ; shp,xi      ; shp,eta
%_________|___________________________________
% node1   |  shp1     ; shp1,xi     ; shp1,eta
% node2   |  shp2     ; shp2,xi     ; shp2,eta 
% Center  |  shp0     ; shp0,xi     ; shp0,eta

shp = zeros(3,3);
xi = xi_gp;
eta= eta_gp;
       
shp(1,1) = eta/2*(1- xi) ; shp(1,2) = -eta/2 ; shp(1,3) = 1/2*(1 - xi);
shp(2,1) = eta/2*(1+ xi) ; shp(2,2) = eta/2  ; shp(2,3) = 1/2*(1 + xi);
shp(3,1) = 1-eta         ; shp(3,2) = 0      ; shp(3,3) = -1          ;
end

% Jacobi matrix
function[detJ,J] = jaco_SBFEM(shp,x_sec)
Nxi = [shp(1,2) 0 shp(2,2) 0 shp(3,2) 0;...
       0 shp(1,2) 0 shp(2,2) 0 shp(3,2)];

Neta = [shp(1,3) 0 shp(2,3) 0 shp(3,3) 0;...
        0 shp(1,3) 0 shp(2,3) 0 shp(3,3)];

J = zeros(2);
J(:,1) = Nxi*x_sec;
J(:,2) = Neta*x_sec;

detJ = det(J);
if abs(detJ) > 1e-10
    invJ = inv(J);
else
    invJ = zeros(2,2);    % To prevent det(J) < 0;
end

end

% B Operator
function[Bmat] = bmat_SBFEM(shp,J)

b1 = [J(1,1) 0;0 J(1,2);...
    J(1,2) 0; 0 J(1,1)];
b2 = [J(2,1) 0;0 J(2,2);...
    J(2,2) 0;0 J(2,1)];

Nxi = [shp(1,2) 0 shp(2,2) 0 shp(3,2) 0;...
       0 shp(1,2) 0 shp(2,2) 0 shp(3,2)];

Neta = [shp(1,3) 0 shp(2,3) 0 shp(3,3) 0;...
        0 shp(1,3) 0 shp(2,3) 0 shp(3,3)];


Bmat = b1*Nxi+b2*Neta;
end


% % Green-Lagrange strain and deformation gradient
function[Eps,F,Fmat] = epsilon01(Bmat,u_mech)
H = Bmat*u_mech;

F = [1;1;0;0];
F = F + H;

Fmat = [F(1), 0, F(3), 0;
        0, F(2), 0, F(4);
        F(3), F(4), F(1), F(2)];

 voigti = [1;1;0];
 Eps = Fmat*F;
  for k2=1:3
     Eps(k2) = 0.5*(Eps(k2)-voigti(k2)); %[E11;E22;2*E12]

end
end

% Geometrical stiffness
function [Smat] = smat01(s)
s = [s(1) s(3);...
    s(3) s(2)];

Smat = [s(1,1) 0 s(1,2) 0;...
    0 s(2,2) 0 s(1,2);...
    s(1,2) 0 s(1,1) 0;...
    0 s(1,2) 0 s(2,2)];
end


function Hmat = hmat(shp)

Hmat = [shp(1,1) 0 shp(2,1) 0 shp(3,1) 0;...
    0 shp(1,1) 0 shp(2,1) 0 shp(3,1)];

end