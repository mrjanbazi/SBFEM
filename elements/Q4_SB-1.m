function [k_elem,r_elem, m_elem, history1,ti_history1] = Q4(x, u_elem, dt, history0, ti_history0, Prob_pars, Mat_props ) 
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
%% Element Routine for SBFEM 2D Element (Task a,b,c)
% Using user-provided subroutines for Shape Functions, Jacobian, and B-Matrix

    
    %% --- 1. Initialization ---
    Prob_pars.ndf = Prob_pars.ndf;          % 2 DOFs (ux, uy)
    Prob_pars.nnodel = Prob_pars.nnodel;   % 4 Nodes (Boundary) <-- Defined as Prob_pars.nnodel
    Prob_pars.ngp = Prob_pars.ngp;          % 4 Gauss Points
    
    % Initialize Final Matrices (Condensed - only boundary nodes)
    k_elem = zeros(Prob_pars.ndf * Prob_pars.nnodel); 
    r_elem = zeros(Prob_pars.ndf * Prob_pars.nnodel, 1); 
    m_elem = zeros(Prob_pars.ndf * Prob_pars.nnodel); 
    
    history1 = history0;
    ti_history1 = ti_history0;

    % Calculate Scaling Center (Geometric Center)
    X0 = mean(x, 1); 
    
    % Initialize Full Stiffness Matrix (Boundary Nodes + Center Node)
    % Size: (Prob_pars.nnodel + 1) * ndf = 10x10
    n_dof_total = Prob_pars.ndf * (Prob_pars.nnodel + 1); 
    K_total = zeros(n_dof_total, n_dof_total);
    
    % Material Matrix (Linear Elasticity for Plane Strain)
    lam = Mat_props.lambda;
    mu  = Mat_props.mu;
    D_mat = [lam+2*mu, lam,      0;
             lam,      lam+2*mu, 0;
             0,        0,        mu];

    %% --- 2. Loop over Sections (SBFEM Logic) ---
    % FIX: Use Prob_pars.nnodel (plural) here
    for i_sec = 1:Prob_pars.nnodel
        
        % Identify nodes of the current section
        node_1 = i_sec;
        node_2 = mod(i_sec, Prob_pars.nnodel) + 1; 
        
        % Section Coordinates: [Center; Node1; Node2]
        % CRITICAL: This order matches shape function rows!
        X_sec = [X0; x(node_1, :); x(node_2, :)]; 
        
        % Assembly Indices [Center DOFs, Node1 DOFs, Node2 DOFs]
        idx_c = Prob_pars.nnodel*Prob_pars.ndf + (1:Prob_pars.ndf);       
        idx_1 = (node_1-1)*Prob_pars.ndf + (1:Prob_pars.ndf);    
        idx_2 = (node_2-1)*Prob_pars.ndf + (1:Prob_pars.ndf);    
        indices = [idx_c, idx_1, idx_2];     
        
        % --- 3. Integration Loop ---
        for n = 1:Prob_pars.ngp
            
            % A) Get Gauss Points (Your Subroutine)
            [gp_pts, gp_wt] = gp01(n, Prob_pars.ngp);
            
            % B) Mapping & Shape Functions
            xi_ref = gp_pts(1);
            eta_gauss = gp_pts(2);
            
            % Mapping [-1,1] -> [0,1] for Radial Direction
            eta = 0.5 * (1 + eta_gauss); 
            Jac_map = 0.5;
            weight = gp_wt * Jac_map;
            
            % Call Shape Function (Your Subroutine)
            shp = SBFEM_shape_functions([xi_ref, eta]);
            
            % C) Jacobian & B-Matrix (Helper Subroutines)
            [detJ, ~, invJ] = SBFEM_jacobian(shp, X_sec);
            
            [Bmat] = SBFEM_B_matrix(shp, invJ);
            
            % D) Stiffness Assembly: K = B' * D * B * detJ * weight
            K_sec = Bmat' * D_mat * Bmat * detJ * weight;
            
            % Add to Global Matrix
            K_total(indices, indices) = K_total(indices, indices) + K_sec;
            
        end 
    end
    
    %% --- 4. Static Condensation (Eliminate Center Node) ---
    num_b_dofs = Prob_pars.ndf * Prob_pars.nnodel; 
    
    K_bb = K_total(1:num_b_dofs, 1:num_b_dofs);         % Boundary-Boundary
    K_b0 = K_total(1:num_b_dofs, num_b_dofs+1:end);     % Boundary-Center
    K_0b = K_total(num_b_dofs+1:end, 1:num_b_dofs);     % Center-Boundary
    K_00 = K_total(num_b_dofs+1:end, num_b_dofs+1:end); % Center-Center
    
    % Schur Complement
    k_elem = K_bb - K_b0 * (K_00 \ K_0b);
    
    % Residual (Linear)
    u_vec = reshape(u_elem', [], 1);
    r_elem = k_elem * u_vec; 

end 


%% --- SUBROUTINES ---

% 1. Gaussian Points (Exactly as you provided)
function [coor,weights] = gp01(n,ngauss)
    switch ngauss
        case 4
            % 4-point Gauss quadrature for 2D element (Square)
            gpcoor = 1/sqrt(3)*[ -1  1  1 -1; -1 -1  1  1];
            gpweights = [1 1 1 1]; 
    end
    coor = gpcoor(:,n); 
    weights = gpweights(n); 
end

% 2. Shape Functions (Exactly as you provided)
function [shp] = SBFEM_shape_functions(coor)
    % shp: [N, dNdxi, dNdeta] | Rows: [center, node1, node2]

shp = zeros(3, 3); 
    xi = coor(1);
    eta = coor(2);
    
    % Circumferential
    N1 = 0.5 * (1 - xi);   N2 = 0.5 * (1 + xi);   
    dN1_dxi = -0.5;        dN2_dxi = 0.5;
    
    % Radial
    N_rad_boundary = eta;  N_rad_center = 1 - eta;      
    dN_rad_boundary_deta = 1; dN_rad_center_deta = -1;
    
    % Combine
    shp(1,1) = N_rad_center;               shp(1,2) = 0;                          shp(1,3) = dN_rad_center_deta;         
    shp(2,1) = N1 * N_rad_boundary;        shp(2,2) = dN1_dxi * N_rad_boundary;   shp(2,3) = N1 * dN_rad_boundary_deta;  
    shp(3,1) = N2 * N_rad_boundary;        shp(3,2) = dN2_dxi * N_rad_boundary;   shp(3,3) = N2 * dN_rad_boundary_deta;  
end

% 3. Jacobian Calculation (Corrected for SBFEM)
function [detJ, J, invJ] = SBFEM_jacobian(shp, node_coords)
    %   J = 2x2 Jacobian matrix
    J = zeros(2,2);
    % Sum over 3 nodes (Center, Node1, Node2)
    for i = 1:3
        J(1,1) = J(1,1) + shp(i,2) * node_coords(i,1); % dx/dxi
        J(1,2) = J(1,2) + shp(i,2) * node_coords(i,2); % dy/dxi
        J(2,1) = J(2,1) + shp(i,3) * node_coords(i,1); % dx/deta
        J(2,2) = J(2,2) + shp(i,3) * node_coords(i,2); % dy/deta
    end
    detJ = det(J);
    %   J_inv = inverse of Jacobian
    if abs(detJ) > 1e-10
        invJ = inv(J);
    else
        invJ = zeros(2,2); 
        % error('Jacobian determinant is zero or too small');
    end
end

% 4. B-Matrix Calculation (Standard Plane Strain)
function [Bmat] = SBFEM_B_matrix(shp, invJ)
    Bmat = zeros(3, 6); % 3 strain components × 6 DOFs (3 nodes × 2 DOF/node)
    
    % Calculating Physical Derivatives dN/dx and dN/dy
    for i = 1:3  % For each node (Center, Node1, Node2)
        
        % Apply Chain Rule
        % dN/dx = (dxi/dx)*dN/dxi + (deta/dx)*dN/deta
        % where [dxi/dx deta/dx] is the first row of invJ
        dN_dx = invJ(1,1) * shp(i,2) + invJ(1,2) * shp(i,3);
        dN_dy = invJ(2,1) * shp(i,2) + invJ(2,2) * shp(i,3);
        
        % Column indices for node i
        % i=1 (Center) -> cols 1,2
        % i=2 (Node1)  -> cols 3,4
        % i=3 (Node2)  -> cols 5,6
        col_start = 2*(i-1) + 1;
        
        Bmat(1, col_start)   = dN_dx;     % eps_xx
        Bmat(2, col_start+1) = dN_dy;     % eps_yy  
        
        Bmat(3, col_start)   = dN_dy;     % gamma_xy
        Bmat(3, col_start+1) = dN_dx;     % gamma_xy
    end
end