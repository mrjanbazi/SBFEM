%% Function to solve a timestep
% INPUT: 
% nodes: 
% U: displacement vector 
% P: external force vector 
% dt: time increment 
% etu: equations to use - array to delete rows and columns of matrix K 



function [U,R,history1,ti_history1] = solve_timestep(System, U, P,dt,history0, ti_history0,Prob_pars, Mat_props, maxiter)


    %%% Start of the Newton-Raphson scheme: 
    steps = 0;
    tol = 1e-8 ;
    convergence = false;

    % Set the initial values for each time step
    [U, System.V(:,end+1), System.A(:,end+1)] = setInitialValues(System, Prob_pars, U, dt);

    while (convergence==false)
        steps = steps + 1;
        
        % initialize displacement increment 
        dU = zeros(Prob_pars.num_nodes*Prob_pars.ndf,1); 
        
        % compute tangent and internal force vector %
        [K,R,M,history1,ti_history1] = get_KandR(System, U, dt,history0, ti_history0, Prob_pars, Mat_props); 
        
        % Compute the effective Stiffness matrix K and the damping matrix C
        C = zeros(size(K));
        K = getKeff(Prob_pars,K,M,C,dt);

        % Compute residual 
        G = P - R- M*System.A(:,end);    % internal force vector - external force vector - displacement boundary conditions for displacement control - Acceleration Force
        
        normG = norm(G(System.etu)); 

        ErrorCriterion = normG/max(norm(R),norm(P));

        disp([num2str(steps), '. Norm of residual/ Error Criterion: ', num2str(normG),'/ ', num2str(ErrorCriterion)]);
        if (steps > maxiter) 
            error(['Solver failed to converge in ', num2str(maxiter), ' timesteps'] );
        elseif (ErrorCriterion<tol)
            disp(['Solver Converged in ' num2str(steps) ' steps']);
            convergence = true;
        else 
            % solve system of equations 
            dU(System.etu) = K(System.etu,System.etu)\G(System.etu);

            % update solution increment
            [U, System.V(:,end), System.A(:,end)] = updateSolution(System, Prob_pars, U, dU, dt);

            % update time-independent history field 
            ti_history0 = ti_history1; 

        end 
            
    end


        
end



%% supplementary function 
function [K,R,M,history1,ti_history1] = get_KandR(System, U, dt,history0, ti_history0,Prob_pars, Mat_props) 
    % Assemble tangential stiffness matrix K and internal force vector R 

    % initialize new history arrays 
    history1 = 0*history0;
    ti_history1 = 0*ti_history0;

    % initialize dof_per_node array 
    dof_per_node = reshape(1:Prob_pars.num_nodes*Prob_pars.ndf,[Prob_pars.ndf,Prob_pars.num_nodes]).';


    % initialize K and R 
    K = zeros(Prob_pars.num_nodes*Prob_pars.ndf,Prob_pars.num_nodes*Prob_pars.ndf );
    M = zeros(Prob_pars.num_nodes*Prob_pars.ndf,Prob_pars.num_nodes*Prob_pars.ndf );
    R = zeros(Prob_pars.num_nodes*Prob_pars.ndf,1);
    
    % element loop 
    for ele=1:Prob_pars.num_ele
        
        % element coordinate matrix x 
        x = System.nodes(System.elements(ele,:),:);

        
        % get element quantity - for one element a bit redundant
        u_matrix = reshape(U,[Prob_pars.ndf,Prob_pars.num_nodes ]).';
        u_elem = u_matrix(System.elements(ele,:),:);
        
        % get history matrices of the current element 
        hist0_elem = history0(:,:,ele); 
        ti_hist0_elem = ti_history0(:,:,ele); 
    
        % call element routine
        [k_elem,r_elem,m_elem,hist1_elem,ti_hist1_elem] = elem_select(x, u_elem, dt, hist0_elem, ti_hist0_elem, Prob_pars, Mat_props); 

        % Assemble matrix
        dof_per_element = reshape(dof_per_node(System.elements(ele,:),:).',[],1);
        
        K(dof_per_element,dof_per_element) = K(dof_per_element,dof_per_element)  + k_elem;
        M(dof_per_element,dof_per_element) = M(dof_per_element,dof_per_element)  + m_elem;
        R(dof_per_element) = R(dof_per_element) + r_elem;
                   
        % write history into history array
        history1(:,:,ele) = hist1_elem; 
        ti_history1(:,:,ele) = ti_hist1_elem; 
    end 

end 

function [U, V, A] = setInitialValues(System, Prob_pars, U, dt)


    if(contains(Prob_pars.StatDyn,'Static'))
    V = System.V(:,end);
    A = System.A(:,end);
    elseif (contains(Prob_pars.StatDyn,'your solver'))
    
        % Here comes the setup of the initial values according to your time integration method
    end
end

function K = getKeff(Prob_pars,K,M,C,dt)
    
    if(contains(Prob_pars.StatDyn,'Static'))
    % In the static case there is no change of the effective Stiffness Matrix    
    K = K;

    elseif (contains(Prob_pars.StatDyn,'your solver'))
    
        % Here comes the effective Stiffness matrix of your time integration method

    end
end

function [U, V, A] = updateSolution(System, Prob_pars, U, dU, dt)

    if(contains(Prob_pars.StatDyn,'Static'))
    V = System.V(:,end);
    A = System.A(:,end);

    elseif (contains(Prob_pars.StatDyn,'your solver'))
    
        % Here comes the update of the quantities according to your time integration method
    end

    U = U + dU;

end

