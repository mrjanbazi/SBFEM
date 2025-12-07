%% select element 
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


function [k_elem,r_elem,m_elem,history1,ti_history1] = elem_select(x, u_elem, dt, history0, ti_history0, Prob_pars, Mat_props )

    switch Prob_pars.element_type
        case 'Q1' 
            [k_elem,r_elem,m_elem,history1,ti_history1] = el_Q1(x, u_elem, dt, history0, ti_history0, Prob_pars, Mat_props );
        case 'Q4' 
            [k_elem,r_elem,m_elem,history1,ti_history1] = Q4(x, u_elem, dt, history0, ti_history0, Prob_pars, Mat_props );
        case 'Q4_SB' 
            [k_elem,r_elem,m_elem,history1,ti_history1] = Q4_SB(x, u_elem, dt, history0, ti_history0, Prob_pars, Mat_props );
    
        otherwise
            k_elem = zeros(Prob_pars.ndf*Prob_pars.nnodel);
            r_elem = zeros(Prob_pars.ndf*Prob_pars.nnodel,1); 
            m_elem = zeros(Prob_pars.ndf*Prob_pars.nnodel);
            history1 = history0; 
            ti_history1 = ti_history0; 
            disp(['Element ', Prob_pars.element_type, ' not implemented yet'])
    end 

end 