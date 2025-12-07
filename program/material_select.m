%% Material select routine 
% INPUT: 
% E: Green Lagrange Strain in Voigt notation [E11;E22;E33;2*E12;2*E13;2*E23]
% dt: time increment 
% h0_gp: history from the beginning of the time step at the current Gauss Point 
% Mat_Props: contains material type and material parameters 

% OUTPUT: 
% S: 2. Piola Kirchhoff stress in Voigt notation [S11,S22,S33,S12,S13,S23] 
% dSdE: Material tangent, derivative of S with respect to E
% hGP1: new trial history variables at the current Gauss Point 


function [S,dSdE,h1_GP] = material_select(E,dt,h0_gp,Mat_props)

    switch Mat_props.type
        case 'St_Venant_Kirchhoff'
            [S,dSdE,h1_GP] = st_venant_Kirchhoff(E,dt,h0_gp,Mat_props);
        case 'St_Venant_Kirchhoff2D'
            [S,dSdE,h1_GP] = st_venant_Kirchhoff2D(E,dt,h0_gp,Mat_props);
        otherwise
            disp(['Material ', Mat_props.type, ' not implemented yet!'])
    end 

end 