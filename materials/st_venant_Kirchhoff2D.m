%% Material routine for St. Venant Kirchhoff Material 
% INPUT: 
% E: Green Lagrange Strain in Voigt notation [E11;E22;2*E12]
% dt: time increment 
% h0_gp: history from the beginning of the time step at the current Gauss Point 
% Mat_Props: contains material type and material parameters 

% OUTPUT: 
% S: 2. Piola Kirchhoff stress in Voigt notation [S11,S22,S12] 
% dSdE: Material tangent, derivative of S with respect to E
% hGP1: new trial history variables at the current Gauss Point 

function [S,dSdE,h1_GP] = st_venant_Kirchhoff2D(E,dt,h0_gp,Mat_props)

    lam = Mat_props.lambda; 
    mu = Mat_props.mu; 
    
    dSdE = [lam+2*mu, lam, 0;
            lam, lam+2*mu, 0;
            0, 0, mu];
    
    S = dSdE*E; 
    
    h1_GP = h0_gp;
end 