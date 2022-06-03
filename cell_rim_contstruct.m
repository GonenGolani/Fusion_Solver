function [rho_vector,z_vector,drho_dz,theta_rho,theta_phi] = cell_rim_contstruct(Cell,phi_res)
% this function recives the cell rim parametres and phi resolution 
% Cell containts R_rim which is the radius of the fusion site and R_curv
% which is the cell curvature

% output is: rho_vector and z_vector as fuction of phi
% mid plane angle at rim: drho_dz
% lipid director rho direction at rim: theta_rho
% lipid director phi direction at rim: theta_phi which must be pi/2!

    rho_vector=Cell.R_rim*ones(1,phi_res);

    z_vector=zeros(1,phi_res);

    if(Cell.R_curv>0)
        drho_dz=-Cell.R_rim./(Cell.R_curv^2-Cell.R_rim^2)^0.5*ones(1,phi_res);      
        theta_rho=acos(Cell.R_rim/abs(Cell.R_curv))*ones(1,phi_res);

    end
    if(Cell.R_curv<0)
        drho_dz=Cell.R_rim./(Cell.R_curv^2-Cell.R_rim^2)^0.5*ones(1,phi_res);
        theta_rho=asin(Cell.R_rim/abs(Cell.R_curv))*ones(1,phi_res)+pi/2;

    end
    
    if(strcmp(Cell.R_curv,'flat'))
        drho_dz=0*ones(1,phi_res);
        theta_rho=pi/2*ones(1,phi_res);
        
    end
    theta_phi=pi/2*ones(1,phi_res);
end




