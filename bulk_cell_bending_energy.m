function [Elastic_energy,resting_energy, proximal_area, distal_area,bulk_splay_energy,bulk_saddel_splay_energy,rest_splay_energy,rest_saddel_splay_energy] = bulk_cell_bending_energy(Cell)
    
    %find cell cyurvature and set radius of distal and proximal monolayer accordingly 
    if (Cell.R_curv>0)
        Rc_proxy=Cell.R_curv+Cell.physical_properties.lipid_length;
        Rc_distal=Cell.R_curv-Cell.physical_properties.lipid_length;
        Rm_proxy=Cell.R_rim*(1+Cell.physical_properties.lipid_length/abs(Cell.R_curv));
        Rm_distal=Cell.R_rim*(1-Cell.physical_properties.lipid_length/abs(Cell.R_curv));
        %curvature
        J_proxy=2/Rc_proxy;
        J_distal=-2/Rc_distal;

    else
        Rc_proxy=abs(Cell.R_curv)-Cell.physical_properties.lipid_length;
        Rc_distal=abs(Cell.R_curv)+Cell.physical_properties.lipid_length;
        Rm_proxy=Cell.R_rim*(1-Cell.physical_properties.lipid_length/abs(Cell.R_curv));
        Rm_distal=Cell.R_rim*(1+Cell.physical_properties.lipid_length/abs(Cell.R_curv));
        
         %curvature
        J_proxy=-2/Rc_proxy;
        J_distal=2/Rc_distal;
    end
    

    distal_area=cell_area(Rc_distal,Rm_distal);
    proximal_area=cell_area(Rc_proxy,Rm_proxy);
    kappa_p=Cell.physical_properties.kappa_up;
    J0_p=Cell.physical_properties.J0_up;
    
    kappa_bar_u=Cell.physical_properties.kappa_bar_up;
    kappa_bar_d=Cell.physical_properties.kappa_bar_down;
    
    kappa_d=Cell.physical_properties.kappa_down;
    J0_d=Cell.physical_properties.J0_down;

    Elastic_energy=(0.5*kappa_p*(J_proxy^2-2*J_proxy*J0_p)+kappa_bar_u./Rc_proxy.^2)*proximal_area...
        +(0.5*kappa_d*(J_distal^2-2*J_distal*J0_d)+kappa_bar_d./Rc_distal.^2).*distal_area;
    
    
    bulk_splay_energy=0.5*kappa_p*(J_proxy.^2-2*J_proxy*J0_p)*proximal_area+0.5*kappa_d*(J_distal.^2-2*J_distal*J0_d).*distal_area;

    bulk_saddel_splay_energy=(kappa_bar_u./Rc_proxy.^2)*proximal_area+(kappa_bar_d./Rc_distal.^2).*distal_area;    

    %### rest
    
    % if volume is indirectly fised
    if Cell.rest_radius_if_volume_fixed~=0
        rest_MP_radius=Cell.rest_radius_if_volume_fixed;
    else % if not fixed
        rest_MP_radius=Cell.R_curv;
    end

    
    if (rest_MP_radius>0)
        Rc_proxy_rest=rest_MP_radius+Cell.physical_properties.lipid_length;
        Rc_distal_rest=rest_MP_radius-Cell.physical_properties.lipid_length;
        %curvature
        J_proxy_rest=2/Rc_proxy_rest;
        J_distal_rest=-2/Rc_distal_rest;

    else
        Rc_proxy_rest=rest_MP_radius-Cell.physical_properties.lipid_length;
        Rc_distal_rest=rest_MP_radius+Cell.physical_properties.lipid_length;
        
         %curvature
        J_proxy_rest=-2/Rc_proxy_rest;
        J_distal_rest=2/Rc_distal_rest;
    end

    
    
    
    
    proximal_area_all=4*pi*Rc_proxy_rest^2;
    distal_area_all=4*pi*Rc_distal_rest^2;
    
    resting_energy=(0.5*kappa_p*(J_proxy_rest^2-2*J_proxy_rest*J0_p)+kappa_bar_u./Rc_proxy_rest.^2).*proximal_area_all...
    +(0.5*kappa_d*(J_distal_rest^2-2*J_distal_rest*J0_d)+kappa_bar_d./Rc_distal_rest.^2).*distal_area_all;
    

    
    rest_splay_energy=0.5*kappa_p*(J_proxy_rest^2-2*J_proxy_rest*J0_p)*proximal_area_all+0.5*kappa_d*(J_distal_rest^2-2*J_distal_rest*J0_d).*distal_area_all;

    rest_saddel_splay_energy=(kappa_bar_u./Rc_proxy_rest.^2)*proximal_area_all+(kappa_bar_d./Rc_distal_rest.^2).*distal_area_all;
end