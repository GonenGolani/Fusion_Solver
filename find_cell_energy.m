function  [cell_energy_cell_area]=find_cell_energy(Cell,Diaphragm,physical_properties)
    
    if Cell.R_curv>0
        up_cell_radius=Cell.R_curv+physical_properties.lipid_length;
        down_cell_radius=Cell.R_curv-physical_properties.lipid_length;     
    end
    
    if Cell.R_curv<0
        up_cell_radius=Cell.R_curv-physical_properties.lipid_length;
        down_cell_radius=Cell.R_curv+physical_properties.lipid_length;    
    end
    
    %total cell area before fusion
    cell_area_up=4*pi*(up_cell_radius).^2;
    cell_area_down=4*pi*(down_cell_radius).^2;
    
    %area in fusion site before fusion
    cell_fusion_site_area_before_up=2*pi*up_cell_radius^2*(1-(1-up_cell_radius^2/Cell.R_curv)^0.5);
    cell_fusion_site_area_before_down=2*pi*down_cell_radius^2*(1-(1-down_cell_radius^2/Cell.R_curv)^0.5);   
   
    %area of fusion site 
    
end