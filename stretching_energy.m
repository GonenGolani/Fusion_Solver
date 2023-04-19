function [stretching_energy,proximal_energy_density,distal_energy_density_cell,distal_energy_density_virus]=...
    stretching_energy(Cell,Virus,Diaphragm,virus_distal_area,virus_proximal_area,Minimazation)

    %find cell cyurvature and set radius of distal and proximal monolayer accordingly 
    if (Cell.R_curv>0)
        Rc_proxy=Cell.R_curv+Cell.physical_properties.lipid_length;
        Rc_distal=Cell.R_curv-Cell.physical_properties.lipid_length;
        Rm_proxy=Cell.R_rim*(1+Cell.physical_properties.lipid_length/abs(Cell.R_curv));
        Rm_distal=Cell.R_rim*(1-Cell.physical_properties.lipid_length/abs(Cell.R_curv));
    else
        Rc_proxy=abs(Cell.R_curv)-Cell.physical_properties.lipid_length;
        Rc_distal=abs(Cell.R_curv)+Cell.physical_properties.lipid_length;
        Rm_proxy=Cell.R_rim*(1-Cell.physical_properties.lipid_length/abs(Cell.R_curv));
        Rm_distal=Cell.R_rim*(1+Cell.physical_properties.lipid_length/abs(Cell.R_curv));
    end
    

    if Minimazation.fix_Cell_volume~=0 %if the cell volume is fixed the refrence radius is the rest one
        Rc_proxy_fixed=(Minimazation.fix_Cell_volume*3/(4*pi)).^(1/3)+Cell.physical_properties.lipid_length;
        proximal_initial_area=4*pi*Rc_proxy_fixed^2+4*pi^2*(Virus.r_sv+Virus.physical_properties.lipid_length)*Virus.r_bv;
    elseif Cell.rest_radius_if_volume_fixed~=0 % if the radius and volume are fixed
        Rc_proxy_fixed=Cell.rest_radius_if_volume_fixed+Cell.physical_properties.lipid_length;
        proximal_initial_area=4*pi*Rc_proxy_fixed^2+4*pi^2*(Virus.r_sv+Virus.physical_properties.lipid_length)*Virus.r_bv;
    else
        proximal_initial_area=4*pi*Rc_proxy^2+4*pi^2*(Virus.r_sv+Virus.physical_properties.lipid_length)*Virus.r_bv;
    end

    [cell_area_outside_of_fusion_site] = cell_area(Rc_proxy,Rm_proxy);
    
    current_proximal_area=Virus.total_area_down+Cell.total_area_up+virus_proximal_area+cell_area_outside_of_fusion_site;


    initial_virus_distal_monolayer=4*pi^2*(Virus.r_sv-Virus.physical_properties.lipid_length)*Virus.r_bv;
    currnet_virus_distal_monolayer=virus_distal_area+Virus.total_area_up+Diaphragm.total_area_up;

    if Minimazation.fix_Cell_volume~=0 %if the cell volume is fixed the refrence radius is the rest one
        Rc_distal_fixed=(Minimazation.fix_Cell_volume*3/(4*pi)).^(1/3)-Cell.physical_properties.lipid_length;
        initial_cell_distal_monolayer=4*pi*Rc_distal_fixed^2;
    elseif Cell.rest_radius_if_volume_fixed~=0 % if the radius and volume are fixed
        Rc_distal_fixed=Cell.rest_radius_if_volume_fixed-Cell.physical_properties.lipid_length;    
        initial_cell_distal_monolayer=4*pi*Rc_distal_fixed^2;
    else
        initial_cell_distal_monolayer=4*pi*Rc_distal^2;
    end
    current_cell_distal_monolayer=cell_area(Rc_distal,Rm_distal)+Cell.total_area_down+Diaphragm.total_area_down;


    if Minimazation.up_down_symmetry==1 
        initial_virus_distal_monolayer=initial_cell_distal_monolayer;
        currnet_virus_distal_monolayer=current_cell_distal_monolayer;
        proximal_initial_area=2*4*pi*Rc_proxy^2;
        current_proximal_area=2*(Cell.total_area_up+cell_area_outside_of_fusion_site);
    end

    if Minimazation.full_lipid_flip_flop==1

    total_initial_area=current_proximal_area+initial_virus_distal_monolayer+initial_cell_distal_monolayer;
    total_current_area=current_proximal_area+currnet_virus_distal_monolayer+current_cell_distal_monolayer;

    energy_density=0.5*Cell.physical_properties.Ks*((total_initial_area-total_current_area)./initial_cell_distal_monolayer)^2;
    stretching_energy=energy_density*total_initial_area;
    proximal_energy_density=energy_density;
    distal_energy_density_cell=energy_density;
    distal_energy_density_virus=energy_density;


    else
    % proximal monolayer

    proximal_energy_density=0.5*Cell.physical_properties.Ks*((proximal_initial_area-current_proximal_area)/proximal_initial_area).^2;
    proximal_energy=proximal_energy_density.*proximal_initial_area;
    
    
    % distal monolayrs
    distal_energy_density_virus=0.5*Virus.physical_properties.Ks*((initial_virus_distal_monolayer-currnet_virus_distal_monolayer)/initial_virus_distal_monolayer).^2;
    virus_distal_energy=distal_energy_density_virus*initial_virus_distal_monolayer;
    
    distal_energy_density_cell=0.5*Cell.physical_properties.Ks*((initial_cell_distal_monolayer-current_cell_distal_monolayer)/initial_cell_distal_monolayer).^2;

    cell_distal_energy=distal_energy_density_cell*initial_cell_distal_monolayer;
  
    stretching_energy=cell_distal_energy+virus_distal_energy+proximal_energy;
    end 



end