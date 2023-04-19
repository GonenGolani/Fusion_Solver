function [proximal_area_change,cell_distal_area_change,virus_distal_area_change,sum_area_change]=...
    area_change_func(Cell,Virus,Diaphragm,virus_distal_area,virus_proximal_area,Minimazation)

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
    


    proximal_initial_area=4*pi*Rc_proxy^2+4*pi^2*(Virus.r_sv+Virus.physical_properties.lipid_length)*Virus.r_bv;
    
    [cell_area_outside_of_fusion_site] = cell_area(Rc_proxy,Rm_proxy);
    
    current_proximal_area=Virus.total_area_down+Cell.total_area_up+virus_proximal_area+cell_area_outside_of_fusion_site;


    initial_virus_distal_monolayer=4*pi^2*(Virus.r_sv-Virus.physical_properties.lipid_length)*Virus.r_bv;
    currnet_virus_distal_monolayer=virus_distal_area+Virus.total_area_up+Diaphragm.total_area_up;


    initial_cell_distal_monolayer=4*pi*Rc_distal^2;
    current_cell_distal_monolayer=cell_area(Rc_distal,Rm_distal)+Cell.total_area_down+Diaphragm.total_area_down;


    if Minimazation.up_down_symmetry==1 
        initial_virus_distal_monolayer=initial_cell_distal_monolayer;
        currnet_virus_distal_monolayer=current_cell_distal_monolayer;
        proximal_initial_area=2*4*pi*Rc_proxy^2;
        current_proximal_area=2*(Cell.total_area_up+cell_area_outside_of_fusion_site);
    end

    
    proximal_area_change=current_proximal_area-proximal_initial_area;
    cell_distal_area_change=current_cell_distal_monolayer-initial_cell_distal_monolayer;
    virus_distal_area_change=currnet_virus_distal_monolayer-initial_virus_distal_monolayer;
    sum_area_change=proximal_area_change+cell_distal_area_change+virus_distal_area_change;

end