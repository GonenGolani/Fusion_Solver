function [Diaphragm,Cell,Virus,Shell,Energy,Total_energy] = find_HD_energy(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc)
    
    %  cell fusion site elastic energy
    [Energy.cell_FS,Cell.total_area_up,Cell.total_area_down,Cell.energy_components]...
    = find_membrane_energy(Cell.tilt,Cell.area_element,Cell.splay,Cell.saddle_splay,Cell.physical_properties);



    if Minimazation.Volume_fixed_based_on_Rc==1
        [Cell,Diaphragm] = recalc_only_cell_mid_plane(Cell,Diaphragm,Minimazation,res_struc);
        Cell.rest_radius_if_volume_fixed=(3*find_Cell_volume(Cell,Diaphragm,res_struc)/(4*pi))^(1/3);
    else
        Cell.rest_radius_if_volume_fixed=0;
    end

    [Energy.cell_bulk, Energy.cell_resting_energy , Cell_proximal_area, Cell_distal_area, Energy.Cell_bulk_splay_energy,...
        Energy.Cell_bulk_saddel_splay_energy,Energy.Cell_rest_splay_energy,Energy.cell_rest_saddel_splay_energy ]...
        = bulk_cell_bending_energy(Cell);

    % if there is up-down symmetry the cell and virus energy are the same
    if Minimazation.up_down_symmetry==1 
        Energy.virus_FS=Energy.cell_FS;
        Virus.total_area_up=Cell.total_area_up;
        Virus.total_area_down=Cell.total_area_down;
        Virus.energy_components=Cell.energy_components;
    
        Energy.virus_bulk=Energy.cell_bulk;
        virus_bulk_distal_area=Cell_distal_area;
        virus_bulk_proximal_area=Cell_proximal_area;
        Energy.virus_resting_energy=Energy.cell_resting_energy;

    else %no up-down symmetry- use 


    %  virus fusion site elastic energy
    [Energy.virus_FS,Virus.total_area_up,Virus.total_area_down,Virus.energy_components]...
    =find_membrane_energy(Virus.tilt,Virus.area_element,Virus.splay,Virus.saddle_splay,Virus.physical_properties);


    %virus bulk elastic energy and area
    [Energy.virus_bulk,virus_bulk_distal_area,virus_bulk_proximal_area,Energy.virus_resting_energy]=bulk_virus_bending_energy_area(Virus,res_struc);
    end

%  diaphragm fusion site elastic energy
    if Minimazation.do_stalk==0
    [Energy.diaphragm_FS,Diaphragm.total_area_up,Diaphragm.total_area_down,Diaphragm.energy_components] ...
    =find_membrane_energy(Diaphragm.tilt,Diaphragm.area_element,Diaphragm.splay,Diaphragm.saddle_splay,Diaphragm.physical_properties);
    else  %if the diaphragm does not exist
        Energy.diaphragm_FS=0;
        Diaphragm.total_area_up=0;
        Diaphragm.total_area_down=0;
        Diaphragm.energy_components.splay_energy=0;
        Diaphragm.energy_components.saddle_splay_energy=0;
        Diaphragm.energy_components.tilt_energy=0;
        Diaphragm.energy_components.splay_grad_square_energy=0;
        Diaphragm.energy_components.splay_grad_tilt_energy=0;
        Diaphragm.energy_components.Energy_density=0;
    end
% virus and cell outside of fusion site energy and area




%stretching energy
if Cell.physical_properties.Ks~=0

    [Energy.stretching,Energy.stretching_proximal_energy_density,Energy.stretching_distal_energy_density_cell,Energy.stretching_distal_energy_density_virus]=...
       stretching_energy(Cell,Virus,Diaphragm,virus_bulk_distal_area,virus_bulk_proximal_area,Minimazation);

else
    Energy.stretching=0;
    Energy.stretching_proximal_energy_density=0;
    Energy.stretching_distal_energy_density_cell=0;
    Energy.stretching_distal_energy_density_virus=0;
end


%overall energy
Total_energy=Energy.cell_FS+Energy.virus_FS+Energy.diaphragm_FS+Energy.virus_bulk+Energy.cell_bulk+Energy.stretching-Energy.cell_resting_energy-Energy.virus_resting_energy;
Energy.reff_energy=Energy.cell_resting_energy+Energy.virus_resting_energy;

Diaphragm.Energy_density=Diaphragm.energy_components.Energy_density+Energy.stretching_distal_energy_density_cell+Energy.stretching_distal_energy_density_virus;
Virus.Energy_density=Virus.energy_components.Energy_density+Energy.stretching_proximal_energy_density+Energy.stretching_distal_energy_density_virus;
Cell.Energy_density=Cell.energy_components.Energy_density+Energy.stretching_proximal_energy_density+Energy.stretching_distal_energy_density_cell;

%energy components

Energy.total_change_splay_energy=Cell.energy_components.splay_energy+Virus.energy_components.splay_energy+Diaphragm.energy_components.splay_energy+...
    Energy.Cell_bulk_splay_energy-Energy.Cell_rest_splay_energy...
    +Energy.virus_bulk-Energy.virus_resting_energy;

Energy.total_change_saddel_splay_energy=Cell.energy_components.saddle_splay_energy+...
    Virus.energy_components.saddle_splay_energy+Diaphragm.energy_components.saddle_splay_energy+...
    Energy.Cell_bulk_saddel_splay_energy-Energy.cell_rest_saddel_splay_energy;

%shell if exist
if Shell.Shell_exist==1
[Shell,Shell.Total_Shell_energy] = Total_shell_energy(Shell);
Total_energy=Total_energy+Shell.Total_Shell_energy;
end

% add global surface tension
if Cell.physical_properties.Surface_tension ~=0 % now the surface tension is equal in all monolayers!
     if Minimazation.do_stalk==1 && Minimazation.up_down_symmetry==1 %simple calcuation
         x_vec_up=Cell.dividing_plane.up.x(1,:);
           z_vec_up=Cell.dividing_plane.up.z(1,:);
           dx_up=diff(x_vec_up);
           dz_up=diff(z_vec_up);
           dx_up=[dx_up dx_up(length(dx_up))];
           dz_up=[dz_up dz_up(length(dz_up))];
           Up_area_change=sum(2*pi*(dz_up.^2+dx_up.^2).^0.5.*x_vec_up-2*pi*dx_up.*x_vec_up);

           x_vec_down=Cell.dividing_plane.down.x(1,:);
           z_vec_down=Cell.dividing_plane.down.z(1,:);
           dx_down=diff(x_vec_down);
           dz_down=diff(z_vec_down);
           dx_down=[dx_down dx_down(length(dx_down))];
           dz_down=[dz_down dz_down(length(dz_down))];

           Down_area_change=sum(2*pi*(dz_down.^2+dx_down.^2).^0.5.*x_vec_down-2*pi*dx_down.*x_vec_down);

           area_change= 2*(Down_area_change+Up_area_change);
     else
        [~,~,~, area_change]= area_change_func(Cell,Virus,Diaphragm,virus_bulk_distal_area,virus_bulk_proximal_area,Minimazation);
     end

     Total_energy=Total_energy+Cell.physical_properties.Surface_tension*area_change;

end

% if added trans-membrane proprties to cell
if Minimazation.Cell_trans_membrane_added_proprties_exist==1
    [Cell.mid_plane.Cell_trans_membrane_energy] = mid_plane_energy(Cell.mid_plane);
    Energy.cell_trans_membrane_energy=Cell.mid_plane.Cell_trans_membrane_energy;
    Total_energy=Total_energy+Cell.mid_plane.Cell_trans_membrane_energy;
end

end

