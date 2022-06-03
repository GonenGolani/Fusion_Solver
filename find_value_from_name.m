function [value] = find_value_from_name(value_to_check,Cell,Virus,Shell,General_physical_properties,Minimazation)
% this function finds the value of a varible in the workspace given its
% name

%% General_physical_properties
if strcmp(value_to_check,'General_physical_properties.tilt_modulus')
    value=General_physical_properties.tilt_modulus;
elseif strcmp(value_to_check,'General_physical_properties.chi')
    value=General_physical_properties.chi;
elseif strcmp(value_to_check,'General_physical_properties.Virus_distal_J0')
    value=General_physical_properties.Virus_distal_J0;
elseif strcmp(value_to_check,'General_physical_properties.Cell_distal_J0')
    value=General_physical_properties.Cell_distal_J0;
elseif strcmp(value_to_check,'General_physical_properties.proximal_J0')
    value=General_physical_properties.proximal_J0;
elseif strcmp(value_to_check,'General_physical_properties.cell_distal_kappa')
    value=General_physical_properties.cell_distal_kappa;
elseif strcmp(value_to_check,'General_physical_properties.virus_distal_kappa')
    value=General_physical_properties.virus_distal_kappa;
elseif strcmp(value_to_check,'General_physical_properties.proximal_kappa')
    value=General_physical_properties.proximal_kappa;
    
  
    
elseif strcmp(value_to_check,'General_physical_properties.Cell_mid_plane_physical_proprties.kappa')
    value=General_physical_properties.Cell_mid_plane_physical_proprties.kappa;
elseif strcmp(value_to_check,'General_physical_properties.Cell_mid_plane_physical_proprties.kappa_bar')
    value=General_physical_properties.Cell_mid_plane_physical_proprties.kappa_bar;    
elseif strcmp(value_to_check,'General_physical_properties.Cell_mid_plane_physical_proprties.J0')
    value=General_physical_properties.Cell_mid_plane_physical_proprties.J0;
%sorting    
elseif strcmp(value_to_check,'General_physical_properties.Cell_mid_plane_physical_proprties.sorting.phi_TM')
    value=General_physical_properties.Cell_mid_plane_physical_proprties.sorting.phi_TM;
elseif strcmp(value_to_check,'General_physical_properties.Cell_mid_plane_physical_proprties.sorting.ratio_p_to_m')
    value=General_physical_properties.Cell_mid_plane_physical_proprties.sorting.ratio_p_to_m;
elseif strcmp(value_to_check,'General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_TM')
    value=General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_TM;    
elseif strcmp(value_to_check,'General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_m')
    value=General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_m;       
elseif strcmp(value_to_check,'General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_p')
    value=General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_p; 
elseif strcmp(value_to_check,'General_physical_properties.Cell_mid_plane_physical_proprties.sorting.kappa_m_kbT_a')
    value=General_physical_properties.Cell_mid_plane_physical_proprties.sorting.kappa_m_kbT_a;   
 
elseif strcmp(value_to_check,'General_physical_properties.Cell_mid_plane_physical_proprties.sorting.affinity')
    value=General_physical_properties.Cell_mid_plane_physical_proprties.sorting.affinity;   
    
%% compartemnts size
elseif strcmp(value_to_check,'Cell.R_curv')
    value=Cell.R_curv;
elseif strcmp(value_to_check,'Cell.R_rim')
    value=Cell.R_rim;
    
    
elseif strcmp(value_to_check,'Virus.r_sv')
    value=Virus.r_sv;    
elseif strcmp(value_to_check,'Virus.r_bv')
    value=Virus.r_bv;    
%%Shell
elseif strcmp(value_to_check,'Shell.Shell_physical_proprties.width')
    value=Shell.Shell_physical_proprties.width;    
    
elseif strcmp(value_to_check,'Shell.Shell_physical_proprties.interaction_distance')
    value=Shell.Shell_physical_proprties.interaction_distance;    
elseif strcmp(value_to_check,'Shell.Shell_physical_proprties.intercation_potential_energy_density')
    value=Shell.Shell_physical_proprties.intercation_potential_energy_density;    
elseif strcmp(value_to_check,'Shell.Shell_physical_proprties.Youngs_modulus')
    value=Shell.Shell_physical_proprties.Youngs_modulus;    
elseif strcmp(value_to_check,'Shell.Shell_physical_proprties.J0')
    Shell.Shell_physical_proprties.J0;    
%%others
elseif strcmp(value_to_check,'Minimazation.fix_Cell_volume')
    value=Minimazation.fix_Cell_volume;  

elseif strcmp(value_to_check,'Minimazation.symmetric_distal_J0')
    value=Minimazation.symmetric_distal_J0;

elseif strcmp(value_to_check,'Minimazation.HD_rim_fixed_cell')
    value=Minimazation.HD_rim_fixed_cell;
      
else
    disp('ERROR in find_value_from_name, no value found')
end



end

