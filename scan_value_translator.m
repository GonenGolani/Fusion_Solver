function [Cell_new,Virus_new,Shell_new,General_physical_properties_new,Minimazation_new] = scan_value_translator(value_to_change,new_value,Cell_old,Virus_old,Shell_old,General_physical_properties_old,Minimazation_old)
% In order to avoid dynamical programing (evalin, eval) this function takes
% a proprtie indicated by string value_to_change and places the new_value

%reset struct
Cell_new=Cell_old;
Virus_new=Virus_old;
Shell_new=Shell_old;
General_physical_properties_new=General_physical_properties_old;
Minimazation_new=Minimazation_old;

%% General_physical_properties
if strcmp(value_to_change,'General_physical_properties.tilt_modulus')
    General_physical_properties_new.tilt_modulus=new_value;
elseif strcmp(value_to_change,'General_physical_properties.chi')
    General_physical_properties_new.chi=new_value;
elseif strcmp(value_to_change,'General_physical_properties.Virus_distal_J0')
    General_physical_properties_new.Virus_distal_J0=new_value;
elseif strcmp(value_to_change,'General_physical_properties.Cell_distal_J0')
    General_physical_properties_new.Cell_distal_J0=new_value;
elseif strcmp(value_to_change,'General_physical_properties.proximal_J0')
    General_physical_properties_new.proximal_J0=new_value;
elseif strcmp(value_to_change,'General_physical_properties.cell_distal_kappa')
    General_physical_properties_new.cell_distal_kappa=new_value;
elseif strcmp(value_to_change,'General_physical_properties.virus_distal_kappa')
    General_physical_properties_new.virus_distal_kappa=new_value;
elseif strcmp(value_to_change,'General_physical_properties.proximal_kappa')
    General_physical_properties_new.proximal_kappa=new_value;
    
elseif strcmp(value_to_change,'General_physical_properties.Cell_mid_plane_physical_proprties.kappa')
    General_physical_properties_new.Cell_mid_plane_physical_proprties.kappa=new_value;
elseif strcmp(value_to_change,'General_physical_properties.Cell_mid_plane_physical_proprties.kappa_bar')
    General_physical_properties_new.Cell_mid_plane_physical_proprties.kappa_bar=new_value;    
elseif strcmp(value_to_change,'General_physical_properties.Cell_mid_plane_physical_proprties.J0')
    General_physical_properties_new.Cell_mid_plane_physical_proprties.J0=new_value;
%sorting    
elseif strcmp(value_to_change,'General_physical_properties.Cell_mid_plane_physical_proprties.sorting.phi_TM')
    General_physical_properties_new.Cell_mid_plane_physical_proprties.sorting.phi_TM=new_value;
elseif strcmp(value_to_change,'General_physical_properties.Cell_mid_plane_physical_proprties.sorting.ratio_p_to_m')
    General_physical_properties_new.Cell_mid_plane_physical_proprties.sorting.ratio_p_to_m=new_value;
elseif strcmp(value_to_change,'General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_TM')
    General_physical_properties_new.Cell_mid_plane_physical_proprties.sorting.j_TM=new_value;    
elseif strcmp(value_to_change,'General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_m')
    General_physical_properties_new.Cell_mid_plane_physical_proprties.sorting.j_m=new_value;       
elseif strcmp(value_to_change,'General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_p')
    General_physical_properties_new.Cell_mid_plane_physical_proprties.sorting.j_p=new_value; 
elseif strcmp(value_to_change,'General_physical_properties.Cell_mid_plane_physical_proprties.sorting.kappa_m_kbT_a')
    General_physical_properties_new.Cell_mid_plane_physical_proprties.sorting.kappa_m_kbT_a=new_value;   
          
%% compartemnts size
elseif strcmp(value_to_change,'Cell.R_curv')
    Cell_new.R_curv=new_value;
elseif strcmp(value_to_change,'Cell.R_rim')
    Cell_new.R_rim=new_value;
    
    
elseif strcmp(value_to_change,'Virus.r_sv')
    Virus_new.r_sv=new_value;    
elseif strcmp(value_to_change,'Virus.r_bv')
    Virus_new.r_bv=new_value;    
%%Shell
elseif strcmp(value_to_change,'Shell.Shell_physical_proprties.width')
    Shell_new.Shell_physical_proprties.width=new_value;    
    
elseif strcmp(value_to_change,'Shell.Shell_physical_proprties.interaction_distance')
    Shell_new.Shell_physical_proprties.interaction_distance=new_value;    
elseif strcmp(value_to_change,'Shell.Shell_physical_proprties.intercation_potential_energy_density')
    Shell_new.Shell_physical_proprties.intercation_potential_energy_density=new_value;    
elseif strcmp(value_to_change,'Shell.Shell_physical_proprties.Youngs_modulus')
    Shell_new.Shell_physical_proprties.Youngs_modulus=new_value;    
elseif strcmp(value_to_change,'Shell.Shell_physical_proprties.J0')
    Shell_new.Shell_physical_proprties.J0=new_value;    
%%others
elseif strcmp(value_to_change,'Minimazation.fix_Cell_volume')
    Minimazation_new.fix_Cell_volume=new_value;  

elseif strcmp(value_to_change,'Minimazation.symmetric_distal_J0')
    Minimazation_new.symmetric_distal_J0=new_value;  

elseif strcmp(value_to_change,'Minimazation.HD_rim_fixed_cell')
    Minimazation_new.HD_rim_fixed_cell=new_value;  
    
else
    disp('ERROR in scan_value_translator, no value assigned')
end






end

