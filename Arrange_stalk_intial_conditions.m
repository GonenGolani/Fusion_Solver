clc
clear all
close all
load('initial_configuration.mat');

Shell.Shell_exist=Shell_exist;
Minimazation.Shell_exist=Shell_exist;
Shell.Shell_physical_proprties=Shell_Shell_physical_proprties;

Cell.R_curv=System_dimensions.Cell_R_curv;
Virus.r_sv=System_dimensions.Virus_r_sv;
Virus.r_bv=System_dimensions.Virus_r_bv;
Diaphragm.r_pore=System_dimensions.Diaphragm_r_pore;
   
[Diaphragm,Cell,Virus,Shell] = Create_intial_coeff_matrix(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc,0);


[Minimazation.DOF_name_and_range] = create_DOF_name_and_range_all_cases(Minimazation,res_struc);

[Diaphragm,Cell,Virus,Shell] = DOF_vector_to_parameters_all_cases(Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,res_struc);



% make cell side (Shell) equal to virus side

Cell.mid_plane_coeff_matrix=-Virus.mid_plane_coeff_matrix;
Cell.up_lipid_director_rho_matrix=-Virus.up_lipid_director_rho_matrix;
Cell.down_lipid_director_rho_matrix=-Virus.down_lipid_director_rho_matrix;
Cell.Gaussian_decay_vector_up.k_vector=Virus.Gaussian_decay_vector_up.k_vector;
Cell.Gaussian_decay_vector_down.k_vector=Virus.Gaussian_decay_vector_down.k_vector;
Cell.R_rim=Virus.Rv;
Shell.Cell.length_from_edge=Virus.Rv;
Minimazation.full_lipid_flip_flop=0;


[Diaphragm,Cell,Virus,Shell] = build_Hemifusion_structure(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc,General_physical_properties);
[Diaphragm,Cell,Virus,Shell,Energy,Total_energy] = find_HD_energy(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);

plot_profile(Diaphragm,Cell,Virus,Shell,1,1,1,Minimazation.up_down_symmetry);


save_config('initial_configuration_STALK_AND_SHELL',Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,res_struc,General_physical_properties);

