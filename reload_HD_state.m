function [Diaphragm,Cell,Virus,Shell,Minimazation,res_struc,General_physical_properties,Energy,Total_energy,DOF_vector,System_dimensions] = reload_HD_state(name,res_struc_old)
%this function takes what is saved in HD_state_thin.mat and reload the configuration

%if res_struc_old exist us it
    

    load(name);
    
    if exist('res_struc_old','var')
      res_struc.phi_res=res_struc_old.phi_res;
      res_struc.rho_res=res_struc_old.rho_res;
        
    end        

    %% backward competability
    
    if isfield(Minimazation,'tension_distance_coupleing')==0
        Minimazation.tension_distance_coupleing.exist=0;
        Minimazation.tension_distance_coupleing.kappa_b=0;
        Minimazation.tension_distance_coupleing.Pext=0;
    end

    if isfield(General_physical_properties,'Surface_tension')==0
        General_physical_properties.Surface_tension=0;
    end
    
    if isfield(General_physical_properties,'Cell_mid_plane_physical_proprties')==0
        General_physical_properties.Cell_mid_plane_physical_proprties.sorting=0;
    else
        if isfield(General_physical_properties.Cell_mid_plane_physical_proprties,'sorting')
            if isfield(General_physical_properties.Cell_mid_plane_physical_proprties.sorting,'j_TM_luminal')==0
                General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_TM_luminal=0; 
            end
        end
    end


    
    if isfield(Minimazation,'Volume_fixed_based_on_Rc')==0
        Minimazation.Volume_fixed_based_on_Rc=0;
    end

    if isfield(Minimazation,'full_lipid_flip_flop')==0
        Minimazation.full_lipid_flip_flop=0; 
    end

    
    if isfield(Minimazation,'Cell_trans_membrane_added_proprties_exist')==0
        Minimazation.Cell_trans_membrane_added_proprties_exist=0;
    end
    if isfield(Minimazation,'number_of_repeats_each_step')==0
        Minimazation.number_of_repeats_each_step=1;
    end
    if isfield(Minimazation,'sorting_protein_in_cell_membrane')==0
        Minimazation.sorting_protein_in_cell_membrane=0;
        
    end
    if isfield(Minimazation,'cell_diaphragm_J0_as_virus_distal')==1
        Minimazation=rmfield(Minimazation,'cell_diaphragm_J0_as_virus_distal');
    end
    if isfield(Minimazation,'fix_Cell_volume')==0
        Minimazation.fix_Cell_volume=0;
    end
    
    if isfield(Minimazation,'fix_inter_membrane_distance')==0
        Minimazation.fix_inter_membrane_distance=0;
    end
    
    if isfield(Minimazation,'do_stalk')==0
        Minimazation.do_stalk=0;
    end

    if isfield(Minimazation,'no_TM_direct_intrinsic_curvature')==0
        Minimazation.no_TM_direct_intrinsic_curvature=0;
    end

    if isfield(Minimazation,'no_sorting')==0
        Minimazation.no_sorting=0;
    end
    

    if Minimazation.sorting_protein_in_cell_membrane==1
        if isfield(General_physical_properties.Cell_mid_plane_physical_proprties.sorting,'full_lipid_flip_flop')==0
                General_physical_properties.Cell_mid_plane_physical_proprties.sorting.full_lipid_flip_flop=0;
        end
    end
    
    if Shell_exist==1
        if isfield(Shell_Shell_physical_proprties,'max_shell_curvature')==0 
            Shell_Shell_physical_proprties.max_shell_curvature=Shell_Shell_physical_proprties.width/2; % max possible value
        end

        if isfield(Shell_Shell_physical_proprties,'P_ratio')==0
            Shell_Shell_physical_proprties.P_ratio=0.5; % like water
        end
    end
    
    if isfield(Minimazation.scan_variable,'delta_value_2')==0
        Minimazation.scan_variable.delta_value_2=0;
        Minimazation.scan_variable.target_value_2=0;
        Minimazation.scan_variable.delta_value_3=0;
        Minimazation.scan_variable.target_value_3=0;
    end
    if isfield(Minimazation,'DOF_name_and_range')==0
        Minimazation.DOF_name_and_range=[];
    end

    %

        
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
    

    [DOF_vector,Minimazation] =  Parameters_to_DOF_all_cases(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);
    
    
    set_back_to_zero=0;
    if General_physical_properties.Splay_grad_square_modulus==0
        General_physical_properties.Splay_grad_square_modulus=1;
        General_physical_properties.Splay_grad_tilt_modulus=1;
        set_back_to_zero=1;
    end

    
    [Diaphragm,Cell,Virus,Shell] = build_Hemifusion_structure(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc,General_physical_properties);
    if set_back_to_zero==1
        
        General_physical_properties.Splay_grad_square_modulus=0;
        General_physical_properties.Splay_grad_tilt_modulus=0;
        Virus.physical_properties.Splay_grad_square_modulus=General_physical_properties.Splay_grad_square_modulus;
        Virus.physical_properties.Splay_grad_tilt_modulus=General_physical_properties.Splay_grad_tilt_modulus;
        Cell.physical_properties.Splay_grad_square_modulus=General_physical_properties.Splay_grad_square_modulus;
        Cell.physical_properties.Splay_grad_tilt_modulus=General_physical_properties.Splay_grad_tilt_modulus;
        Diaphragm.physical_properties.Splay_grad_square_modulus=General_physical_properties.Splay_grad_square_modulus;
        Diaphragm.physical_properties.Splay_grad_tilt_modulus=General_physical_properties.Splay_grad_tilt_modulus;
    end
    [Diaphragm,Cell,Virus,Shell,Energy,Total_energy] = find_HD_energy(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);

    

end

