function [Diaphragm,Cell,Virus] = DOF_vector_to_parameters_with_smoothing_polynom(DOF_vector,Diaphragm,Cell,Virus,Minimazation,res_struc)

    DOF_name_and_range=Minimazation.DOF_name_and_range;

    
    
    %we strat with the rims 
    %  rim position and general proprties
    Diaphragm.rim.x0=DOF_vector(DOF_name_and_range.Diaphragm_rim_x0);  % ellipse size in x direction
    Diaphragm.rim.y0=DOF_vector(DOF_name_and_range.Diaphragm_rim_y0);  % ellipse size in y direction
    Diaphragm.z0=DOF_vector(DOF_name_and_range.Diaphragm_rim_z0);      % set Diaphragm hight at rho=0 
    
    
    Virus.hight_above_plane=DOF_vector(DOF_name_and_range.Virus_hight_above_plane);    % hight of virus surface above cell surface
    
    Virus.Rv=DOF_vector(DOF_name_and_range.Virus_Rv); % fusion site size in torus         

    Cell.R_rim=DOF_vector(DOF_name_and_range.Cell_R_rim); 
    
   
    %% rims
    % diaphragm diaphragm rim hight
    Diaphragm.rim.hl_vector=DOF_vector(DOF_name_and_range.Diaphragm_rim_hl_vector);      
    % diaphragm diaphragm rim hight dirivitive
    Diaphragm.rim.drho_dz_vector=DOF_vector(DOF_name_and_range.Diaphragm_rim_drho_dz_vector); 
    
    %virus angle at diphragm rim edge
    Virus.rim.drho_dz_vector=DOF_vector(DOF_name_and_range.Virus_rim_drho_dz_vector);     
    
    %cell angle at diphragm rim edge
    Cell.rim.drho_dz_vector=DOF_vector(DOF_name_and_range.Cell_rim_drho_dz_vector);
    %%
    if Minimazation.exponential_decay_tilt==1
         Diaphragm.Gaussian_decay_vector_up.k_vector_smooth_poly=DOF_vector(DOF_name_and_range.Diaphragm_Gaussian_decay_vector_up_k_vector_smooth_poly);

         Diaphragm.Gaussian_decay_vector_down.k_vector_smooth_poly=DOF_vector(DOF_name_and_range.Diaphragm_Gaussian_decay_vector_down_k_vector_smooth_poly);

         Virus.Gaussian_decay_vector_up.k_vector_smooth_poly=DOF_vector(DOF_name_and_range.Virus_Gaussian_decay_vector_up_k_vector_smooth_poly);

         Virus.Gaussian_decay_vector_down.k_vector_smooth_poly=DOF_vector(DOF_name_and_range.Virus_Gaussian_decay_vector_down_k_vector_smooth_poly);

         Cell.Gaussian_decay_vector_up.k_vector_smooth_poly=DOF_vector(DOF_name_and_range.Cell_Gaussian_decay_vector_up_k_vector_smooth_poly);

         Cell.Gaussian_decay_vector_down.k_vector_smooth_poly=DOF_vector(DOF_name_and_range.Cell_Gaussian_decay_vector_down_k_vector_smooth_poly);

    end    

    %the DOF intial vector of all the matricies
    if res_struc.poly_degree_z-4>0
        Diaphragm.mid_plane_smooth_polynom_matrix=DOF_vector(DOF_name_and_range.MATRIX_Diaphragm_mid_plane_smooth_polynom_matrix);
        Cell.mid_plane_smooth_polynom_matrix=DOF_vector(DOF_name_and_range.MATRIX_Cell_mid_plane_smooth_polynom_matrix);
        Virus.mid_plane_smooth_polynom_matrix=DOF_vector(DOF_name_and_range.MATRIX_Virus_mid_plane_smooth_polynom_matrix);
    end
   
    %rho lipid director
    if res_struc.poly_degree_LD_rho-4>0
        Diaphragm.up_lipid_director_rho_smooth_polynom_matrix=DOF_vector(DOF_name_and_range.MATRIX_Diaphragm_up_lipid_director_rho_smooth_polynom_matrix);
        Diaphragm.down_lipid_director_rho_smooth_polynom_matrix=DOF_vector(DOF_name_and_range.MATRIX_Diaphragm_down_lipid_director_rho_smooth_polynom_matrix);
        Cell.up_lipid_director_rho_smooth_polynom_matrix=DOF_vector(DOF_name_and_range.MATRIX_Cell_up_lipid_director_rho_smooth_polynom_matrix);
        Cell.down_lipid_director_rho_smooth_polynom_matrix=DOF_vector(DOF_name_and_range.MATRIX_Cell_down_lipid_director_rho_smooth_polynom_matrix);
        Virus.up_lipid_director_rho_smooth_polynom_matrix=DOF_vector(DOF_name_and_range.MATRIX_Virus_up_lipid_director_rho_smooth_polynom_matrix);
        Virus.down_lipid_director_rho_smooth_polynom_matrix=DOF_vector(DOF_name_and_range.MATRIX_Virus_down_lipid_director_rho_smooth_polynom_matrix);
     end
        
        %phi lipid director     
     if res_struc.poly_degree_LD_phi-4>0
        Diaphragm.up_lipid_director_phi_smooth_polynom_matrix=DOF_vector(DOF_name_and_range.MATRIX_Diaphragm_up_lipid_director_phi_smooth_polynom_matrix);
        Diaphragm.down_lipid_director_phi_smooth_polynom_matrix=DOF_vector(DOF_name_and_range.MATRIX_Diaphragm_down_lipid_director_phi_smooth_polynom_matrix);
        Cell.up_lipid_director_phi_smooth_polynom_matrix=DOF_vector(DOF_name_and_range.MATRIX_Cell_up_lipid_director_phi_smooth_polynom_matrix);
        Cell.down_lipid_director_phi_smooth_polynom_matrix=DOF_vector(DOF_name_and_range.MATRIX_Cell_down_lipid_director_phi_smooth_polynom_matrix);
        Virus.up_lipid_director_phi_smooth_polynom_matrix=DOF_vector(DOF_name_and_range.MATRIX_Virus_up_lipid_director_phi_smooth_polynom_matrix);
        Virus.down_lipid_director_phi_smooth_polynom_matrix=DOF_vector(DOF_name_and_range.MATRIX_Virus_down_lipid_director_phi_smooth_polynom_matrix);
     end

    

end

