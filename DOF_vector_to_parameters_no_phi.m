function [Diaphragm,Cell,Virus,Shell] = DOF_vector_to_parameters_no_phi(Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,res_struc)
    %trnaslate the degrees of freedom (DOF) vector to the values
    % this is first step of minimaztion so polynum coefficiants are the
    % smae for all phi!
    %minimaztion parametrs in of the structure 
    DOF_name_and_range=Minimazation.DOF_name_and_range;

    int=1; %position vector
    

    Diaphragm.rim.x0=DOF_vector(DOF_name_and_range.Diaphragm_rim_x0);   % ellipse size in x direction
    Diaphragm.rim.y0=DOF_vector(DOF_name_and_range.Diaphragm_rim_y0);   % ellipse size in y direction
    Diaphragm.z0=DOF_vector(DOF_name_and_range.Diaphragm_rim_z0);       % set Diaphragm hight at rho=0 
    
    
    Virus.hight_above_plane=DOF_vector(DOF_name_and_range.Virus_hight_above_plane);    % hight of virus surface above cell surface
    
    Virus.Rv=DOF_vector(DOF_name_and_range.Virus_Rv); % fusion site size in torus         

    Cell.R_rim=DOF_vector(DOF_name_and_range.Cell_R_rim);      
    
    % max rim int
   
    %% rims
    % diaphragm diaphragm rim hight
    Diaphragm.rim.hl_vector=DOF_vector(DOF_name_and_range.Diaphragm_rim_hl_vector);      
    
    % diaphragm diaphragm rim hight dirivitive
    Diaphragm.rim.drho_dz_vector=DOF_vector(DOF_name_and_range.Diaphragm_rim_drho_dz_vector); 
    
    %virus angle at diphragm rim edge
    Virus.rim.drho_dz_vector=DOF_vector(DOF_name_and_range.Virus_rim_drho_dz_vector);     
    
    %cell angle at diphragm rim edge
    Cell.rim.drho_dz_vector=DOF_vector(DOF_name_and_range.Cell_rim_drho_dz_vector);

    
        %% exponential decay of lipid director rho
    if Minimazation.exponential_decay_tilt==1
        Diaphragm.Gaussian_decay_vector_up.k_vector(1:res_struc.phi_res)=ones(1,res_struc.phi_res)*DOF_vector(DOF_name_and_range.Diaphragm_Gaussian_decay_vector_up_k_vector); 
        
        Diaphragm.Gaussian_decay_vector_down.k_vector(1:res_struc.phi_res)=ones(1,res_struc.phi_res)*DOF_vector(DOF_name_and_range.Diaphragm_Gaussian_decay_vector_down_k_vector);
        
        Virus.Gaussian_decay_vector_up.k_vector(1:res_struc.phi_res)=ones(1,res_struc.phi_res)*DOF_vector(DOF_name_and_range.Virus_Gaussian_decay_vector_up_k_vector);        

        Virus.Gaussian_decay_vector_down.k_vector(1:res_struc.phi_res)=ones(1,res_struc.phi_res)*DOF_vector(DOF_name_and_range.Virus_Gaussian_decay_vector_down_k_vector);      

        Cell.Gaussian_decay_vector_up.k_vector(1:res_struc.phi_res)=ones(1,res_struc.phi_res)*DOF_vector(DOF_name_and_range.Cell_Gaussian_decay_vector_up_k_vector);        

        Cell.Gaussian_decay_vector_down.k_vector(1:res_struc.phi_res)=ones(1,res_struc.phi_res)*DOF_vector(DOF_name_and_range.Cell_Gaussian_decay_vector_down_k_vector);       

    end
    
    %% mid plane
    % set first 4 columbs to zero
    Diaphragm.mid_plane_coeff_matrix=zeros(res_struc.phi_res,4); 
    Virus.mid_plane_coeff_matrix=zeros(res_struc.phi_res,4);
    Cell.mid_plane_coeff_matrix=zeros(res_struc.phi_res,4);

    
    Diaphragm.mid_plane_coeff_matrix=insert_value_to_matrix_col(DOF_vector(DOF_name_and_range.Diaphragm_mid_plane_coeff_matrix), res_struc.poly_degree_z-4,res_struc.phi_res,4); 
    
    Virus.mid_plane_coeff_matrix=insert_value_to_matrix_col(DOF_vector(DOF_name_and_range.Virus_mid_plane_coeff_matrix), res_struc.poly_degree_z-4,res_struc.phi_res,4);     
    
    Cell.mid_plane_coeff_matrix=insert_value_to_matrix_col(DOF_vector(DOF_name_and_range.Cell_mid_plane_coeff_matrix), res_struc.poly_degree_z-4,res_struc.phi_res,4);

    %% lipid director rho
    Diaphragm.up_lipid_director_rho_matrix=zeros(res_struc.phi_res,4); 
    Diaphragm.down_lipid_director_rho_matrix=zeros(res_struc.phi_res,4); 
    Virus.up_lipid_director_rho_matrix=zeros(res_struc.phi_res,2); 
    Virus.down_lipid_director_rho_matrix=zeros(res_struc.phi_res,2); 
    Cell.up_lipid_director_rho_matrix=zeros(res_struc.phi_res,2); 
    Cell.down_lipid_director_rho_matrix=zeros(res_struc.phi_res,2); 
    
    Diaphragm.up_lipid_director_rho_matrix=insert_value_to_matrix_col(DOF_vector(DOF_name_and_range.Diaphragm_up_lipid_director_rho_matrix), res_struc.poly_degree_LD_rho-4,res_struc.phi_res,4); 
   
    Diaphragm.down_lipid_director_rho_matrix=insert_value_to_matrix_col(DOF_vector(DOF_name_and_range.Diaphragm_down_lipid_director_rho_matrix), res_struc.poly_degree_LD_rho-4,res_struc.phi_res,4); 
    Virus.up_lipid_director_rho_matrix=insert_value_to_matrix_col(DOF_vector(DOF_name_and_range.Virus_up_lipid_director_rho_matrix), res_struc.poly_degree_LD_rho-2,res_struc.phi_res,2); 
    Virus.down_lipid_director_rho_matrix=insert_value_to_matrix_col(DOF_vector(DOF_name_and_range.Virus_down_lipid_director_rho_matrix), res_struc.poly_degree_LD_rho-2,res_struc.phi_res,2); 
    Cell.up_lipid_director_rho_matrix=insert_value_to_matrix_col(DOF_vector(DOF_name_and_range.Cell_up_lipid_director_rho_matrix), res_struc.poly_degree_LD_rho-2,res_struc.phi_res,2); 
    Cell.down_lipid_director_rho_matrix=insert_value_to_matrix_col(DOF_vector(DOF_name_and_range.Cell_down_lipid_director_rho_matrix), res_struc.poly_degree_LD_rho-2,res_struc.phi_res,2); 
    
    %% lipid director phi
    Diaphragm.up_lipid_director_phi_matrix=zeros(res_struc.phi_res,4); 
    Diaphragm.down_lipid_director_phi_matrix=zeros(res_struc.phi_res,4); 
    Virus.up_lipid_director_phi_matrix=zeros(res_struc.phi_res,4); 
    Virus.down_lipid_director_phi_matrix=zeros(res_struc.phi_res,4); 
    Cell.up_lipid_director_phi_matrix=zeros(res_struc.phi_res,4); 
    Cell.down_lipid_director_phi_matrix=zeros(res_struc.phi_res,4); 
    


    
    %% shell
    if Shell.Shell_exist==1
    % the diaphragm shell
        Shell.Diaphragm.z_i=DOF_vector(DOF_name_and_range.Shell_Diaphragm_z_i);
        Shell.Diaphragm.z_f=DOF_vector(DOF_name_and_range.Shell_Diaphragm_z_f);
        Shell.Diaphragm.dzdrho_edge=DOF_vector(DOF_name_and_range.Shell_Diaphragm_dzdrho_edge);
        Shell.Diaphragm.length_from_edge=DOF_vector(DOF_name_and_range.Shell_Diaphragm_length_from_edge);

        % the Cell shell
        Shell.Cell.z_i=DOF_vector(DOF_name_and_range.Shell_Cell_z_i);
        Shell.Cell.z_f=DOF_vector(DOF_name_and_range.Shell_Cell_z_f); 
        Shell.Cell.dzdrho_edge=DOF_vector(DOF_name_and_range.Shell_Cell_dzdrho_edge);
        Shell.Cell.length_from_edge=DOF_vector(DOF_name_and_range.Shell_Cell_length_from_edge); 


        %% mid plane to mid-plane coefficiant matrix 

        Shell.Diaphragm.MP_to_MP_distance_matrix=insert_value_to_matrix_col(DOF_vector(DOF_name_and_range.Shell_Diaphragm_MP_to_MP_distance_matrix), res_struc.poly_degree_shell-4,res_struc.phi_res,4); 
        int=int+res_struc.poly_degree_shell-4;         
 
        Shell.Cell.MP_to_MP_distance_matrix=insert_value_to_matrix_col(DOF_vector(DOF_name_and_range.Shell_Cell_MP_to_MP_distance_matrix), res_struc.poly_degree_shell-4,res_struc.phi_res,4); 

    end

end

