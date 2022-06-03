function [Diaphragm,Cell,Virus,Shell] = Create_intial_coeff_matrix(Diaphragm_old,Cell_old,Virus_old,Shell_old,Minimazation,res_struc,intital_tilt_decay_length)

    Diaphragm=Diaphragm_old;
    Cell=Cell_old;
    Virus=Virus_old;
    Shell=Shell_old;    


    Diaphragm.mid_plane_coeff_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_z); 
    Cell.mid_plane_coeff_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_z); 
    Virus.mid_plane_coeff_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_z);

    Diaphragm.up_lipid_director_rho_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_LD_rho); 
    Diaphragm.down_lipid_director_rho_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_LD_rho); 
    Diaphragm.up_lipid_director_phi_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_LD_phi);
    Diaphragm.down_lipid_director_phi_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_LD_phi);

    Cell.up_lipid_director_rho_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_LD_rho); 
    Cell.down_lipid_director_rho_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_LD_rho); 
    Cell.up_lipid_director_phi_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_LD_phi);
    Cell.down_lipid_director_phi_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_LD_phi);

    Virus.up_lipid_director_rho_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_LD_rho); 
    Virus.down_lipid_director_rho_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_LD_rho); 
    Virus.up_lipid_director_phi_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_LD_phi);
    Virus.down_lipid_director_phi_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_LD_phi);
    
    if Shell.Shell_exist==1
        Shell.Diaphragm.MP_to_MP_distance_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_shell); 
        Shell.Cell.MP_to_MP_distance_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_shell); 
        Shell.Junction.mid_plane_coeff_matrix=zeros(res_struc.phi_res,4); 

    end
    if Minimazation.non_axial_symmetric_polynom==1
        if res_struc.poly_degree_z-4>0
            Diaphragm.mid_plane_smooth_polynom_matrix=zeros(res_struc.poly_degree_z-4,res_struc.polar_angle_polynom);
            Cell.mid_plane_smooth_polynom_matrix=zeros(res_struc.poly_degree_z-4,res_struc.polar_angle_polynom);
            Virus.mid_plane_smooth_polynom_matrix=zeros(res_struc.poly_degree_z-4,res_struc.polar_angle_polynom);
        else
            Diaphragm.mid_plane_smooth_polynom_matrix=[];
            Cell.mid_plane_smooth_polynom_matrix=[];
            Virus.mid_plane_smooth_polynom_matrix=[];

        end
        
        if res_struc.poly_degree_LD_rho-4>0
            Diaphragm.up_lipid_director_rho_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_rho-4,res_struc.polar_angle_polynom);
            Diaphragm.down_lipid_director_rho_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_rho-4,res_struc.polar_angle_polynom);   
            Cell.up_lipid_director_rho_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_rho-4,res_struc.polar_angle_polynom);
            Cell.down_lipid_director_rho_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_rho-4,res_struc.polar_angle_polynom);
            Virus.up_lipid_director_rho_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_rho-4,res_struc.polar_angle_polynom);
            Virus.down_lipid_director_rho_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_rho-4,res_struc.polar_angle_polynom); 
        else
            Diaphragm.up_lipid_director_rho_smooth_polynom_matrix=[];
            Diaphragm.down_lipid_director_rho_smooth_polynom_matrix=[];  
            Cell.up_lipid_director_rho_smooth_polynom_matrix=[];
            Cell.down_lipid_director_rho_smooth_polynom_matrix=[];
            Virus.up_lipid_director_rho_smooth_polynom_matrix=[];
            Virus.down_lipid_director_rho_smooth_polynom_matrix=[];           
        end
        
        if res_struc.poly_degree_LD_phi-4>0
            Diaphragm.up_lipid_director_phi_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_phi-4,res_struc.polar_angle_polynom);
            Diaphragm.down_lipid_director_phi_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_phi-4,res_struc.polar_angle_polynom); 
            Cell.up_lipid_director_phi_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_phi-4,res_struc.polar_angle_polynom);
            Cell.down_lipid_director_phi_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_phi-4,res_struc.polar_angle_polynom); 
            Virus.up_lipid_director_phi_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_phi-4,res_struc.polar_angle_polynom);
            Virus.down_lipid_director_phi_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_phi-4,res_struc.polar_angle_polynom); 
        else
            Diaphragm.up_lipid_director_phi_smooth_polynom_matrix=[]; 
            Diaphragm.down_lipid_director_phi_smooth_polynom_matrix=[];
            Cell.up_lipid_director_phi_smooth_polynom_matrix=[];
            Cell.down_lipid_director_phi_smooth_polynom_matrix=[]; 
            Virus.up_lipid_director_phi_smooth_polynom_matrix=[];
            Virus.down_lipid_director_phi_smooth_polynom_matrix=[];
        end
    end
    if Minimazation.exponential_decay_tilt==1
        Diaphragm.Gaussian_decay_vector_up.k_vector=ones(res_struc.phi_res,1)*intital_tilt_decay_length;
        Diaphragm.Gaussian_decay_vector_down.k_vector=ones(res_struc.phi_res,1)*intital_tilt_decay_length;
        Virus.Gaussian_decay_vector_up.k_vector=ones(res_struc.phi_res,1)*intital_tilt_decay_length;
        Virus.Gaussian_decay_vector_down.k_vector=ones(res_struc.phi_res,1)*intital_tilt_decay_length;
        Cell.Gaussian_decay_vector_up.k_vector=ones(res_struc.phi_res,1)*intital_tilt_decay_length;
        Cell.Gaussian_decay_vector_down.k_vector=ones(res_struc.phi_res,1)*intital_tilt_decay_length;
        
        if Minimazation.non_axial_symmetric_polynom==1
            Diaphragm.Gaussian_decay_vector_up.k_vector_smooth_poly=zeros(1,res_struc.polar_angle_polynom);
            Diaphragm.Gaussian_decay_vector_down.k_vector_smooth_poly=zeros(1,res_struc.polar_angle_polynom);
            Virus.Gaussian_decay_vector_up.k_vector_smooth_poly=zeros(1,res_struc.polar_angle_polynom);
            Virus.Gaussian_decay_vector_down.k_vector_smooth_poly=zeros(1,res_struc.polar_angle_polynom);
            Cell.Gaussian_decay_vector_up.k_vector_smooth_poly=zeros(1,res_struc.polar_angle_polynom);
            Cell.Gaussian_decay_vector_down.k_vector_smooth_poly=zeros(1,res_struc.polar_angle_polynom);
            
            Diaphragm.Gaussian_decay_vector_up.k_vector_smooth_poly(1)=intital_tilt_decay_length;
            Diaphragm.Gaussian_decay_vector_down.k_vector_smooth_poly(1)=intital_tilt_decay_length;
            Virus.Gaussian_decay_vector_up.k_vector_smooth_poly(1)=intital_tilt_decay_length;
            Virus.Gaussian_decay_vector_down.k_vector_smooth_poly(1)=intital_tilt_decay_length;
            Cell.Gaussian_decay_vector_up.k_vector_smooth_poly(1)=intital_tilt_decay_length;
            Cell.Gaussian_decay_vector_down.k_vector_smooth_poly(1)=intital_tilt_decay_length;
            
        end
        
    end
end

