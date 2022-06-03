function [Diaphragm,Cell,Virus] = save_initial_smoothing_polynom_to_matrix(Diaphragm,Cell,Virus,res_struc)
%this function saves the zero order coefficient to the smoothing matrix. the highr order can be minimaized 
if res_struc.poly_degree_z-4>0
    Diaphragm.mid_plane_smooth_polynom_matrix(1:res_struc.poly_degree_z-4,1)=Diaphragm.mid_plane_coeff_matrix(1,5:res_struc.poly_degree_z)';
    Cell.mid_plane_smooth_polynom_matrix(1:res_struc.poly_degree_z-4,1)=Cell.mid_plane_coeff_matrix(1,5:res_struc.poly_degree_z)';
    Virus.mid_plane_smooth_polynom_matrix(1:res_struc.poly_degree_z-4,1)=Virus.mid_plane_coeff_matrix(1,5:res_struc.poly_degree_z)';
end

if res_struc.poly_degree_LD_rho-4>0
    Diaphragm.up_lipid_director_rho_smooth_polynom_matrix(1:res_struc.poly_degree_LD_rho-4,1)=Diaphragm.up_lipid_director_rho_matrix(1,5:res_struc.poly_degree_LD_rho)';
    Diaphragm.down_lipid_director_rho_smooth_polynom_matrix(1:res_struc.poly_degree_LD_rho-4,1)=Diaphragm.down_lipid_director_rho_matrix(1,5:res_struc.poly_degree_LD_rho)';
    Cell.up_lipid_director_rho_smooth_polynom_matrix(1:res_struc.poly_degree_LD_rho-4,1)=Cell.up_lipid_director_rho_matrix(1,5:res_struc.poly_degree_LD_rho)';
    Cell.down_lipid_director_rho_smooth_polynom_matrix(1:res_struc.poly_degree_LD_rho-4,1)=Cell.down_lipid_director_rho_matrix(1,5:res_struc.poly_degree_LD_rho)';
    Virus.up_lipid_director_rho_smooth_polynom_matrix(1:res_struc.poly_degree_LD_rho-4,1)=Virus.up_lipid_director_rho_matrix(1,5:res_struc.poly_degree_LD_rho)';
    Virus.down_lipid_director_rho_smooth_polynom_matrix(1:res_struc.poly_degree_LD_rho-4,1)=Virus.down_lipid_director_rho_matrix(1,5:res_struc.poly_degree_LD_rho)';
end

if res_struc.poly_degree_LD_phi-4>0
    Diaphragm.up_lipid_director_phi_smooth_polynom_matrix(1:res_struc.poly_degree_LD_phi-4,1)=Diaphragm.up_lipid_director_phi_matrix(1,5:res_struc.poly_degree_LD_phi)';
    Diaphragm.down_lipid_director_phi_smooth_polynom_matrix(1:res_struc.poly_degree_LD_phi-4,1)=Diaphragm.down_lipid_director_phi_matrix(1,5:res_struc.poly_degree_LD_phi)';
    Cell.up_lipid_director_phi_smooth_polynom_matrix(1:res_struc.poly_degree_LD_phi-4,1)=Cell.up_lipid_director_phi_matrix(1,5:res_struc.poly_degree_LD_phi)';
    Cell.down_lipid_director_phi_smooth_polynom_matrix(1:res_struc.poly_degree_LD_phi-4,1)=Cell.down_lipid_director_phi_matrix(1,5:res_struc.poly_degree_LD_phi)';
    Virus.up_lipid_director_phi_smooth_polynom_matrix(1:res_struc.poly_degree_LD_phi-4,1)=Virus.up_lipid_director_phi_matrix(1,5:res_struc.poly_degree_LD_phi)';
    Virus.down_lipid_director_phi_smooth_polynom_matrix(1:res_struc.poly_degree_LD_phi-4,1)=Virus.down_lipid_director_phi_matrix(1,5:res_struc.poly_degree_LD_phi)';
end

% k vector
Diaphragm.Gaussian_decay_vector_up.k_vector_smooth_poly(1)=Diaphragm.Gaussian_decay_vector_up.k_vector(1);
Diaphragm.Gaussian_decay_vector_down.k_vector_smooth_poly(1)=Diaphragm.Gaussian_decay_vector_down.k_vector(1);
Cell.Gaussian_decay_vector_up.k_vector_smooth_poly(1)=Cell.Gaussian_decay_vector_up.k_vector(1);
Cell.Gaussian_decay_vector_down.k_vector_smooth_poly(1)=Cell.Gaussian_decay_vector_down.k_vector(1);
Virus.Gaussian_decay_vector_up.k_vector_smooth_poly(1)=Virus.Gaussian_decay_vector_up.k_vector(1);
Virus.Gaussian_decay_vector_down.k_vector_smooth_poly(1)=Virus.Gaussian_decay_vector_down.k_vector(1);

end

