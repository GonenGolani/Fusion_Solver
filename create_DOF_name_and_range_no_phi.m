function [DOF_name_and_range] = create_DOF_name_and_range_no_phi(Minimazation,res_struc)
    %trnaslate the degrees of freedom (DOF) vector to the values
    %minimaztion parametrs in of the structure 
%  rim position and general proprties
%ellipse size in x direction
int=1;
DOF_name_and_range.Diaphragm_rim_x0=int;
int=int+1; 
%ellipse size in y direction
DOF_name_and_range.Diaphragm_rim_y0=int; 
int=int+1;
% set Diaphragm hight at rho=0
DOF_name_and_range.Diaphragm_rim_z0=int;
int=int+1;

% hight of virus surface above cell surface
DOF_name_and_range.Virus_hight_above_plane=int;
int=int+1;

% fusion site size in virus side 
DOF_name_and_range.Virus_Rv=int;
int=int+1;  

% fusion site size in virus side     
DOF_name_and_range.Cell_R_rim=int;
int=int+1;


% rims
% diaphragm diaphragm rim hight
DOF_name_and_range.Diaphragm_rim_hl_vector=int:int+res_struc.poly_degree_rim-1;
int=int+res_struc.poly_degree_rim; 
% diaphragm diaphragm rim hight dirivitive
DOF_name_and_range.Diaphragm_rim_drho_dz_vector=int:int+res_struc.poly_degree_rim-1;
int=int+res_struc.poly_degree_rim; 

%virus angle at diphragm rim edge
DOF_name_and_range.Virus_rim_drho_dz_vector=int:int+res_struc.poly_degree_rim-1;
int=int+res_struc.poly_degree_rim;

%cell angle at diphragm rim edge
DOF_name_and_range.Cell_rim_drho_dz_vector=int:int+res_struc.poly_degree_rim-1;
int=int+res_struc.poly_degree_rim; 



 % tilt exponential decay
if Minimazation.exponential_decay_tilt==1
    
    DOF_name_and_range.Diaphragm_Gaussian_decay_vector_up_k_vector=int;
    int=int+1;
    
    DOF_name_and_range.Diaphragm_Gaussian_decay_vector_down_k_vector=int;
    int=int+1;
    
    DOF_name_and_range.Virus_Gaussian_decay_vector_up_k_vector=int;
    int=int+1;
    
    DOF_name_and_range.Virus_Gaussian_decay_vector_down_k_vector=int;
    int=int+1;
    
    DOF_name_and_range.Cell_Gaussian_decay_vector_up_k_vector=int;
    int=int+1;
    
    DOF_name_and_range.Cell_Gaussian_decay_vector_down_k_vector=int;
    int=int+1;
end

% mid plane

DOF_name_and_range.Diaphragm_mid_plane_coeff_matrix=int:int+res_struc.poly_degree_z-5;
int=int+res_struc.poly_degree_z-4;  

DOF_name_and_range.Virus_mid_plane_coeff_matrix=int:int+res_struc.poly_degree_z-5;
int=int+res_struc.poly_degree_z-4;

DOF_name_and_range.Cell_mid_plane_coeff_matrix=int:int+res_struc.poly_degree_z-5;
int=int+res_struc.poly_degree_z-4;

% lipid director rho

% the diaphragm has less DOF becuse axis symmetry!
DOF_name_and_range.Diaphragm_up_lipid_director_rho_matrix=int:int+res_struc.poly_degree_LD_rho-5;
int=int+res_struc.poly_degree_LD_rho-4;

DOF_name_and_range.Diaphragm_down_lipid_director_rho_matrix=int:int+res_struc.poly_degree_LD_rho-5;
int=int+res_struc.poly_degree_LD_rho-4;

DOF_name_and_range.Virus_up_lipid_director_rho_matrix=int:int+res_struc.poly_degree_LD_rho-3;
int=int+res_struc.poly_degree_LD_rho-2;

DOF_name_and_range.Virus_down_lipid_director_rho_matrix=int:int+res_struc.poly_degree_LD_rho-3;
int=int+res_struc.poly_degree_LD_rho-2;

DOF_name_and_range.Cell_up_lipid_director_rho_matrix=int:int+res_struc.poly_degree_LD_rho-3;
int=int+res_struc.poly_degree_LD_rho-2;

DOF_name_and_range.Cell_down_lipid_director_rho_matrix=int:int+res_struc.poly_degree_LD_rho-3;
int=int+res_struc.poly_degree_LD_rho-2; 

if Minimazation.Shell_exist==1

    DOF_name_and_range.Shell_Diaphragm_z_i=int;
    int=int+1;

    DOF_name_and_range.Shell_Diaphragm_z_f=int;
    int=int+1;

    DOF_name_and_range.Shell_Diaphragm_dzdrho_edge=int;
    int=int+1;

    DOF_name_and_range.Shell_Diaphragm_length_from_edge=int;
    int=int+1;
    
    % the Cell shell
    DOF_name_and_range.Shell_Cell_z_i=int;
    int=int+1;

    DOF_name_and_range.Shell_Cell_z_f=int;
    int=int+1;

    DOF_name_and_range.Shell_Cell_dzdrho_edge=int;
    int=int+1;

    DOF_name_and_range.Shell_Cell_length_from_edge=int;
    int=int+1;
    
    % max rim int
   
    %% mid plane to mid-plane coefficiant matrix 

    DOF_name_and_range.Shell_Diaphragm_MP_to_MP_distance_matrix=(int:int+res_struc.poly_degree_shell-5);    
    int=int+res_struc.poly_degree_shell-4;  
    
    DOF_name_and_range.Shell_Cell_MP_to_MP_distance_matrix=(int:int+res_struc.poly_degree_shell-5);    

end

end

