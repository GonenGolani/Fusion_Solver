function [DOF_vector,DOF_name_and_range] = parameters_to_DOF_vector_no_phi(Diaphragm,Cell,Virus,Minimazation,res_struc)
    %trnaslate the degrees of freedom (DOF) vector to the values
    %minimaztion parametrs in of the structure 
%  rim position and general proprties
%ellipse size in x direction
int=1;
DOF_name_and_range.Diaphragm_rim_x0=int;
DOF_vector(int)=Diaphragm.rim.x0;  
int=int+1; 
%ellipse size in y direction
DOF_name_and_range.Diaphragm_rim_y0=int; 
DOF_vector(int)=Diaphragm.rim.y0;
int=int+1;
% set Diaphragm hight at rho=0
DOF_name_and_range.Diaphragm_rim_z0=int;
DOF_vector(int)=Diaphragm.z0;   
int=int+1;

% hight of virus surface above cell surface
DOF_name_and_range.Virus_hight_above_plane=int;
DOF_vector(int)=Virus.hight_above_plane;    
int=int+1;

% fusion site size in virus side 
DOF_name_and_range.Virus_Rv=int;
DOF_vector(int)=Virus.Rv;
int=int+1;  

% fusion site size in virus side     
DOF_name_and_range.Cell_R_rim=int;
DOF_vector(int)=Cell.R_rim; 
int=int+1;


% rims
% diaphragm diaphragm rim hight
DOF_name_and_range.Diaphragm_rim_hl_vector=int:int+res_struc.poly_degree_rim-1;
DOF_vector(int:int+res_struc.poly_degree_rim-1)=Diaphragm.rim.hl_vector;      
int=int+res_struc.poly_degree_rim; 
% diaphragm diaphragm rim hight dirivitive
DOF_name_and_range.Diaphragm_rim_drho_dz_vector=int:int+res_struc.poly_degree_rim-1;
DOF_vector(int:int+res_struc.poly_degree_rim-1)=Diaphragm.rim.drho_dz_vector; 
int=int+res_struc.poly_degree_rim; 

%virus angle at diphragm rim edge
DOF_name_and_range.Virus_rim_drho_dz_vector=int:int+res_struc.poly_degree_rim-1;
DOF_vector(int:int+res_struc.poly_degree_rim-1)=Virus.rim.drho_dz_vector;     
int=int+res_struc.poly_degree_rim;

%cell angle at diphragm rim edge
DOF_name_and_range.Cell_rim_drho_dz_vector=int:int+res_struc.poly_degree_rim-1;
DOF_vector(int:int+res_struc.poly_degree_rim-1)=Cell.rim.drho_dz_vector;
int=int+res_struc.poly_degree_rim; 



 % tilt exponential decay
if Minimazation.exponential_decay_tilt==1
    
    DOF_name_and_range.Diaphragm_Gaussian_decay_vector_up_k_vector=int;
    DOF_vector(int)=Diaphragm.Gaussian_decay_vector_up.k_vector(1);    
    int=int+1;
    
    DOF_name_and_range.Diaphragm_Gaussian_decay_vector_down_k_vector=int;
    DOF_vector(int)=Diaphragm.Gaussian_decay_vector_down.k_vector(1);  
    int=int+1;
    
    DOF_name_and_range.Virus_Gaussian_decay_vector_up_k_vector=int;
    DOF_vector(int)=Virus.Gaussian_decay_vector_up.k_vector(1);        
    int=int+1;
    
    DOF_name_and_range.Virus_Gaussian_decay_vector_down_k_vector=int;
    DOF_vector(int)=Virus.Gaussian_decay_vector_down.k_vector(1);      
    int=int+1;
    
    DOF_name_and_range.Cell_Gaussian_decay_vector_up_k_vector=int;
    DOF_vector(int)=Cell.Gaussian_decay_vector_up.k_vector(1);         
    int=int+1;
    
    DOF_name_and_range.Cell_Gaussian_decay_vector_down_k_vector=int;
    DOF_vector(int)=Cell.Gaussian_decay_vector_down.k_vector(1);       
    int=int+1;
end

% mid plane

DOF_name_and_range.Diaphragm_mid_plane_coeff_matrix=int:int+res_struc.poly_degree_z-5;
DOF_vector(int:int+res_struc.poly_degree_z-5)=Diaphragm.mid_plane_coeff_matrix(1,5:res_struc.poly_degree_z); 
int=int+res_struc.poly_degree_z-4;  

DOF_name_and_range.Virus_mid_plane_coeff_matrix=int:int+res_struc.poly_degree_z-5;
DOF_vector(int:int+res_struc.poly_degree_z-5)=Virus.mid_plane_coeff_matrix(1,5:res_struc.poly_degree_z);      
int=int+res_struc.poly_degree_z-4;

DOF_name_and_range.Cell_mid_plane_coeff_matrix=int:int+res_struc.poly_degree_z-5;
DOF_vector(int:int+res_struc.poly_degree_z-5)=Cell.mid_plane_coeff_matrix(1,5:res_struc.poly_degree_z); 
int=int+res_struc.poly_degree_z-4;

% lipid director rho

% the diaphragm has less DOF becuse axis symmetry!
DOF_name_and_range.Diaphragm_up_lipid_director_rho_matrix=int:int+res_struc.poly_degree_LD_rho-5;
DOF_vector(int:int+res_struc.poly_degree_LD_rho-5)=Diaphragm.up_lipid_director_rho_matrix(1,5:res_struc.poly_degree_LD_rho);
int=int+res_struc.poly_degree_LD_rho-4;

DOF_name_and_range.Diaphragm_down_lipid_director_rho_matrix=int:int+res_struc.poly_degree_LD_rho-5;
DOF_vector(int:int+res_struc.poly_degree_LD_rho-5)=Diaphragm.down_lipid_director_rho_matrix(1,5:res_struc.poly_degree_LD_rho); 
int=int+res_struc.poly_degree_LD_rho-4;

DOF_name_and_range.Virus_up_lipid_director_rho_matrix=int:int+res_struc.poly_degree_LD_rho-3;
DOF_vector(int:int+res_struc.poly_degree_LD_rho-3)=Virus.up_lipid_director_rho_matrix(1,3:res_struc.poly_degree_LD_rho);
int=int+res_struc.poly_degree_LD_rho-2;

DOF_name_and_range.Virus_down_lipid_director_rho_matrix=int:int+res_struc.poly_degree_LD_rho-3;
DOF_vector(int:int+res_struc.poly_degree_LD_rho-3)=Virus.down_lipid_director_rho_matrix(1,3:res_struc.poly_degree_LD_rho); 
int=int+res_struc.poly_degree_LD_rho-2;

DOF_name_and_range.Cell_up_lipid_director_rho_matrix=int:int+res_struc.poly_degree_LD_rho-3;
DOF_vector(int:int+res_struc.poly_degree_LD_rho-3)=Cell.up_lipid_director_rho_matrix(1,3:res_struc.poly_degree_LD_rho); 
int=int+res_struc.poly_degree_LD_rho-2;

DOF_name_and_range.Cell_down_lipid_director_rho_matrix=int:int+res_struc.poly_degree_LD_rho-3;
DOF_vector(int:int+res_struc.poly_degree_LD_rho-3)=Cell.down_lipid_director_rho_matrix(1,3:res_struc.poly_degree_LD_rho); 
int=int+res_struc.poly_degree_LD_rho-2; 


end

