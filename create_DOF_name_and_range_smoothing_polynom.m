function [DOF_name_and_range] = create_DOF_name_and_range_smoothing_polynom(Minimazation,res_struc)

% just initialize to nothing
if res_struc.poly_degree_z-4>0
    DOF_name_and_range.Diaphragm_mid_plane_smooth_polynom_matrix=[];
    DOF_name_and_range.Cell_mid_plane_smooth_polynom_matrix=[];
    DOF_name_and_range.Virus_mid_plane_smooth_polynom_matrix=[];
    DOF_name_and_range.MATRIX_Diaphragm_mid_plane_smooth_polynom_matrix=zeros(res_struc.poly_degree_z-4,res_struc.polar_angle_polynom);
    DOF_name_and_range.MATRIX_Cell_mid_plane_smooth_polynom_matrix=zeros(res_struc.poly_degree_z-4,res_struc.polar_angle_polynom);
    DOF_name_and_range.MATRIX_Virus_mid_plane_smooth_polynom_matrix=zeros(res_struc.poly_degree_z-4,res_struc.polar_angle_polynom);
end

if res_struc.poly_degree_LD_rho-4>0
    DOF_name_and_range.Diaphragm_up_lipid_director_rho_smooth_polynom_matrix=[];
    DOF_name_and_range.Diaphragm_down_lipid_director_rho_smooth_polynom_matrix=[];
    DOF_name_and_range.Cell_up_lipid_director_rho_smooth_polynom_matrix=[];
    DOF_name_and_range.Cell_down_lipid_director_rho_smooth_polynom_matrix=[];
    DOF_name_and_range.Virus_up_lipid_director_rho_smooth_polynom_matrix=[];
    DOF_name_and_range.Virus_down_lipid_director_rho_smooth_polynom_matrix=[];

    DOF_name_and_range.MATRIX_Diaphragm_up_lipid_director_rho_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_rho-4,res_struc.polar_angle_polynom);
    DOF_name_and_range.MATRIX_Diaphragm_down_lipid_director_rho_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_rho-4,res_struc.polar_angle_polynom);
    DOF_name_and_range.MATRIX_Cell_up_lipid_director_rho_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_rho-4,res_struc.polar_angle_polynom);
    DOF_name_and_range.MATRIX_Cell_down_lipid_director_rho_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_rho-4,res_struc.polar_angle_polynom);
    DOF_name_and_range.MATRIX_Virus_up_lipid_director_rho_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_rho-4,res_struc.polar_angle_polynom);
    DOF_name_and_range.MATRIX_Virus_down_lipid_director_rho_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_rho-4,res_struc.polar_angle_polynom);
end

if res_struc.poly_degree_LD_phi-4>0
    DOF_name_and_range.Diaphragm_up_lipid_director_phi_smooth_polynom_matrix=[];
    DOF_name_and_range.Diaphragm_down_lipid_director_phi_smooth_polynom_matrix=[];
    DOF_name_and_range.Cell_up_lipid_director_phi_smooth_polynom_matrix=[];
    DOF_name_and_range.Cell_down_lipid_director_phi_smooth_polynom_matrix=[];
    DOF_name_and_range.Virus_up_lipid_director_phi_smooth_polynom_matrix=[];
    DOF_name_and_range.Virus_down_lipid_director_phi_smooth_polynom_matrix=[];
    DOF_name_and_range.MATRIX_Diaphragm_up_lipid_director_phi_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_phi-4,res_struc.polar_angle_polynom);
    DOF_name_and_range.MATRIX_Diaphragm_down_lipid_director_phi_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_phi-4,res_struc.polar_angle_polynom);
    DOF_name_and_range.MATRIX_Cell_up_lipid_director_phi_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_phi-4,res_struc.polar_angle_polynom);
    DOF_name_and_range.MATRIX_Cell_down_lipid_director_phi_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_phi-4,res_struc.polar_angle_polynom);
    DOF_name_and_range.MATRIX_Virus_up_lipid_director_phi_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_phi-4,res_struc.polar_angle_polynom);
    DOF_name_and_range.MATRIX_Virus_down_lipid_director_phi_smooth_polynom_matrix=zeros(res_struc.poly_degree_LD_phi-4,res_struc.polar_angle_polynom);
end

int=1; %position vector

%  rim position and general proprties
%ellipse size in x direction
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


%% rims
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

%% 

if Minimazation.exponential_decay_tilt==1
     DOF_name_and_range.Diaphragm_Gaussian_decay_vector_up_k_vector_smooth_poly=int:int+res_struc.polar_angle_polynom-1;
     int=int+res_struc.polar_angle_polynom;

     DOF_name_and_range.Diaphragm_Gaussian_decay_vector_down_k_vector_smooth_poly=int:int+res_struc.polar_angle_polynom-1;
     int=int+res_struc.polar_angle_polynom;

     DOF_name_and_range.Virus_Gaussian_decay_vector_up_k_vector_smooth_poly=int:int+res_struc.polar_angle_polynom-1;
     int=int+res_struc.polar_angle_polynom;

     DOF_name_and_range.Virus_Gaussian_decay_vector_down_k_vector_smooth_poly=int:int+res_struc.polar_angle_polynom-1;         
     int=int+res_struc.polar_angle_polynom;

     DOF_name_and_range.Cell_Gaussian_decay_vector_up_k_vector_smooth_poly=int:int+res_struc.polar_angle_polynom-1;                 
     int=int+res_struc.polar_angle_polynom;

     DOF_name_and_range.Cell_Gaussian_decay_vector_down_k_vector_smooth_poly=int:int+res_struc.polar_angle_polynom-1;                 
     int=int+res_struc.polar_angle_polynom;
end


%% do each order sepratly order term of the smoothing polynom
int_poly_order=1;
while int_poly_order<=res_struc.polar_angle_polynom
    %mid plane
    if res_struc.poly_degree_z-4>0

        DOF_name_and_range.Diaphragm_mid_plane_smooth_polynom_matrix=[DOF_name_and_range.Diaphragm_mid_plane_smooth_polynom_matrix , int:int+res_struc.poly_degree_z-5];
        DOF_name_and_range.MATRIX_Diaphragm_mid_plane_smooth_polynom_matrix(:,int_poly_order)= int:int+res_struc.poly_degree_z-5;
        int=int+res_struc.poly_degree_z-4;

        DOF_name_and_range.Cell_mid_plane_smooth_polynom_matrix=[DOF_name_and_range.Cell_mid_plane_smooth_polynom_matrix , int:int+res_struc.poly_degree_z-5];
        DOF_name_and_range.MATRIX_Cell_mid_plane_smooth_polynom_matrix(:,int_poly_order)= int:int+res_struc.poly_degree_z-5;
        int=int+res_struc.poly_degree_z-4;

        DOF_name_and_range.Virus_mid_plane_smooth_polynom_matrix=[DOF_name_and_range.Virus_mid_plane_smooth_polynom_matrix , int:int+res_struc.poly_degree_z-5];
        DOF_name_and_range.MATRIX_Virus_mid_plane_smooth_polynom_matrix(:,int_poly_order)= int:int+res_struc.poly_degree_z-5;
        int=int+res_struc.poly_degree_z-4;
    end     
    %rho lipid director

    if res_struc.poly_degree_LD_rho-4>0

        DOF_name_and_range.Diaphragm_up_lipid_director_rho_smooth_polynom_matrix=[DOF_name_and_range.Diaphragm_up_lipid_director_rho_smooth_polynom_matrix,int:int+res_struc.poly_degree_LD_rho-5];            
        DOF_name_and_range.MATRIX_Diaphragm_up_lipid_director_rho_smooth_polynom_matrix(:,int_poly_order)= int:int+res_struc.poly_degree_LD_rho-5;
        int=int+res_struc.poly_degree_LD_rho-4;

        DOF_name_and_range.Diaphragm_down_lipid_director_rho_smooth_polynom_matrix=[DOF_name_and_range.Diaphragm_down_lipid_director_rho_smooth_polynom_matrix,int:int+res_struc.poly_degree_LD_rho-5];            
        DOF_name_and_range.MATRIX_Diaphragm_down_lipid_director_rho_smooth_polynom_matrix(:,int_poly_order)= int:int+res_struc.poly_degree_LD_rho-5;
        int=int+res_struc.poly_degree_LD_rho-4;

        DOF_name_and_range.Cell_up_lipid_director_rho_smooth_polynom_matrix=[DOF_name_and_range.Cell_up_lipid_director_rho_smooth_polynom_matrix,int:int+res_struc.poly_degree_LD_rho-5];            
        DOF_name_and_range.MATRIX_Cell_up_lipid_director_rho_smooth_polynom_matrix(:,int_poly_order)= int:int+res_struc.poly_degree_LD_rho-5;
        int=int+res_struc.poly_degree_LD_rho-4;

        DOF_name_and_range.Cell_down_lipid_director_rho_smooth_polynom_matrix=[DOF_name_and_range.Cell_down_lipid_director_rho_smooth_polynom_matrix,int:int+res_struc.poly_degree_LD_rho-5];                        
        DOF_name_and_range.MATRIX_Cell_down_lipid_director_rho_smooth_polynom_matrix(:,int_poly_order)= int:int+res_struc.poly_degree_LD_rho-5;
        int=int+res_struc.poly_degree_LD_rho-4;

        DOF_name_and_range.Virus_up_lipid_director_rho_smooth_polynom_matrix=[DOF_name_and_range.Virus_up_lipid_director_rho_smooth_polynom_matrix,int:int+res_struc.poly_degree_LD_rho-5];                        
        DOF_name_and_range.MATRIX_Virus_up_lipid_director_rho_smooth_polynom_matrix(:,int_poly_order)= int:int+res_struc.poly_degree_LD_rho-5;
        int=int+res_struc.poly_degree_LD_rho-4;

        DOF_name_and_range.Virus_down_lipid_director_rho_smooth_polynom_matrix=[DOF_name_and_range.Virus_down_lipid_director_rho_smooth_polynom_matrix,int:int+res_struc.poly_degree_LD_rho-5];                        
        DOF_name_and_range.MATRIX_Virus_down_lipid_director_rho_smooth_polynom_matrix(:,int_poly_order)= int:int+res_struc.poly_degree_LD_rho-5;
        int=int+res_struc.poly_degree_LD_rho-4; 
    end


    
    %phi lipid director
    if res_struc.poly_degree_LD_phi-4>0

        DOF_name_and_range.Diaphragm_up_lipid_director_phi_smooth_polynom_matrix=[DOF_name_and_range.Diaphragm_up_lipid_director_phi_smooth_polynom_matrix,int:int+res_struc.poly_degree_LD_phi-5];
        DOF_name_and_range.MATRIX_Diaphragm_up_lipid_director_phi_smooth_polynom_matrix(:,int_poly_order)= int:int+res_struc.poly_degree_LD_phi-5;
        int=int+res_struc.poly_degree_LD_phi-4;

        DOF_name_and_range.Diaphragm_down_lipid_director_phi_smooth_polynom_matrix=[DOF_name_and_range.Diaphragm_down_lipid_director_phi_smooth_polynom_matrix,int:int+res_struc.poly_degree_LD_phi-5];
        DOF_name_and_range.MATRIX_Diaphragm_down_lipid_director_phi_smooth_polynom_matrix(:,int_poly_order)= int:int+res_struc.poly_degree_LD_phi-5;
        int=int+res_struc.poly_degree_LD_phi-4;
        
        DOF_name_and_range.Cell_up_lipid_director_phi_smooth_polynom_matrix=[DOF_name_and_range.Cell_up_lipid_director_phi_smooth_polynom_matrix,int:int+res_struc.poly_degree_LD_phi-5];
        DOF_name_and_range.MATRIX_Cell_up_lipid_director_phi_smooth_polynom_matrix(:,int_poly_order)= int:int+res_struc.poly_degree_LD_phi-5;
        int=int+res_struc.poly_degree_LD_phi-4;

        DOF_name_and_range.Cell_down_lipid_director_phi_smooth_polynom_matrix=[DOF_name_and_range.Cell_down_lipid_director_phi_smooth_polynom_matrix,int:int+res_struc.poly_degree_LD_phi-5];
        DOF_name_and_range.MATRIX_Cell_down_lipid_director_phi_smooth_polynom_matrix(:,int_poly_order)= int:int+res_struc.poly_degree_LD_phi-5;
        int=int+res_struc.poly_degree_LD_phi-4; 

        DOF_name_and_range.Virus_up_lipid_director_phi_smooth_polynom_matrix=[DOF_name_and_range.Virus_up_lipid_director_phi_smooth_polynom_matrix,int:int+res_struc.poly_degree_LD_phi-5];
        DOF_name_and_range.MATRIX_Virus_up_lipid_director_phi_smooth_polynom_matrix(:,int_poly_order)= int:int+res_struc.poly_degree_LD_phi-5;
        int=int+res_struc.poly_degree_LD_phi-4;

        DOF_name_and_range.Virus_down_lipid_director_phi_smooth_polynom_matrix=[DOF_name_and_range.Virus_down_lipid_director_phi_smooth_polynom_matrix,int:int+res_struc.poly_degree_LD_phi-5];        
        DOF_name_and_range.MATRIX_Virus_down_lipid_director_phi_smooth_polynom_matrix(:,int_poly_order)= int:int+res_struc.poly_degree_LD_phi-5;
        int=int+res_struc.poly_degree_LD_phi-4;
    end


    int_poly_order=int_poly_order+1;
end


end

