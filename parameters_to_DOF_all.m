function [DOF_vector_all] = parameters_to_DOF_all(Diaphragm,Cell,Virus,res_struc)

    int=1; %position vector
    
    
    %we strat with the rims 
    %  rim position and general proprties
    DOF_vector_all(int)=Diaphragm.rim.x0; int=int+1; %ellipse size in x direction
    DOF_vector_all(int)=Diaphragm.rim.y0; int=int+1;%ellipse size in y direction
    DOF_vector_all(int)=Diaphragm.z0;     int=int+1;% set Diaphragm hight at rho=0 
    
    
    DOF_vector_all(int)=Virus.hight_above_plane;    % hight of virus surface above cell surface
    int=int+1;
    
    DOF_vector_all(int)=Virus.Rv; % fusion site size in torus         
    int=int+1;   

    DOF_vector_all(int)=Cell.R_rim; 
    DOF_stages.general=int;       
    int=int+1;
    
    % max rim int
   
    %% rims
    % diaphragm diaphragm rim hight
    DOF_vector_all(int:int+res_struc.poly_degree_rim-1)=Diaphragm.rim.hl_vector;      
    int=int+res_struc.poly_degree_rim; 
    % diaphragm diaphragm rim hight dirivitive
    DOF_vector_all(int:int+res_struc.poly_degree_rim-1)=Diaphragm.rim.drho_dz_vector; 
    int=int+res_struc.poly_degree_rim; 
    
    %virus angle at diphragm rim edge
    DOF_vector_all(int:int+res_struc.poly_degree_rim-1)=Virus.rim.drho_dz_vector;     
    int=int+res_struc.poly_degree_rim;
    
    %cell angle at diphragm rim edge
    DOF_vector_all(int:int+res_struc.poly_degree_rim-1)=Cell.rim.drho_dz_vector;
    int=int+res_struc.poly_degree_rim;  
    
    %% mid plane
    add_DOF_int=(res_struc.poly_degree_z-4)*res_struc.phi_res;
    
    DOF_vector_all(int:int+add_DOF_int-1)=array_to_vector(Diaphragm.mid_plane_coeff_matrix(:,5:res_struc.poly_degree_z),res_struc.poly_degree_z-4,res_struc.phi_res);
    int=int+add_DOF_int;  
 
    DOF_vector_all(int:int+add_DOF_int-1)=array_to_vector(Virus.mid_plane_coeff_matrix(:,5:res_struc.poly_degree_z),res_struc.poly_degree_z-4,res_struc.phi_res);
    int=int+add_DOF_int;  
    
    DOF_vector_all(int:int+add_DOF_int-1)=array_to_vector(Cell.mid_plane_coeff_matrix(:,5:res_struc.poly_degree_z),res_struc.poly_degree_z-4,res_struc.phi_res);
    int=int+add_DOF_int;  
    
    %% lipid director rho 
    add_DOF_int=(res_struc.poly_degree_LD_rho-4)*res_struc.phi_res;
    
    DOF_vector_all(int:int+add_DOF_int-1)=array_to_vector(Diaphragm.up_lipid_director_rho_matrix(:,5:res_struc.poly_degree_LD_rho),res_struc.poly_degree_LD_rho-4,res_struc.phi_res);
    int=int+add_DOF_int;  
    
    DOF_vector_all(int:int+add_DOF_int-1)=array_to_vector(Diaphragm.down_lipid_director_rho_matrix(:,5:res_struc.poly_degree_LD_rho),res_struc.poly_degree_LD_rho-4,res_struc.phi_res);
    int=int+add_DOF_int;  

    add_DOF_int=(res_struc.poly_degree_LD_rho-2)*res_struc.phi_res;
    DOF_vector_all(int:int+add_DOF_int-1)=array_to_vector(Virus.up_lipid_director_rho_matrix(:,3:res_struc.poly_degree_LD_rho),res_struc.poly_degree_LD_rho-2,res_struc.phi_res);
    int=int+add_DOF_int;  
    
    DOF_vector_all(int:int+add_DOF_int-1)=array_to_vector(Virus.down_lipid_director_rho_matrix(:,3:res_struc.poly_degree_LD_rho),res_struc.poly_degree_LD_rho-2,res_struc.phi_res);
    int=int+add_DOF_int; 
    
    DOF_vector_all(int:int+add_DOF_int-1)=array_to_vector(Cell.up_lipid_director_rho_matrix(:,3:res_struc.poly_degree_LD_rho),res_struc.poly_degree_LD_rho-2,res_struc.phi_res);
    int=int+add_DOF_int;  
    
    DOF_vector_all(int:int+add_DOF_int-1)=array_to_vector(Cell.down_lipid_director_rho_matrix(:,3:res_struc.poly_degree_LD_rho),res_struc.poly_degree_LD_rho-2,res_struc.phi_res);
    int=int+add_DOF_int; 
    
    
    %% lipid director phi
    add_DOF_int=(res_struc.poly_degree_LD_phi-4)*res_struc.phi_res;
    
    DOF_vector_all(int:int+add_DOF_int-1)=array_to_vector(Diaphragm.up_lipid_director_phi_matrix(:,5:res_struc.poly_degree_LD_phi),res_struc.poly_degree_LD_phi-4,res_struc.phi_res);
    int=int+add_DOF_int;  

    DOF_vector_all(int:int+add_DOF_int-1)=array_to_vector(Diaphragm.down_lipid_director_phi_matrix(:,5:res_struc.poly_degree_LD_phi),res_struc.poly_degree_LD_phi-4,res_struc.phi_res);
    int=int+add_DOF_int;      
    
    DOF_vector_all(int:int+add_DOF_int-1)=array_to_vector(Virus.up_lipid_director_phi_matrix(:,5:res_struc.poly_degree_LD_phi),res_struc.poly_degree_LD_phi-4,res_struc.phi_res);
    int=int+add_DOF_int;  

    DOF_vector_all(int:int+add_DOF_int-1)=array_to_vector(Virus.down_lipid_director_phi_matrix(:,5:res_struc.poly_degree_LD_phi),res_struc.poly_degree_LD_phi-4,res_struc.phi_res);
    int=int+add_DOF_int;
    
    DOF_vector_all(int:int+add_DOF_int-1)=array_to_vector(Cell.up_lipid_director_phi_matrix(:,5:res_struc.poly_degree_LD_phi),res_struc.poly_degree_LD_phi-4,res_struc.phi_res);
    int=int+add_DOF_int;  

    DOF_vector_all(int:int+add_DOF_int-1)=array_to_vector(Cell.down_lipid_director_phi_matrix(:,5:res_struc.poly_degree_LD_phi),res_struc.poly_degree_LD_phi-4,res_struc.phi_res);
   
    
    
end

