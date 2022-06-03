function [Diaphragm,Cell,Virus] = DOF_to_parameters(Diaphragm,Cell,Virus,DOF_vector_all,res_struc)
int=1; %position vector
    
    
    %we strat with the rims 
    %  rim position and general proprties
    Diaphragm.rim.x0=DOF_vector_all(int); int=int+1; %ellipse size in x direction
    Diaphragm.rim.y0=DOF_vector_all(int); int=int+1;%ellipse size in y direction
    Diaphragm.z0=DOF_vector_all(int);     int=int+1;% set Diaphragm hight at rho=0 
    
    
    Virus.hight_above_plane=DOF_vector_all(int);    % hight of virus surface above cell surface
    int=int+1;
    
    Virus.Rv=DOF_vector_all(int); % fusion site size in torus         
    int=int+1;   

    Cell.R_rim=DOF_vector_all(int); 
    int=int+1;
    
    % max rim int
   
    %% rims
    % diaphragm diaphragm rim hight
    Diaphragm.rim.hl_vector = vector_to_array(DOF_vector_all(int:int+res_struc.poly_degree_rim-1),res_struc.poly_degree_rim,1);
    int=int+res_struc.poly_degree_rim; 
    % diaphragm diaphragm rim hight dirivitive
    Diaphragm.rim.drho_dz_vector = vector_to_array(DOF_vector_all(int:int+res_struc.poly_degree_rim-1),res_struc.poly_degree_rim,1);
    int=int+res_struc.poly_degree_rim; 
    
    %virus angle at diphragm rim edge
    Virus.rim.drho_dz_vector = vector_to_array(DOF_vector_all(int:int+res_struc.poly_degree_rim-1),res_struc.poly_degree_rim,1);
    int=int+res_struc.poly_degree_rim;
    
    %cell angle at diphragm rim edge
    Cell.rim.drho_dz_vector = vector_to_array(DOF_vector_all(int:int+res_struc.poly_degree_rim-1),res_struc.poly_degree_rim,1);
    int=int+res_struc.poly_degree_rim;  
    
    %% mid plane
    add_DOF_int=(res_struc.poly_degree_z-4)*res_struc.phi_res;
    Diaphragm.mid_plane_coeff_matrix(:,5:res_struc.poly_degree_z) =...
        vector_to_array(DOF_vector_all(int:int+add_DOF_int-1),res_struc.poly_degree_z-4,res_struc.phi_res);
    int=int+add_DOF_int;  
 
    Virus.mid_plane_coeff_matrix(:,5:res_struc.poly_degree_z) =...
        vector_to_array(DOF_vector_all(int:int+add_DOF_int-1),res_struc.poly_degree_z-4,res_struc.phi_res);
    int=int+add_DOF_int;  

    Cell.mid_plane_coeff_matrix(:,5:res_struc.poly_degree_z) =...
        vector_to_array(DOF_vector_all(int:int+add_DOF_int-1),res_struc.poly_degree_z-4,res_struc.phi_res);
    int=int+add_DOF_int;  
    
    %% lipid director rho 
    add_DOF_int=(res_struc.poly_degree_LD_rho-4)*res_struc.phi_res;
    
    Diaphragm.up_lipid_director_rho_matrix(:,5:res_struc.poly_degree_LD_rho) =...
        vector_to_array(DOF_vector_all(int:int+add_DOF_int-1),res_struc.poly_degree_LD_rho-4,res_struc.phi_res);  
    int=int+add_DOF_int;  
    
    Diaphragm.down_lipid_director_rho_matrix(:,5:res_struc.poly_degree_LD_rho) =...
        vector_to_array(DOF_vector_all(int:int+add_DOF_int-1),res_struc.poly_degree_LD_rho-4,res_struc.phi_res);  
    int=int+add_DOF_int; 

    add_DOF_int=(res_struc.poly_degree_LD_rho-2)*res_struc.phi_res;

    Virus.up_lipid_director_rho_matrix(:,3:res_struc.poly_degree_LD_rho) =...
        vector_to_array(DOF_vector_all(int:int+add_DOF_int-1),res_struc.poly_degree_LD_rho-2,res_struc.phi_res);  
    int=int+add_DOF_int;  
    
    Virus.down_lipid_director_rho_matrix(:,3:res_struc.poly_degree_LD_rho) =...
        vector_to_array(DOF_vector_all(int:int+add_DOF_int-1),res_struc.poly_degree_LD_rho-2,res_struc.phi_res);  
    int=int+add_DOF_int;
    
    Cell.up_lipid_director_rho_matrix(:,3:res_struc.poly_degree_LD_rho) =...
        vector_to_array(DOF_vector_all(int:int+add_DOF_int-1),res_struc.poly_degree_LD_rho-2,res_struc.phi_res);  
    int=int+add_DOF_int;  
    
    Cell.down_lipid_director_rho_matrix(:,3:res_struc.poly_degree_LD_rho) =...
        vector_to_array(DOF_vector_all(int:int+add_DOF_int-1),res_struc.poly_degree_LD_rho-2,res_struc.phi_res);  
    int=int+add_DOF_int; 
    
    
    %% lipid director phi
    add_DOF_int=(res_struc.poly_degree_LD_phi-4)*res_struc.phi_res;

    Diaphragm.up_lipid_director_phi_matrix(:,5:res_struc.poly_degree_LD_phi) =...
        vector_to_array(DOF_vector_all(int:int+add_DOF_int-1),res_struc.poly_degree_LD_phi-4,res_struc.phi_res);  
    int=int+add_DOF_int; 

    Diaphragm.down_lipid_director_phi_matrix(:,5:res_struc.poly_degree_LD_phi) =...
        vector_to_array(DOF_vector_all(int:int+add_DOF_int-1),res_struc.poly_degree_LD_phi-4,res_struc.phi_res);  
    int=int+add_DOF_int;      

    Virus.up_lipid_director_phi_matrix(:,5:res_struc.poly_degree_LD_phi) =...
        vector_to_array(DOF_vector_all(int:int+add_DOF_int-1),res_struc.poly_degree_LD_phi-4,res_struc.phi_res);  
    int=int+add_DOF_int; 

    Virus.down_lipid_director_phi_matrix(:,5:res_struc.poly_degree_LD_phi) =...
        vector_to_array(DOF_vector_all(int:int+add_DOF_int-1),res_struc.poly_degree_LD_phi-4,res_struc.phi_res);  
    int=int+add_DOF_int;    

    Cell.up_lipid_director_phi_matrix(:,5:res_struc.poly_degree_LD_phi) =...
        vector_to_array(DOF_vector_all(int:int+add_DOF_int-1),res_struc.poly_degree_LD_phi-4,res_struc.phi_res);  
    int=int+add_DOF_int; 

    Cell.down_lipid_director_phi_matrix(:,5:res_struc.poly_degree_LD_phi) =...
        vector_to_array(DOF_vector_all(int:int+add_DOF_int-1),res_struc.poly_degree_LD_phi-4,res_struc.phi_res);  

end

