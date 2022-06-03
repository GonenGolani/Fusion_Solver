function [DOF_vector,DOF_name_and_range] = Shell_parameters_to_DOF_vector(old_DOF_vector,old_DOF_name_and_range,Shell,res_struc)
    %trnaslate the degrees of freedom (DOF) vector to the values
    %minimaztion parametrs in of the structure 
    
    %old_DOF_vector - the DOF vector of the HD without the shell
    
    DOF_vector=old_DOF_vector;
    DOF_name_and_range=old_DOF_name_and_range;
    
    
    int=length(DOF_vector)+1; %position vector
    
    
   
    % the diaphragm shell
    DOF_vector(int)=Shell.Diaphragm.z_i;                 
    DOF_name_and_range.Shell_Diaphragm_z_i=int;
    int=int+1;

    DOF_vector(int)=Shell.Diaphragm.z_f;                 
    DOF_name_and_range.Shell_Diaphragm_z_f=int;
    int=int+1;

    DOF_vector(int)=Shell.Diaphragm.dzdrho_edge;
    DOF_name_and_range.Shell_Diaphragm_dzdrho_edge=int;
    int=int+1;

    DOF_vector(int)=Shell.Diaphragm.length_from_edge;
    DOF_name_and_range.Shell_Diaphragm_length_from_edge=int;
    int=int+1;
    
    % the Cell shell
    DOF_vector(int)=Shell.Cell.z_i;
    DOF_name_and_range.Shell_Cell_z_i=int;
    int=int+1;

    DOF_vector(int)=Shell.Cell.z_f;
    DOF_name_and_range.Shell_Cell_z_f=int;
    int=int+1;

    DOF_vector(int)=Shell.Cell.dzdrho_edge;
    DOF_name_and_range.Shell_Cell_dzdrho_edge=int;
    int=int+1;

    DOF_vector(int)=Shell.Cell.length_from_edge;
    DOF_name_and_range.Shell_Cell_length_from_edge=int;
    int=int+1;
    
    % max rim int
   
    %% mid plane to mid-plane coefficiant matrix 

    DOF_vector(int:int+res_struc.poly_degree_shell-5)=Shell.Diaphragm.MP_to_MP_distance_matrix(1,5:res_struc.poly_degree_shell); 
    DOF_name_and_range.Shell_Diaphragm_MP_to_MP_distance_matrix=(int:int+res_struc.poly_degree_shell-5);    
    int=int+res_struc.poly_degree_shell-4;  
    
    DOF_vector(int:int+res_struc.poly_degree_shell-5)=Shell.Cell.MP_to_MP_distance_matrix(1,5:res_struc.poly_degree_shell); 
    DOF_name_and_range.Shell_Cell_MP_to_MP_distance_matrix=(int:int+res_struc.poly_degree_shell-5);    


    
end

