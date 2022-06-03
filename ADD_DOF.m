function [new_DOF_vector,new_Minimazation] = ADD_DOF(old_res_struc,new_res_struc,old_DOF_vector,Minimazation)
% this function update the DOF by adding new arrays to the old_DOF_vector
% the values that were already minimized are kept

    new_DOF_vector=old_DOF_vector;
    new_Minimazation=Minimazation;
    %first 6 cannot change
    new_DOF_vector(1:6)=old_DOF_vector(1:6);
    new_Minimazation.steps_size_vector(1:6)=Minimazation.steps_size_vector(1:6);
    % max rim old_int
    int=7;
    %% rims
    % diaphragm diaphragm rim hight
    rim_diff=new_res_struc.poly_degree_rim-old_res_struc.poly_degree_rim;
    
    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,rim_diff,int,0);
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (Minimazation.steps_size_vector,rim_diff,int,new_Minimazation.steps_size_vector(int-1));
    
    int=int+new_res_struc.poly_degree_rim; 
    % diaphragm diaphragm rim hight dirivitive
    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,rim_diff,int,0);
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,rim_diff,int,new_Minimazation.steps_size_vector(int-1));
    int=int+new_res_struc.poly_degree_rim; 
    
    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,rim_diff,int,0);
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,rim_diff,int,new_Minimazation.steps_size_vector(int-1));
    int=int+new_res_struc.poly_degree_rim;
    
    %cell angle at diphragm rim edge
    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,rim_diff,int,0);
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,rim_diff,int,new_Minimazation.steps_size_vector(int-1));
    int=int+new_res_struc.poly_degree_rim;  

    %% mid plane
    z_diff=new_res_struc.poly_degree_z-old_res_struc.poly_degree_z;

    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,z_diff,int,0);
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,z_diff,int,new_Minimazation.steps_size_vector(int-1));
    int=int+new_res_struc.poly_degree_z-4;  
    
    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,z_diff,int,0);
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,z_diff,int,new_Minimazation.steps_size_vector(int-1));
    int=int+new_res_struc.poly_degree_z-4;
    
    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,z_diff,int,0); 
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,z_diff,int,new_Minimazation.steps_size_vector(int-1));    
    int=int+new_res_struc.poly_degree_z-4;

    %% lipid director rho
    LD_rho_diff=new_res_struc.poly_degree_LD_rho-old_res_struc.poly_degree_LD_rho;

    
    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,LD_rho_diff,int,0);
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,LD_rho_diff,int,new_Minimazation.steps_size_vector(int-1));    
    int=int+new_res_struc.poly_degree_LD_rho-4;
    
    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,LD_rho_diff,int,0); 
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,LD_rho_diff,int,new_Minimazation.steps_size_vector(int-1));  
    int=int+new_res_struc.poly_degree_LD_rho-4;
    
    [new_DOF_vector] =Add_value_to_vector(new_DOF_vector,LD_rho_diff,int,0); 
     [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,LD_rho_diff,int,new_Minimazation.steps_size_vector(int-1));  
    int=int+new_res_struc.poly_degree_LD_rho-2;
    
    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,LD_rho_diff,int,0);
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,LD_rho_diff,int,new_Minimazation.steps_size_vector(int-1));      
    int=int+new_res_struc.poly_degree_LD_rho-2;
    
    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,LD_rho_diff,int,0);
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,LD_rho_diff,int,new_Minimazation.steps_size_vector(int-1));      
    int=int+new_res_struc.poly_degree_LD_rho-2; 
    
    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,LD_rho_diff,int,0);
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,LD_rho_diff,int,new_Minimazation.steps_size_vector(int-1));      
    int=int+new_res_struc.poly_degree_LD_rho-2;
    
    
    %% lipid director phi
    LD_phi_diff=new_res_struc.poly_degree_LD_phi-old_res_struc.poly_degree_LD_phi;

    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,LD_phi_diff,int,0);
    
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,LD_phi_diff,int,new_Minimazation.steps_size_vector(int-1));      
    int=int+new_res_struc.poly_degree_LD_phi-4;
    
    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,LD_phi_diff,int,0);
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,LD_phi_diff,int,new_Minimazation.steps_size_vector(int-1));
    int=int+new_res_struc.poly_degree_LD_phi-4;
    
    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,LD_phi_diff,int,0);
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,LD_phi_diff,int,new_Minimazation.steps_size_vector(int-1));
    int=int+new_res_struc.poly_degree_LD_phi-4;
    
    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,LD_phi_diff,int,0);
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,LD_phi_diff,int,new_Minimazation.steps_size_vector(int-1));    
    int=int+new_res_struc.poly_degree_LD_phi-4;
    
    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,LD_phi_diff,int,0);
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,LD_phi_diff,int,new_Minimazation.steps_size_vector(int-1));       
    int=int+new_res_struc.poly_degree_LD_phi-4;
    
    [new_DOF_vector] = Add_value_to_vector(new_DOF_vector,LD_phi_diff,int,0);  
    [new_Minimazation.steps_size_vector] = Add_value_to_vector...
        (new_Minimazation.steps_size_vector,LD_phi_diff,int,new_Minimazation.steps_size_vector(int-1));
    

end

