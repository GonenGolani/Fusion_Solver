function [fix_DOF,DOF_vector] = Force_DOF_constreints(DOF_vector,Diaphragm,Cell,Virus,Minimazation,res_struc)
% this function fixes some of the DOF given the constrints specified by the user
% the output is the DOF_vector after constreints fixing and fix_DOF: a
% vector conteining 1 in each DOF fixed and zero where it is not fixed

DOF_name_and_range=Minimazation.DOF_name_and_range;

fix_DOF=DOF_vector*0;


if Minimazation.axial_symmetry==1   % if there is axial symmetry x0=y0
    DOF_vector(DOF_name_and_range.Diaphragm_rim_x0)=Diaphragm.rim.y0;
    fix_DOF(DOF_name_and_range.Diaphragm_rim_x0)=1;
    % rim is treated becuse res_struc.poly_degree_rim=1
    if res_struc.poly_degree_rim~=1
        fprintf('ERROR!!!! res_struc.poly_degree_rim~=1 in axial symmetry');        
    end
    % all the higher pllynom smotthing minimaztion in non axial symmetric
    % cases does not enter the DOF_vector  
    
    
    
% if there is up_down -> fix the dipharagm in the middle
if Minimazation.up_down_symmetry==1
    
 
    DOF_vector(DOF_name_and_range.Diaphragm_rim_z0)=Virus.hight_above_plane/2;
    fix_DOF(DOF_name_and_range.Diaphragm_rim_z0)=1; 
   

    if Minimazation.HD_rim_fixed==0
        DOF_vector(DOF_name_and_range.Virus_Rv)=Cell.R_rim;      
        fix_DOF(DOF_name_and_range.Virus_Rv)=1;
    end
    
    %new part
    % make all polynoms symmetric around the diaphragm mid plane
    
    DOF_vector(DOF_name_and_range.Diaphragm_Gaussian_decay_vector_down_k_vector)=Diaphragm.Gaussian_decay_vector_up.k_vector(1);
    fix_DOF(DOF_name_and_range.Diaphragm_Gaussian_decay_vector_down_k_vector)=1;
    
    DOF_vector(DOF_name_and_range.Virus_Gaussian_decay_vector_up_k_vector)=Cell.Gaussian_decay_vector_down.k_vector(1);
    fix_DOF(DOF_name_and_range.Virus_Gaussian_decay_vector_up_k_vector)=1;
   
    DOF_vector(DOF_name_and_range.Virus_Gaussian_decay_vector_down_k_vector)=Cell.Gaussian_decay_vector_up.k_vector(1);
    fix_DOF(DOF_name_and_range.Virus_Gaussian_decay_vector_down_k_vector)=1;
    

    
    % cell and virus mid plane
    %z0_vec=zeros(1,res_struc.poly_degree_z)+Diaphragm.z0;
    DOF_vector(DOF_name_and_range.Virus_mid_plane_coeff_matrix)=-Cell.mid_plane_coeff_matrix(1,5:res_struc.poly_degree_z);
    fix_DOF(DOF_name_and_range.Virus_mid_plane_coeff_matrix)=1;
    
    % lipid director
    DOF_vector(DOF_name_and_range.Diaphragm_down_lipid_director_rho_matrix)=-Diaphragm.up_lipid_director_rho_matrix(1,5:res_struc.poly_degree_LD_rho);
    fix_DOF(DOF_name_and_range.Diaphragm_down_lipid_director_rho_matrix)=1;
    
    DOF_vector(DOF_name_and_range.Virus_up_lipid_director_rho_matrix)=-Cell.down_lipid_director_rho_matrix(1,3:res_struc.poly_degree_LD_rho);
    fix_DOF(DOF_name_and_range.Virus_up_lipid_director_rho_matrix)=1;

    DOF_vector(DOF_name_and_range.Virus_down_lipid_director_rho_matrix)=-Cell.up_lipid_director_rho_matrix(1,3:res_struc.poly_degree_LD_rho);
    fix_DOF(DOF_name_and_range.Virus_down_lipid_director_rho_matrix)=1;
   
    
    %set the junction tangent angle of the (Virus.rim.drho_dz_vector=Cell.rim.drho_dz_vector)
    DOF_vector(DOF_name_and_range.Virus_rim_drho_dz_vector)=-Cell.rim.drho_dz_vector;
    fix_DOF(DOF_name_and_range.Virus_rim_drho_dz_vector)=1;
    
end


    
end

%fixed distance between cell and virus
if Minimazation.fix_inter_membrane_distance~=0
    DOF_vector(DOF_name_and_range.Virus_hight_above_plane)=Minimazation.fix_inter_membrane_distance; 
    fix_DOF(DOF_name_and_range.Virus_hight_above_plane)=1;
end

if Minimazation.HD_rim_fixed>0
    DOF_vector(DOF_name_and_range.Virus_Rv)=Minimazation.HD_rim_fixed;
    DOF_vector(DOF_name_and_range.Cell_R_rim)=Minimazation.HD_rim_fixed;

    fix_DOF(DOF_name_and_range.Virus_Rv)=1;
    fix_DOF(DOF_name_and_range.Cell_R_rim)=1;
end

if Minimazation.HD_rim_fixed_virus>0
    DOF_vector(DOF_name_and_range.Virus_Rv)=Minimazation.HD_rim_fixed_virus;
    fix_DOF(DOF_name_and_range.Virus_Rv)=1; 
end

if Minimazation.HD_rim_fixed_cell>0
    DOF_vector(6)=Minimazation.HD_rim_fixed_cell;
    fix_DOF(DOF_name_and_range.Cell_R_rim)=1; 
end


%flat diaphragm    
diaphragm_radius=(Diaphragm.rim.x0^2+Diaphragm.rim.y0^2)^0.5/2^0.5;
% if the diaphragm is fixed to flat or very small
if Minimazation.flat_diaphragm==1 || diaphragm_radius<Minimazation.diaphragm_flat_below 
    % diaphragm rim 
    DOF_vector(DOF_name_and_range.Diaphragm_rim_hl_vector)=0; %set al hl_vector to 0   
    DOF_vector(DOF_name_and_range.Diaphragm_rim_hl_vector(1))=Diaphragm.z0; %the rim hight: hl_vector(1)
    fix_DOF(DOF_name_and_range.Diaphragm_rim_hl_vector)=1; 
    
    %tangent of Diaphragm rim
    DOF_vector(DOF_name_and_range.Diaphragm_rim_drho_dz_vector)=0; %the rim hight: hl_vector(1)
    fix_DOF(DOF_name_and_range.Diaphragm_rim_drho_dz_vector)=1; 
    

end


if Minimazation.do_stalk==1
    fix_DOF(DOF_name_and_range.Diaphragm_rim_x0)=1; %x0
    fix_DOF(DOF_name_and_range.Diaphragm_rim_y0)=1; %y0
    DOF_vector(DOF_name_and_range.Diaphragm_rim_x0)=2*res_struc.pore_open_res;
    DOF_vector(DOF_name_and_range.Diaphragm_rim_y0)=2*res_struc.pore_open_res;
    
    %if stalk exist fix angle to be Cell.rim.drho_dz_vector and Virus.rim.drho_dz_vector +-1
    DOF_vector(DOF_name_and_range.Virus_rim_drho_dz_vector)=0;
    DOF_vector(DOF_name_and_range.Virus_rim_drho_dz_vector(1))=1;
    fix_DOF(DOF_name_and_range.Virus_rim_drho_dz_vector)=1;

    DOF_vector(DOF_name_and_range.Cell_rim_drho_dz_vector)=0;
    DOF_vector(DOF_name_and_range.Cell_rim_drho_dz_vector(1))=-1;
    fix_DOF(DOF_name_and_range.Cell_rim_drho_dz_vector)=1;
    
    % fix diaphragm hight
    DOF_vector(DOF_name_and_range.Diaphragm_rim_hl_vector)=0; %set al hl_vector to 0   
    DOF_vector(DOF_name_and_range.Diaphragm_rim_hl_vector(1))=Diaphragm.z0; %the rim hight: hl_vector(1)
    fix_DOF(DOF_name_and_range.Diaphragm_rim_hl_vector)=1; 
    
    %tangent of Diaphragm rim
    DOF_vector(DOF_name_and_range.Diaphragm_rim_drho_dz_vector)=0; %the rim hight: hl_vector(1)
    fix_DOF(DOF_name_and_range.Diaphragm_rim_drho_dz_vector)=1; 
    
    if  Minimazation.axial_symmetry==1
        
    end
    if strcmp(Minimazation.step,'first') || Minimazation.axial_symmetry==1 
    %Mid plane coeff
        fix_DOF(DOF_name_and_range.Diaphragm_mid_plane_coeff_matrix)=1;
        fix_DOF(DOF_name_and_range.Diaphragm_up_lipid_director_rho_matrix)=1;
        fix_DOF(DOF_name_and_range.Diaphragm_down_lipid_director_rho_matrix)=1;
        %fix_DOF(DOF_name_and_range.Diaphragm_up_lipid_director_phi_matrix)=1;
        %fix_DOF(DOF_name_and_range.Diaphragm_down_lipid_director_phi_matrix)=1;
        
    end
    if strcmp(Minimazation.step,'second') &&  Minimazation.axial_symmetry==0
        fix_DOF(DOF_name_and_range.Diaphragm_Gaussian_decay_vector_up_k_vector_smooth_poly)=1;
        fix_DOF(DOF_name_and_range.Diaphragm_Gaussian_decay_vector_down_k_vector_smooth_poly)=1;
        fix_DOF(DOF_name_and_range.Diaphragm_mid_plane_smooth_polynom_matrix)=1;
        fix_DOF(DOF_name_and_range.Diaphragm_up_lipid_director_rho_smooth_polynom_matrix)=1;
        fix_DOF(DOF_name_and_range.Diaphragm_down_lipid_director_rho_smooth_polynom_matrix)=1;
        %fix_DOF(DOF_name_and_range.Diaphragm_up_lipid_director_phi_smooth_polynom_matrix)=1;
        %fix_DOF(DOF_name_and_range.Diaphragm_down_lipid_director_phi_smooth_polynom_matrix)=1;
    end
    
    if Minimazation.Shell_exist==1
    
    %set the shell-diaphragm to be flat 
    DOF_vector(DOF_name_and_range.Shell_Diaphragm_z_f)=DOF_vector(DOF_name_and_range.Shell_Diaphragm_z_i);
    fix_DOF(DOF_name_and_range.Shell_Diaphragm_z_f)=1;

    DOF_vector(DOF_name_and_range.Shell_Diaphragm_dzdrho_edge)=0;
    fix_DOF(DOF_name_and_range.Shell_Diaphragm_dzdrho_edge)=1;

    DOF_vector(DOF_name_and_range.Shell_Diaphragm_length_from_edge)=0; %the stalk is non-intercating
    fix_DOF(DOF_name_and_range.Shell_Diaphragm_length_from_edge)=1;

    DOF_vector(DOF_name_and_range.Shell_Diaphragm_MP_to_MP_distance_matrix)=0; %the stalk is flat
    fix_DOF(DOF_name_and_range.Shell_Diaphragm_MP_to_MP_distance_matrix)=1;
    
    % set the junction area to be also flat and equal to the diaphragm

    DOF_vector(DOF_name_and_range.Shell_Cell_dzdrho_edge)=0;
    fix_DOF(DOF_name_and_range.Shell_Cell_dzdrho_edge)=1;

    end 

    
end





end

