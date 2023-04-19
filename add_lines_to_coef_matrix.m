function [Diaphragm,Cell,Virus,Shell] = add_lines_to_coef_matrix(Diaphragm_old,Cell_old,Virus_old,Shell_old,Minimazation,res_struc)

%this function corects the coefficiant matrix if the number of dgress of
%freedon is changed. If more are added then it the new DOF are zero if less
%then it just delete the last one

Diaphragm=Diaphragm_old;
Cell=Cell_old;
Virus=Virus_old;
Shell=Shell_old;    


old_z_mid_plane_poly=length(Diaphragm.mid_plane_coeff_matrix(1,:));
old_poly_degree_LD_rho=length(Diaphragm.up_lipid_director_rho_matrix(1,:));
old_poly_degree_LD_phi=length(Diaphragm.up_lipid_director_phi_matrix(1,:));


if old_z_mid_plane_poly<res_struc.poly_degree_z
    Diaphragm.mid_plane_coeff_matrix(:,old_z_mid_plane_poly+1:res_struc.poly_degree_z)=zeros(res_struc.phi_res,res_struc.poly_degree_z-old_z_mid_plane_poly); 
    Cell.mid_plane_coeff_matrix(:,old_z_mid_plane_poly+1:res_struc.poly_degree_z)=zeros(res_struc.phi_res,res_struc.poly_degree_z-old_z_mid_plane_poly);
    Virus.mid_plane_coeff_matrix(:,old_z_mid_plane_poly+1:res_struc.poly_degree_z)=zeros(res_struc.phi_res,res_struc.poly_degree_z-old_z_mid_plane_poly);
end
if old_z_mid_plane_poly>res_struc.poly_degree_z
    Diaphragm.mid_plane_coeff_matrix(:,res_struc.poly_degree_z+1:old_z_mid_plane_poly)=[];
    Cell.mid_plane_coeff_matrix(:,res_struc.poly_degree_z+1:old_z_mid_plane_poly)=[];
    Virus.mid_plane_coeff_matrix(:,res_struc.poly_degree_z+1:old_z_mid_plane_poly)=[];

end



if old_poly_degree_LD_rho<res_struc.poly_degree_LD_rho

    Diaphragm.up_lipid_director_rho_matrix(:,old_poly_degree_LD_rho+1:res_struc.poly_degree_LD_rho)=zeros(res_struc.phi_res,res_struc.poly_degree_LD_rho-old_poly_degree_LD_rho);
    Diaphragm.down_lipid_director_rho_matrix(:,old_poly_degree_LD_rho+1:res_struc.poly_degree_LD_rho)=zeros(res_struc.phi_res,res_struc.poly_degree_LD_rho-old_poly_degree_LD_rho);

    Cell.up_lipid_director_rho_matrix(:,old_poly_degree_LD_rho+1:res_struc.poly_degree_LD_rho)=zeros(res_struc.phi_res,res_struc.poly_degree_LD_rho-old_poly_degree_LD_rho);
    Cell.down_lipid_director_rho_matrix(:,old_poly_degree_LD_rho+1:res_struc.poly_degree_LD_rho)=zeros(res_struc.phi_res,res_struc.poly_degree_LD_rho-old_poly_degree_LD_rho);

    Virus.up_lipid_director_rho_matrix(:,old_poly_degree_LD_rho+1:res_struc.poly_degree_LD_rho)=zeros(res_struc.phi_res,res_struc.poly_degree_LD_rho-old_poly_degree_LD_rho);
    Virus.down_lipid_director_rho_matrix(:,old_poly_degree_LD_rho+1:res_struc.poly_degree_LD_rho)=zeros(res_struc.phi_res,res_struc.poly_degree_LD_rho-old_poly_degree_LD_rho);

end

if old_poly_degree_LD_rho>res_struc.poly_degree_LD_rho

    Diaphragm.up_lipid_director_rho_matrix(:,res_struc.poly_degree_LD_rho+1:old_poly_degree_LD_rho)=[];
    Diaphragm.down_lipid_director_rho_matrix(:,res_struc.poly_degree_LD_rho+1:old_poly_degree_LD_rho)=[];

    Cell.up_lipid_director_rho_matrix(:,res_struc.poly_degree_LD_rho+1:old_poly_degree_LD_rho)=[];
    Cell.down_lipid_director_rho_matrix(:,res_struc.poly_degree_LD_rho+1:old_poly_degree_LD_rho)=[];

    Virus.up_lipid_director_rho_matrix(:,res_struc.poly_degree_LD_rho+1:old_poly_degree_LD_rho)=[];
    Virus.down_lipid_director_rho_matrix(:,res_struc.poly_degree_LD_rho+1:old_poly_degree_LD_rho)=[];

end

if old_poly_degree_LD_phi<res_struc.poly_degree_LD_phi

    Diaphragm.up_lipid_director_phi_matrix(:,old_poly_degree_LD_phi+1:res_struc.poly_degree_LD_phi)=zeros(res_struc.phi_res,res_struc.poly_degree_LD_phi-old_poly_degree_LD_phi);
    Diaphragm.down_lipid_director_phi_matrix(:,old_poly_degree_LD_phi+1:res_struc.poly_degree_LD_phi)=zeros(res_struc.phi_res,res_struc.poly_degree_LD_phi-old_poly_degree_LD_phi);

    Cell.up_lipid_director_phi_matrix(:,old_poly_degree_LD_phi+1:res_struc.poly_degree_LD_phi)=zeros(res_struc.phi_res,res_struc.poly_degree_LD_phi-old_poly_degree_LD_phi);
    Cell.down_lipid_director_phi_matrix(:,old_poly_degree_LD_phi+1:res_struc.poly_degree_LD_phi)=zeros(res_struc.phi_res,res_struc.poly_degree_LD_phi-old_poly_degree_LD_phi);

    Virus.up_lipid_director_phi_matrix(:,old_poly_degree_LD_phi+1:res_struc.poly_degree_LD_phi)=zeros(res_struc.phi_res,res_struc.poly_degree_LD_phi-old_poly_degree_LD_phi);
    Virus.down_lipid_director_phi_matrix(:,old_poly_degree_LD_phi+1:res_struc.poly_degree_LD_phi)=zeros(res_struc.phi_res,res_struc.poly_degree_LD_phi-old_poly_degree_LD_phi);

end

if old_poly_degree_LD_phi>res_struc.poly_degree_LD_phi

    Diaphragm.up_lipid_director_phi_matrix(:,res_struc.poly_degree_LD_phi+1:old_poly_degree_LD_phi)=[];
    Diaphragm.down_lipid_director_phi_matrix(:,res_struc.poly_degree_LD_phi+1:old_poly_degree_LD_phi)=[];

    Cell.up_lipid_director_phi_matrix(:,res_struc.poly_degree_LD_phi+1:old_poly_degree_LD_phi)=[];
    Cell.down_lipid_director_phi_matrix(:,res_struc.poly_degree_LD_phi+1:old_poly_degree_LD_phi)=[];

    Virus.up_lipid_director_phi_matrix(:,res_struc.poly_degree_LD_phi+1:old_poly_degree_LD_phi)=[];
    Virus.down_lipid_director_phi_matrix(:,res_struc.poly_degree_LD_phi+1:old_poly_degree_LD_phi)=[];

end


%%
if Shell.Shell_exist==1
    
    old_shell_poly=length(Shell.Diaphragm.MP_to_MP_distance_matrix(1,:));

    if old_shell_poly<res_struc.poly_degree_shell
        Shell.Diaphragm.MP_to_MP_distance_matrix(:,old_shell_poly+1:res_struc.poly_degree_shell)=zeros(res_struc.phi_res,res_struc.poly_degree_shell-old_shell_poly); 
        Shell.Cell.MP_to_MP_distance_matrix(:,old_shell_poly+1:res_struc.poly_degree_shell)=zeros(res_struc.phi_res,res_struc.poly_degree_shell-old_shell_poly);
    end
    if old_z_mid_plane_poly>res_struc.poly_degree_z
        Shell.Diaphragm.MP_to_MP_distance_matrix(:,res_struc.poly_degree_shell+1:old_shell_poly)=[];
        Shell.Cell.MP_to_MP_distance_matrix(:,res_struc.poly_degree_shell+1:old_shell_poly)=[];

    end

end

%%
if Minimazation.axial_symmetry==0
    old_polar_angle_polynom=length(Diaphragm.mid_plane_smooth_polynom_matrix(1,:));


    %create empty matrix in maximal size
    Diaphragm.mid_plane_smooth_polynom_matrix=zeros(max(res_struc.poly_degree_z,old_z_mid_plane_poly)-4,max(res_struc.polar_angle_polynom,old_polar_angle_polynom));
    Cell.mid_plane_smooth_polynom_matrix=zeros(max(res_struc.poly_degree_z,old_z_mid_plane_poly)-4,max(res_struc.polar_angle_polynom,old_polar_angle_polynom));
    Virus.mid_plane_smooth_polynom_matrix=zeros(max(res_struc.poly_degree_z,old_z_mid_plane_poly)-4,max(res_struc.polar_angle_polynom,old_polar_angle_polynom));
    
    Diaphragm.up_lipid_director_rho_smooth_polynom_matrix=zeros(max(res_struc.poly_degree_LD_rho,old_poly_degree_LD_rho)-4,max(res_struc.polar_angle_polynom,old_polar_angle_polynom));
    Diaphragm.down_lipid_director_rho_smooth_polynom_matrix=zeros(max(res_struc.poly_degree_LD_rho,old_poly_degree_LD_rho)-4,max(res_struc.polar_angle_polynom,old_polar_angle_polynom));
    Cell.up_lipid_director_rho_smooth_polynom_matrix=zeros(max(res_struc.poly_degree_LD_rho,old_poly_degree_LD_rho)-4,max(res_struc.polar_angle_polynom,old_polar_angle_polynom));
    Cell.down_lipid_director_rho_smooth_polynom_matrix=zeros(max(res_struc.poly_degree_LD_rho,old_poly_degree_LD_rho)-4,max(res_struc.polar_angle_polynom,old_polar_angle_polynom));
    Virus.up_lipid_director_rho_smooth_polynom_matrix=zeros(max(res_struc.poly_degree_LD_rho,old_poly_degree_LD_rho)-4,max(res_struc.polar_angle_polynom,old_polar_angle_polynom));
    Virus.down_lipid_director_rho_smooth_polynom_matrix=zeros(max(res_struc.poly_degree_LD_rho,old_poly_degree_LD_rho)-4,max(res_struc.polar_angle_polynom,old_polar_angle_polynom));

    Diaphragm.Gaussian_decay_vector_up.k_vector_smooth_poly=zeros(1,max(res_struc.polar_angle_polynom,old_polar_angle_polynom));
    Diaphragm.Gaussian_decay_vector_down.k_vector_smooth_poly=zeros(1,max(res_struc.polar_angle_polynom,old_polar_angle_polynom));
    Cell.Gaussian_decay_vector_up.k_vector_smooth_poly=zeros(1,max(res_struc.polar_angle_polynom,old_polar_angle_polynom));
    Cell.Gaussian_decay_vector_down.k_vector_smooth_poly=zeros(1,max(res_struc.polar_angle_polynom,old_polar_angle_polynom));
    Virus.Gaussian_decay_vector_up.k_vector_smooth_poly=zeros(1,max(res_struc.polar_angle_polynom,old_polar_angle_polynom));
    Virus.Gaussian_decay_vector_down.k_vector_smooth_poly=zeros(1,max(res_struc.polar_angle_polynom,old_polar_angle_polynom));

    %plug in the old matrix layer
    Diaphragm.mid_plane_smooth_polynom_matrix(1:old_z_mid_plane_poly-4,1:old_polar_angle_polynom)=Diaphragm_old.mid_plane_smooth_polynom_matrix;
    Cell.mid_plane_smooth_polynom_matrix(1:old_z_mid_plane_poly-4,1:old_polar_angle_polynom)=Cell_old.mid_plane_smooth_polynom_matrix;
    Virus.mid_plane_smooth_polynom_matrix(1:old_z_mid_plane_poly-4,1:old_polar_angle_polynom)=Virus_old.mid_plane_smooth_polynom_matrix;

    Diaphragm.up_lipid_director_rho_smooth_polynom_matrix(1:old_poly_degree_LD_rho-4,1:old_polar_angle_polynom)=Diaphragm_old.up_lipid_director_rho_smooth_polynom_matrix;
    Diaphragm.down_lipid_director_rho_smooth_polynom_matrix(1:old_poly_degree_LD_rho-4,1:old_polar_angle_polynom)=Diaphragm_old.down_lipid_director_rho_smooth_polynom_matrix;
    Cell.up_lipid_director_rho_smooth_polynom_matrix(1:old_poly_degree_LD_rho-4,1:old_polar_angle_polynom)=Cell_old.up_lipid_director_rho_smooth_polynom_matrix;
    Cell.down_lipid_director_rho_smooth_polynom_matrix(1:old_poly_degree_LD_rho-4,1:old_polar_angle_polynom)=Cell_old.down_lipid_director_rho_smooth_polynom_matrix;
    Virus.up_lipid_director_rho_smooth_polynom_matrix(1:old_poly_degree_LD_rho-4,1:old_polar_angle_polynom)=Virus_old.up_lipid_director_rho_smooth_polynom_matrix;
    Virus.down_lipid_director_rho_smooth_polynom_matrix(1:old_poly_degree_LD_rho-4,1:old_polar_angle_polynom)=Virus_old.down_lipid_director_rho_smooth_polynom_matrix;

    Diaphragm.Gaussian_decay_vector_up.k_vector_smooth_poly(1:old_polar_angle_polynom)=Diaphragm_old.Gaussian_decay_vector_up.k_vector_smooth_poly;
    Diaphragm.Gaussian_decay_vector_down.k_vector_smooth_poly(1:old_polar_angle_polynom)=Diaphragm_old.Gaussian_decay_vector_down.k_vector_smooth_poly;
    Cell.Gaussian_decay_vector_up.k_vector_smooth_poly(1:old_polar_angle_polynom)=Cell_old.Gaussian_decay_vector_up.k_vector_smooth_poly;
    Cell.Gaussian_decay_vector_down.k_vector_smooth_poly(1:old_polar_angle_polynom)=Cell_old.Gaussian_decay_vector_down.k_vector_smooth_poly;
    Virus.Gaussian_decay_vector_up.k_vector_smooth_poly(1:old_polar_angle_polynom)=Virus_old.Gaussian_decay_vector_up.k_vector_smooth_poly;
    Virus.Gaussian_decay_vector_down.k_vector_smooth_poly(1:old_polar_angle_polynom)=Virus_old.Gaussian_decay_vector_down.k_vector_smooth_poly;

    if old_polar_angle_polynom>res_struc.polar_angle_polynom
        Diaphragm_old.mid_plane_smooth_polynom_matrix(:,res_struc.polar_angle_polynom+1:old_polar_angle_polynom)=[];
        Cell.mid_plane_smooth_polynom_matrix(:,res_struc.polar_angle_polynom+1:old_polar_angle_polynom)=[];
        Virus.mid_plane_smooth_polynom_matrix(:,res_struc.polar_angle_polynom+1:old_polar_angle_polynom)=[];

        Diaphragm.down_lipid_director_rho_smooth_polynom_matrix(:,res_struc.polar_angle_polynom+1:old_polar_angle_polynom)=[];
        Diaphragm.up_lipid_director_rho_smooth_polynom_matrix(:,res_struc.polar_angle_polynom+1:old_polar_angle_polynom)=[];
        Cell.up_lipid_director_rho_smooth_polynom_matrix(:,res_struc.polar_angle_polynom+1:old_polar_angle_polynom)=[];
        Cell.down_lipid_director_rho_smooth_polynom_matrix(:,res_struc.polar_angle_polynom+1:old_polar_angle_polynom)=[];
        Virus.up_lipid_director_rho_smooth_polynom_matrix(:,res_struc.polar_angle_polynom+1:old_polar_angle_polynom)=[];
        Virus.down_lipid_director_rho_smooth_polynom_matrix(:,res_struc.polar_angle_polynom+1:old_polar_angle_polynom)=[];
    
        Virus.Gaussian_decay_vector_up.k_vector_smooth_poly(res_struc.polar_angle_polynom+1:old_polar_angle_polynom)=[];
        Diaphragm.Gaussian_decay_vector_up.k_vector_smooth_poly(res_struc.polar_angle_polynom+1:old_polar_angle_polynom)=[];
        Cell.Gaussian_decay_vector_up.k_vector_smooth_poly(res_struc.polar_angle_polynom+1:old_polar_angle_polynom)=[];
        Diaphragm.Gaussian_decay_vector_down.k_vector_smooth_poly(res_struc.polar_angle_polynom+1:old_polar_angle_polynom)=[];
        Cell.Gaussian_decay_vector_down.k_vector_smooth_poly(res_struc.polar_angle_polynom+1:old_polar_angle_polynom)=[];
        Virus.Gaussian_decay_vector_down.k_vector_smooth_poly(res_struc.polar_angle_polynom+1:old_polar_angle_polynom)=[];
    end




    if old_z_mid_plane_poly>res_struc.poly_degree_z
        Diaphragm.mid_plane_smooth_polynom_matrix(res_struc.poly_degree_z+1-4:old_z_mid_plane_poly-4,:)=[];
        Cell.mid_plane_smooth_polynom_matrix(res_struc.poly_degree_z+1-4:old_z_mid_plane_poly-4,:)=[];
        Virus.mid_plane_smooth_polynom_matrix(res_struc.poly_degree_z+1-4:old_z_mid_plane_poly-4,:)=[];
    end

    if old_poly_degree_LD_rho>res_struc.poly_degree_LD_rho

        Diaphragm.down_lipid_director_rho_smooth_polynom_matrix(res_struc.poly_degree_LD_rho+1-4:old_poly_degree_LD_rho-4,:)=[];
        Diaphragm.up_lipid_director_rho_smooth_polynom_matrix(res_struc.poly_degree_LD_rho+1-4:old_poly_degree_LD_rho-4,:)=[];
        Cell.up_lipid_director_rho_smooth_polynom_matrix(res_struc.poly_degree_LD_rho+1-4:old_poly_degree_LD_rho-4,:)=[];
        Cell.down_lipid_director_rho_smooth_polynom_matrix(res_struc.poly_degree_LD_rho+1-4:old_poly_degree_LD_rho-4,:)=[];
        Virus.up_lipid_director_rho_smooth_polynom_matrix(res_struc.poly_degree_LD_rho+1-4:old_poly_degree_LD_rho-4,:)=[];
        Virus.down_lipid_director_rho_smooth_polynom_matrix(res_struc.poly_degree_LD_rho+1-4:old_poly_degree_LD_rho-4,:)=[];

    end

end



end

