function [Diaphragm,Cell,Virus,Shell, exit_status_global,Minimazation] = ...
    build_Hemifusion_structure(Diaphragm_old,Cell_old,Virus_old,Shell_old,Minimazation,res_struc,General_physical_properties)

Diaphragm=Diaphragm_old;
Cell=Cell_old;
Virus=Virus_old;
Shell=Shell_old;

% exit_status_membrane 1 if successful 0 if error 
exit_status_global=1;
exit_status_membrane=1;
exit_status_shell=1;


lipid_length=General_physical_properties.lipid_length;
streatching_modulus=General_physical_properties.streatching_modulus;
tilt_modulus=General_physical_properties.tilt_modulus;
chi=General_physical_properties.chi;

% matrix containg cos(2*int*phi)  col is int (starting from 0) rows is phi 0->pi/2
[cosine_matrix] = cosine_matrix_creatore(res_struc.phi_res,res_struc.poly_degree_rim); 

if Shell.Shell_exist==1
    Shell.Shell_physical_proprties.kappa=Shell.Shell_physical_proprties.Youngs_modulus*Shell.Shell_physical_proprties.width^3/(12*(1-Shell.Shell_physical_proprties.P_ratio^2));
    Shell.Shell_physical_proprties.kappa_bar=-Shell.Shell_physical_proprties.Youngs_modulus*Shell.Shell_physical_proprties.width^3/(12*(1-Shell.Shell_physical_proprties.P_ratio));
    Shell.Shell_physical_proprties.streatch_modulus=Shell.Shell_physical_proprties.Youngs_modulus*Shell.Shell_physical_proprties.width/(2*(1-Shell.Shell_physical_proprties.P_ratio));
    Shell.Shell_physical_proprties.monolayer_width=General_physical_properties.lipid_length;
end

%build matrix polynom if it is given by smoothing function
if Minimazation.non_axial_symmetric_polynom==1 && Minimazation.axial_symmetry==0 && strcmp(Minimazation.step,'second') 
    
	[Diaphragm.mid_plane_coeff_matrix] = angle_coefficiants_polynom...
        (Diaphragm.mid_plane_coeff_matrix,Diaphragm.mid_plane_smooth_polynom_matrix,5,res_struc.poly_degree_z,res_struc,'hamonic');
	[Cell.mid_plane_coeff_matrix] = angle_coefficiants_polynom...
        (Cell.mid_plane_coeff_matrix,Cell.mid_plane_smooth_polynom_matrix,5,res_struc.poly_degree_z,res_struc,'hamonic');
    [Virus.mid_plane_coeff_matrix] = angle_coefficiants_polynom...
        (Virus.mid_plane_coeff_matrix,Virus.mid_plane_smooth_polynom_matrix,5,res_struc.poly_degree_z,res_struc,'hamonic');
    
    [Diaphragm.up_lipid_director_rho_matrix] = angle_coefficiants_polynom...
        (Diaphragm.up_lipid_director_rho_matrix,Diaphragm.up_lipid_director_rho_smooth_polynom_matrix,5,res_struc.poly_degree_LD_rho,res_struc,'hamonic');
    [Diaphragm.down_lipid_director_rho_matrix] = angle_coefficiants_polynom...
        (Diaphragm.down_lipid_director_rho_matrix,Diaphragm.down_lipid_director_rho_smooth_polynom_matrix,5,res_struc.poly_degree_LD_rho,res_struc,'hamonic');
    
    [Cell.up_lipid_director_rho_matrix] = angle_coefficiants_polynom...
        (Cell.up_lipid_director_rho_matrix,Cell.up_lipid_director_rho_smooth_polynom_matrix,5,res_struc.poly_degree_LD_rho,res_struc,'hamonic');
    [Cell.down_lipid_director_rho_matrix] = angle_coefficiants_polynom...
        (Cell.down_lipid_director_rho_matrix,Cell.down_lipid_director_rho_smooth_polynom_matrix,5,res_struc.poly_degree_LD_rho,res_struc,'hamonic'); 

    [Virus.up_lipid_director_rho_matrix] = angle_coefficiants_polynom...
        (Virus.up_lipid_director_rho_matrix,Virus.up_lipid_director_rho_smooth_polynom_matrix,5,res_struc.poly_degree_LD_rho,res_struc,'hamonic');
    [Virus.down_lipid_director_rho_matrix] = angle_coefficiants_polynom...
        (Virus.down_lipid_director_rho_matrix,Virus.down_lipid_director_rho_smooth_polynom_matrix,5,res_struc.poly_degree_LD_rho,res_struc,'hamonic'); 
    
    
    if Minimazation.exponential_decay_tilt==1
        
        [Diaphragm.Gaussian_decay_vector_up.k_vector] = angle_coefficiants_polynom...
            (Diaphragm.Gaussian_decay_vector_up.k_vector,Diaphragm.Gaussian_decay_vector_up.k_vector_smooth_poly,1,1,res_struc,'hamonic');       
        [Diaphragm.Gaussian_decay_vector_down.k_vector] = angle_coefficiants_polynom...
            (Diaphragm.Gaussian_decay_vector_down.k_vector,Diaphragm.Gaussian_decay_vector_down.k_vector_smooth_poly,1,1,res_struc,'hamonic');           
        
        [Cell.Gaussian_decay_vector_up.k_vector] = angle_coefficiants_polynom...
            (Cell.Gaussian_decay_vector_up.k_vector,Cell.Gaussian_decay_vector_up.k_vector_smooth_poly,1,1,res_struc,'hamonic');       
        [Cell.Gaussian_decay_vector_down.k_vector] = angle_coefficiants_polynom...
            (Cell.Gaussian_decay_vector_down.k_vector,Cell.Gaussian_decay_vector_down.k_vector_smooth_poly,1,1,res_struc,'hamonic');           
        
        [Virus.Gaussian_decay_vector_up.k_vector] = angle_coefficiants_polynom...
            (Virus.Gaussian_decay_vector_up.k_vector,Virus.Gaussian_decay_vector_up.k_vector_smooth_poly,1,1,res_struc,'hamonic');       
        [Virus.Gaussian_decay_vector_down.k_vector] = angle_coefficiants_polynom...
            (Virus.Gaussian_decay_vector_down.k_vector,Virus.Gaussian_decay_vector_down.k_vector_smooth_poly,1,1,res_struc,'hamonic');    
        
        
    end
end

if Minimazation.exponential_decay_tilt==1
    Virus.Gaussian_decay_vector_up.exponential_decay_tilt=1;
    Cell.Gaussian_decay_vector_up.exponential_decay_tilt=1;
    Diaphragm.Gaussian_decay_vector_up.exponential_decay_tilt=1;
    Virus.Gaussian_decay_vector_down.exponential_decay_tilt=1;
    Cell.Gaussian_decay_vector_down.exponential_decay_tilt=1;
    Diaphragm.Gaussian_decay_vector_down.exponential_decay_tilt=1;
else
    Virus.Gaussian_decay_vector_up.exponential_decay_tilt=0;
    Cell.Gaussian_decay_vector_up.exponential_decay_tilt=0;
    Diaphragm.Gaussian_decay_vector_up.exponential_decay_tilt=0;
    Virus.Gaussian_decay_vector_down.exponential_decay_tilt=0;
    Cell.Gaussian_decay_vector_down.exponential_decay_tilt=0;
    Diaphragm.Gaussian_decay_vector_down.exponential_decay_tilt=0;  
end
if Minimazation.HD_rim_fixed>0 
    Virus.Rv=Minimazation.HD_rim_fixed;
    Cell.R_rim=Minimazation.HD_rim_fixed;
end

%calcualte the shape of the Diaphragm rim
[rho_vector_Diaphragm_rim,z_vector_Diaphragm_rim,drho_dz_Diaphragm_rim] = diaphragm_rim_contstruct...
    (Diaphragm.rim,cosine_matrix,res_struc.phi_res);

%create Diaphragm mid plane
% constraint Diaphragm 
Diaphragm.constraint.mid_plane.zi=ones(1,res_struc.phi_res)*Diaphragm.z0; % set hight of center diaphragm
Diaphragm.constraint.mid_plane.zf=z_vector_Diaphragm_rim; % hight at edge of Diaphragm
Diaphragm.constraint.mid_plane.dzdrhoi=zeros(1,res_struc.phi_res); %set zero angle at rho=0
Diaphragm.constraint.mid_plane.dzdrhof=drho_dz_Diaphragm_rim; %dervitive of hight at edge
Diaphragm.constraint.mid_plane.rhoi=ones(1,res_struc.phi_res)*Diaphragm.r_pore; %center of Diaphragm is rho=Diaphragm.r_pore
Diaphragm.constraint.mid_plane.rhof=rho_vector_Diaphragm_rim; %find position of Diaphragm rho

%coeffciant matrix
% force boundary conditions
[Diaphragm.mid_plane_coeff_matrix,exit_status,Minimazation] = Force_poly_phi...
    (Diaphragm.mid_plane_coeff_matrix,Diaphragm.constraint.mid_plane,res_struc.phi_res,res_struc.poly_degree_z,Minimazation); 
if exit_status==0
    exit_status_membrane=0;
end

%create virues fusion site mid-plane

Virus.z_v_center=Virus.r_bv+Virus.r_sv+Virus.hight_above_plane; % torus center

%tagent angle between virus mid plan and diphragm


% find the angle and z at the torus where rho=Rv
[rho_vector_virues,z_vector_virues,drho_dz_virues] = virus_rim_contstruct(Virus,res_struc.phi_res);

Virus.constraint.mid_plane.zi=z_vector_Diaphragm_rim; %z of Diaphragm rim is inner boudary
Virus.constraint.mid_plane.zf=z_vector_virues; %z of torus rim is outer boudary
Virus.constraint.mid_plane.dzdrhoi=Virus.rim.drho_dz_vector*cosine_matrix; %angle at edge
Virus.constraint.mid_plane.dzdrhof=drho_dz_virues;  %angle at torus rim is outer boudary angle
Virus.constraint.mid_plane.rhoi=rho_vector_Diaphragm_rim; % Diaphragm rim
Virus.constraint.mid_plane.rhof=rho_vector_virues;   %fusion site size at torus

[Virus.mid_plane_coeff_matrix,exit_status,Minimazation] = Force_poly_phi...
    (Virus.mid_plane_coeff_matrix,Virus.constraint.mid_plane,res_struc.phi_res,res_struc.poly_degree_z,Minimazation); % force boundary conditions
if exit_status==0
    exit_status_membrane=0;
end



% create mid plane of cell membrane
[rho_vector_cell,z_vector_cell,drho_dz_cell,~,~] = cell_rim_contstruct(Cell,res_struc.phi_res);

Cell.constraint.mid_plane.zi=z_vector_Diaphragm_rim; %inner boundary is Diaphragm rim
Cell.constraint.mid_plane.zf=z_vector_cell;
Cell.constraint.mid_plane.dzdrhoi=Cell.rim.drho_dz_vector*cosine_matrix;
Cell.constraint.mid_plane.dzdrhof=drho_dz_cell;
Cell.constraint.mid_plane.rhoi=rho_vector_Diaphragm_rim;
Cell.constraint.mid_plane.rhof=rho_vector_cell;


[Cell.mid_plane_coeff_matrix,exit_status,Minimazation] = Force_poly_phi...
    (Cell.mid_plane_coeff_matrix,Cell.constraint.mid_plane,res_struc.phi_res,res_struc.poly_degree_z,Minimazation);
if exit_status==0
    exit_status_membrane=0;
end


%% Tilt angle

%boundary conditions at the rim
alpha_d=atan(drho_dz_Diaphragm_rim); % angle at edge of diaphragm from rho axis
alpha_c=atan(Cell.constraint.mid_plane.dzdrhoi); % angle at edge of cell from rho axis
alpha_v=atan(Virus.constraint.mid_plane.dzdrhoi); % angle at edge of virus from rho axis

theta_rho_diaphragm_virues=(alpha_v-alpha_d)/2; % angle between lipid director and normal at diaphragm virues interface
theta_rho_diaphragm_cell=(alpha_c-alpha_d)/2; % angle between lipid director and normal at diaphragm virues interface
theta_rho_virues_cell=(pi+alpha_c-alpha_v)/2; % angle between lipid director and normal at cell virues interface



%% DIAPHRAGM

% lipid director in rho direction diaphragm
% constraint Diaphragm for rho direction up and down
Diaphragm.constraint.rho_lipid_angle_up.zi=zeros(1,res_struc.phi_res); % parallel to normal at rho=0
Diaphragm.constraint.rho_lipid_angle_up.zf=theta_rho_diaphragm_virues; % tilt angle at edge of Diaphragm and virus
Diaphragm.constraint.rho_lipid_angle_up.dzdrhoi=ones(1,res_struc.phi_res)*0; % drivitive of angle at edge of Diaphragm and virus
Diaphragm.constraint.rho_lipid_angle_up.dzdrhof=ones(1,res_struc.phi_res)*0; % drivitive of angle at edge of Diaphragm and virus
Diaphragm.constraint.rho_lipid_angle_up.rhoi=ones(1,res_struc.phi_res)*Diaphragm.r_pore; %center of Diaphragm is rho=0
Diaphragm.constraint.rho_lipid_angle_up.rhof=rho_vector_Diaphragm_rim; %find position of Diaphragm rho
if Minimazation.do_stalk==1
    Diaphragm.constraint.rho_lipid_angle_up.zf=zeros(1,res_struc.phi_res); % no tilt if stalk exist
end
% find coeff of lipid director in rho direction in diaphragm
[Diaphragm.up_lipid_director_rho_matrix,exit_status,Minimazation] = Force_poly_phi...
    (Diaphragm.up_lipid_director_rho_matrix,Diaphragm.constraint.rho_lipid_angle_up,res_struc.phi_res,res_struc.poly_degree_LD_rho,Minimazation);
if exit_status==0
    exit_status_membrane=0;
end

% constraint Diaphragm for rho direction down 
Diaphragm.constraint.rho_lipid_angle_down.zi=ones(1,res_struc.phi_res)*0; % set angle of center diaphragm to -pi/2
Diaphragm.constraint.rho_lipid_angle_down.zf=theta_rho_diaphragm_cell; % tilt angle at edge of Diaphragm and virus
Diaphragm.constraint.rho_lipid_angle_down.dzdrhoi=ones(1,res_struc.phi_res)*0; % drivitive of angle at edge of Diaphragm and virus
Diaphragm.constraint.rho_lipid_angle_down.dzdrhof=ones(1,res_struc.phi_res)*0; % drivitive of angle at edge of Diaphragm and virus
Diaphragm.constraint.rho_lipid_angle_down.rhoi=ones(1,res_struc.phi_res)*Diaphragm.r_pore; %center of Diaphragm is rho=*Diaphragm.r_pore
Diaphragm.constraint.rho_lipid_angle_down.rhof=rho_vector_Diaphragm_rim; %find position of Diaphragm rho
if Minimazation.do_stalk==1
    Diaphragm.constraint.rho_lipid_angle_down.zf=zeros(1,res_struc.phi_res); % no tilt if stalk exist
end
% find coeff of lipid director in rho direction in diaphragm

[Diaphragm.down_lipid_director_rho_matrix,exit_status,Minimazation] = Force_poly_phi...
    (Diaphragm.down_lipid_director_rho_matrix,Diaphragm.constraint.rho_lipid_angle_down,res_struc.phi_res,res_struc.poly_degree_LD_rho,Minimazation);
if exit_status==0
    exit_status_membrane=0;
end



% constraint Diaphragm for phi direction up 
Diaphragm.constraint.phi_lipid_angle_up.zi=ones(1,res_struc.phi_res)*0; %director phi_direction at phi=0
Diaphragm.constraint.phi_lipid_angle_up.zf=ones(1,res_struc.phi_res)*0; %director phi_direction at phi=0
Diaphragm.constraint.phi_lipid_angle_up.dzdrhoi=zeros(1,res_struc.phi_res); %mirror symmetry at phi=0
Diaphragm.constraint.phi_lipid_angle_up.dzdrhof=zeros(1,res_struc.phi_res);  %mirror symmetry at phi=pi/2
Diaphragm.constraint.phi_lipid_angle_up.rhoi=ones(1,res_struc.phi_res)*Diaphragm.r_pore;   % start is phi=*Diaphragm.r_pore
Diaphragm.constraint.phi_lipid_angle_up.rhof=rho_vector_Diaphragm_rim;   %end is phi=pi/2
% find coeff of lipid director in rho direction in diaphragm

[Diaphragm.up_lipid_director_phi_matrix,exit_status,Minimazation] = Force_poly_phi...
    (Diaphragm.up_lipid_director_phi_matrix,Diaphragm.constraint.phi_lipid_angle_up,res_struc.phi_res,res_struc.poly_degree_LD_phi,Minimazation); 
if exit_status==0
    exit_status_membrane=0;
end



% constraint Diaphragm for phi direction down 
Diaphragm.constraint.phi_lipid_angle_down.zi=ones(1,res_struc.phi_res)*0; %director phi_direction at phi=0
Diaphragm.constraint.phi_lipid_angle_down.zf=ones(1,res_struc.phi_res)*0; %director phi_direction at phi=pi/2
Diaphragm.constraint.phi_lipid_angle_down.dzdrhoi=zeros(1,res_struc.phi_res); %mirror symmetry at phi=0
Diaphragm.constraint.phi_lipid_angle_down.dzdrhof=zeros(1,res_struc.phi_res);  %mirror symmetry at phi=pi/2
Diaphragm.constraint.phi_lipid_angle_down.rhoi=ones(1,res_struc.phi_res)*Diaphragm.r_pore;   % start is phi=0
Diaphragm.constraint.phi_lipid_angle_down.rhof=rho_vector_Diaphragm_rim;   %end is phi=pi/2

% find coeff of lipid director in rho direction in diaphragm


[Diaphragm.down_lipid_director_phi_matrix,exit_status,Minimazation] = Force_poly_phi....
    (Diaphragm.down_lipid_director_phi_matrix,Diaphragm.constraint.phi_lipid_angle_down,res_struc.phi_res,res_struc.poly_degree_LD_phi,Minimazation); 
if exit_status==0
    exit_status_membrane=0;
end




%% virues membrane
% lipid director in rho direction diaphragm
% constraint Diaphragm for rho direction up and down

Virus.constraint.rho_lipid_angle_up.zi=-theta_rho_diaphragm_virues; % tilt angle at edge of Diaphragm and virus
Virus.constraint.rho_lipid_angle_up.zf=ones(1,res_struc.phi_res)*0; % at edge must be perelel to virus normal
Virus.constraint.rho_lipid_angle_up.dzdrhoi=ones(1,res_struc.phi_res)*0; % drivitive of angle at edge of Diaphragm and virus
Virus.constraint.rho_lipid_angle_up.dzdrhof=ones(1,res_struc.phi_res)*0; % drivitive of angle at edge of Diaphragm and virus
Virus.constraint.rho_lipid_angle_up.rhoi=rho_vector_Diaphragm_rim; %center of Diaphragm is rho=0
Virus.constraint.rho_lipid_angle_up.rhof=rho_vector_virues; %find position of Diaphragm rho
if Minimazation.do_stalk==1
    Virus.constraint.rho_lipid_angle_up.zi=-pi/4*ones(1,res_struc.phi_res); % 45deg tilt angle
end
% find coeff of lipid director in rho direction in diaphragm
[Virus.up_lipid_director_rho_matrix,exit_status,Minimazation] = Force_poly_phi...
    (Virus.up_lipid_director_rho_matrix,Virus.constraint.rho_lipid_angle_up,res_struc.phi_res,res_struc.poly_degree_LD_rho,Minimazation);
if exit_status==0
    exit_status_membrane=0;
end



% constraint Diaphragm for phi direction up 
Virus.constraint.phi_lipid_angle_up.zi=ones(1,res_struc.phi_res)*0;  %director phi_direction at n=N
Virus.constraint.phi_lipid_angle_up.zf=ones(1,res_struc.phi_res)*0;  %director phi_direction at n=N
Virus.constraint.phi_lipid_angle_up.dzdrhoi=zeros(1,res_struc.phi_res); %mirror symmetry at phi=0
Virus.constraint.phi_lipid_angle_up.dzdrhof=zeros(1,res_struc.phi_res); %mirror symmetry at phi=pi/2
Virus.constraint.phi_lipid_angle_up.rhoi=rho_vector_Diaphragm_rim;    % start is phi=0
Virus.constraint.phi_lipid_angle_up.rhof=rho_vector_virues;   %end is phi=pi/2
% find coeff of lipid director in rho direction in diaphragm


[Virus.up_lipid_director_phi_matrix,exit_status,Minimazation] = ...
    Force_poly_phi(Virus.up_lipid_director_phi_matrix,Virus.constraint.phi_lipid_angle_up,res_struc.phi_res,res_struc.poly_degree_LD_phi,Minimazation); 
if exit_status==0
    exit_status_membrane=0;
end



%down lipid director
Virus.constraint.rho_lipid_angle_down.zi=theta_rho_virues_cell; % tilt angle at edge of Diaphragm and virus
Virus.constraint.rho_lipid_angle_down.zf=ones(1,res_struc.phi_res)*0; % at edge must be perelel to virus normal
Virus.constraint.rho_lipid_angle_down.dzdrhoi=ones(1,res_struc.phi_res)*0; % drivitive of angle at edge of Diaphragm and virus
Virus.constraint.rho_lipid_angle_down.dzdrhof=ones(1,res_struc.phi_res)*0; % drivitive of angle at edge of Diaphragm and virus
Virus.constraint.rho_lipid_angle_down.rhoi=rho_vector_Diaphragm_rim; %center of Diaphragm is rho=0
Virus.constraint.rho_lipid_angle_down.rhof=rho_vector_virues; %find position of Diaphragm rho
% find coeff of lipid director in rho direction in diaphragm
[Virus.down_lipid_director_rho_matrix,exit_status,Minimazation] = Force_poly_phi...
    (Virus.down_lipid_director_rho_matrix,Virus.constraint.rho_lipid_angle_down,res_struc.phi_res,res_struc.poly_degree_LD_rho,Minimazation);
if exit_status==0
    exit_status_membrane=0;
end

Virus.constraint.phi_lipid_angle_down.zi=ones(1,res_struc.phi_res)*0;  %director phi_direction at phi=0
Virus.constraint.phi_lipid_angle_down.zf=ones(1,res_struc.phi_res)*0;  %director phi_direction at phi=pi/2
Virus.constraint.phi_lipid_angle_down.dzdrhoi=zeros(1,res_struc.phi_res); %mirror symmetry at phi=0
Virus.constraint.phi_lipid_angle_down.dzdrhof=zeros(1,res_struc.phi_res); %mirror symmetry at phi=pi/2
Virus.constraint.phi_lipid_angle_down.rhoi=rho_vector_Diaphragm_rim;    % start is phi=0
Virus.constraint.phi_lipid_angle_down.rhof=rho_vector_virues;   %end is phi=pi/2
% find coeff of lipid director in rho direction in diaphragm

[Virus.down_lipid_director_phi_matrix,exit_status,Minimazation] = Force_poly_phi...
    (Virus.down_lipid_director_phi_matrix,Virus.constraint.phi_lipid_angle_down,res_struc.phi_res,res_struc.poly_degree_LD_phi,Minimazation); 
if exit_status==0
    exit_status_membrane=0;
end



%% Cell membrane

% lipid director in rho direction diaphragm
% constraint Diaphragm for rho direction up and down
Cell.constraint.rho_lipid_angle_up.zi=-theta_rho_virues_cell;           % tilt angle at edge of Diaphragm and virus
Cell.constraint.rho_lipid_angle_up.zf=ones(1,res_struc.phi_res)*0;      % at edge must be perelel to virus normal
Cell.constraint.rho_lipid_angle_up.dzdrhoi=ones(1,res_struc.phi_res)*0; % drivitive of angle at edge of Diaphragm and virus
Cell.constraint.rho_lipid_angle_up.dzdrhof=ones(1,res_struc.phi_res)*0; % drivitive of angle at edge of Diaphragm and virus
Cell.constraint.rho_lipid_angle_up.rhoi=rho_vector_Diaphragm_rim;       % edge of diaphragm rim
Cell.constraint.rho_lipid_angle_up.rhof=rho_vector_cell;                % cell edge
% find coeff of lipid director in rho direction in diaphragm

[Cell.up_lipid_director_rho_matrix,exit_status,Minimazation] = Force_poly_phi...
    (Cell.up_lipid_director_rho_matrix,Cell.constraint.rho_lipid_angle_up,res_struc.phi_res,res_struc.poly_degree_LD_rho,Minimazation);
if exit_status==0
    exit_status_membrane=0;
end

% constraint Diaphragm for phi direction up 
Cell.constraint.phi_lipid_angle_up.zi=ones(1,res_struc.phi_res)*0;     % tilt angle at junction =0
Cell.constraint.phi_lipid_angle_up.zf=ones(1,res_struc.phi_res)*0;     % tilt angle at edge =0
Cell.constraint.phi_lipid_angle_up.dzdrhoi=zeros(1,res_struc.phi_res); % mirror symmetry at phi=0
Cell.constraint.phi_lipid_angle_up.dzdrhof=zeros(1,res_struc.phi_res); % mirror symmetry at phi=pi/2
Cell.constraint.phi_lipid_angle_up.rhoi=rho_vector_Diaphragm_rim;      % edge of diaphragm rim
Cell.constraint.phi_lipid_angle_up.rhof=rho_vector_cell;               % cell edge
% Cell coeff of lipid director in rho direction in diaphragm


[Cell.up_lipid_director_phi_matrix,exit_status,Minimazation] = Force_poly_phi...
    (Cell.up_lipid_director_phi_matrix,Cell.constraint.phi_lipid_angle_up,res_struc.phi_res,res_struc.poly_degree_LD_phi,Minimazation); 
if exit_status==0
    exit_status_membrane=0;
end


%down lipid director
Cell.constraint.rho_lipid_angle_down.zi=-theta_rho_diaphragm_cell;           % tilt angle at edge of Diaphragm and virus
Cell.constraint.rho_lipid_angle_down.zf=ones(1,res_struc.phi_res)*0;      % at edge must be perelel to virus normal
Cell.constraint.rho_lipid_angle_down.dzdrhoi=ones(1,res_struc.phi_res)*0; % drivitive of angle at edge of Diaphragm and virus
Cell.constraint.rho_lipid_angle_down.dzdrhof=ones(1,res_struc.phi_res)*0; % drivitive of angle at edge of Diaphragm and virus
Cell.constraint.rho_lipid_angle_down.rhoi=rho_vector_Diaphragm_rim;       % edge of diaphragm rim
Cell.constraint.rho_lipid_angle_down.rhof=rho_vector_cell;                % cell edge
% find coeff of lipid director in rho direction in diaphragm
if Minimazation.do_stalk==1
    Cell.constraint.rho_lipid_angle_down.zi=pi/4*ones(1,res_struc.phi_res); % 45deg tilt angle
end
[Cell.down_lipid_director_rho_matrix,exit_status,Minimazation] = Force_poly_phi...
    (Cell.down_lipid_director_rho_matrix,Cell.constraint.rho_lipid_angle_down,res_struc.phi_res,res_struc.poly_degree_LD_rho,Minimazation);
if exit_status==0
    exit_status_membrane=0;
end


Cell.constraint.phi_lipid_angle_down.zi=ones(1,res_struc.phi_res)*0;  %director phi_direction at phi=0
Cell.constraint.phi_lipid_angle_down.zf=ones(1,res_struc.phi_res)*0;  %director phi_direction at phi=pi/2
Cell.constraint.phi_lipid_angle_down.dzdrhoi=zeros(1,res_struc.phi_res); %mirror symmetry at phi=0
Cell.constraint.phi_lipid_angle_down.dzdrhof=zeros(1,res_struc.phi_res); %mirror symmetry at phi=pi/2
Cell.constraint.phi_lipid_angle_down.rhoi=rho_vector_Diaphragm_rim;    % start is phi=0
Cell.constraint.phi_lipid_angle_down.rhof=rho_vector_cell;   %end is phi=pi/2
% find coeff of lipid director in rho direction in diaphragm

[Cell.down_lipid_director_phi_matrix,exit_status,Minimazation] = Force_poly_phi...
    (Cell.down_lipid_director_phi_matrix,Cell.constraint.phi_lipid_angle_down,res_struc.phi_res,res_struc.poly_degree_LD_phi,Minimazation); 
if exit_status==0
    exit_status_membrane=0;
end


%% physical properties
Virus.physical_properties.lipid_length=lipid_length; %[nm]
Virus.physical_properties.kappa_up=General_physical_properties.virus_distal_kappa; 
Virus.physical_properties.kappa_down=General_physical_properties.proximal_kappa;
Virus.physical_properties.kappa_bar_up=Virus.physical_properties.kappa_up*chi; %kappa units
Virus.physical_properties.kappa_bar_down=Virus.physical_properties.kappa_down*chi; %kappa units
Virus.physical_properties.Ks=streatching_modulus; %kappa/nm^2
Virus.physical_properties.kappa_tilt=tilt_modulus; %kappa/nm^2
Virus.physical_properties.J0_up=General_physical_properties.Virus_distal_J0;
Virus.physical_properties.J0_down=General_physical_properties.proximal_J0;
Virus.physical_properties.Splay_grad_square_modulus=General_physical_properties.Splay_grad_square_modulus;
Virus.physical_properties.Splay_grad_tilt_modulus=General_physical_properties.Splay_grad_tilt_modulus;
Virus.physical_properties.Surface_tension=General_physical_properties.Surface_tension;


Cell.physical_properties.lipid_length=lipid_length; %[nm]
Cell.physical_properties.kappa_up=General_physical_properties.proximal_kappa; 
Cell.physical_properties.kappa_down=General_physical_properties.cell_distal_kappa; 
Cell.physical_properties.kappa_bar_up=Cell.physical_properties.kappa_up*chi; %kappa units
Cell.physical_properties.kappa_bar_down=Cell.physical_properties.kappa_down*chi; %kappa unitscell_physical_properties.Ks=1; %kappa/nm^2
Cell.physical_properties.Ks=streatching_modulus; %kappa/nm^2
Cell.physical_properties.kappa_tilt=tilt_modulus; %kappa/nm^2
Cell.physical_properties.J0_up=General_physical_properties.proximal_J0;
Cell.physical_properties.J0_down=General_physical_properties.Cell_distal_J0;
Cell.physical_properties.Splay_grad_square_modulus=General_physical_properties.Splay_grad_square_modulus;
Cell.physical_properties.Splay_grad_tilt_modulus=General_physical_properties.Splay_grad_tilt_modulus;
Cell.physical_properties.Surface_tension=General_physical_properties.Surface_tension;


Diaphragm.physical_properties.lipid_length=lipid_length; %[nm]
Diaphragm.physical_properties.kappa_up=Virus.physical_properties.kappa_up; %[nm]
Diaphragm.physical_properties.kappa_down=Cell.physical_properties.kappa_down; %[nm]
Diaphragm.physical_properties.kappa_bar_up=Diaphragm.physical_properties.kappa_up*chi; %kappa units
Diaphragm.physical_properties.kappa_bar_down=Diaphragm.physical_properties.kappa_down*chi; %kappa unitsVirus_physical_properties.kappa_bar_down=kappa_down*chi; %kappa unitsDiaphragm_physical_properties.Ks=1; %kappa/nm^2
Diaphragm.physical_properties.Ks=streatching_modulus; %kappa/nm^2
Diaphragm.physical_properties.kappa_tilt=tilt_modulus; %kappa/nm^2
Diaphragm.physical_properties.J0_up=Virus.physical_properties.J0_up;
Diaphragm.physical_properties.J0_down=Cell.physical_properties.J0_down;
Diaphragm.physical_properties.Splay_grad_square_modulus=General_physical_properties.Splay_grad_square_modulus;
Diaphragm.physical_properties.Splay_grad_tilt_modulus=General_physical_properties.Splay_grad_tilt_modulus;
Diaphragm.physical_properties.Surface_tension=General_physical_properties.Surface_tension;




%% save shape



%cell
[Cell.mid_plane,Cell.lipid_director]=Cell_shape(Diaphragm,Cell,res_struc);

[Cell.dividing_plane,Cell.tilt,Cell.lipid_length,Cell.area_element,Cell.splay,Cell.saddle_splay]=...
    find_membrane_geometry(Cell.mid_plane,Cell.lipid_director,Cell.physical_properties,res_struc.phi_res,res_struc.rho_res);

%Virus
[Virus.mid_plane,Virus.lipid_director]=Virus_shape(Diaphragm,Virus,res_struc);

[Virus.dividing_plane,Virus.tilt,Virus.lipid_length,Virus.area_element,Virus.splay,Virus.saddle_splay]=...
    find_membrane_geometry(Virus.mid_plane,Virus.lipid_director,Virus.physical_properties,res_struc.phi_res,res_struc.rho_res);

%Diaphragm
[Diaphragm.mid_plane,Diaphragm.lipid_director]=diaphragm_shape(Diaphragm,res_struc);

[Diaphragm.dividing_plane,Diaphragm.tilt,Diaphragm.lipid_length,Diaphragm.area_element,Diaphragm.splay,Diaphragm.saddle_splay]=...
    find_membrane_geometry(Diaphragm.mid_plane,Diaphragm.lipid_director,Diaphragm.physical_properties,res_struc.phi_res,res_struc.rho_res);


%% Shell
if Shell.Shell_exist==1
 [Shell.Diaphragm,Shell.Cell,Shell.Junction,exit_status_shell,Minimazation] = build_shell(Shell.Diaphragm,Shell.Cell,Diaphragm,Cell,res_struc,Minimazation);
end

%% cell trans-membrane proprties
if Minimazation.Cell_trans_membrane_added_proprties_exist==1
    Cell.mid_plane.mid_plane_physical_proprties.kappa=General_physical_properties.Cell_mid_plane_physical_proprties.kappa;
    Cell.mid_plane.mid_plane_physical_proprties.kappa_bar=General_physical_properties.Cell_mid_plane_physical_proprties.kappa_bar;
    Cell.mid_plane.mid_plane_physical_proprties.J0=General_physical_properties.Cell_mid_plane_physical_proprties.J0;
    
    [Cell.mid_plane.Total_curvature,Cell.mid_plane.Gaussian_curvature]=mid_plane_geomtry(Cell.mid_plane,res_struc);
end

% sorting between diaphragm and distal cell membrane
if strcmp(Minimazation.sorting_protein_in_cell_membrane,'splay_based') || strcmp(Minimazation.sorting_protein_in_cell_membrane,'affinity_based')

    j_s=General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_s;
    j_0=General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_0;
    j_TM_cytosolic=General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_TM;
    j_TM_luminal=General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_TM_luminal;
    ratio_s_to_0=General_physical_properties.Cell_mid_plane_physical_proprties.sorting.ratio_s_to_0;
    phi_TM=General_physical_properties.Cell_mid_plane_physical_proprties.sorting.phi_TM;
    
    phi_0_R=(1-phi_TM)/(1+ratio_s_to_0);
    phi_s_R=(1-phi_TM)*ratio_s_to_0/(1+ratio_s_to_0);

    J_0_cell_cytosolic=j_s*phi_s_R+phi_0_R*j_0+phi_TM*j_TM_cytosolic;
    J_0_cell_luminal=j_s*phi_s_R+phi_0_R*j_0+phi_TM*j_TM_luminal;
    J0_R_no_TM=(j_0+j_s*ratio_s_to_0)/(1+ratio_s_to_0);

    % set all curvature based on lipid compostion
    Diaphragm.physical_properties.J0_up=J0_R_no_TM;
    Diaphragm.physical_properties.J0_down=J0_R_no_TM;

    Cell.physical_properties.J0_up=J_0_cell_luminal;
    Cell.physical_properties.J0_down=J_0_cell_cytosolic;

    Virus.physical_properties.J0_up=J0_R_no_TM;
    Virus.physical_properties.J0_down=J0_R_no_TM;
   
end
 
if strcmp(Minimazation.sorting_protein_in_cell_membrane,'splay_based') && Minimazation.no_sorting==0
    kappa_m_kbT_a=General_physical_properties.Cell_mid_plane_physical_proprties.sorting.kappa_m_kbT_a;
    % implemet the new curvature to the directly contacting monolayers
    [~,Virus.physical_properties.J0_down,~] = Diaphargm_J_m_trans_membrane_out_2componenets(kappa_m_kbT_a,j_s,j_0,j_TM_luminal,phi_s_R,phi_0_R,phi_TM,Minimazation.MAX_ITR);    
    [~,Diaphragm.physical_properties.J0_down,~] = Diaphargm_J_m_trans_membrane_out_2componenets(kappa_m_kbT_a,j_s,j_0,j_TM_cytosolic,phi_s_R,phi_0_R,phi_TM,Minimazation.MAX_ITR);
    %if there is flipflop virus inner as cytsolic side!
    if Minimazation.full_lipid_flip_flop==1
        Diaphragm.physical_properties.J0_up=Diaphragm.physical_properties.J0_down;
        Virus.physical_properties.J0_up=Diaphragm.physical_properties.J0_down;
    end


end

if strcmp(Minimazation.sorting_protein_in_cell_membrane,'affinity_based') && Minimazation.no_sorting==0

    Affinity=General_physical_properties.Cell_mid_plane_physical_proprties.sorting.affinity;
    bound_lipids=General_physical_properties.Cell_mid_plane_physical_proprties.sorting.number_of_intercating_lipids;
    phi_chol_reservior=ratio_s_to_0/(1+ratio_s_to_0);

    [phi_s_free] = Special_lipid_mole_fraction_at_IFITM3_free_part(Affinity,phi_TM,bound_lipids,phi_chol_reservior);

    J0_virus_diaphragm=phi_s_free*j_s+(1-phi_s_free)*j_0;
    Diaphragm.physical_properties.J0_up=J0_virus_diaphragm;
    Diaphragm.physical_properties.J0_down=J0_virus_diaphragm;
    Virus.physical_properties.J0_up=J0_virus_diaphragm;
    Virus.physical_properties.J0_down=J0_virus_diaphragm;

    if Minimazation.no_TM_direct_intrinsic_curvature==1
        Cell.physical_properties.J0_up=J0_R_no_TM;
        Cell.physical_properties.J0_down=J0_R_no_TM;
    end

end


if Minimazation.relaxed_spherical_VLP==1
    J_sm=General_physical_properties.proximal_J0;
    Cell.physical_properties.J0_up=General_physical_properties.proximal_J0+(1-lipid_length*J_sm)/Cell.R_curv;
    Cell.physical_properties.J0_down=General_physical_properties.Cell_distal_J0-(1-lipid_length*J_sm)/Cell.R_curv;
end

if Minimazation.relaxed_cylindrical_VLP==1
    rbv=Virus.r_bv;
    rsv_proxy=Virus.r_sv+Virus.physical_properties.lipid_length;
    rsv_distal=Virus.r_sv-Virus.physical_properties.lipid_length;
    J2A_proximal=(rbv^4/(rbv^2*rsv_proxy^2-rsv_proxy^4))^0.5;
    J2A_distal=(rbv^4/(rbv^2*rsv_distal^2-rsv_distal^4))^0.5;
    Delta_J_cylinder=(J2A_proximal+J2A_distal)/(2*rbv);
    Virus.physical_properties.J0_down=General_physical_properties.proximal_J0+Delta_J_cylinder;

    Virus.physical_properties.J0_up=General_physical_properties.Virus_distal_J0-Delta_J_cylinder;
end



if exit_status_membrane==0 || exit_status_shell==0 % if error 
    exit_status_global=0;
end


end

