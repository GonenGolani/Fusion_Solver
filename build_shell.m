function [Shell_Diaphragm,Shell_Cell,Shell_Junction,exit_status_global,Minimazation] =...
    build_shell(Shell_Diaphragm,Shell_Cell,Diaphragm,Cell,res_struc,Minimazation)

exit_status_global=1;
[cosine_matrix] = cosine_matrix_creatore(res_struc.phi_res,res_struc.poly_degree_rim); 
Shell_Junction.mid_plane_coeff_matrix=zeros(res_struc.phi_res,4); 
%calcualte the shape of the Diaphragm rim
[rho_vector_Diaphragm_rim,~,~] = diaphragm_rim_contstruct(Diaphragm.rim,cosine_matrix,res_struc.phi_res);

%[rho_vector_cell,~,~,theta_rho_cell,~] = cell_rim_contstruct(Cell,res_struc.phi_res);


%boundary conditions at the rim


%% make sure that the postion of the edges are ok

%conditions
if min(rho_vector_Diaphragm_rim)-abs(Shell_Diaphragm.length_from_edge)-2*Diaphragm.r_pore<=0 % if  edge reach the middele of the diaphragm  
    Shell_Diaphragm.length_from_edge=0;
end

if abs(Shell_Diaphragm.length_from_edge)<=0 % if diaphargm edge reach edge of the diaphragm
    Shell_Diaphragm.length_from_edge=0;
end

if abs(Shell_Cell.length_from_edge)>=Cell.R_rim-max(rho_vector_Diaphragm_rim)-Diaphragm.r_pore % if cell edge reach edge of fusion sit
    Shell_Cell.length_from_edge=Cell.R_rim-max(rho_vector_Diaphragm_rim)-Diaphragm.r_pore;
end

if abs(Shell_Cell.length_from_edge)<=0 % if cell edge reach edge of the diaphragm
    Shell_Cell.length_from_edge=Diaphragm.r_pore;
end

%% create Diaphragm mid plane
% constraint Diaphragm 
Shell_Diaphragm.constraint.zi=ones(1,res_struc.phi_res)*Shell_Diaphragm.z_i;     % set hight of center diaphragm
Shell_Diaphragm.constraint.zf=ones(1,res_struc.phi_res)*Shell_Diaphragm.z_f;     % hight at edge of Diaphragm
Shell_Diaphragm.constraint.dzdrhoi=zeros(1,res_struc.phi_res);                            %set zero angle at rho=0
Shell_Diaphragm.constraint.dzdrhof=ones(1,res_struc.phi_res)*Shell_Diaphragm.dzdrho_edge; %dervitive of hight at edge
Shell_Diaphragm.constraint.rhoi=ones(1,res_struc.phi_res)*Diaphragm.r_pore; %center of Diaphragm is rho=Diaphragm.r_pore
Shell_Diaphragm.constraint.rhof=rho_vector_Diaphragm_rim-abs(Shell_Diaphragm.length_from_edge); %find position of Diaphragm rho

%coeffciant matrix
% force boundary conditions
[Shell_Diaphragm.MP_to_MP_distance_matrix,exit_status,Minimazation] = ...
    Force_poly_phi(Shell_Diaphragm.MP_to_MP_distance_matrix,Shell_Diaphragm.constraint,res_struc.phi_res,res_struc.poly_degree_shell,Minimazation); 
if exit_status==0
    exit_status_global=0;
end

% create diaphragm shell based on the diphargm mid plane
Diphargm_inner_rim.x0=Diaphragm.r_pore;
Diphargm_inner_rim.y0=Diaphragm.r_pore;
Diphargm_outer_rim.x0=Diaphragm.rim.x0-Shell_Diaphragm.length_from_edge;
Diphargm_outer_rim.y0=Diaphragm.rim.x0-Shell_Diaphragm.length_from_edge;


%create the mid to mid plane based on the matrix
[~,~,Shell_Diaphragm.distance_from_membrane]...
    = mid_plane_surface_creator(Shell_Diaphragm.MP_to_MP_distance_matrix,Diphargm_inner_rim,Diphargm_outer_rim,res_struc);

%create the diaphragm mid plane to be used for bulding the shell
[diaphragm_mid_plane_for_shell.x,diaphragm_mid_plane_for_shell.y,diaphragm_mid_plane_for_shell.z]...
    = mid_plane_surface_creator(Diaphragm.mid_plane_coeff_matrix,Diphargm_inner_rim,Diphargm_outer_rim,res_struc);

%find normal to shell 
[Shell_Diaphragm.normal.x,Shell_Diaphragm.normal.y,Shell_Diaphragm.normal.z]...
    =surfnorm(diaphragm_mid_plane_for_shell.x,diaphragm_mid_plane_for_shell.y,diaphragm_mid_plane_for_shell.z);

%flip normal if needed needs to point down!
if Shell_Diaphragm.normal.z(1,1)>0
    Shell_Diaphragm.normal.x=Shell_Diaphragm.normal.x*-1;
    Shell_Diaphragm.normal.y=Shell_Diaphragm.normal.y*-1;
    Shell_Diaphragm.normal.z=Shell_Diaphragm.normal.z*-1;   
end

% Position of the Shell mid plane 
Shell_Diaphragm.mid_plane.x=diaphragm_mid_plane_for_shell.x+Shell_Diaphragm.normal.x.*Shell_Diaphragm.distance_from_membrane; 
Shell_Diaphragm.mid_plane.y=diaphragm_mid_plane_for_shell.y+Shell_Diaphragm.normal.y.*Shell_Diaphragm.distance_from_membrane; 
Shell_Diaphragm.mid_plane.z=diaphragm_mid_plane_for_shell.z+Shell_Diaphragm.normal.z.*Shell_Diaphragm.distance_from_membrane; 


% geometry of the shell
[Shell_Diaphragm.normal,Shell_Diaphragm.total_curvature,Shell_Diaphragm.Gaussian_curvature,Shell_Diaphragm.mid_plane]...
    =find_shell_geometry(Shell_Diaphragm.mid_plane,res_struc.phi_res,res_struc.rho_res);

%% create mid plane of cell membrane

Shell_Cell.constraint.zi=ones(1,res_struc.phi_res)*Shell_Cell.z_i; %inner boundary is Diaphragm rim
Shell_Cell.constraint.zf=ones(1,res_struc.phi_res)*Shell_Cell.z_f;
Shell_Cell.constraint.dzdrhoi=ones(1,res_struc.phi_res)*Shell_Cell.dzdrho_edge;
Shell_Cell.constraint.dzdrhof=ones(1,res_struc.phi_res)*0;
Shell_Cell.constraint.rhoi=rho_vector_Diaphragm_rim+Shell_Cell.length_from_edge;
Shell_Cell.constraint.rhof=Cell.R_rim*ones(1,res_struc.phi_res);


[Shell_Cell.MP_to_MP_distance_matrix,exit_status,Minimazation] = Force_poly_phi...
    (Shell_Cell.MP_to_MP_distance_matrix,Shell_Cell.constraint,res_struc.phi_res,res_struc.poly_degree_shell,Minimazation);
if exit_status==0
    exit_status_global=0;
end

Cell_inner_rim.x0=Diaphragm.rim.x0+Shell_Cell.length_from_edge;
Cell_inner_rim.y0=Diaphragm.rim.y0+Shell_Cell.length_from_edge;
Cell_outer_rim.x0=Cell.R_rim;
Cell_outer_rim.y0=Cell.R_rim;

%create the mid to mid plane based on the matrix
[~,~,Shell_Cell.distance_from_membrane]...
    = mid_plane_surface_creator(Shell_Cell.MP_to_MP_distance_matrix,Cell_inner_rim,Cell_outer_rim,res_struc);


% create the mid plane for the shell construct
[Cell_mid_plane_for_shell.x,Cell_mid_plane_for_shell.y,Cell_mid_plane_for_shell.z]...
    = mid_plane_surface_creator(Cell.mid_plane_coeff_matrix,Cell_inner_rim,Cell_outer_rim,res_struc);

%find normal to shell 
[Shell_Cell.normal.x,Shell_Cell.normal.y,Shell_Cell.normal.z]...
    =surfnorm(Cell_mid_plane_for_shell.x,Cell_mid_plane_for_shell.y,Cell_mid_plane_for_shell.z);

%flip normal if needed needs to point down!
if Shell_Cell.normal.z(1,1)>0
    Shell_Cell.normal.x=Shell_Cell.normal.x*-1;
    Shell_Cell.normal.y=Shell_Cell.normal.y*-1;
    Shell_Cell.normal.z=Shell_Cell.normal.z*-1;   
end

% Position of the Shell mid plane 
Shell_Cell.mid_plane.x=Cell_mid_plane_for_shell.x+Shell_Cell.normal.x.*Shell_Cell.distance_from_membrane; 
Shell_Cell.mid_plane.y=Cell_mid_plane_for_shell.y+Shell_Cell.normal.y.*Shell_Cell.distance_from_membrane; 
Shell_Cell.mid_plane.z=Cell_mid_plane_for_shell.z+Shell_Cell.normal.z.*Shell_Cell.distance_from_membrane; 


% geometry of the shell
[Shell_Cell.normal,Shell_Cell.total_curvature,Shell_Cell.Gaussian_curvature,Shell_Cell.mid_plane]=find_shell_geometry...
        (Shell_Cell.mid_plane,res_struc.phi_res,res_struc.rho_res);


%% Shell shape at junction area

% find derivitive at edge of shell on diaphragm side
Diaphragm_normal_rho_sign=(Shell_Diaphragm.normal.x(:,res_struc.rho_res)>=0 & Shell_Diaphragm.normal.y(:,res_struc.rho_res)>=0)*2-1;
Diaphragm_normal_rho_edge=Diaphragm_normal_rho_sign.*(Shell_Diaphragm.normal.x(:,res_struc.rho_res).^2+Shell_Diaphragm.normal.y(:,res_struc.rho_res).^2).^0.5;

Shell_junction_inner_rim=(Shell_Diaphragm.mid_plane.x(:,res_struc.rho_res).^2+Shell_Diaphragm.mid_plane.y(:,res_struc.rho_res).^2).^0.5;

% find derivitive at edge of shell on cell side
Cell_normal_rho_sign=(Shell_Cell.normal.x(:,1)<=0 & Shell_Cell.normal.y(:,1)<=0)*2-1;
Cell_normal_rho_edge=Cell_normal_rho_sign.*(Shell_Cell.normal.x(:,1).^2+Shell_Cell.normal.y(:,1).^2).^0.5;

Shell_junction_outer_rim=(Shell_Cell.mid_plane.x(:,1).^2+Shell_Cell.mid_plane.y(:,1).^2).^0.5;

% build the position matrix of the junction area
Shell_Junction.constraint.zi=Shell_Diaphragm.mid_plane.z(:,res_struc.rho_res); 
Shell_Junction.constraint.zf=Shell_Cell.mid_plane.z(:,1);
Shell_Junction.constraint.dzdrhoi=-Diaphragm_normal_rho_edge./Shell_Diaphragm.normal.z(:,res_struc.rho_res);
Shell_Junction.constraint.dzdrhof=Cell_normal_rho_edge./Shell_Cell.normal.z(:,1);
Shell_Junction.constraint.rhoi=Shell_junction_inner_rim;
Shell_Junction.constraint.rhof=Shell_junction_outer_rim;

[Shell_Junction.mid_plane_coeff_matrix,exit_status,Minimazation] = Force_poly_phi...
    (Shell_Junction.mid_plane_coeff_matrix,Shell_Junction.constraint,res_struc.phi_res,4,Minimazation);
if exit_status==0
    exit_status_global=0;
end

%build shell mid-plane
[Shell_Junction.mid_plane.x,Shell_Junction.mid_plane.y,Shell_Junction.mid_plane.z]...
    = mid_plane_surface_creator_genral_rim(Shell_Junction.mid_plane_coeff_matrix,Shell_junction_inner_rim,Shell_junction_outer_rim,res_struc.phi_res,res_struc.rho_res_Junction);

[Shell_Junction.normal,Shell_Junction.total_curvature,Shell_Junction.Gaussian_curvature,Shell_Junction.mid_plane]=find_shell_geometry...
        (Shell_Junction.mid_plane,res_struc.phi_res,res_struc.rho_res_Junction);


end

