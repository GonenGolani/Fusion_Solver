function [rho_vector,z_vector,dz_drho] = virus_rim_contstruct(Virus,phi_res)
%calcaulte the rim proprties of the virus 
%rim radial posotion: rho_vector
% hight: z_vector
% tangent in ib rho direction: drho_dz
% lipid angle in rho direction theta_rho
% lipid angle in phi direction theta_phi
    

% Fusion site rim in the rho phi system
phi=linspace(0,pi/2,phi_res);
z_v_center_vec=Virus.z_v_center*ones(1,phi_res);
r_bv_vec=Virus.r_bv*ones(1,phi_res);
r_sv_vec=Virus.r_sv*ones(1,phi_res);


% from the parametric description
% helping functions
sinV=Virus.Rv*sin(phi)/Virus.r_sv;
cosV=(1-sinV.^2).^0.5; % can be +- for inner or outer intersection;
y_star=(r_bv_vec+r_sv_vec.*cosV); %just to short things
cosU=Virus.Rv*cos(phi)./y_star;
sinU=(1-cosU.^2).^0.5;


% the contact rim between Fusion site and virus
xt=(Virus.r_bv*cosU+Virus.r_sv*cosU.*cosV);
yt=Virus.r_sv*sinV;
z_vector =z_v_center_vec-(Virus.r_bv*sinU+Virus.r_sv*sinU.*cosV);  % can be +- for lower or upper;
rho_vector=(xt.^2+yt.^2).^0.5;



%tangent
Rb=Virus.r_bv;
Rs=Virus.r_sv;
Rr=Virus.Rv;

%tangent at edge
dsinVdrho=sin(phi)/Rs;
dcosVdrho=-sinV./cosV.*dsinVdrho;
dcosUdrho=(y_star-Rr*Rs.*dcosVdrho).*cos(phi)./y_star.^2;
dsinUdrho=-cosU./sinU.*dcosUdrho;

dz_drho=-(dsinUdrho.*y_star+dcosVdrho.*Rs.*sinU);


end