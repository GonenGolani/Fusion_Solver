function  [mid_plane,lipid_director]=diaphragm_shape(Diaphragm,res_struct)
%this function is a helper to diaphragm_energy and plot_diaphragm
%it calcualtes the diaphragm mid-plane and save it as a strcture of x y and z
% Normal vector to mid-plane N 
% up and down lipid director n

    
%mid plane
inner_rim.x0=Diaphragm.r_pore;
inner_rim.y0=Diaphragm.r_pore;

Diaphragm.Gaussian_decay_vector_up.x0=Diaphragm.rim.x0;
Diaphragm.Gaussian_decay_vector_up.y0=Diaphragm.rim.y0;
Diaphragm.Gaussian_decay_vector_down.x0=Diaphragm.rim.x0;
Diaphragm.Gaussian_decay_vector_down.y0=Diaphragm.rim.y0;


[mid_plane.x,mid_plane.y,mid_plane.z] = mid_plane_surface_creator(Diaphragm.mid_plane_coeff_matrix,inner_rim,Diaphragm.rim,res_struct);
[mid_plane.normal.x,mid_plane.normal.y,mid_plane.normal.z]=surfnorm(mid_plane.x,mid_plane.y,mid_plane.z);
RHO=(mid_plane.x.^2+mid_plane.y.^2).^0.5;
sinphi=mid_plane.y./RHO;
cosphi=mid_plane.x./RHO;
mid_plane.normal.rho=cosphi.*mid_plane.normal.x+sinphi.*mid_plane.normal.y;
mid_plane.normal.phi=-sinphi.*mid_plane.normal.x+cosphi.*mid_plane.normal.y;

%lipid director in up direction
[lipid_director.up.x,lipid_director.up.y,lipid_director.up.z,lipid_director.up.rho,lipid_director.up.phi] =lipid_director_creator...
    (Diaphragm.up_lipid_director_rho_matrix,Diaphragm.up_lipid_director_phi_matrix,Diaphragm.Gaussian_decay_vector_up,mid_plane,inner_rim,Diaphragm.rim,res_struct,1);
%lipid director in down direction
[lipid_director.down.x,lipid_director.down.y,lipid_director.down.z,lipid_director.up.rho,lipid_director.up.phi] =lipid_director_creator...
    (Diaphragm.down_lipid_director_rho_matrix,Diaphragm.down_lipid_director_phi_matrix,Diaphragm.Gaussian_decay_vector_down,mid_plane,inner_rim,Diaphragm.rim,res_struct,-1);






end
