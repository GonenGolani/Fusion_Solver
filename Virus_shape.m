function  [mid_plane,lipid_director]=Virus_shape(Diaphragm,Virus,res_struct)
    
    
inner_rim.x0=Diaphragm.rim.x0;
inner_rim.y0=Diaphragm.rim.y0;
outer_rim.x0=Virus.Rv;
outer_rim.y0=Virus.Rv;

Virus.Gaussian_decay_vector_up.x0=Diaphragm.rim.x0;
Virus.Gaussian_decay_vector_up.y0=Diaphragm.rim.y0;
Virus.Gaussian_decay_vector_down.x0=Diaphragm.rim.x0;
Virus.Gaussian_decay_vector_down.y0=Diaphragm.rim.y0;

[mid_plane.x,mid_plane.y,mid_plane.z] = mid_plane_surface_creator(Virus.mid_plane_coeff_matrix,inner_rim,outer_rim,res_struct);

[mid_plane.normal.x,mid_plane.normal.y,mid_plane.normal.z]=surfnorm(mid_plane.x,mid_plane.y,mid_plane.z);
RHO=(mid_plane.x.^2+mid_plane.y.^2).^0.5;
sinphi=mid_plane.y./RHO;
cosphi=mid_plane.x./RHO;
mid_plane.normal.rho=cosphi.*mid_plane.normal.x+sinphi.*mid_plane.normal.y;
mid_plane.normal.phi=-sinphi.*mid_plane.normal.x+cosphi.*mid_plane.normal.y;

%lipid director in up direction
[lipid_director.up.x,lipid_director.up.y,lipid_director.up.z,lipid_director.up.rho,lipid_director.up.phi] =lipid_director_creator...
    (Virus.up_lipid_director_rho_matrix,Virus.up_lipid_director_phi_matrix,Virus.Gaussian_decay_vector_up,mid_plane,inner_rim,outer_rim,res_struct,1);
%lipid director in down direction
[lipid_director.down.x,lipid_director.down.y,lipid_director.down.z,lipid_director.up.rho,lipid_director.up.phi] =lipid_director_creator...
    (Virus.down_lipid_director_rho_matrix,Virus.down_lipid_director_phi_matrix,Virus.Gaussian_decay_vector_down,mid_plane,inner_rim,outer_rim,res_struct,-1);


end