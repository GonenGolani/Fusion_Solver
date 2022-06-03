function [Cell,Diaphragm] = recalc_only_cell_mid_plane(Cell_old,Diaphragm_old,Minimazation,res_struc)

Cell=Cell_old;
Diaphragm=Diaphragm_old;


% matrix containg cos(2*int*phi)  col is int (starting from 0) rows is phi 0->pi/2
[cosine_matrix] = cosine_matrix_creatore(res_struc.phi_res,res_struc.poly_degree_rim); 

%calcualte the shape of the Diaphragm rim
[rho_vector_Diaphragm_rim,z_vector_Diaphragm_rim,~] = diaphragm_rim_contstruct...
    (Diaphragm.rim,cosine_matrix,res_struc.phi_res);


% create mid plane of cell membrane
[rho_vector_cell,z_vector_cell,drho_dz_cell,~,~] = cell_rim_contstruct(Cell,res_struc.phi_res);

Cell.constraint.mid_plane.zi=z_vector_Diaphragm_rim; %inner boundary is Diaphragm rim
Cell.constraint.mid_plane.zf=z_vector_cell;
Cell.constraint.mid_plane.dzdrhoi=Cell.rim.drho_dz_vector*cosine_matrix;
Cell.constraint.mid_plane.dzdrhof=drho_dz_cell;
Cell.constraint.mid_plane.rhoi=rho_vector_Diaphragm_rim;
Cell.constraint.mid_plane.rhof=rho_vector_cell;


[Cell.mid_plane_coeff_matrix,~,~] = Force_poly_phi...
    (Cell.mid_plane_coeff_matrix,Cell.constraint.mid_plane,res_struc.phi_res,res_struc.poly_degree_z,Minimazation);


[Cell.mid_plane,~]=Cell_shape(Diaphragm,Cell,res_struc);

Cell.mid_plane.rho=(Cell.mid_plane.x.^2+Cell.mid_plane.y.^2).^0.5;
[Cell.mid_plane.drho,~]=gradient(Cell.mid_plane.rho);
Diaphragm.mid_plane.rho=(Diaphragm.mid_plane.x.^2+Diaphragm.mid_plane.y.^2).^0.5;
[Diaphragm.mid_plane.drho,~]=gradient(Diaphragm.mid_plane.rho);


end

