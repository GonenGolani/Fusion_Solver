function [total_Volume] = find_Cell_volume(Cell,Diaphragm,res_struc)
%find cell area out side of Fusion site
volume_bulk=2*pi/3*(1+(1-(Cell.R_rim/Cell.R_curv)^2)^0.5)*Cell.R_curv^3;

% find area between z=0 to diaphragm
dphi=pi/2/res_struc.phi_res;
volume_above_diaphragm=4*sum(Diaphragm.mid_plane.rho.*Diaphragm.mid_plane.drho.*Diaphragm.mid_plane.z*dphi,'all');
volume_above_cell=4*sum(Cell.mid_plane.rho.*Cell.mid_plane.drho.*Cell.mid_plane.z*dphi,'all');

total_Volume=volume_bulk+volume_above_diaphragm+volume_above_cell;

end

