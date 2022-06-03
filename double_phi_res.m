function [New_Diaphragm,New_Cell,New_Virus,New_res_struc] = double_phi_res(old_Diaphragm,old_Cell,old_Virus,old_res_struc)

New_Virus=old_Virus;
New_Cell=old_Cell;
New_Diaphragm=old_Diaphragm;
New_res_struc=old_res_struc;

New_res_struc.phi_res=old_res_struc.phi_res*2;

[New_Virus.mid_plane_coeff_matrix] = double_array(old_Virus.mid_plane_coeff_matrix);
[New_Virus.up_lipid_director_rho_matrix] = double_array(old_Virus.up_lipid_director_rho_matrix);
[New_Virus.down_lipid_director_rho_matrix] = double_array(old_Virus.down_lipid_director_rho_matrix);
[New_Virus.up_lipid_director_phi_matrix] = double_array(old_Virus.up_lipid_director_phi_matrix);
[New_Virus.down_lipid_director_phi_matrix] = double_array(old_Virus.down_lipid_director_phi_matrix);


[New_Diaphragm.mid_plane_coeff_matrix] = double_array(old_Diaphragm.mid_plane_coeff_matrix);
[New_Diaphragm.up_lipid_director_rho_matrix] = double_array(old_Diaphragm.up_lipid_director_rho_matrix);
[New_Diaphragm.down_lipid_director_rho_matrix] = double_array(old_Diaphragm.down_lipid_director_rho_matrix);
[New_Diaphragm.up_lipid_director_phi_matrix] = double_array(old_Diaphragm.up_lipid_director_phi_matrix);
[New_Diaphragm.down_lipid_director_phi_matrix] = double_array(old_Diaphragm.down_lipid_director_phi_matrix);


[New_Cell.mid_plane_coeff_matrix] = double_array(old_Cell.mid_plane_coeff_matrix);
[New_Cell.up_lipid_director_rho_matrix] = double_array(old_Cell.up_lipid_director_rho_matrix);
[New_Cell.down_lipid_director_rho_matrix] = double_array(old_Cell.down_lipid_director_rho_matrix);
[New_Cell.up_lipid_director_phi_matrix] = double_array(old_Cell.up_lipid_director_phi_matrix);
[New_Cell.down_lipid_director_phi_matrix] = double_array(old_Cell.down_lipid_director_phi_matrix);


end
