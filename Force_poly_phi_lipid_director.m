function [out_matrix_coeff] = Force_poly_phi_lipid_director(coeff_matrix,constraint,phi_res,poly_res)
% This function recives the coefficient matrix (columms - hm coeff and rows
% the phi angle)
%%%%%%%
%  constraint is structure with 4 vectors: rhoi, rhof, tanalphai and tanalphaf 
% the cells in each vector corresponds to phi 
%%%%%%%
% the output is the same matrix as coeff_matrix but with the first columms
% replaced so the boundary condtions are forced.

out_matrix_coeff=zeros(phi_res,poly_res);


int=1;
while int<=phi_res
    angle_i=constraint.angle_i(int);
    angle_f=constraint.angle_f(int);
    rhoi=constraint.rhoi(int);
    rhof=constraint.rhof(int);
    out_matrix_coeff(int,:)= poly_force_boundary_no_div(angle_i,angle_f,rhoi,rhof,coeff_matrix(int,:));
    int=int+1;
end


   


end

