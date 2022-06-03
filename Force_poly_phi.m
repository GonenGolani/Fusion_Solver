function [out_matrix_coeff,exit_status,Minimazation] = Force_poly_phi(coeff_matrix,constraint,res,poly_degree,Minimazation)
% This function recives the coefficient matrix (columms - hm coeff and rows
% the phi angle)
%%%%%%%
%  constraint is structure with 4 vectors: rhoi, rhof, tanalphai and tanalphaf 
% the cells in each vector corresponds to phi 
%%%%%%%
% the output is the same matrix as coeff_matrix but with the first columms
% replaced so the boundary condtions are forced.
% exit_status 1 if successful 0 if error 

out_matrix_coeff=zeros(res,poly_degree);
int=1;
while int<=res
    zi=constraint.zi(int);
    zf=constraint.zf(int);
    dzdrhoi=constraint.dzdrhoi(int);
    dzdrhof=constraint.dzdrhof(int);
    rhoi=constraint.rhoi(int);
    rhof=constraint.rhof(int);
    [out_matrix_coeff(int,:),exit_status,Minimazation]= poly_force_boundary(zi,zf,dzdrhoi,dzdrhof,rhoi,rhof,coeff_matrix(int,:),Minimazation);
    if exit_status==0
        return;
    end
    int=int+1;
end


   


end

