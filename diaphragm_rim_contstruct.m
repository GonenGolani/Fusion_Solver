function [rho_vector,z_vector,drho_dz] = diaphragm_rim_contstruct(Diphragm_rim,cosine_matrix,phi_res)
% this function recives the diaphragm rim parametres and phi resolution 
% the cosine_matrix is a matrix containing the expenston of cos(2*l*phi)
% for different phi so the formation of the vec
% Diphragm_rim containts x0 and y0 for the rho vector
% Diphragm_rim.hl_vector and Diphragm_rim.drho_dz_vector 
% containes the plynomial coefficciants of the hight and derivitive at the edge of the rim 

% output is: rho_vector and z_vector as fuction of phi
% input: stracture called Diphragm_rim containes:
  
phi_range=linspace(0,pi/2,phi_res);
rho_vector=(Diphragm_rim.x0.^2.*cos(phi_range).^2+Diphragm_rim.y0.^2.*sin(phi_range).^2).^0.5;

z_vector=Diphragm_rim.hl_vector*cosine_matrix;

drho_dz=Diphragm_rim.drho_dz_vector*cosine_matrix;


end

