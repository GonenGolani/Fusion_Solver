function [a_n_out] = poly_force_boundary_no_div(angle_i,angle_f,rhoi,rhof,a_n_in)

% Input= vector a_n_in with ploynum coeffients with length N.
% first 4 does not matter becuse found from boundary conditions
% boundary condtions:
% outer edge radius: Rad
% angle at rho=rhof: Ang_0 at 0 in deg
% angle at=rho=Rad: Ang_R at 0 in deg
% hight at rho=rhof: z0 
% hight at rho=Rad: z1 
% output= assign to 4 first places values using. give vector a_n_out 
% solved by inversing the matrix



A_mat=[1 rhoi ; 1 rhof ];


arr_size=size(a_n_in);
an_length=arr_size(2);

thorw_ot_first_2_vector=ones(1,an_length);
thorw_ot_first_2_vector(1:2)=0;
% RHS side solution vector
RHS_vector(:,1)=angle_i-mypolyval(a_n_in.*thorw_ot_first_2_vector,rhoi);
RHS_vector(:,2)=angle_f-mypolyval(a_n_in.*thorw_ot_first_2_vector,rhof);

%out an vector 1-2
A_inv=inv(A_mat);
a_1_to_2=A_inv*RHS_vector';

%a_1_to_4=A_mat/RHS_vector'; %solve the equations (better then inv function!)

a_n_out=a_n_in;
a_n_out(:,1:2)=a_1_to_2';


end

