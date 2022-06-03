function [a_n_out,exit_status,Minimazation] = poly_force_boundary(zi,zf,dzdrhoi,dzdrhof,rhoi,rhof,a_n_in,Minimazation)

% Input= vector a_n_in with ploynum coeffients with length N.
% first 4 does not matter becuse found from boundary conditions
% boundary condtions:
% outer edge radius: Rad
% angle at rho=rhof: Ang_0 at 0 in deg
% angle at=rho=Rad: Ang_R at 0 in deg
% hight at rho=rhof: z0 
% hight at rho=Rad: z1 
% output= assign to 4 first places values using. give vector a_n_out 
% exit_status 1 if successful 0 if error 
% solved by inversing the matrix

exit_status=1;

A_mat=[1 rhoi rhoi^2 rhoi^3 ; 1 rhof rhof^2 rhof^3 ; 0 1 2*rhoi 3*rhoi^2 ; 0 1 2*rhof 3*rhof^2];


arr_size=size(a_n_in);
an_length=arr_size(2);

thorw_ot_first_4_vector=ones(1,an_length);
thorw_ot_first_4_vector(1:4)=0;
% RHS side solution vector
RHS_vector(:,1)=zi-mypolyval(a_n_in.*thorw_ot_first_4_vector,rhoi);
RHS_vector(:,2)=zf-mypolyval(a_n_in.*thorw_ot_first_4_vector,rhof);
RHS_vector(:,3)=dzdrhoi-mypolyder(a_n_in.*thorw_ot_first_4_vector,rhoi);
RHS_vector(:,4)=dzdrhof-mypolyder(a_n_in.*thorw_ot_first_4_vector,rhof);

%clear previous warnings
clear warnMsg warnId
lastwarn('')

%out an vector 1-4
A_inv=inv(A_mat);
%break if inversion is not ok
[~, warnId] = lastwarn ();
if isempty (warnId)==0
    disp('ERROR! BREAK!');
    if Minimazation.break_error==20 % break if error is hapening too many times
        save('last_cofig_before_break.mat');
        fid=fopen('Run terminated.txt','wt'); fprintf (fid,''); fclose (fid);
        exit;
    end
    exit_status=0;
    clear warnMsg warnId
    Minimazation.break_error=Minimazation.break_error+1;
end

if exit_status==1 %everything is ok

    a_1_to_4=A_inv*RHS_vector';

    %a_1_to_4=A_mat/RHS_vector'; %solve the equations (better then inv function!)

    a_n_out=a_n_in;
    a_n_out(:,1:4)=a_1_to_4';
else % problem
    a_n_out=0;
end

end

