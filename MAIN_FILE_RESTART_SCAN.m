
% Fusion intermidate solver
clc
clear all
close all
tic;
this_folder=pwd;


%use if first step is diffrent
if exist('first step number.txt','file')==2
    first_step_number=importdata('first step number.txt');
else
    first_step_number=1; 
end


%use if first step is diffrent
if exist('repeats before start.txt','file')==2
    repeats_before_step_1=importdata('repeats before start.txt');
else
    repeats_before_step_1=1; 
end

if isfile('Add Shell.txt')
    [~,~,~,Shell,General_physical_properties,Minimazation,res_struc,~] = read_conditions;
    Minimazation.Shell_exist=1;
    cd(this_folder);
    

    
    [Diaphragm,Cell,Virus,~,~,res_struc,~,~,~,~,~] = reload_HD_state('initial_configuration.mat',res_struc);

    %create the new shell
    Shell.Diaphragm.MP_to_MP_distance_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_shell); 
    Shell.Cell.MP_to_MP_distance_matrix=zeros(res_struc.phi_res,res_struc.poly_degree_shell); 
    Shell.Junction.mid_plane_coeff_matrix=zeros(res_struc.phi_res,4); 

    if Minimazation.do_stalk==1
        Shell.Diaphragm.z_i=Shell.Cell.z_f+Diaphragm.z0;
        Shell.Diaphragm.z_f=Shell.Cell.z_f+Diaphragm.z0;
    end

else 
    [~,~,~,~,General_physical_properties,Minimazation,res_struc,~] = read_conditions;
    
    cd(this_folder);
    
    [Diaphragm,Cell,Virus,Shell,~,~,~,~,~,~,~] = reload_HD_state('initial_configuration.mat',res_struc);

end

%add or remove lines to matrix if needed
[Diaphragm,Cell,Virus,Shell] = add_lines_to_coef_matrix(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);


%cell proprties
cd('Parameters'); cd('System dimensions');
Cell.R_curv=importdata('Cell radius.txt');  % curvature of cell
Virus.r_sv=importdata('Virus small radius.txt'); % small radius of torus
Virus.r_bv=importdata('Virus big radius.txt'); % large radius of torus
cd(this_folder);

Minimazation.step='second';


Minimazation.break_error=0; %number of times the simulation was breaked and returned to last satble configuration
clear DOF_vector

if Minimazation.axial_symmetry==0 && Minimazation.non_axial_symmetric_polynom==1 && length(Cell.rim.drho_dz_vector)==1
        % reset the rim
        Cell.rim.drho_dz_vector(2:res_struc.poly_degree_rim)=0;
        Virus.rim.drho_dz_vector(2:res_struc.poly_degree_rim)=0;
        Diaphragm.rim.hl_vector(2:res_struc.poly_degree_rim)=0;
        Diaphragm.rim.drho_dz_vector(2:res_struc.poly_degree_rim)=0;
        
        [Diaphragm,Cell,Virus] = save_initial_smoothing_polynom_to_matrix(Diaphragm,Cell,Virus,res_struc);      
end

[DOF_vector,Minimazation] = Parameters_to_DOF_all_cases(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);

[Diaphragm,Cell,Virus,Shell] = build_Hemifusion_structure(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc,General_physical_properties);
[Diaphragm,Cell,Virus,Shell,~,TOTAL_ENERGY_AFTER] = find_HD_energy(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);
Minimazation.energy_each_step=TOTAL_ENERGY_AFTER;

%[Diaphragm,Cell,Virus,Shell] = build_Hemifusion_structure(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc,General_physical_properties); %for test can be deleted
%plot_profile(Diaphragm,Cell,Virus,Shell,1,1,1); %for test can be deleted
%[Diaphragm,Cell,Virus,Shell,Energy,TOTAL_ENERGY_AFTER] = find_HD_energy(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc); %for test can be deleted

%% first min  to get the system to right position


if repeats_before_step_1>0
[Diaphragm,Cell,Virus,Shell,~,DOF_vector,Minimazation] = ...
    My_converger_all_DOF(Diaphragm,Cell,Virus,Shell,DOF_vector,res_struc,Minimazation,General_physical_properties,'only first');
end

int_repeats=1;
while int_repeats<=repeats_before_step_1

    [Diaphragm,Cell,Virus,Shell,~,DOF_vector,Minimazation] = ...
        My_converger_all_DOF(Diaphragm,Cell,Virus,Shell,DOF_vector,res_struc,Minimazation,General_physical_properties,'all');

    int_repeats=int_repeats+1;
end

save_name=['step_' num2str(first_step_number) '.mat'];

save_config(save_name,Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,res_struc,General_physical_properties);


%% iterative varible change

% epoch 1
[Diaphragm,Cell,Virus,Shell,General_physical_properties,Minimazation,DOF_vector,last_step_epoch1] = ...
    interate_value(Minimazation.scan_variable.variable_name,Minimazation.scan_variable.target_value,Minimazation.scan_variable.delta_value,first_step_number+1 ...
    ,Diaphragm,Cell,Virus,Shell,General_physical_properties,Minimazation,res_struc,DOF_vector);
if Minimazation.scan_variable.delta_value_2~=0
    % epoch 2
    [Diaphragm,Cell,Virus,Shell,General_physical_properties,Minimazation,DOF_vector,last_step_epoch2] = ...
        interate_value(Minimazation.scan_variable.variable_name,Minimazation.scan_variable.target_value_2,Minimazation.scan_variable.delta_value_2,last_step_epoch1 ...
        ,Diaphragm,Cell,Virus,Shell,General_physical_properties,Minimazation,res_struc,DOF_vector);

    if Minimazation.scan_variable.delta_value_3~=0
        % epoch 3
        [Diaphragm,Cell,Virus,Shell,General_physical_properties,Minimazation,DOF_vector] = ...
            interate_value(Minimazation.scan_variable.variable_name,Minimazation.scan_variable.target_value_3,Minimazation.scan_variable.delta_value_3,last_step_epoch2 ...
            ,Diaphragm,Cell,Virus,Shell,General_physical_properties,Minimazation,res_struc,DOF_vector);
    end
end


%% end time measure and save output
elapsedTime = toc; % in seconds
fid=fopen('AAA run time(sec).txt','wt'); fprintf (fid,'%f',elapsedTime); fclose (fid);
exit;

