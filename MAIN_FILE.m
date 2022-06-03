

clc
clear all
close all
this_folder=pwd;

[Diaphragm,Virus,Cell,Shell,General_physical_properties,Minimazation,res_struc,intital_tilt_decay_length] = read_conditions;


Minimazation.break_error=0; %number of times the simulation was breaked and returned to last satble configuration
cd(this_folder);


% start measure cose time
tic;
%% 


%coefficiant matrix initialazation
[Diaphragm,Cell,Virus,Shell] = Create_intial_coeff_matrix(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc,intital_tilt_decay_length);



Minimazation.step='first'; 

%% first step of minimaztion - do not allow change in coefficiantans in the phi direction 
[DOF_vector,Minimazation] = Parameters_to_DOF_all_cases(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);

[~,DOF_vector] = Force_DOF_constreints(DOF_vector,Diaphragm,Cell,Virus,Minimazation,res_struc); % fix constreints
[Diaphragm,Cell,Virus,Shell] = DOF_vector_to_parameters_all_cases(Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,res_struc); % correct accordinly
[Diaphragm,Cell,Virus,Shell,~,Minimazation] = build_Hemifusion_structure(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc,General_physical_properties);

if Minimazation.fix_Cell_volume~=0
   [Cell.R_curv] = force_cell_volume_constraint(Cell,Diaphragm,Minimazation,res_struc);
   [Diaphragm,Cell,Virus,Shell,~,Minimazation] = build_Hemifusion_structure(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc,General_physical_properties);
end


[Diaphragm,Cell,Virus,Shell,~,Last_energy] = find_HD_energy(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);
Minimazation.energy_each_step=Last_energy;
[DOF_vector,Minimazation] = Parameters_to_DOF_all_cases(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);


%First minimaztion in low res


[Diaphragm,Cell,Virus,Shell,~,DOF_vector,Minimazation] = ...
    My_converger_all_DOF(Diaphragm,Cell,Virus,Shell,DOF_vector,res_struc,Minimazation,General_physical_properties,'only first');

Minimazation.temperature=Minimazation.temperature/2;


[Diaphragm,Cell,Virus,Shell,~,DOF_vector,Minimazation] = ...
    My_converger_all_DOF(Diaphragm,Cell,Virus,Shell,DOF_vector,res_struc,Minimazation,General_physical_properties,'all');

%% second round increasing tolerance
Minimazation.temperature=0;

[Diaphragm,Cell,Virus,Shell,~,DOF_vector,Minimazation] = ...
    My_converger_all_DOF(Diaphragm,Cell,Virus,Shell,DOF_vector,res_struc,Minimazation,General_physical_properties,'all');




if Minimazation.axial_symmetry==1 % if this is last step do 
    
    [Diaphragm,Cell,Virus,Shell,~,DOF_vector,Minimazation] = ...
        My_converger_all_DOF(Diaphragm,Cell,Virus,Shell,DOF_vector,res_struc,Minimazation,General_physical_properties,'all');
    
end

%save('After_first_min.mat')
fprintf('First step done\n');
save_config('HD_state_thin.mat',Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,res_struc,General_physical_properties);
Minimazation.step='second';


%% second step  - minimize each phi seprately no DOF adding
if Minimazation.axial_symmetry==0

    clear DOF_vector;
    if Minimazation.non_axial_symmetric_polynom==1
        [Diaphragm,Cell,Virus] = save_initial_smoothing_polynom_to_matrix(Diaphragm,Cell,Virus,res_struc); 
    end
    [DOF_vector,Minimazation] = Parameters_to_DOF_all_cases(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);
  
    [Diaphragm,Cell,Virus,Shell,~,Minimazation] = build_Hemifusion_structure(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc,General_physical_properties);

[Diaphragm,Cell,Virus,Shell,~,~] = find_HD_energy(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);
end

%last minimaztion

[Diaphragm,Cell,Virus,Shell,~,DOF_vector,Minimazation] = ...
    My_converger_all_DOF(Diaphragm,Cell,Virus,Shell,DOF_vector,res_struc,Minimazation,General_physical_properties,'all');

fprintf('HD minimaztion done\n');



%% 

if Minimazation.do_iterative==1
save_config('step_1.mat',Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,res_struc,General_physical_properties);
% epoch 1
[Diaphragm,Cell,Virus,Shell,General_physical_properties,Minimazation,DOF_vector,last_step_epoch1] = ...
    interate_value(Minimazation.scan_variable.variable_name,Minimazation.scan_variable.target_value,Minimazation.scan_variable.delta_value,2 ...
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

end





%end time measure and save output
elapsedTime = toc; % in seconds
fid=fopen('AAA run time(sec).txt','wt'); fprintf (fid,'%f',elapsedTime); fclose (fid);
exit;


