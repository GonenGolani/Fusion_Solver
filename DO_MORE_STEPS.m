

clc
clear all
close all
this_folder=pwd;

% start measure cose time
tic;

[Diaphragm,Virus,Cell,Shell,General_physical_properties,Minimazation,res_struc,intital_tilt_decay_length] = read_conditions;
cd(this_folder);

%load new staff
scan_variable_temp=Minimazation.scan_variable;
Target_Tolerance=Minimazation.Tolerance;
res_struc_temp=res_struc;
last_step_id=importdata('last_step_id.txt');

[Diaphragm,Cell,Virus,Shell,Minimazation,res_struc,General_physical_properties,Energy,Total_energy,DOF_vector,System_dimensions]...
    = reload_HD_state(['step_' num2str(last_step_id) '.mat']);

%
Minimazation.scan_variable=scan_variable_temp;
res_struc=res_struc_temp;
Minimazation.Tolerance=Target_Tolerance;

%find if need to increase or decrease value
clear eval
value_now=eval(Minimazation.scan_variable.variable_name);
Minimazation.scan_variable.delta_value=abs(Minimazation.scan_variable.delta_value); % make sure the delta is positive
n=Minimazation.scan_variable.value_n;

%number of steps
max_steps=ceil(abs(Minimazation.scan_variable.target_value^n-value_now^n)/Minimazation.scan_variable.delta_value^n)+last_step_id;
%steps direction
sign=1;
if Minimazation.scan_variable.target_value-value_now<0
    sign=-1;
end


int_step=last_step_id+1;
value=(value_now^n+sign*Minimazation.scan_variable.delta_value^n)^(1/n);
while int_step<=max_steps

    evalin('base',[Minimazation.scan_variable.variable_name,'=',num2str(value)]); % reassign new value to desired 'variable_name'

    if Minimazation.relaxed_spherical_VLP==1
        alpha=General_physical_properties.cell_distal_kappa/General_physical_properties.proximal_kappa;
        General_physical_properties.Cell_distal_J0=(General_physical_properties.proximal_J0-2*(1+alpha)/Cell.R_curv)/alpha;
    end

    if Minimazation.relaxed_cylindrical_VLP==1
        alpha=General_physical_properties.virus_distal_kappa/General_physical_properties.proximal_kappa;
        General_physical_properties.Virus_distal_J0=(General_physical_properties.proximal_J0-(1+alpha)/Virus.r_sv)/alpha;
    end

    if Minimazation.symmetric_J0==1
        General_physical_properties.Virus_distal_J0=General_physical_properties.proximal_J0;
        General_physical_properties.Cell_distal_J0=General_physical_properties.proximal_J0;
    end

    if Minimazation.symmetric_distal_J0~=0
        General_physical_properties.Virus_distal_J0=Minimazation.symmetric_distal_J0;
        General_physical_properties.Cell_distal_J0=Minimazation.symmetric_distal_J0;  
    end


    Minimazation.energy_each_step=0;
    %rebuild HD
    [Diaphragm,Cell,Virus,Shell,~,Minimazation] = build_Hemifusion_structure(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc,General_physical_properties);
    [Diaphragm,Cell,Virus,Shell,~,~] = find_HD_energy(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);



    %re minimize the energy
    [Diaphragm,Cell,Virus,Shell,~,DOF_vector,Minimazation] = My_converger_all_DOF...
        (Diaphragm,Cell,Virus,Shell,DOF_vector,res_struc,Minimazation,General_physical_properties,'only first');

    [Diaphragm,Cell,Virus,Shell,~,DOF_vector,Minimazation] = My_converger_all_DOF...
        (Diaphragm,Cell,Virus,Shell,DOF_vector,res_struc,Minimazation,General_physical_properties,'all'); 

    Minimazation.Tolerance=Target_Tolerance;


    %save results
    save_name=['step_' num2str(int_step) '.mat'];
    save_config(save_name,Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,res_struc,General_physical_properties);

    %next step
    int_step=int_step+1;
    value=(value^n+sign*Minimazation.scan_variable.delta_value^n)^(1/n);
end






%end time measure and save output
elapsedTime = toc; % in seconds
fid=fopen('AAA run time second run(sec).txt','wt'); fprintf (fid,'%f',elapsedTime); fclose (fid);
exit;


