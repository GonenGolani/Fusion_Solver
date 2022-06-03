function [Diaphragm,Cell,Virus,Shell,General_physical_properties,Minimazation,DOF_vector,last_step] = ...
    interate_value(value_to_change,target_value,delta_value,initial_step,Diaphragm_old,Cell_old,Virus_old,Shell_old,General_physical_properties_old,Minimazation_old,res_struc,DOF_vector_old)

Diaphragm=Diaphragm_old;
Cell=Cell_old;
Virus=Virus_old;
Shell=Shell_old;
General_physical_properties=General_physical_properties_old;
Minimazation=Minimazation_old;
DOF_vector=DOF_vector_old;

%find if need to increase or decrease value
value_now=find_value_from_name(value_to_change,Cell,Virus,Shell,General_physical_properties,Minimazation);
delta_value=abs(delta_value); % make sure the delta is positive
n=Minimazation.scan_variable.value_n;
%number of steps
max_steps=ceil(abs(target_value^n-value_now^n)/delta_value^n);
%steps direction
sign=1;
if target_value-value_now<0
sign=-1;
end



int_step=initial_step;
last_step=int_step;
while int_step<max_steps+initial_step
    
    value_now=(value_now^n+sign*delta_value^n)^(1/n);
    [Cell,Virus,Shell,General_physical_properties,Minimazation] = ...
        scan_value_translator(value_to_change,value_now,Cell,Virus,Shell,General_physical_properties,Minimazation);

    
    if Minimazation.fix_Cell_volume~=0
        [Cell.R_curv] = force_cell_volume_constraint(Cell,Diaphragm,Minimazation,res_struc);
    end
    
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



    %save results
    save_name=['step_' num2str(int_step) '.mat'];
    save_config(save_name,Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,res_struc,General_physical_properties);
    [~,~,~,~,~,energy_now] = find_HD_energy(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);

    fid=fopen('Energy in step.txt','a+'); fprintf (fid,'%.0f %f\n',int_step,energy_now); fclose (fid);
    %next step
    int_step=int_step+1;
    last_step=int_step;

end

end

