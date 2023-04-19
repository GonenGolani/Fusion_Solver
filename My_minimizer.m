function [Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,final_Energy,exit_status] = My_minimizer...
    (Diaphragm,Cell,Virus,Shell,DOF_vector,int_DOF,res_struc,Minimazation,General_physical_properties)
% this function minimize the energy using the grdient decent method
% If the energy is deacreasing continue in the same direction and double
% next step size
% if energy is deacreasing go back and reduce size by 1/3



% it takes on paramter in DOF_vector(int_DOF) and minimize the energy up tp Tolerance
% the function recives the inital step size but not its direction, which is
% chosen randomly.
% the output is the new stracture, final energy and the new DOF_vector
% vector
%initilize vectors
exit_status=1; % exit_status_membrane 1 if successful 0 if error 

energy_vec=NaN(1,Minimazation.MAX_ITR);
steps_size=NaN(1,Minimazation.MAX_ITR);
count_increase_steps=0;
% the energy vector record the energy at each step of the iteration
[Diaphragm,Cell,Virus,Shell,~,Minimazation] = build_Hemifusion_structure(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc,General_physical_properties);
[Diaphragm,Cell,Virus,Shell,~,energy_vec(1)] = find_HD_energy(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);

negative_cut_off=-100;
if energy_vec(1)<negative_cut_off
    fprintf('Bad intitial condition in "my_minimizer" (<-100)\n');
end

initial_energy=energy_vec(1);
final_Energy=initial_energy+1;
% the position vector record the DOF value at each step
steps_size(1:2)=Minimazation.steps_size_vector(int_DOF);

save_Diaphragm0=Diaphragm;
save_Cell0=Cell;
save_Virus0=Virus;
save_DOF_vector0=DOF_vector;
Minimazation_0=Minimazation;
Minimazation_break_error_save0=Minimazation.break_error;


energy_decrease=1;


int=2; % start from second step
while int<=Minimazation.MAX_ITR %run until hit MAX_ITR which is global paramter

    if energy_decrease==1
        save_Diaphragm_int=Diaphragm;
        save_Cell_int=Cell;
        save_Virus_int=Virus;
        save_DOF_vector_int=DOF_vector;

        Minimazation_break_error_save=Minimazation.break_error;
    end
    %save last configuration in crgy increases

    %change DOF

    DOF_vector(int_DOF)=DOF_vector(int_DOF)+steps_size(int); % change paramter

    [~,DOF_vector] = Force_DOF_constreints(DOF_vector,Diaphragm,Cell,Virus,Minimazation,res_struc); % if there is another DOF change accordinly
    %avoid neagtive values for x0 and y0
    if int_DOF==1 ||  int_DOF==2
        if DOF_vector(1)<=0.05
            if Shell.Shell_exist==0
            DOF_vector(1)=0.1;
            else
                DOF_vector(1)=0.2;
            end
        end
        if DOF_vector(2)<=0.05
            if Shell.Shell_exist==0
                DOF_vector(2)=0.1;
            else
                DOF_vector(2)=0.2;
            end
        end
    end        

    %avoid over lap between virus and cell
    if int_DOF==4 %virus hight above cell
        if DOF_vector(4)<=General_physical_properties.lipid_length*3
            DOF_vector(4)=General_physical_properties.lipid_length*3;
        end
    end
    if int_DOF==6 %force cell rim smaller than radius
        if DOF_vector(6)>Cell.R_curv-0.05
            DOF_vector(6)=Cell.R_curv-0.05;
        end
    end
    if Minimazation.flat_diaphragm==1 && int_DOF==3 % if the diaphragm is fixed to flat and the z0 is minimized
       DOF_vector(7)=DOF_vector(3);
    end

    if Minimazation.fix_Cell_volume~=0
        [Cell.R_curv] = force_cell_volume_constraint(Cell,Diaphragm,Minimazation,res_struc);
    end

    % re calc shape and energy
    [~,DOF_vector] = Force_DOF_constreints(DOF_vector,Diaphragm,Cell,Virus,Minimazation,res_struc); % if there is another DOF change accordinly
    [Diaphragm,Cell,Virus,Shell] = DOF_vector_to_parameters_all_cases(Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,res_struc);

    [Diaphragm,Cell,Virus,Shell,exit_status,Minimazation] = build_Hemifusion_structure(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc,General_physical_properties);

    [DOF_vector,Minimazation] = Parameters_to_DOF_all_cases(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc); % to apply corrections if made to DOF_vector by build_Hemifusion_structure

    if exit_status==0 % break if error
        return;
    end
    [Diaphragm,Cell,Virus,Shell,~,energy_vec(int)] = find_HD_energy(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);
    final_Energy=energy_vec(int);
    Minimazation.steps_size_vector(int_DOF)=steps_size(int);

    %%

     if (energy_vec(int)<energy_vec(int-1)) % if the energy decreases keep direction
        count_increase_steps=0;
        energy_decrease=1;
        steps_size(int+1)=steps_size(int)*(rand/3+1);
        number_of_min_step=length(Minimazation.energy_each_step);
        Minimazation.energy_each_step(number_of_min_step+1)=final_Energy;
     end

    %stop condition
    if abs(energy_vec(int)-energy_vec(int-1))<Minimazation.Tolerance
        number_of_min_step=length(Minimazation.energy_each_step);
        Minimazation.energy_each_step(number_of_min_step+1)=final_Energy;
        break;
    end

     % if the energy increases flip direction and do half step next time (to avoid same position)
     % if the energy has very negative value consider it as error and
     % do the same thing, erese error energy
     do_reverse=0;
     if energy_vec(int)<negative_cut_off && abs(energy_vec(int)-energy_vec(int-1))>100
         do_reverse=1;
         energy_decrease=0;
         final_Energy=energy_vec(int-1);
         energy_vec(int)=energy_vec(int-1);
     end

     if (energy_vec(int)>energy_vec(int-1)) || do_reverse==1
         count_increase_steps=count_increase_steps+1;
         energy_increase=energy_vec(int)-energy_vec(int-1);
         energy_decrease=0;
         % do simulated anealing convergence
         if Minimazation.temperature==0 % break if temprature is zero
            %move one step back
            Diaphragm=save_Diaphragm_int;
            Cell=save_Cell_int;
            Virus=save_Virus_int;
            DOF_vector=save_DOF_vector_int;
            Minimazation.break_error=Minimazation_break_error_save;
            if energy_vec(int)/energy_vec(int-1)>1.5 
                % if the change is very big just return to first position
                steps_size(int+1)=-steps_size(int)/20; 
            else
                steps_size(int+1)=-steps_size(int)/3;
            end
         else  % break if temprature is not zero use Simulated Annealing algorithem
             if exp(-energy_increase/Minimazation.temperature)>rand
                Diaphragm=save_Diaphragm_int;
                Cell=save_Cell_int;
                Virus=save_Virus_int;
                DOF_vector=save_DOF_vector_int;
                Minimazation.break_error=Minimazation_break_error_save;
                steps_size(int+1)=-steps_size(int)/3;
             else
                steps_size(int+1)=steps_size(int)*(2*rand-1);  
             end
         end

         if count_increase_steps==10
             break;
         end

     end



    if int==Minimazation.MAX_ITR
        %Diaphragm=save_Diaphragm0;
        %Cell=save_Cell0;
        %Virus=save_Virus0;
        %DOF_vector=save_DOF_vector0;
        %Minimazation=Minimazation_0;
        %final_Energy=initial_energy;
        %Minimazation.steps_size_vector(int_DOF)=Minimazation.steps_size_vector(int_DOF)/5;
        %if exist('crushed point.mat','file')==0
        %    save('crushed point.mat');
        %end
        fprintf('Minimazation not complete\n');
        %break;
    end

    int=int+1;
end

%check if procedure worked, if not restore intial configuration
Minimazation.steps_size_vector(int_DOF)=steps_size(int);
if initial_energy<final_Energy ||  final_Energy<negative_cut_off
    Diaphragm=save_Diaphragm0;
    Cell=save_Cell0;
    Virus=save_Virus0;
    DOF_vector=save_DOF_vector0;
    [Diaphragm,Cell,Virus,Shell,exit_status,Minimazation] = build_Hemifusion_structure(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc,General_physical_properties);
    [Diaphragm,Cell,Virus,Shell,~,final_Energy] = find_HD_energy(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);
    Minimazation.break_error=Minimazation_break_error_save0;
    Minimazation.steps_size_vector(int_DOF)=Minimazation.steps_size_vector(int_DOF)/10;
    %fprintf('Restoring old configuration DOF: %.0f\n',int_DOF);  
else
    number_of_steps=length(Minimazation.energy_each_step);
    fprintf('Energy now: %f step: %.0f DOF: %.0f\n',final_Energy,number_of_steps,int_DOF);   
end   
    
end
