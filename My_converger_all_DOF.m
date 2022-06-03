function [Diaphragm,Cell,Virus,Shell,Final_Energy,DOF_vector,Minimazation] = ...
    My_converger_all_DOF(Diaphragm_old,Cell_old,Virus_old,Shell_old,DOF_vector_old,res_struc,Minimazation_old,General_physical_properties,which)
    % minimize the energy by iterative search.
    % at each step change one parameter in DOF_vector. 
    
    


    
    Diaphragm=Diaphragm_old;
    Cell=Cell_old;
    Virus=Virus_old;
    Shell=Shell_old;
    DOF_vector=DOF_vector_old;
    Minimazation=Minimazation_old;
    



    
    energy_vec_all=NaN(1,Minimazation.MAX_ITR);

    
    [fix_DOF,DOF_vector] = Force_DOF_constreints(DOF_vector,Diaphragm,Cell,Virus,Minimazation,res_struc); % fix constreints
    


 
    [Diaphragm,Cell,Virus,Shell] = DOF_vector_to_parameters_all_cases(Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,res_struc); % correct accordinly
    


    [Diaphragm,Cell,Virus,Shell] = build_Hemifusion_structure(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc,General_physical_properties); 
    [Diaphragm,Cell,Virus,Shell,~,initial_energy] = find_HD_energy(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);

    
    %save initial state
    save_Diaphragm_int=Diaphragm;
    save_Cell_int=Cell;
    save_Virus_int=Virus;
    save_DOF_vector_int=DOF_vector;
    Minimazation_break_error_save=Minimazation.break_error;
    energy_decrease=0; % assume that energy is increasing
    
    energy_vec_all(1)=initial_energy;
    
    %create a random vector of 1 and -1 and multiply it with initial_steps_size
    Minimazation.steps_size_vector=((randi(2,[1,length(DOF_vector)])-1)*2-1).*rand(1,length(DOF_vector))*Minimazation.max_step;
    
    

    %do recursive minimaztion loops until global minima is achived
    int_loops=2;
    while  int_loops<=Minimazation.MAX_ITR
        
         % save last stable config if lowered the energy
         
        if energy_decrease==1
        save_Diaphragm_int=Diaphragm;
        save_Cell_int=Cell;
        save_Virus_int=Virus;
        save_DOF_vector_int=DOF_vector;
        Minimazation_break_error_save=Minimazation.break_error;
        end
         %we next run over all DOF until global minima achived      

        last_energy=energy_vec_all(int_loops-1);
        

        
        %generate random order to minimize DOF vector
         if strcmp(which,'only first')
             if  Minimazation.exponential_decay_tilt==1
                DOF_number=(7+2*res_struc.poly_degree_rim+6);        
             else
                DOF_number=(7+2*res_struc.poly_degree_rim);
             end
         else
            DOF_number=length(DOF_vector);
         end
        DOV_min_vector=randperm(DOF_number);   
        
        
        
        int_DOF=1;      
        while int_DOF<=DOF_number %run over all DOF 
            
            minimize_int=DOV_min_vector(int_DOF);
            if fix_DOF(minimize_int)~=1 % check if the number is on the fixed list
                if strcmp(which,'only 0') % minimize degreess of freedom that are zero
                    if DOF_vector(minimize_int)==0
                    [Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,last_energy,exit_status]...
                        =My_minimizer(Diaphragm,Cell,Virus,Shell,DOF_vector,minimize_int,res_struc,Minimazation,General_physical_properties);
                    end
                else
                [Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,last_energy,exit_status]...
                    =My_minimizer(Diaphragm,Cell,Virus,Shell,DOF_vector,minimize_int,res_struc,Minimazation,General_physical_properties);                
                end
                
                
                if exit_status==0 ||  last_energy>initial_energy % if there is a problem
                    disp('old state reload!');
                    Diaphragm=save_Diaphragm_int;
                    Cell=save_Cell_int;
                    Virus=save_Virus_int;
                    DOF_vector=save_DOF_vector_int; 
                    Minimazation.break_error=Minimazation_break_error_save;
                    energy_decrease=0;
                else
                    energy_decrease=1;
                end
                
            end
            int_DOF=int_DOF+1;
            
        end
        
        energy_vec_all(int_loops)=last_energy;
        Final_Energy=energy_vec_all(int_loops);
        
        if abs(energy_vec_all(int_loops)-energy_vec_all(int_loops-1))<Minimazation.Tolerance
            break;
        end
        

        if int_loops==Minimazation.MAX_ITR
            fprintf('Minimazation not complete in My_converger_all_DOF\n');
            break;
        end
        %[Minimazation.steps_size_vector] = reduce_step_site_if_DOF_is0(Minimazation.steps_size_vector,DOF_vector);

        int_loops=int_loops+1;
    end
    

end

