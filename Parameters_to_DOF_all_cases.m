function [DOF_vector,Minimazation] = Parameters_to_DOF_all_cases(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc)


if strcmp(Minimazation.step,'first')
    [DOF_vector,Minimazation.DOF_name_and_range] = parameters_to_DOF_vector_no_phi(Diaphragm,Cell,Virus,Minimazation,res_struc);
    if Shell.Shell_exist==1
        [DOF_vector,Minimazation.DOF_name_and_range] = Shell_parameters_to_DOF_vector(DOF_vector,Minimazation.DOF_name_and_range,Shell,res_struc);
    end
end

if strcmp(Minimazation.step,'second')
    if Minimazation.axial_symmetry==0
        if Minimazation.non_axial_symmetric_polynom==1
                [DOF_vector,Minimazation.DOF_name_and_range] = parameters_to_DOF_vector_with_smoothing_polynom(Diaphragm,Cell,Virus,Minimazation,res_struc);
            if Shell.Shell_exist==1  
               disp('Shell not updated for non symetric cases!');
            end                
        else  
            [DOF_vector] = parameters_to_DOF_all(Diaphragm,Cell,Virus,res_struc);
            if Shell.Shell_exist==1            
                disp('Shell not updated for non symetric cases!');
                [DOF_vector,Minimazation.DOF_name_and_range] = Shell_parameters_to_DOF_vector_no_axis_symmetry(DOF_vector,Minimazation.DOF_name_and_range,Shell,res_struc);
            end
        end    
    
    else % with axial symmetry
        [DOF_vector,Minimazation.DOF_name_and_range] = parameters_to_DOF_vector_no_phi(Diaphragm,Cell,Virus,Minimazation,res_struc);
        if Shell.Shell_exist==1
            [DOF_vector,Minimazation.DOF_name_and_range] = Shell_parameters_to_DOF_vector(DOF_vector,Minimazation.DOF_name_and_range,Shell,res_struc);
        end        
    end   
end

end

