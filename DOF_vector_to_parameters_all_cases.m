function [Diaphragm,Cell,Virus,Shell] = DOF_vector_to_parameters_all_cases(Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,res_struc)

        %#### any update in here must be also be implemented in "create_DOF_name_and_range_all_cases"
        if strcmp(Minimazation.step,'first')
            
            [Diaphragm,Cell,Virus,Shell] = DOF_vector_to_parameters_no_phi(Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,res_struc);
            
        end
        
        
        if strcmp(Minimazation.step,'second')  
            %do separte minimaztion to each polar angle phi
            if Minimazation.non_axial_symmetric_polynom==0 && Minimazation.axial_symmetry==0 
                 disp('not updated for non symetric cases with general phi polynom!');
            end
            % use smoothing polynom to minimize the matrix coefficiants
            if Minimazation.non_axial_symmetric_polynom==1 && Minimazation.axial_symmetry==0
                [Diaphragm,Cell,Virus] = DOF_vector_to_parameters_with_smoothing_polynom(DOF_vector,Diaphragm,Cell,Virus,Minimazation,res_struc);
            end
            %smae as first step
            if   Minimazation.axial_symmetry==1 
                [Diaphragm,Cell,Virus,Shell] = DOF_vector_to_parameters_no_phi(Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,res_struc);
            end
        end  


end

