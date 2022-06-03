function [DOF_name_and_range] = create_DOF_name_and_range_all_cases(Minimazation,res_struc)


    if strcmp(Minimazation.step,'first') || (strcmp(Minimazation.step,'second') && Minimazation.axial_symmetry==1 )
        [DOF_name_and_range] = create_DOF_name_and_range_no_phi(Minimazation,res_struc);
    end
    
    if strcmp(Minimazation.step,'second') && Minimazation.axial_symmetry==0 
       [DOF_name_and_range] = create_DOF_name_and_range_smoothing_polynom(Minimazation,res_struc);  
    end

end

