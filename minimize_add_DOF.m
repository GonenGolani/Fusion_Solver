function [Diaphragm,Cell,Virus,DOF_vector,res_struc,Minimazation,General_physical_properties] =...
    minimize_add_DOF(Diaphragm,Cell,Virus,DOF_vector,res_struc,Minimazation,General_physical_properties,initial_energy)
% this function gradualy add DOF until minimal value is achived


energy_vec=NaN(1,Minimazation.MAX_ITR);

energy_vec(1)=initial_energy;

%poly z increase
int=2;
while res_struc.poly_degree_z<Minimazation.MAX_DOF_MP
    old_res_struc=res_struc;
    save_Diaphragm=Diaphragm;
    save_Cell=Cell;
    save_Virus=Virus;
    save_DOF_vector=DOF_vector;
    save_Minimazation=Minimazation;
    
    res_struc.poly_degree_z=res_struc.poly_degree_z+1;
    
    [DOF_vector,Minimazation] = ADD_DOF(old_res_struc,res_struc,DOF_vector,Minimazation);
    [Diaphragm,Cell,Virus] = DOF_vector_to_parameters_no_phi(Diaphragm,Cell,Virus,Minimazation,DOF_vector,res_struc);
    [Diaphragm,Cell,Virus,energy_vec(int),DOF_vector,Minimazation] = ...
        My_converger_all_DOF(Diaphragm,Cell,Virus,DOF_vector,res_struc,Minimazation,energy_vec(int-1),General_physical_properties,'only 0');
    [Diaphragm,Cell,Virus,energy_vec(int),DOF_vector,Minimazation] = ...
        My_converger_all_DOF(Diaphragm,Cell,Virus,DOF_vector,res_struc,Minimazation,energy_vec(int),General_physical_properties,'all');

    if  abs(energy_vec(int)-energy_vec(int-1))<Minimazation.Tolerance || energy_vec(int)>energy_vec(int-1)
        Diaphragm=save_Diaphragm;
        Cell=save_Cell;
        Virus=save_Virus;
        DOF_vector=save_DOF_vector;
        res_struc=old_res_struc;  
        Minimazation=save_Minimazation;
        break;
    end
    int=int+1;
end

%poly rim increase
while res_struc.poly_degree_rim<Minimazation.MAX_DOF_rim
    old_res_struc=res_struc;
    save_Diaphragm=Diaphragm;
    save_Cell=Cell;
    save_Virus=Virus;
    save_DOF_vector=DOF_vector;
    save_Minimazation=Minimazation;
    
    res_struc.poly_degree_rim=res_struc.poly_degree_rim+1;
    
    [DOF_vector,Minimazation] = ADD_DOF(old_res_struc,res_struc,DOF_vector,Minimazation);
    
    [Diaphragm,Cell,Virus] = DOF_vector_to_parameters_no_phi(Diaphragm,Cell,Virus,Minimazation,DOF_vector,res_struc);
    
        [Diaphragm,Cell,Virus,energy_vec(int),DOF_vector,Minimazation] = ...
    My_converger_all_DOF(Diaphragm,Cell,Virus,DOF_vector,res_struc,Minimazation,energy_vec(int-1),General_physical_properties,'only 0');
    [Diaphragm,Cell,Virus,energy_vec(int),DOF_vector,Minimazation] = ...
        My_converger_all_DOF(Diaphragm,Cell,Virus,DOF_vector,res_struc,Minimazation,energy_vec(int),General_physical_properties,'all');
    
    if abs(energy_vec(int)-energy_vec(int-1))<Minimazation.Tolerance  || energy_vec(int)>energy_vec(int-1)
        Diaphragm=save_Diaphragm;
        Cell=save_Cell;
        Virus=save_Virus;
        DOF_vector=save_DOF_vector;
        res_struc=old_res_struc;
        Minimazation=save_Minimazation;
        break;
    end
    int=int+1;
end

%poly poly_degree_LD_rho increase
while res_struc.poly_degree_LD_rho<Minimazation.MAX_DOF_angle_RHO
    old_res_struc=res_struc;
    save_Diaphragm=Diaphragm;
    save_Cell=Cell;
    save_Virus=Virus;
    save_DOF_vector=DOF_vector;
    save_Minimazation=Minimazation;
    
    res_struc.poly_degree_LD_rho=res_struc.poly_degree_LD_rho+1;
    [DOF_vector,Minimazation] = ADD_DOF(old_res_struc,res_struc,DOF_vector,Minimazation);
    [Diaphragm,Cell,Virus] = DOF_vector_to_parameters_no_phi(Diaphragm,Cell,Virus,Minimazation,DOF_vector,res_struc);

    [Diaphragm,Cell,Virus,energy_vec(int),DOF_vector,Minimazation] = ...
        My_converger_all_DOF(Diaphragm,Cell,Virus,DOF_vector,res_struc,Minimazation,energy_vec(int-1),General_physical_properties,'only 0');
    [Diaphragm,Cell,Virus,energy_vec(int),DOF_vector,Minimazation] = ...
        My_converger_all_DOF(Diaphragm,Cell,Virus,DOF_vector,res_struc,Minimazation,energy_vec(int),General_physical_properties,'all');
    
    if abs(energy_vec(int)-energy_vec(int-1))<Minimazation.Tolerance  || energy_vec(int)>energy_vec(int-1)
        Diaphragm=save_Diaphragm;
        Cell=save_Cell;
        Virus=save_Virus;
        DOF_vector=save_DOF_vector;
        res_struc=old_res_struc;
        Minimazation=save_Minimazation;
        break;
    end
    int=int+1;
end

%poly poly_degree_LD_phi increase
while res_struc.poly_degree_LD_phi<Minimazation.MAX_DOF_angle_PHI
    old_res_struc=res_struc;
    save_Diaphragm=Diaphragm;
    save_Cell=Cell;
    save_Virus=Virus;
    save_DOF_vector=DOF_vector;
    save_Minimazation=Minimazation;
    
    res_struc.poly_degree_LD_phi=res_struc.poly_degree_LD_phi+1;
    [DOF_vector,Minimazation] = ADD_DOF(old_res_struc,res_struc,DOF_vector,Minimazation);
    [Diaphragm,Cell,Virus] = DOF_vector_to_parameters_no_phi(Diaphragm,Cell,Virus,Minimazation,DOF_vector,res_struc);
    
    [Diaphragm,Cell,Virus,energy_vec(int),DOF_vector,Minimazation] = ...
        My_converger_all_DOF(Diaphragm,Cell,Virus,DOF_vector,res_struc,Minimazation,energy_vec(int-1),General_physical_properties,'all');
    [Diaphragm,Cell,Virus,energy_vec(int),DOF_vector,Minimazation] = ...
        My_converger_all_DOF(Diaphragm,Cell,Virus,DOF_vector,res_struc,Minimazation,energy_vec(int),General_physical_properties,'all');
    
    if abs(energy_vec(int)-energy_vec(int-1))<Minimazation.Tolerance  || energy_vec(int)>energy_vec(int-1)
        Diaphragm=save_Diaphragm;
        Cell=save_Cell;
        Virus=save_Virus;
        DOF_vector=save_DOF_vector;
        res_struc=old_res_struc;
        Minimazation=save_Minimazation;
        break;
    end
    int=int+1;
end



end

