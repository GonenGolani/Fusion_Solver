function [surface_energy_vec,pore_radius_vec] = Pore_surface_energy...
    (Diaphragm,Cell,Virus,DOF_vector,res_struc,Minimazation,initial_energy,General_physical_properties)

    %create a folder to save the pores structures if it does not exist
    this_folder=pwd;
    if exist('Pore structure files')~=7
        mkdir 'Pore structure files'
    end
    cd 'Pore structure files'

    fileID1 = fopen('pore_radius_and_energy.txt','w');
    fileID2 = fopen('Diaphragm x0 y0.txt','w');
    
    %Diaphragm.r_pore=General_physical_properties.lipid_length; % initial pore size is lipid length
    surface_energy_vec(1)=initial_energy;
    pore_radius_vec(1)=Diaphragm.r_pore; 
    
    %first line
    fprintf(fileID1,'%0.4f %0.4f\n',pore_radius_vec(1),surface_energy_vec(1));
    fprintf(fileID2,'%0.4f %0.4f\n',Diaphragm.rim.x0,Diaphragm.rim.x0);
    
    int=2;
    
    while Diaphragm.r_pore<Diaphragm.rim.x0 && Diaphragm.r_pore<Diaphragm.rim.y0
        Diaphragm.r_pore=Diaphragm.r_pore+res_struc.pore_open_res;
        
        [Diaphragm,Cell,Virus] = build_Hemifusion_structure(Diaphragm,Cell,Virus,res_struc,General_physical_properties);
        [Diaphragm,Cell,Virus,~,initial_energy] = find_HD_energy(Diaphragm,Cell,Virus,res_struc);

        [Diaphragm,Cell,Virus,surface_energy_vec(int),DOF_vector,Minimazation] = My_converger_all_DOF...
            (Diaphragm,Cell,Virus,DOF_vector,res_struc,Minimazation,initial_energy,General_physical_properties,'all');
        
        pore_radius_vec(int)=Diaphragm.r_pore;
        file_name=['Pore size ' num2str(pore_radius_vec(int), 3) '.mat'];
        save(file_name) % save results
        fprintf(fileID1,'%0.4f %0.4f\n',pore_radius_vec(int),surface_energy_vec(int));
        fprintf(fileID2,'%0.4f %0.4f\n',Diaphragm.rim.x0,Diaphragm.rim.x0);

        int=int+1;
    end
    fclose(fileID1);
    fclose(fileID2);

    cd(this_folder);

end

