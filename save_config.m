function [nothing]=save_config(name,Diaphragm,Cell,Virus,Shell,Minimazation,DOF_vector,res_struc,General_physical_properties)
System_dimensions.Cell_R_curv=Cell.R_curv;
System_dimensions.Virus_r_sv=Virus.r_sv;
System_dimensions.Virus_r_bv=Virus.r_bv;
System_dimensions.Diaphragm_r_pore=Diaphragm.r_pore;
Shell_exist=Shell.Shell_exist;
Shell_Shell_physical_proprties=Shell.Shell_physical_proprties;

save(name,'Minimazation','DOF_vector','res_struc','System_dimensions','General_physical_properties','Shell_exist','Shell_Shell_physical_proprties');
nothing=0;
end

