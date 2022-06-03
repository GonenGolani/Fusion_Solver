function [Shell,Total_energy] = Total_shell_energy(Shell)

    [Shell.Diaphragm,Shell.Diaphragm.membrane_interaction_energy_density]=find_local_shell_energy(Shell.Diaphragm,Shell.Shell_physical_proprties,1);
    [Shell.Cell,Shell.Cell.membrane_interaction_energy_density]=find_local_shell_energy(Shell.Cell,Shell.Shell_physical_proprties,1);
    [Shell.Junction,Shell.Junction.membrane_interaction_energy_density]=find_local_shell_energy(Shell.Junction,Shell.Shell_physical_proprties,0);
    
    Total_energy=Shell.Diaphragm.total_energy+Shell.Cell.total_energy+Shell.Junction.total_energy;
    
end

