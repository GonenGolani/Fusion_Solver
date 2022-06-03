function  [Local_shell,membrane_interaction_energy_density]=find_local_shell_energy(Local_shell,Shell_physical_proprties,include_membrane_intercation)

    
    %this function finds the local shell energy in the diffrent regions
    % abbreviations
    area_element=Local_shell.mid_plane.area_element;
    total_curvature=Local_shell.total_curvature;
    Gaussian_curvature=Local_shell.Gaussian_curvature;
    det_metric=Local_shell.mid_plane.metric_det;
    undeformed_area_element=Local_shell.mid_plane.area_element./det_metric.^0.5;

    kappa=Shell_physical_proprties.kappa;
    kappa_bar=Shell_physical_proprties.kappa_bar;
    streatch_modulus=Shell_physical_proprties.streatch_modulus;
    J0=Shell_physical_proprties.J0;
    shell_half_width=Shell_physical_proprties.width/2;
    monolayer_width=Shell_physical_proprties.monolayer_width;
    
    %the 4 factor is becuse we calc only qurter of system!
    
    Local_shell.total_area=4*sum(area_element,'all');
    
    Local_shell.stretching_energy=4*sum(0.5*streatch_modulus.*(det_metric.^0.5-1).^2.*undeformed_area_element,'all');
    
    Local_shell.bending_energy=4*sum(0.5*kappa.*(total_curvature-J0).^2.*area_element,'all');
    
    
    Local_shell.Saddle_splay_energy=4*sum(kappa_bar.*Gaussian_curvature.*area_element,'all');
    
    %add non-linear effects to prevent too curvature larger then width
    detJK=(total_curvature.^2-4*Gaussian_curvature).^0.5;
    c1=(total_curvature+detJK)/2;
    c2=(total_curvature-detJK)/2;
    if max(abs(c1),[],'all')<Shell_physical_proprties.max_shell_curvature && max(abs(c2),[],'all')<Shell_physical_proprties.max_shell_curvature
        prohibition_energy_curvature=0;
    else
        prohibition_energy_curvature=1000;
    end
    

    
    prohibition_energy_overlap=0;
    u0=Shell_physical_proprties.intercation_potential_energy_density;    
    if include_membrane_intercation==1
        z_tag=Local_shell.distance_from_membrane;
        z0=Shell_physical_proprties.interaction_distance;

        if strcmp(Shell_physical_proprties.Shell_membrane_interaction,'Lennard_Jones')
            membrane_interaction_energy_density=4*u0*((z0./z_tag).^12-(z0./z_tag).^6+1/4);
        end
        if strcmp(Shell_physical_proprties.Shell_membrane_interaction,'Lennard_Jones_modified')
            membrane_interaction_energy_density=u0*((z0./z_tag).^12-2*(z0./z_tag).^6+1); %2 order expension as in spring
        end
        if strcmp(Shell_physical_proprties.Shell_membrane_interaction,'Spring')
            membrane_interaction_energy_density=u0*(((z_tag-z0)./z0).^2);
        end
        %steric repulsion term regrdless of interaction type
        if min(z_tag-monolayer_width-shell_half_width,[],'all')<0
            prohibition_energy_overlap=1000;
        end
        
    end
    if include_membrane_intercation==0 %the shell is disconected
        membrane_interaction_energy_density=u0;
    end
    

    Local_shell.membrane_interaction_energy=4*sum(membrane_interaction_energy_density.*area_element,'all');  
    
    Local_shell.total_energy=Local_shell.stretching_energy+Local_shell.bending_energy+Local_shell.Saddle_splay_energy...
        +Local_shell.membrane_interaction_energy+prohibition_energy_curvature+prohibition_energy_overlap;
    
    Local_shell.energy_density=0.5*streatch_modulus.*(det_metric.^0.5-1).^2+0.5*kappa.*(total_curvature-J0).^2 ...
        +kappa_bar.*Gaussian_curvature+membrane_interaction_energy_density;
    
end