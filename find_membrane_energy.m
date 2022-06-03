function  [total_energy,total_area_up,total_area_down,energy_components]...
    =find_membrane_energy(tilt,area_element,splay,saddle_splay,physical_properties)
    
    %the 4 factor is becuse we calc only qurter of system!
    total_area_up=4*sum(area_element.up,'all');
    total_area_down=4*sum(area_element.down,'all');
    
    tilt_energy_up=4*sum(0.5*physical_properties.kappa_tilt*(tilt.up.x.^2+tilt.up.y.^2+tilt.up.z.^2).*area_element.up,'all');
    tilt_energy_down=4*sum(0.5*physical_properties.kappa_tilt*(tilt.down.x.^2+tilt.down.y.^2+tilt.down.z.^2).*area_element.down,'all');
    
    splay_energy_up=4*sum(0.5*physical_properties.kappa_up.*(splay.up.^2-2*splay.up*physical_properties.J0_up).*area_element.up,'all');
    splay_energy_down=4*sum(0.5*physical_properties.kappa_down.*(splay.down.^2-2*splay.down*physical_properties.J0_down).*area_element.down,'all');
    
    saddle_splay_energy_up=4*sum(physical_properties.kappa_bar_up.*saddle_splay.up.*area_element.up,'all');
    saddle_splay_energy_down=4*sum(physical_properties.kappa_bar_down.*saddle_splay.down.*area_element.down,'all');

    if physical_properties.Splay_grad_square_modulus~=0
        splay_grad_square_energy_density_up=0.5*physical_properties.Splay_grad_square_modulus*(splay.splay_grad_up).^2;
        splay_grad_square_energy_density_down=0.5*physical_properties.Splay_grad_square_modulus*(splay.splay_grad_down).^2;
        splay_grad_square_energy_up=4*sum(splay_grad_square_energy_density_up.*area_element.up,'all');
        splay_grad_square_energy_down=4*sum(splay_grad_square_energy_density_down.*area_element.down,'all');
    else
        splay_grad_square_energy_density_up=0;
        splay_grad_square_energy_density_down=0;
        splay_grad_square_energy_up=0;
        splay_grad_square_energy_down=0;

    end

    if physical_properties.Splay_grad_tilt_modulus~=0
        splay_grad_tilt_energy_density_up=0.5*physical_properties.Splay_grad_tilt_modulus*splay.splay_grad_up.*tilt.up.rho;
        splay_grad_tilt_energy_density_down=0.5*physical_properties.Splay_grad_tilt_modulus*splay.splay_grad_down.*tilt.down.rho;
        splay_grad_tilt_energy_up=4*sum(splay_grad_tilt_energy_density_up.*area_element.up,'all');
        splay_grad_tilt_energy_down=4*sum(splay_grad_tilt_energy_density_down.*area_element.down,'all');
    else
        splay_grad_tilt_energy_up=0;
        splay_grad_tilt_energy_down=0;
    end


    total_energy=tilt_energy_up+tilt_energy_down+splay_energy_up+splay_energy_down+saddle_splay_energy_up...
            +saddle_splay_energy_down+splay_grad_square_energy_up+splay_grad_square_energy_down...
            +splay_grad_tilt_energy_up+splay_grad_tilt_energy_down;
    
    energy_components.Energy_density=0.5*physical_properties.kappa_tilt*(tilt.up.x.^2+tilt.up.y.^2+tilt.up.z.^2)+...
        0.5*physical_properties.kappa_tilt*(tilt.down.x.^2+tilt.down.y.^2+tilt.down.z.^2)...
        +0.5*physical_properties.kappa_up.*(splay.up.^2-2*splay.up*physical_properties.J0_up) ...
        +0.5*physical_properties.kappa_down.*(splay.down.^2-2*splay.down*physical_properties.J0_down) ...
        +physical_properties.kappa_bar_up.*saddle_splay.up + physical_properties.kappa_bar_down.*saddle_splay.down+...
        +splay_grad_square_energy_density_up+splay_grad_square_energy_density_down;
    
    
    energy_components.splay_energy=splay_energy_up+splay_energy_down;
    energy_components.saddle_splay_energy=saddle_splay_energy_up+saddle_splay_energy_down;
    energy_components.tilt_energy=tilt_energy_up+tilt_energy_down;
    energy_components.splay_grad_square_energy=splay_grad_square_energy_up+splay_grad_square_energy_down;
    energy_components.splay_grad_tilt_energy=splay_grad_tilt_energy_up+splay_grad_tilt_energy_down;
    
end