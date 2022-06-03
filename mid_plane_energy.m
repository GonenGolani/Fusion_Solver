function [energy] = mid_plane_energy(mid_plane)
    
    %this function finds the energy assoacited with the mid plane 
    % abbreviations
    area_element=mid_plane.area_element;
    total_curvature=mid_plane.Total_curvature;
    Gaussian_curvature=mid_plane.Gaussian_curvature;
    
    kappa=mid_plane.mid_plane_physical_proprties.kappa;
    kappa_bar=mid_plane.mid_plane_physical_proprties.kappa_bar;
    J0=mid_plane.mid_plane_physical_proprties.J0;
  
    
    %the 4 factor is becuse we calc only qurter of system!
    
    mid_plane.total_area=4*sum(area_element,'all');
    
    
    mid_plane.bending_energy=4*sum(0.5*kappa.*(total_curvature-J0).^2.*area_element,'all');
    
    mid_plane.Saddle_splay_energy=4*sum(kappa_bar.*Gaussian_curvature.*area_element,'all');
 
    energy=mid_plane.bending_energy+mid_plane.Saddle_splay_energy;
    
end

