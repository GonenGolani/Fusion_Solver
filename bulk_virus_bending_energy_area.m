function [Elastic_energy,distal_area,proximal_area,virus_resting_energy] = bulk_virus_bending_energy_area(Virus,res_struc)


%find the size of the fusion site boundary

% Fusion site in the rho phi system
phi=linspace(0,pi/2,res_struc.phi_res);
rho=linspace(0.0001,Virus.Rv,res_struc.rho_res);

[RHO,PHI]=meshgrid(rho,phi);


% from the parametric description
% helping functions
sinV=RHO.*sin(PHI)./Virus.r_sv;
cosV=(1-sinV.^2).^0.5; % can be +- for inner or outer intersection;
y_star=(Virus.r_bv+Virus.r_sv.*cosV); %just to short things
cosU=RHO.*cos(PHI)./y_star;
sinU=-(1-cosU.^2).^0.5;

%the surface
mid_plane.x=RHO.*cos(PHI);
mid_plane.y=RHO.*sin(PHI);
mid_plane.z=(Virus.r_bv+Virus.r_sv.*cosV).*sinU;

[mid_plane.normal.x,mid_plane.normal.y,mid_plane.normal.z] = surfnorm(mid_plane.x,mid_plane.y,mid_plane.z);

lipid_director.up.x=mid_plane.normal.x;
lipid_director.up.y=mid_plane.normal.y;
lipid_director.up.z=mid_plane.normal.z;

lipid_director.down.x=-mid_plane.normal.x;
lipid_director.down.y=-mid_plane.normal.y;
lipid_director.down.z=-mid_plane.normal.z;


% geometry
[~,~,~,area_element,splay,saddle_splay]=...
    find_membrane_geometry(mid_plane,lipid_director,Virus.physical_properties,res_struc.phi_res,res_struc.rho_res);



%elastic energy and area in the missing fusion site
Elastic_energy_density_distal=0.5*Virus.physical_properties.kappa_up*(splay.up.^2-2*splay.up*Virus.physical_properties.J0_up)+... 
    Virus.physical_properties.kappa_bar_up*saddle_splay.up;
Elastic_energy_density_proxy=0.5*Virus.physical_properties.kappa_down*(splay.down.^2-2*splay.down*Virus.physical_properties.J0_down)+... 
    Virus.physical_properties.kappa_bar_down*saddle_splay.down;



distal_elastic_energy_in=4*sum(Elastic_energy_density_distal.*area_element.up,'all');
proxy_elastic_energy_in=4*sum(Elastic_energy_density_proxy.*area_element.down,'all');

distal_area_in=4*sum(area_element.up,'all');
proximal_area_in=4*sum(area_element.down,'all');

%elastic energy and area in the rest of the torus

rbv=Virus.r_bv;




rsv_distal=Virus.r_sv-Virus.physical_properties.lipid_length;
elastic_energy_distal_all=2*pi^2*Virus.physical_properties.kappa_up*...
    ((rbv^4/(rbv^2*rsv_distal^2-rsv_distal^4))^0.5+rbv*Virus.physical_properties.J0_up);

area_distal_all=4*pi^2*rsv_distal*rbv;


rsv_proxy=Virus.r_sv+Virus.physical_properties.lipid_length;
elastic_energy_proxy_all=2*pi^2*Virus.physical_properties.kappa_down*...
    ((rbv^4/(rbv^2*rsv_proxy^2-rsv_proxy^4))^0.5-rbv*Virus.physical_properties.J0_down);

area_proxy_all=4*pi^2*rsv_proxy*rbv;

% all inus what is in

Elastic_energy=elastic_energy_proxy_all+elastic_energy_distal_all-proxy_elastic_energy_in-distal_elastic_energy_in;

distal_area=area_distal_all-distal_area_in;

proximal_area=area_proxy_all-proximal_area_in;


virus_resting_energy=elastic_energy_proxy_all+elastic_energy_distal_all;
end



