function  [dividing_plane,tilt,lipid_length,area_element,splay,saddle_splay]=find_membrane_geometry...
        (mid_plane,lipid_director,physical_properties,phi_res,rho_res)
% this function finds the membrane normal tilt splay saddle splay and twist 
% input: mid_plane - structure containing array of x y and z position of  mid-plane
%        lipid_director structure containing array of nx ny and nz in up
%        and down direction
%        normal- normal to mid plane pointing up

% output - 
%          dividing_plane position up and down, area_element, metric
%          normal - normal defined as ex X ey / abs(ex X ey) of the mid-plane
%          tilt - structure lick lipid_director but of tilt vector
%          lipid tail length -lipid_length in up and down
%          area_element - area calcualted using detrminant of metric
%          dA=dx*dy*(1+dzdx^2+dzdy^2)^0.5
%          splay - trace of splay tensor
%          saddle_splay - det of splay tensor

%find normal
%check if flat




%normelize lipid director
lipid_director_up_length=(lipid_director.up.x.^2+lipid_director.up.y.^2+lipid_director.up.z.^2).^0.5;
lipid_director_down_length=(lipid_director.down.x.^2+lipid_director.down.y.^2+lipid_director.down.z.^2).^0.5;
lipid_director.up.x=lipid_director.up.x./lipid_director_up_length;
lipid_director.up.y=lipid_director.up.y./lipid_director_up_length;
lipid_director.up.z=lipid_director.up.z./lipid_director_up_length;
lipid_director.down.x=lipid_director.down.x./lipid_director_down_length;
lipid_director.down.y=lipid_director.down.y./lipid_director_down_length;
lipid_director.down.z=lipid_director.down.z./lipid_director_down_length;

%tilt
%product of lipid director and normal
Nn_up=mid_plane.normal.x.*lipid_director.up.x+mid_plane.normal.y.*lipid_director.up.y+mid_plane.normal.z.*lipid_director.up.z;
Nn_down=-mid_plane.normal.x.*lipid_director.down.x-mid_plane.normal.y.*lipid_director.down.y-mid_plane.normal.z.*lipid_director.down.z;


tilt.up.x=lipid_director.up.x./Nn_up-mid_plane.normal.x; %tilt in x up monolayer
tilt.up.y=lipid_director.up.y./Nn_up-mid_plane.normal.y; %tilt in y up monolayer
tilt.up.z=lipid_director.up.z./Nn_up-mid_plane.normal.z; %tilt in z up monolayer



tilt.down.x=lipid_director.down.x./Nn_down+mid_plane.normal.x; %tilt in x down monolayer
tilt.down.y=lipid_director.down.y./Nn_down+mid_plane.normal.y; %tilt in y down monolayer
tilt.down.z=lipid_director.down.z./Nn_down+mid_plane.normal.z; %tilt in z down monolayer



lipid_length.up=physical_properties.lipid_length*(1+tilt.up.x.^2+tilt.up.y.^2+tilt.up.z.^2).^0.5; %length of lipid tail up monolayer
lipid_length.down=physical_properties.lipid_length*(1+tilt.down.x.^2+tilt.down.y.^2+tilt.down.z.^2).^0.5; %length of lipid tail up monolayer


%% dividing plane diivding plane (DP) is MD+lipid length*lipid director
dividing_plane.up.x=mid_plane.x+lipid_director.up.x.*lipid_length.up; 
dividing_plane.up.y=mid_plane.y+lipid_director.up.y.*lipid_length.up; 
dividing_plane.up.z=mid_plane.z+lipid_director.up.z.*lipid_length.up; 

dividing_plane.down.x=mid_plane.x+lipid_director.down.x.*lipid_length.down; 
dividing_plane.down.y=mid_plane.y+lipid_director.down.y.*lipid_length.down; 
dividing_plane.down.z=mid_plane.z+lipid_director.down.z.*lipid_length.down; 


%% find length elements
dphi=(pi/2)/phi_res;



dividing_plane.up.rho=(dividing_plane.up.x.^2+dividing_plane.up.y.^2).^0.5;
dividing_plane.up.cosphi=dividing_plane.up.x./dividing_plane.up.rho;
dividing_plane.up.sinphi=dividing_plane.up.y./dividing_plane.up.rho;

dividing_plane.down.rho=(dividing_plane.down.x.^2+dividing_plane.down.y.^2).^0.5;
dividing_plane.down.cosphi=dividing_plane.down.x./dividing_plane.down.rho;
dividing_plane.down.sinphi=dividing_plane.down.y./dividing_plane.down.rho;

dividing_plane.up.drho=my_diff(dividing_plane.up.rho,'x',phi_res,rho_res);
dividing_plane.up.rho_dphi=dividing_plane.up.rho*dphi;

dividing_plane.down.drho=my_diff(dividing_plane.down.rho,'x',phi_res,rho_res);
dividing_plane.down.rho_dphi=dividing_plane.down.rho*dphi;

[dividing_plane.up.dz_direction_rho,dividing_plane.up.dz_direction_phi]=gradient(dividing_plane.up.z);
[dividing_plane.down.dz_direction_rho,dividing_plane.down.dz_direction_phi]=gradient(dividing_plane.down.z);


dividing_plane.up.dhdrho=dividing_plane.up.dz_direction_rho./dividing_plane.up.drho;
dividing_plane.up.dhdrhophi=dividing_plane.up.dz_direction_phi./dividing_plane.up.rho_dphi;

dividing_plane.down.dhdrho=dividing_plane.down.dz_direction_rho./dividing_plane.down.drho;
dividing_plane.down.dhdrhophi=dividing_plane.down.dz_direction_phi./dividing_plane.down.rho_dphi;


dividing_plane.up.dhdx=dividing_plane.up.dhdrho.*dividing_plane.up.cosphi-dividing_plane.up.dhdrhophi.*dividing_plane.up.sinphi;
dividing_plane.up.dhdy=dividing_plane.up.dhdrho.*dividing_plane.up.sinphi+dividing_plane.up.dhdrhophi.*dividing_plane.up.cosphi;

dividing_plane.down.dhdx=dividing_plane.down.dhdrho.*dividing_plane.down.cosphi-dividing_plane.down.dhdrhophi.*dividing_plane.down.sinphi;
dividing_plane.down.dhdy=dividing_plane.down.dhdrho.*dividing_plane.down.sinphi+dividing_plane.down.dhdrhophi.*dividing_plane.down.cosphi;



%% metric
%detrminanat
dividing_plane.up.metric_det=(1+dividing_plane.up.dhdx.^2+dividing_plane.up.dhdy.^2);
dividing_plane.down.metric_det=(1+dividing_plane.down.dhdx.^2+dividing_plane.down.dhdy.^2);

% covarient metric elements
dividing_plane.up.metric.xx=1+dividing_plane.up.dhdx.^2;
dividing_plane.up.metric.xy=dividing_plane.up.dhdx.*dividing_plane.up.dhdy;
dividing_plane.up.metric.yy=1+dividing_plane.up.dhdy.^2;

dividing_plane.down.metric.xx=1+dividing_plane.down.dhdx.^2;
dividing_plane.down.metric.xy=dividing_plane.down.dhdx.*dividing_plane.down.dhdy;
dividing_plane.down.metric.yy=1+dividing_plane.down.dhdy.^2;


%contra-varient elements
dividing_plane.up.inv_metric.xx=dividing_plane.up.metric.yy./dividing_plane.up.metric_det;
dividing_plane.up.inv_metric.yy=dividing_plane.up.metric.xx./dividing_plane.up.metric_det;
dividing_plane.up.inv_metric.xy=-dividing_plane.up.dhdx.*dividing_plane.up.dhdy./dividing_plane.up.metric_det;

dividing_plane.down.inv_metric.xx=dividing_plane.down.metric.yy./dividing_plane.down.metric_det;
dividing_plane.down.inv_metric.yy=dividing_plane.down.metric.xx./dividing_plane.down.metric_det;
dividing_plane.down.inv_metric.xy=-dividing_plane.down.dhdx.*dividing_plane.up.dhdy./dividing_plane.down.metric_det;


%% 

[lipid_director.up.dnx_direction_rho,lipid_director.up.dnx_direction_phi]=gradient(lipid_director.up.x);
[lipid_director.up.dny_direction_rho,lipid_director.up.dny_direction_phi]=gradient(lipid_director.up.y);
[lipid_director.up.dnz_direction_rho,lipid_director.up.dnz_direction_phi]=gradient(lipid_director.up.z);

[lipid_director.down.dnx_direction_rho,lipid_director.down.dnx_direction_phi]=gradient(lipid_director.down.x);
[lipid_director.down.dny_direction_rho,lipid_director.down.dny_direction_phi]=gradient(lipid_director.down.y);
[lipid_director.down.dnz_direction_rho,lipid_director.down.dnz_direction_phi]=gradient(lipid_director.down.z);

lipid_director.up.dnxdrho=lipid_director.up.dnx_direction_rho./dividing_plane.up.drho;
lipid_director.up.dnxdrhophi=lipid_director.up.dnx_direction_phi./dividing_plane.up.rho_dphi;
lipid_director.up.dnydrho=lipid_director.up.dny_direction_rho./dividing_plane.up.drho;
lipid_director.up.dnydrhophi=lipid_director.up.dny_direction_phi./dividing_plane.up.rho_dphi;
lipid_director.up.dnzdrho=lipid_director.up.dnz_direction_rho./dividing_plane.up.drho;
lipid_director.up.dnzdrhophi=lipid_director.up.dnz_direction_phi./dividing_plane.up.rho_dphi;

lipid_director.down.dnxdrho=lipid_director.down.dnx_direction_rho./dividing_plane.down.drho;
lipid_director.down.dnxdrhophi=lipid_director.down.dnx_direction_phi./dividing_plane.down.rho_dphi;
lipid_director.down.dnydrho=lipid_director.down.dny_direction_rho./dividing_plane.down.drho;
lipid_director.down.dnydrhophi=lipid_director.down.dny_direction_phi./dividing_plane.down.rho_dphi;
lipid_director.down.dnzdrho=lipid_director.down.dnz_direction_rho./dividing_plane.down.drho;
lipid_director.down.dnzdrhophi=lipid_director.down.dnz_direction_phi./dividing_plane.down.rho_dphi;

lipid_director.up.dnxdx=lipid_director.up.dnxdrho.*dividing_plane.up.cosphi-lipid_director.up.dnxdrhophi.*dividing_plane.up.sinphi;
lipid_director.up.dnydy=lipid_director.up.dnydrho.*dividing_plane.up.sinphi+lipid_director.up.dnydrhophi.*dividing_plane.up.cosphi;
lipid_director.up.dnzdx=lipid_director.up.dnzdrho.*dividing_plane.up.cosphi-lipid_director.up.dnzdrhophi.*dividing_plane.up.sinphi;
lipid_director.up.dnzdy=lipid_director.up.dnzdrho.*dividing_plane.up.sinphi+lipid_director.up.dnzdrhophi.*dividing_plane.up.cosphi;


lipid_director.down.dnxdx=lipid_director.down.dnxdrho.*dividing_plane.down.cosphi-lipid_director.down.dnxdrhophi.*dividing_plane.down.sinphi;
lipid_director.down.dnydy=lipid_director.down.dnydrho.*dividing_plane.down.sinphi+lipid_director.down.dnydrhophi.*dividing_plane.down.cosphi;
lipid_director.down.dnzdx=lipid_director.down.dnzdrho.*dividing_plane.down.cosphi-lipid_director.down.dnzdrhophi.*dividing_plane.down.sinphi;
lipid_director.down.dnzdy=lipid_director.down.dnzdrho.*dividing_plane.down.sinphi+lipid_director.down.dnzdrhophi.*dividing_plane.down.cosphi;

%%
%lipid splay tensor
splay_tensor.up.xx=lipid_director.up.dnxdx+dividing_plane.up.dhdx.*lipid_director.up.dnzdx;
splay_tensor.up.xy=dividing_plane.up.dhdx.*lipid_director.up.dnzdy;
splay_tensor.up.yy=lipid_director.up.dnydy+dividing_plane.up.dhdy.*lipid_director.up.dnzdy;
splay_tensor.up.yx=dividing_plane.up.dhdy.*lipid_director.up.dnzdx;

splay_tensor.down.xx=lipid_director.down.dnxdx+dividing_plane.down.dhdx.*lipid_director.down.dnzdx;
splay_tensor.down.xy=dividing_plane.down.dhdx.*lipid_director.down.dnzdy;
splay_tensor.down.yy=lipid_director.down.dnydy+dividing_plane.down.dhdy.*lipid_director.down.dnzdy;
splay_tensor.down.yx=dividing_plane.down.dhdy.*lipid_director.down.dnzdx;
%splay and saddle_splay
splay.up=splay_tensor.up.xx.*dividing_plane.up.inv_metric.xx+splay_tensor.up.xy.*dividing_plane.up.inv_metric.xy...
    +splay_tensor.up.yx.*dividing_plane.up.inv_metric.xy+splay_tensor.up.yy.*dividing_plane.up.inv_metric.yy;

splay.down=splay_tensor.down.xx.*dividing_plane.down.inv_metric.xx+splay_tensor.down.xy.*dividing_plane.down.inv_metric.xy...
    +splay_tensor.down.yx.*dividing_plane.down.inv_metric.xy+splay_tensor.down.yy.*dividing_plane.down.inv_metric.yy;

saddle_splay.up=(splay_tensor.up.xx.*splay_tensor.up.yy-splay_tensor.up.yx.*splay_tensor.up.yx)./dividing_plane.up.metric_det;

saddle_splay.down=(splay_tensor.down.xx.*splay_tensor.down.yy-splay_tensor.down.yx.*splay_tensor.down.yx)./dividing_plane.down.metric_det;

[~,~,~,area_element.up] = my_surfnorm(dividing_plane.up.x,dividing_plane.up.y,dividing_plane.up.z);
[~,~,~,area_element.down] = my_surfnorm(dividing_plane.down.x,dividing_plane.down.y,dividing_plane.down.z);


%%
if physical_properties.Splay_grad_square_modulus~=0 || physical_properties.Splay_grad_tilt_modulus~=0
    
    [splay_change_rho_up,~]=gradient(splay.up);
    [splay_change_rho_down,~]=gradient(splay.down);

    splay.splay_grad_up=(1+dividing_plane.up.dhdrho.^2).^-0.5.*splay_change_rho_up./dividing_plane.up.drho;
    splay.splay_grad_down=(1+dividing_plane.down.dhdrho.^2).^-0.5.*splay_change_rho_down./dividing_plane.down.drho;

end
if physical_properties.Splay_grad_tilt_modulus~=0

    tilt.up.rho=dividing_plane.up.cosphi.*tilt.up.x+dividing_plane.up.sinphi.*tilt.up.y;    %tilt in rho up monolayer
    tilt.up.phi=-dividing_plane.up.sinphi.*tilt.up.x+dividing_plane.up.cosphi.*tilt.up.y;   %tilt in phi up monolayer

    tilt.down.rho=dividing_plane.down.cosphi.*tilt.down.x+dividing_plane.down.sinphi.*tilt.down.y;  %tilt in rho down monolayer
    tilt.down.phi=-dividing_plane.down.sinphi.*tilt.down.x+dividing_plane.down.cosphi.*tilt.down.y; %tilt in phi down monolayer

    
end




end
