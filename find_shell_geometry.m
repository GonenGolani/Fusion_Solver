function  [normal,total_curvature,Gaussian_curvature,shell_mid_plane]=find_shell_geometry...
        (shell_mid_plane,phi_res,rho_res)
% input: membrane_mid_plane - structure containing array of x y and z position of the original membrane mid-plane
%        distance_from_membrane - structure cotnaining the distnace of the shell from the mid-plane


% output - 
%          shell_mid_plane position up and down, area_element, metric
%          normal - normal defined as ex X ey / abs(ex X ey) of the mid-plane
%          area_element - area calcualted using detrminant of metric
%          dA=dx*dy*(1+dzdx^2+dzdy^2)^0.5
%          total_curvature - trace of curvature tensor
%          Gaussian_curvature - det of curvature tensor


[normal.x,normal.y,normal.z,shell_mid_plane.area_element] = my_surfnorm(shell_mid_plane.x,shell_mid_plane.y,shell_mid_plane.z);

norm_length=(normal.x.^2+normal.y.^2+normal.z.^2).^0.5;
normal.x=normal.x./norm_length;
normal.y=normal.y./norm_length;
normal.z=normal.z./norm_length;

%flip normal if needed needs to point down!
if normal.z(1,1)>0
    normal.x=normal.x*-1;
    normal.y=normal.y*-1;
    normal.z=normal.z*-1;   
end

%% find length elements
dphi=(pi/2)/phi_res;



shell_mid_plane.rho=(shell_mid_plane.x.^2+shell_mid_plane.y.^2).^0.5;
shell_mid_plane.cosphi=shell_mid_plane.x./shell_mid_plane.rho;
shell_mid_plane.sinphi=shell_mid_plane.y./shell_mid_plane.rho;


shell_mid_plane.drho=my_diff(shell_mid_plane.rho,'x',phi_res,rho_res);
shell_mid_plane.rho_dphi=shell_mid_plane.rho*dphi;


[shell_mid_plane.dz_direction_rho,shell_mid_plane.dz_direction_phi]=gradient(shell_mid_plane.z);

shell_mid_plane.dhdrho=shell_mid_plane.dz_direction_rho./shell_mid_plane.drho;
shell_mid_plane.dhdrhophi=shell_mid_plane.dz_direction_phi./shell_mid_plane.rho_dphi;


shell_mid_plane.dhdx=shell_mid_plane.dhdrho.*shell_mid_plane.cosphi-shell_mid_plane.dhdrhophi.*shell_mid_plane.sinphi;
shell_mid_plane.dhdy=shell_mid_plane.dhdrho.*shell_mid_plane.sinphi+shell_mid_plane.dhdrhophi.*shell_mid_plane.cosphi;


%% metric
%detrminanat
shell_mid_plane.metric_det=(1+shell_mid_plane.dhdx.^2+shell_mid_plane.dhdy.^2);

% covarient metric elements
shell_mid_plane.metric.xx=1+shell_mid_plane.dhdx.^2;
shell_mid_plane.metric.xy=shell_mid_plane.dhdx.*shell_mid_plane.dhdy;
shell_mid_plane.metric.yy=1+shell_mid_plane.dhdy.^2;


%contra-varient elements
shell_mid_plane.inv_metric.xx=shell_mid_plane.metric.yy./shell_mid_plane.metric_det;
shell_mid_plane.inv_metric.yy=shell_mid_plane.metric.xx./shell_mid_plane.metric_det;
shell_mid_plane.inv_metric.xy=-shell_mid_plane.dhdx.*shell_mid_plane.dhdy./shell_mid_plane.metric_det;



[normal.dnx_direction_rho,normal.dnx_direction_phi]=gradient(normal.x);
[normal.dny_direction_rho,normal.dny_direction_phi]=gradient(normal.y);
[normal.dnz_direction_rho,normal.dnz_direction_phi]=gradient(normal.z);

normal.dnxdrho=normal.dnx_direction_rho./shell_mid_plane.drho;
normal.dnxdrhophi=normal.dnx_direction_phi./shell_mid_plane.rho_dphi;
normal.dnydrho=normal.dny_direction_rho./shell_mid_plane.drho;
normal.dnydrhophi=normal.dny_direction_phi./shell_mid_plane.rho_dphi;
normal.dnzdrho=normal.dnz_direction_rho./shell_mid_plane.drho;
normal.dnzdrhophi=normal.dnz_direction_phi./shell_mid_plane.rho_dphi;

normal.dnxdx=normal.dnxdrho.*shell_mid_plane.cosphi-normal.dnxdrhophi.*shell_mid_plane.sinphi;
normal.dnydy=normal.dnydrho.*shell_mid_plane.sinphi+normal.dnydrhophi.*shell_mid_plane.cosphi;
normal.dnzdx=normal.dnzdrho.*shell_mid_plane.cosphi-normal.dnzdrhophi.*shell_mid_plane.sinphi;
normal.dnzdy=normal.dnzdrho.*shell_mid_plane.sinphi+normal.dnzdrhophi.*shell_mid_plane.cosphi;



%%
%% Curvature tensor
curvature_tensor.xx=normal.dnxdx+shell_mid_plane.dhdx.*normal.dnzdx;
curvature_tensor.xy=shell_mid_plane.dhdx.*normal.dnzdy;
curvature_tensor.yy=normal.dnydy+shell_mid_plane.dhdy.*normal.dnzdy;
curvature_tensor.yx=shell_mid_plane.dhdy.*normal.dnzdx;


%Total and Gaussian
total_curvature=curvature_tensor.xx.*shell_mid_plane.inv_metric.xx+curvature_tensor.xy.*shell_mid_plane.inv_metric.xy...
    +curvature_tensor.yx.*shell_mid_plane.inv_metric.xy+curvature_tensor.yy.*shell_mid_plane.inv_metric.yy;


Gaussian_curvature=(curvature_tensor.xx.*curvature_tensor.yy-curvature_tensor.yx.*curvature_tensor.yx)./shell_mid_plane.metric_det;



end
