function [total_curvature,Gaussian_curvature] = mid_plane_geomtry(mid_plane,res_struc)

% input: membrane_mid_plane - structure containing array of x y and z position of the original membrane mid-plane
%        distance_from_membrane - structure cotnaining the distnace of the shell from the mid-plane


% output - 
%          mid_plane position up and down, area_element, metric
%          normal - normal defined as ex X ey / abs(ex X ey) of the mid-plane
%          area_element - area calcualted using detrminant of metric
%          dA=dx*dy*(1+dzdx^2+dzdy^2)^0.5
%          total_curvature - trace of curvature tensor
%          Gaussian_curvature - det of curvature tensor

    normal=mid_plane.normal;
    norm_length=(normal.x.^2+normal.y.^2+normal.z.^2).^0.5;
    normal.x=normal.x./norm_length;
    normal.y=normal.y./norm_length;
    normal.z=normal.z./norm_length;

    %flip normal if needed needs to point down!
    if normal.z(1,1)<0
        normal.x=normal.x*-1;
        normal.y=normal.y*-1;
        normal.z=normal.z*-1;   
    end

    %% find length elements
    dphi=(pi/2)/res_struc.phi_res;



    mid_plane.rho=(mid_plane.x.^2+mid_plane.y.^2).^0.5;
    mid_plane.cosphi=mid_plane.x./mid_plane.rho;
    mid_plane.sinphi=mid_plane.y./mid_plane.rho;


    mid_plane.drho=my_diff(mid_plane.rho,'x',res_struc.phi_res,res_struc.rho_res);
    mid_plane.rho_dphi=mid_plane.rho*dphi;


    [mid_plane.dz_direction_rho,mid_plane.dz_direction_phi]=gradient(mid_plane.z);

    mid_plane.dhdrho=mid_plane.dz_direction_rho./mid_plane.drho;
    mid_plane.dhdrhophi=mid_plane.dz_direction_phi./mid_plane.rho_dphi;


    mid_plane.dhdx=mid_plane.dhdrho.*mid_plane.cosphi-mid_plane.dhdrhophi.*mid_plane.sinphi;
    mid_plane.dhdy=mid_plane.dhdrho.*mid_plane.sinphi+mid_plane.dhdrhophi.*mid_plane.cosphi;


    %% metric
    %detrminanat
    mid_plane.metric_det=(1+mid_plane.dhdx.^2+mid_plane.dhdy.^2);

    % covarient metric elements
    mid_plane.metric.xx=1+mid_plane.dhdx.^2;
    mid_plane.metric.xy=mid_plane.dhdx.*mid_plane.dhdy;
    mid_plane.metric.yy=1+mid_plane.dhdy.^2;


    %contra-varient elements
    mid_plane.inv_metric.xx=mid_plane.metric.yy./mid_plane.metric_det;
    mid_plane.inv_metric.yy=mid_plane.metric.xx./mid_plane.metric_det;
    mid_plane.inv_metric.xy=-mid_plane.dhdx.*mid_plane.dhdy./mid_plane.metric_det;



    [normal.dnx_direction_rho,normal.dnx_direction_phi]=gradient(normal.x);
    [normal.dny_direction_rho,normal.dny_direction_phi]=gradient(normal.y);
    [normal.dnz_direction_rho,normal.dnz_direction_phi]=gradient(normal.z);

    normal.dnxdrho=normal.dnx_direction_rho./mid_plane.drho;
    normal.dnxdrhophi=normal.dnx_direction_phi./mid_plane.rho_dphi;
    normal.dnydrho=normal.dny_direction_rho./mid_plane.drho;
    normal.dnydrhophi=normal.dny_direction_phi./mid_plane.rho_dphi;
    normal.dnzdrho=normal.dnz_direction_rho./mid_plane.drho;
    normal.dnzdrhophi=normal.dnz_direction_phi./mid_plane.rho_dphi;

    normal.dnxdx=normal.dnxdrho.*mid_plane.cosphi-normal.dnxdrhophi.*mid_plane.sinphi;
    normal.dnydy=normal.dnydrho.*mid_plane.sinphi+normal.dnydrhophi.*mid_plane.cosphi;
    normal.dnzdx=normal.dnzdrho.*mid_plane.cosphi-normal.dnzdrhophi.*mid_plane.sinphi;
    normal.dnzdy=normal.dnzdrho.*mid_plane.sinphi+normal.dnzdrhophi.*mid_plane.cosphi;



    %%
    %% Curvature tensor
    curvature_tensor.xx=normal.dnxdx+mid_plane.dhdx.*normal.dnzdx;
    curvature_tensor.xy=mid_plane.dhdx.*normal.dnzdy;
    curvature_tensor.yy=normal.dnydy+mid_plane.dhdy.*normal.dnzdy;
    curvature_tensor.yx=mid_plane.dhdy.*normal.dnzdx;


    %Total and Gaussian
    total_curvature=curvature_tensor.xx.*mid_plane.inv_metric.xx+curvature_tensor.xy.*mid_plane.inv_metric.xy...
        +curvature_tensor.yx.*mid_plane.inv_metric.xy+curvature_tensor.yy.*mid_plane.inv_metric.yy;


    Gaussian_curvature=(curvature_tensor.xx.*curvature_tensor.yy-curvature_tensor.yx.*curvature_tensor.yx)./mid_plane.metric_det;



end

