function [x_arr,y_arr,z_arr] = mid_plane_surface_creator_genral_rim(coeff_matrix,inner_rim,outer_rim,phi_res,rho_res)
% coeff_matrix is the coefficiant matrix for that describe z=sum[h_l(phi)*rho^l]
% inner_rim and outer_rim is described by the genral function
% inner_rim(phi) and outer_rim(phi)
% phi_res is the angular resolution
% rho_res is the radial resolution

% output: x,y,z coordinates of the midplane


    
    phi_vec=linspace(0,pi/2,phi_res);
    
    % set size of arries
    x_arr=zeros(phi_res,rho_res);
    y_arr=zeros(phi_res,rho_res);
    z_arr=zeros(phi_res,rho_res);

  
    int=1;
    while int<=phi_res %run over all phi
        local_phi=phi_vec(int); %find local phi
        
        
        
        rho_vec=linspace(inner_rim(int),outer_rim(int),rho_res); % create vector containing rho
        x_arr(int,:)=rho_vec*cos(local_phi); % project on x axis
        y_arr(int,:)=rho_vec*sin(local_phi); % project on y axis
        
        if int==phi_res %becuse of rounding error check if phi=pi/2 set cos to 0
            x_arr(int,:)=rho_vec*0; % project on x axis
            y_arr(int,:)=rho_vec*1; % project on y axis
        end
        z_arr(int,:)=vector_my_poly_val(coeff_matrix(int,:),rho_vec); %find z vlaue and each rho for given phi
        

        int=int+1;
        
    end

    
end

