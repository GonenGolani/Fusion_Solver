function [x_arr,y_arr,z_arr] = mid_plane_surface_creator(coeff_matrix,inner_rim,outer_rim,res_struc)
% coeff_matrix is the coefficiant matrix for that describe z=sum[h_l(phi)*rho^l]
% inner_rim and outer_rim is structure containing rho at the x and y
% direction x0 and y0
% output: x,y,z coordinates of the midplane

    
    phi_res=res_struc.phi_res;
    rho_res=res_struc.rho_res;
    
    phi_vec=linspace(0,pi/2,phi_res);
    
    % set size of arries
    x_arr=zeros(phi_res,rho_res);
    y_arr=zeros(phi_res,rho_res);
    z_arr=zeros(phi_res,rho_res);

  
    int=1;
    while int<=phi_res %run over all phi
        local_phi=phi_vec(int); %find local phi
        %find  minimal and maximal rho for each phi
        %inner loop dimensions
        if (inner_rim.x0==0 || inner_rim.y0==0) % if we are at rho=0 
            rho_min=0;  % set rho to 0
        else  %if else find it based on ellipse 
            mahana_min=(inner_rim.x0^2.*sin(local_phi).^2+inner_rim.y0^2.*cos(local_phi).^2).^0.5;
            rho_min=inner_rim.x0*inner_rim.y0./mahana_min;
        end
        %outer loop rho for each phi
        mahana_max=(outer_rim.x0^2.*sin(local_phi).^2+outer_rim.y0^2.*cos(local_phi).^2).^0.5;
        rho_max=outer_rim.x0*outer_rim.y0./mahana_max;
        
        rho_vec=linspace(rho_min,rho_max,rho_res); % create vector containing rho
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

