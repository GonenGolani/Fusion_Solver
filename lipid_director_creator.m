function [n_x,n_y,n_z,n_rho,n_phi] = lipid_director_creator(rho_angle_matrix,phi_angle_matrix,Gaussian_decay_vector,mid_plane,inner_rim,outer_rim,res_struct,up_down)
% rho_angle_matrix is the coefficiant matrix for that describe theta_rho=sum[theta_l(phi)*rho^l]
% phi_angle_matrix is the coefficiant matrix for that describe theta_phi=sum[theta_l(rho)*cos(2*l*phi)]
% inner_rim and outer_rim is structure containing rho at the x and y
% up_down indicate the direction the direction in z 
% Gaussian_decay_vector is a structure contening 1) k_vector in the length phi_res with coeffciants that will be used for the Gaussian.
% 2) vector conteining the x0 and y0 of the diphragm rim

% direction x0 and y0 are the ellipse coordinets in the x and y directions
% output: x,y,z coordinates of the midplane

    
    phi_res=res_struct.phi_res;
    rho_res=res_struct.rho_res;
    
    phi_vec=linspace(0,pi/2,phi_res);
    
    % set size of arries
    theta_rho=zeros(phi_res,rho_res);
    theta_phi=zeros(phi_res,rho_res);
    n_rho=zeros(phi_res,rho_res);
    n_phi=zeros(phi_res,rho_res);
    
    
    % surface normal is polar coordintes
    RHO=(mid_plane.x.^2+mid_plane.y.^2).^0.5;
    sinphi=mid_plane.y./RHO;
    cosphi=mid_plane.x./RHO;
    
    N_rho=mid_plane.normal.rho;
    N_phi=mid_plane.normal.phi;
    N_z=mid_plane.normal.z;
    
    if N_z(1,1)<0 && up_down==1
        N_rho=-N_rho;
        N_phi=-N_phi;
        N_z=-N_z;
    end
    if N_z(1,1)>0 && up_down==-1
        N_rho=-N_rho;
        N_phi=-N_phi;
        N_z=-N_z;
    end    
    
    %% create theta_rho
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
        theta_rho(int,:)=vector_my_poly_val(rho_angle_matrix(int,:),rho_vec); %find angle between normal to rho direction
        theta_phi(int,:)=vector_my_poly_val(phi_angle_matrix(int,:),rho_vec); %find angle between normal to phi direction
        
        if Gaussian_decay_vector.exponential_decay_tilt==1 % Apply exponential decay to theta_rho
            
            mahana_diaphragm_rim=(Gaussian_decay_vector.x0^2.*sin(local_phi).^2+Gaussian_decay_vector.y0^2.*cos(local_phi).^2).^0.5;
            rho_rim=Gaussian_decay_vector.x0*Gaussian_decay_vector.y0./mahana_diaphragm_rim;
            theta_rho(int,:)=exp(-abs(Gaussian_decay_vector.k_vector(int)).*(rho_vec-rho_rim).^2).*theta_rho(int,:);
        end
        int=int+1;
        
    end
    

    n_rho=cos(theta_rho).*N_rho-sin(theta_rho).*N_z;        %lipid director projection on rho direction cos(theta+acos(N_rho))
    n_z=sin(theta_rho).*N_rho+cos(theta_rho).*N_z;          %lipid director projection on z direction sin(theta+acos(N_rho))
    n_phi=cos(theta_phi).*N_phi-sin(theta_phi).*N_z;        %lipid director projection on phi direction
    
    %normolize
    n_abs=(n_rho.^2+n_phi.^2+n_z.^2).^0.5;
    n_rho=n_rho./n_abs;
    n_phi=n_phi./n_abs;
    n_z=n_z./n_abs;
    % convert to x-y-z system
    n_x=n_rho.*cosphi-n_phi.*sinphi;  %lipid director projection on x direction
    n_y=n_rho.*sinphi+n_phi.*cosphi;  %lipid director projection on y direction
    

    
end

