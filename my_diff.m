function [diff_array] = my_diff(array,dim,phi_res,rho_res)
    
    
    %finds diff in dim directions and fills row or column with linear extrpolation 
    %  

    
    if(dim=='y')
        [~,diff_array]=gradient(array);
        %diff_array=diff(array,1,1);
        diff_array(phi_res,:)=add_diff_last_line(diff_array(phi_res-1,:),diff_array(phi_res-2,:),phi_res-1,phi_res-2,phi_res); 
    end
    
    if(dim=='x')
        [diff_array,~]=gradient(array);
        %diff_array=diff(array,1,2);
        diff_array(:,rho_res)=add_diff_last_line(diff_array(:,rho_res-1),diff_array(:,rho_res-2),rho_res-1,rho_res-2,rho_res); 
    end


end
