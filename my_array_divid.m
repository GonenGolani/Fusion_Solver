function [output] = my_array_divid(arr1,arr2,res_struct)
    output=arr1./arr2;
    %find if nan exist  - if so change to 0
    output(isnan(output))=0;
    
    %correct last line of x
    output(res_struct.phi_res,:)=add_diff_last_line(output(res_struct.phi_res-1,:),output(res_struct.phi_res-2,:)...
        ,res_struct.phi_res-1,res_struct.phi_res-2,res_struct.phi_res); 
end

