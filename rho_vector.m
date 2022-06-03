function [out_vector] = rho_vector(rho,length)
% input rho value and length
% out: [rho^0, rho^1, ... , rho^legnth]
    int=0;
    out_vector=ones(1,length);
    while int<length
        out_vector(int+1)=rho^int;
        int=int+1;
    end
   
   
end

