function [out_vec] = my_harmonic_expanstion(an,x)
% f(x)= a1+a2*cos(2*x)+a3**cos(4*x)+...+aN**cos(2*N*x)
    out_vec=zeros(1,length(x));
    max_order=length(an);
    int=1;
    while int<=max_order
        out_vec=out_vec+an(int)*cos(2*x*(int-1));
        int=int+1;
    end

end