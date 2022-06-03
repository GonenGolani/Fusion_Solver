function [cosine_matrix] = cosine_matrix_creatore(phi_res,rim_coeff_number)

cosine_matrix=zeros(rim_coeff_number,phi_res);
int1=1;
phi_vector=linspace(0,pi/2,phi_res);
while int1<=rim_coeff_number
    int2=1;
    while int2<=phi_res
        cosine_matrix(int1,int2)=cos(2*(int1-1)*phi_vector(int2));
        int2=int2+1;
    end
    int1=int1+1;
end
   
end

