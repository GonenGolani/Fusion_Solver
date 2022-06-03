function [z_vector] = vector_my_poly_val(coeff_vector,rho_vector)
    % this function recives the rho vector and polynomial expanstion
    % coeffciantns h_l and returns the z=h0+h1*rho+...hl*rho^l value at each point
    poly_order=length(coeff_vector);
    int=1;
    z_vector=0;
    while int<=poly_order
        rho_vector_poly=rho_vector.^(int-1);
        z_vector=z_vector+rho_vector_poly*coeff_vector(int);
        
        int=int+1;
    end
    
end

