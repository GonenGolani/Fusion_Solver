function [phi_s_free] = Special_lipid_mole_fraction_at_IFITM3_free_part(Affinity,phi_IFITM3,bound_lipids,phi_reservior)

E=(1-exp(-Affinity));
A=(1-phi_IFITM3*bound_lipids)*E;
B=-(1-phi_IFITM3*bound_lipids*E+phi_reservior*E);
C=phi_reservior;

phi_s_free=(-B-(B^2-4*A*C)^0.5)/(2*A);

end