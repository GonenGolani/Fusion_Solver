function [phi_p_D,J_0D,Error] = Diaphargm_J_m_trans_membrane_out_2componenets(kappa_m_kbT_a,j_p,j_m,j_TM,phi_p_R,phi_m_R,phi_TM,RES)
% this function gets the bending rigidty times lipid area diveded by kbT (kappa_m./kbT.*a) kappa_m_kbT_a,
% the curvature of three lipid componenets jp jm and
%  and the molar ration of the three componenets at the reserviour: phi_p_R
% , phi_m_R and  And the protien that canot diffuse to the diaphragm
% phi_TM nad induced curvature J_TM


    phi_D_vec=linspace(1/RES,1-1/RES,RES);
    J0R=phi_p_R*j_p+phi_m_R*j_m+j_TM*phi_TM;
    J0D_vec=j_p*phi_D_vec+j_m*(1-phi_D_vec);

    %eq=phi_D_vec./(1-phi_D_vec)-phi_p_R./phi_m_R.*exp(-kappa_m_kbT_a*(j_p-j_m)*(phi_D_vec*j_p+(1-phi_D_vec)*j_m-J0R));
    exp_term=exp(kappa_m_kbT_a.*(j_p-j_m).*(J0D_vec-J0R));
    
    eq=phi_D_vec-phi_p_R./(exp_term.*phi_m_R+phi_p_R);

    [Error,index]=min(abs(eq));
    phi_p_D=phi_D_vec(index);

    J_0D=j_p*phi_p_D+j_m*(1-phi_p_D);

end

