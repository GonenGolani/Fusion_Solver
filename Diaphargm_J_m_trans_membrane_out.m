function [phi_0_D_sol,phi_p_D_sol,phi_m_D_sol,J_0D,error] = Diaphargm_J_m_trans_membrane_out(kappa_m_kbT_a,j_0,j_p,j_m,j_TM,phi_0_R,phi_p_R,phi_m_R,phi_TM,RES)
% this function gets the bending rigidty times lipid area diveded by kbT (kappa_m./kbT.*a) kappa_m_kbT_a,
% the curvature of three lipid componenets jp jm and
% j_0. and the molar ration of the three componenets at the reserviour: phi_p_R
% , phi_m_R and phi_0_R. And the protien that canot diffuse to the diaphragm
% phi_TM nad induced curvature J_TM

% numerical solve three equations


int=1;
delta=1/RES;
while int<=RES
    phi_p_D_vec(int,1:RES)=linspace(delta,delta*(int-1),RES);
    phi_m_D_vec(int,1:RES)=linspace(2*delta,1-delta*(int),RES);

    int=int+1;
end

if phi_0_R==0
       eq1=log(phi_p_D_vec./phi_p_R)-kappa_m_kbT_a.*(j_p).*(j_p.*(phi_p_D_vec-phi_p_R)+j_m.*(phi_m_D_vec-phi_m_R)-phi_TM.*j_TM);
       eq2=log(phi_m_D_vec./phi_m_R)-kappa_m_kbT_a.*(j_m).*(j_p.*(phi_p_D_vec-phi_p_R)+j_m.*(phi_m_D_vec-phi_m_R)-phi_TM.*j_TM);
else
    eq1=log(phi_p_D_vec./phi_p_R.*phi_0_R./(1-phi_p_D_vec-phi_m_D_vec))-kappa_m_kbT_a.*(j_p-j_0).*(j_0.*((1-phi_p_D_vec-phi_m_D_vec)-phi_0_R)+j_p.*(phi_p_D_vec-phi_p_R)+j_m.*(phi_m_D_vec-phi_m_R)-phi_TM.*j_TM);
    eq2=log(phi_m_D_vec./phi_m_R.*phi_0_R./(1-phi_p_D_vec-phi_m_D_vec))-kappa_m_kbT_a.*(j_m-j_0).*(j_0.*((1-phi_p_D_vec-phi_m_D_vec)-phi_0_R)+j_p.*(phi_p_D_vec-phi_p_R)+j_m.*(phi_m_D_vec-phi_m_R)-phi_TM.*j_TM);
end
%throw out non real elements
[min_y,index_y_temp]=min(abs(eq1)+abs(eq2),[],2);
[error,index_x]=min(min_y);
index_y=index_y_temp(index_x);

phi_p_D_sol=phi_p_D_vec(index_x,index_y);
phi_m_D_sol=phi_m_D_vec(index_x,index_y);
phi_0_D_sol=1-phi_p_D_sol-phi_m_D_sol;

J_0D=j_0*phi_0_D_sol+phi_m_D_sol*j_m+phi_p_D_sol*j_p;


end

