function  Plot_fusing_membrane(Diaphragm,Virus,Cell,res)
% outline: Plot the virus fusion site and cell membranes
% Virus structure
%input: Virus stact:
% torus small (inner) radius r_sv 
% torus larger radius r_bv
% torus center in z direction is z_v_center
% radius of fusion site: Rv



%% torus in self coordinets
u=linspace(0,2*pi,2*res); 
v=linspace(0,2*pi,2*res);
[U,V]=meshgrid(u,v);
X_Virues_full=(Virus.r_bv+Virus.r_sv*cos(V)).*cos(U);
Y_Virues_full=Virus.r_sv*sin(V);
Z_Virues_full=Virus.z_v_center+(Virus.r_bv+Virus.r_sv*cos(V)).*sin(U);
%Normal vector to torus


%% Fusion site rim in the rho phi system
phi=linspace(0,2*pi,res);
rho_vec=linspace(Virus.Rv,Virus.r_sv,res);
z_v_center_vec=Virus.z_v_center*ones(1,res);
r_bv_vec=Virus.r_bv*ones(1,res);
r_sv_vec=Virus.r_sv*ones(1,res);


%plot torus using phi a and rho corrdinates
[PHI,RHO]=meshgrid(phi,rho_vec);

sinV=RHO.*sin(PHI)./Virus.r_sv;
cosV=(1-sinV.^2).^0.5; % can be +- for inner or outer intersection;
y_star=(r_bv_vec+r_sv_vec.*cosV); %just to short things
cosU=RHO.*cos(PHI)./y_star;
sinU=-(1-cosU.^2).^0.5;


X_Virues=(Virus.r_bv+Virus.r_sv*cosV).*cosU;
Y_Virues=Virus.r_sv*sinV;
Z_Virues=Virus.z_v_center+(Virus.r_bv+Virus.r_sv*cosV).*sinU;

N_Torus_x=cosU.*cosV; 
N_Torus_y=sinV; 
N_Torus_z=sinU.*cosV; 

%plot rim
sinV_rim=Virus.Rv*sin(phi)/Virus.r_sv;
cosV_rim=(1-sinV_rim.^2).^0.5; % can be +- for inner or outer intersection;
y_star_rim=(r_bv_vec+r_sv_vec.*cosV_rim); %just to short things
cosU_rim=Virus.Rv*cos(phi)./y_star_rim;
sinU_rim=-(1-cosU_rim.^2).^0.5;




dsinVdrho_rim=sin(phi)/Virus.r_sv; 
dcosVdrho_rim=-dsinVdrho_rim.*sinV_rim./cosV_rim; %can be +-
dcosUdrho_rim=(y_star.^2.*cos(phi)+Virus.Rv^2.*cos(phi).^2)./y_star_rim.^3;
dsinUdrho_rim=-dcosUdrho_rim.*sinU_rim./cosU_rim; %can be +-

% the contact rim between Fusion site and virus
xt_virues=(Virus.r_bv*cosU_rim+Virus.r_sv*cosU_rim.*cosV_rim);
yt_virues=Virus.r_sv*sinV_rim;
zt_virues=z_v_center_vec+(Virus.r_bv*sinU_rim+Virus.r_sv*sinU_rim.*cosV_rim);  % can be +- for lower or upper;

% tangent angle
dzdrho_rim=y_star_rim.*dsinUdrho_rim+Virus.r_sv.*sinU_rim.*dcosVdrho_rim;



% normal to torus at the rim
Nx_rim=cosV_rim.*cosU_rim;
Ny_rim=sinV_rim;
Nz_rim=cosV_rim.*sinU_rim;

%% cell membrane plot

cell_rho_vector=linspace(Cell.R_rim,Cell.R_curv,res);
cell_phi_vector=linspace(0,2*pi,res);
X_cell=cell_rho_vector'*cos(cell_phi_vector);
Y_cell=cell_rho_vector'*sin(cell_phi_vector);
if (Cell.R_curv>0)
    Z_cell=(-(Cell.R_curv^2-Cell.R_rim^2)^0.5+(Cell.R_curv^2-cell_rho_vector.^2).^0.5)'*ones(1,res);
end

if (Cell.R_curv<0)
    Z_cell=((Cell.R_curv^2-Cell.R_rim^2)^0.5-(Cell.R_curv^2-cell_rho_vector.^2).^0.5)'*ones(1,res);

end


if(strcmp(Cell.R_curv,'flat'))
        Z_cell=0*ones(res,res);
end

% cell rim
xt_cell=Cell.R_rim*cos(cell_phi_vector);
yt_cell=Cell.R_rim*sin(cell_phi_vector);
zt_cell=zeros(1,res);
%% plot

lim_plot=1.5*max(Cell.R_rim,Virus.Rv);

figure(4)
hold on
% viures
%s1=surf(X_Virues,Y_Virues,Z_cell,'FaceColor','none'); %plot using rho phi coordinates
s2=surf(X_Virues_full,Y_Virues_full,Z_Virues_full,'FaceColor','none');  %plot using U V coordinates
s2.EdgeColor='k';

P3=plot3(xt_virues,yt_virues,zt_virues,'r');
P3.LineWidth=2;

%cell
s2=surf(X_cell,Y_cell,Z_cell); 
s2.FaceColor='none';
P4=plot3(xt_cell,yt_cell,zt_cell,'r');
P4.LineWidth=2;
axis equal
xlim([-lim_plot lim_plot])
ylim([-lim_plot lim_plot])
zlim([-lim_plot lim_plot])
view(90,0)
xlabel('x')
ylabel('y')
zlabel('z')



end

