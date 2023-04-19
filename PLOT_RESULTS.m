
clc
clear all
close all
addpath('C:\Users\gonen\Google Drive\Fusion viral project\Fusion_Solver\source files\')
cd('C:\Fusion intermidate Solver\Results\Stalk_only\Change_l_J0_-0.22_STALK\1')
[Diaphragm,Cell,Virus,Shell,Minimazation,res_struc,General_physical_properties,Energy,Total_energy,DOF_vector,System_dimensions]...
    = reload_HD_state('step_16.mat');

[Diaphragm,Cell,Virus,Shell,Energy,Total_energy] = find_HD_energy(Diaphragm,Cell,Virus,Shell,Minimazation,res_struc);



Plot_fusing_membrane(Diaphragm,Virus,Cell,100);

plot_profile(Diaphragm,Cell,Virus,Shell,1,1,0,Minimazation.up_down_symmetry);



plot_HD(Diaphragm,Cell,Virus,Shell,res_struc);
    Kappa_tilt=Diaphragm.physical_properties.kappa_tilt;
    chi=Diaphragm.physical_properties.kappa_bar_up;
    
    plot_gemotric_properties(Diaphragm.dividing_plane.up,Diaphragm.tilt.up,Diaphragm.splay.up,Diaphragm.saddle_splay.up,Diaphragm.splay.splay_grad_up,...
        Kappa_tilt,chi,Diaphragm.physical_properties.J0_up,Minimazation,res_struc,'Distal','b');
        
    plot_gemotric_properties(Virus.dividing_plane.up,Virus.tilt.up,Virus.splay.up,Virus.saddle_splay.up,Virus.splay.splay_grad_up,...
        Kappa_tilt,chi,Virus.physical_properties.J0_up,Minimazation,res_struc,'none','b');
    
    plot_gemotric_properties(Diaphragm.dividing_plane.down,Diaphragm.tilt.down,Diaphragm.splay.down,Diaphragm.saddle_splay.down,Diaphragm.splay.splay_grad_down,...
        Kappa_tilt,chi,Diaphragm.physical_properties.J0_down,Minimazation,res_struc,'Spherical Inner monolayer','m');
    
    plot_gemotric_properties(Cell.dividing_plane.down,Cell.tilt.down,Cell.splay.down,Cell.saddle_splay.down,Cell.splay.splay_grad_down,...
        Kappa_tilt,chi,Cell.physical_properties.J0_down,Minimazation,res_struc,'none','m');


    plot_gemotric_properties(Virus.dividing_plane.down,Virus.tilt.down,Virus.splay.down,Virus.saddle_splay.down,Virus.splay.splay_grad_down,...
        Kappa_tilt,chi,Virus.physical_properties.J0_down,Minimazation,res_struc,'Proximal','r');
    plot_gemotric_properties(Cell.dividing_plane.up,Cell.tilt.up,Cell.splay.up,Cell.saddle_splay.up,Cell.splay.splay_grad_up,...
        Kappa_tilt,chi,Cell.physical_properties.J0_up,Minimazation,res_struc,'Spherical outer monolayer','b');


%membraen stress
figure(15)
hold on
plot(Diaphragm.mid_plane.x(1,:),Diaphragm.Energy_density(1,:),'k','LineWidth',2);
plot(Cell.mid_plane.x(1,:),Cell.Energy_density(1,:),'k','LineWidth',2);
xL=xlabel('{\rho} [nm]');
yL=ylabel('Membrane stress [\kappa_m{\cdot}nm^-^2]');
set(gca,'FontSize',18);
set(gca,'FontWeight','bold');
xL.FontSize=18;
xL.FontWeight='bold';
yL.FontSize=18;
yL.FontWeight='bold';
set(gcf,'color','white')
legend;


