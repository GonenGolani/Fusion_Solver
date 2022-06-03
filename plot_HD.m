function  plot_HD(Diaphragm,Cell,Virus,Shell,res_struc)

    %% plot the diaphragm rim
    

    close all
       
    % rim coordinets are the last row 
    xt_MP=Diaphragm.mid_plane.x(:,res_struc.rho_res);
    yt_MP=Diaphragm.mid_plane.y(:,res_struc.rho_res);
    zt_MP=Diaphragm.mid_plane.z(:,res_struc.rho_res);
    
    xt_DP_up=Diaphragm.dividing_plane.up.x(:,res_struc.rho_res);
    yt_DP_up=Diaphragm.dividing_plane.up.y(:,res_struc.rho_res);
    zt_DP_up=Diaphragm.dividing_plane.up.z(:,res_struc.rho_res);    
    
    xt_DP_down=Diaphragm.dividing_plane.down.x(:,res_struc.rho_res);
    yt_DP_down=Diaphragm.dividing_plane.down.y(:,res_struc.rho_res);
    zt_DP_down=Diaphragm.dividing_plane.down.z(:,res_struc.rho_res);   

    % rim plot
    figure(5)
    hold on;
    %mid plan rim
    plot3(xt_MP,yt_MP,zt_MP,'Color','c','LineWidth',2);
    plot3(-xt_MP,yt_MP,zt_MP,'Color','c','LineWidth',2);
    plot3(xt_MP,-yt_MP,zt_MP,'Color','c','LineWidth',2);
    plot3(-xt_MP,-yt_MP,zt_MP,'Color','c','LineWidth',2);
    
    %up dividing plane rim
    plot3(xt_DP_up,yt_DP_up,zt_DP_up,'Color','c','LineWidth',2);
    plot3(-xt_DP_up,yt_DP_up,zt_DP_up,'Color','c','LineWidth',2);
    plot3(xt_DP_up,-yt_DP_up,zt_DP_up,'Color','c','LineWidth',2);
    plot3(-xt_DP_up,-yt_DP_up,zt_DP_up,'Color','c','LineWidth',2);

    %down dividing plane rim
    plot3(xt_DP_down,yt_DP_down,zt_DP_down,'Color','c','LineWidth',2);
    plot3(-xt_DP_down,yt_DP_down,zt_DP_down,'Color','c','LineWidth',2);
    plot3(xt_DP_down,-yt_DP_down,zt_DP_down,'Color','c','LineWidth',2);
    plot3(-xt_DP_down,-yt_DP_down,zt_DP_down,'Color','c','LineWidth',2);
    axis equal;

    max_rho_dia=max(abs(Diaphragm.rim.x0),abs(Diaphragm.rim.y0));
    [max_rho,index]=max([max_rho_dia ; Cell.R_rim-max_rho_dia ; Virus.Rv-max_rho_dia]);
     factor_diaphragm=2*ceil(max_rho/max_rho_dia); 
     factor_cell=ceil(2*max_rho/(Cell.R_rim-max_rho_dia));
     factor_virus=ceil(2*max_rho/(Virus.Rv-max_rho_dia));
    
    %% plot the membrane
    plot_membrane_shape(Diaphragm.mid_plane,Diaphragm.lipid_director,Diaphragm.lipid_length,Diaphragm.dividing_plane...
        ,Diaphragm.tilt,Diaphragm.splay,Diaphragm.saddle_splay,res_struc,'m','b',factor_diaphragm);
    
    plot_membrane_shape(Cell.mid_plane,Cell.lipid_director,Cell.lipid_length,Cell.dividing_plane,Cell.tilt...
        ,Cell.splay,Cell.saddle_splay,res_struc,'r','b',factor_cell);

    plot_membrane_shape(Virus.mid_plane,Virus.lipid_director,Virus.lipid_length,Virus.dividing_plane,Virus.tilt...
        ,Virus.splay,Virus.saddle_splay,res_struc,'m','r',factor_virus);

    %%plot only diaphragm
    figure(9)
    hold on
    %mid plane surface
    surf( Diaphragm.mid_plane.x, Diaphragm.mid_plane.y, Diaphragm.mid_plane.z, Diaphragm.Energy_density,'EdgeColor','none');
    surf(-Diaphragm.mid_plane.x, Diaphragm.mid_plane.y, Diaphragm.mid_plane.z, Diaphragm.Energy_density,'EdgeColor','none');
    surf( Diaphragm.mid_plane.x,-Diaphragm.mid_plane.y, Diaphragm.mid_plane.z, Diaphragm.Energy_density,'EdgeColor','none');
    surf(-Diaphragm.mid_plane.x,-Diaphragm.mid_plane.y, Diaphragm.mid_plane.z, Diaphragm.Energy_density,'EdgeColor','none');
    colorbar
    title('energy density diaphragm');
    xL=xlabel('nm');
    yL=ylabel('nm');
    set(gca,'FontSize',18); set(gca,'FontWeight','bold'); set(gca,'LineWidth',2);
    xL.FontSize=28; xL.FontWeight='bold';  yL.FontSize=28;  yL.FontWeight='bold';
    axis equal;
    set(gcf,'color','white');

    if Shell.Shell_exist==1
       plot_Shell(Shell.Diaphragm,Shell.Cell,Shell.Junction,Shell.Shell_physical_proprties,res_struc);
    end
    

end

