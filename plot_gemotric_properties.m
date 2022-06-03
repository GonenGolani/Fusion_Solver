function  plot_gemotric_properties(plane,tilt,splay,saddle_splay,splay_grad,kappa_tilt,chi,J_sm,Minimazation,res_struc,legend_name,color_plot)


    
if Minimazation.axial_symmetry==1
    rho_min=(plane.x(1,1)^2+plane.y(1,1)^2)^0.5;
    rho_max=(plane.x(1,res_struc.rho_res)^2+plane.y(1,res_struc.rho_res)^2)^0.5;
    rho_vec=linspace(rho_min,rho_max,res_struc.rho_res);
    
    
    tilt_square=(tilt.x(1,:).^2+tilt.y(1,:).^2+tilt.z(1,:).^2);
    
    energy_density=0.5*(splay(1,:).^2-2*splay(1,:)*J_sm)+saddle_splay(1,:)*chi+0.5*tilt_square*kappa_tilt;

    
    figure(10)
    hold on
    if strcmp(legend_name,'none')==0
        plot(rho_vec,tilt_square.^0.5,'LineWidth',2,'DisplayName',legend_name,'Color',color_plot);
    else
        p1=plot(rho_vec,tilt_square.^0.5,'LineWidth',2,'Color',color_plot);
        set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    xlabel('{\rho} [nm]')
    ylabel('Tilt')
    set(gca,'FontSize',18); set(gca,'FontWeight','bold'); set(gca,'LineWidth',2);
    xL.FontSize=28; xL.FontWeight='bold';  yL.FontSize=28;  yL.FontWeight='bold';
    set(gcf,'color','white')
    legend;

    
    figure(11)
    hold on
    if strcmp(legend_name,'none')==0
        plot(rho_vec,saddle_splay(1,:),'LineWidth',2,'DisplayName',legend_name,'Color',color_plot);
    else
        p1=plot(rho_vec,saddle_splay(1,:),'LineWidth',2,'Color',color_plot);
        set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    xlabel('{\rho} [nm]')
    ylabel('Saddle Splay [nm^-^2]')
    set(gca,'FontSize',18); set(gca,'FontWeight','bold'); set(gca,'LineWidth',2);
    xL.FontSize=28; xL.FontWeight='bold';  yL.FontSize=28;  yL.FontWeight='bold';
    set(gcf,'color','white')
    legend;

    
    
    figure(12)
    hold on
    if strcmp(legend_name,'none')==0
        plot(rho_vec,splay(1,:),'LineWidth',2,'DisplayName',legend_name,'Color',color_plot);
    else
        p1=plot(rho_vec,splay(1,:),'LineWidth',2,'Color',color_plot);
        set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    xlabel('{\rho} [nm]')
    ylabel('Splay [nm^-^1]')
    set(gca,'FontSize',18); set(gca,'FontWeight','bold'); set(gca,'LineWidth',2);
    xL.FontSize=28; xL.FontWeight='bold';  yL.FontSize=28;  yL.FontWeight='bold';
    set(gcf,'color','white')
    legend;
    
    figure(13)
    hold on
    if strcmp(legend_name,'none')==0
        plot(rho_vec,splay_grad(1,:),'LineWidth',2,'DisplayName',legend_name,'Color',color_plot);
    else
        p1=plot(rho_vec,splay_grad(1,:),'LineWidth',2,'Color',color_plot);
        set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    xlabel('{\rho} [nm]')
    ylabel('Grad Splay [nm^-^2]')
    set(gca,'FontSize',18); set(gca,'FontWeight','bold'); set(gca,'LineWidth',2);
    xL.FontSize=28; xL.FontWeight='bold';  yL.FontSize=28;  yL.FontWeight='bold';
    set(gcf,'color','white')
    legend;
    
    figure(14)
    hold on
    if strcmp(legend_name,'none')==0
        plot(rho_vec,energy_density(1,:),'LineWidth',2,'DisplayName',legend_name,'Color',color_plot);
    else
        p1=plot(rho_vec,energy_density(1,:),'LineWidth',2,'Color',color_plot);
        set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    xlabel('{\rho} [nm]')
    ylabel('Energy density [{\kappa}_m*nm^-^2]')
    set(gca,'FontSize',18); set(gca,'FontWeight','bold'); set(gca,'LineWidth',2);
    xL.FontSize=28; xL.FontWeight='bold';  yL.FontSize=28;  yL.FontWeight='bold';
    set(gcf,'color','white')
    legend;
    
end
    
if Minimazation.axial_symmetry==0
    
    int=1;
    while int<=res_struc.phi_res
        rho_min(int)=(plane.x(int,1)^2+plane.y(int,1)^2)^0.5;
        rho_max(int)=(plane.x(int,res_struc.rho_res)^2+plane.y(int,res_struc.rho_res)^2)^0.5;
        rho_vec(int,:)=linspace(rho_min(int),rho_max(int),res_struc.rho_res);
        phi_vec=linspace(0,pi/2,res_struc.phi_res);
        tilt_square(int,:)=(tilt.x(int,:).^2+tilt.y(int,:).^2+tilt.z(int,:).^2);

        figure(10)
        hold on
        plot3(rho_vec(int,:),phi_vec(int)*180/pi*ones(1,res_struc.rho_res),tilt_square(int,:),'LineWidth',2);
        xlabel('{\rho} [nm]')
        ylabel('Azimuth angle [deg]')
        zlabel('Tilt square')
        set(gca,'FontSize',18); set(gca,'FontWeight','bold'); set(gca,'LineWidth',2);
        xL.FontSize=28; xL.FontWeight='bold';  yL.FontSize=28;  yL.FontWeight='bold';
        set(gcf,'color','white');

        figure(11)
        hold on
        plot3(rho_vec(int,:),phi_vec(int)*180/pi*ones(1,res_struc.rho_res),saddle_splay(int,:),'LineWidth',2);

        xlabel('{\rho} [nm]')
        ylabel('Azimuth angle [deg]')
        zlabel('Saddle Splay [nm^-^2]')
        set(gca,'FontSize',18); set(gca,'FontWeight','bold'); set(gca,'LineWidth',2);
        xL.FontSize=28; xL.FontWeight='bold';  yL.FontSize=28;  yL.FontWeight='bold';
        set(gcf,'color','white');

        figure(12)
         hold on
        plot3(rho_vec(int,:),phi_vec(int)*180/pi*ones(1,res_struc.rho_res),splay(int,:),'LineWidth',2);
        xlabel('{\rho} [nm]')
        ylabel('Azimuth angle [deg]')
        zlabel('Splay [nm^-^1]')
        set(gca,'FontSize',18); set(gca,'FontWeight','bold'); set(gca,'LineWidth',2);
        xL.FontSize=28; xL.FontWeight='bold';  yL.FontSize=28;  yL.FontWeight='bold';
        set(gcf,'color','white');
        
        
        figure(13)
         hold on
        plot3(rho_vec(int,:),phi_vec(int)*180/pi*ones(1,res_struc.rho_res),splay_grad(int,:),'LineWidth',2);
        xlabel('{\rho} [nm]')
        ylabel('Azimuth angle [deg]')
        zlabel('Grad Splay [nm^-^2]')
        set(gca,'FontSize',18); set(gca,'FontWeight','bold'); set(gca,'LineWidth',2);
        xL.FontSize=28; xL.FontWeight='bold';  yL.FontSize=28;  yL.FontWeight='bold';
        set(gcf,'color','white');
        
        clear rho_vec
        int=int+1;
    end
end


end







