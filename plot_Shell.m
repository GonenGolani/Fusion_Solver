function  plot_Shell(Shell_Diaphragm,Shell_Cell,Shell_Junction,Shell_physical_proprties,res_struc)

% plot the diaphragm rim

Shell_color=[0 0.6 0.3];




%upper bound of line
Shell_Diaphragm.upper_plane.x=Shell_Diaphragm.mid_plane.x+Shell_physical_proprties.width*Shell_Diaphragm.normal.x/2;
Shell_Diaphragm.upper_plane.y=Shell_Diaphragm.mid_plane.y+Shell_physical_proprties.width*Shell_Diaphragm.normal.y/2;
Shell_Diaphragm.upper_plane.z=Shell_Diaphragm.mid_plane.z+Shell_physical_proprties.width*Shell_Diaphragm.normal.z/2;

Shell_Diaphragm.lower_plane.x=Shell_Diaphragm.mid_plane.x-Shell_physical_proprties.width*Shell_Diaphragm.normal.x/2;
Shell_Diaphragm.lower_plane.y=Shell_Diaphragm.mid_plane.y-Shell_physical_proprties.width*Shell_Diaphragm.normal.y/2;
Shell_Diaphragm.lower_plane.z=Shell_Diaphragm.mid_plane.z-Shell_physical_proprties.width*Shell_Diaphragm.normal.z/2;

Shell_Cell.upper_plane.x=Shell_Cell.mid_plane.x+Shell_physical_proprties.width*Shell_Cell.normal.x/2;
Shell_Cell.upper_plane.y=Shell_Cell.mid_plane.y+Shell_physical_proprties.width*Shell_Cell.normal.y/2;
Shell_Cell.upper_plane.z=Shell_Cell.mid_plane.z+Shell_physical_proprties.width*Shell_Cell.normal.z/2;

Shell_Cell.lower_plane.x=Shell_Cell.mid_plane.x-Shell_physical_proprties.width*Shell_Cell.normal.x/2;
Shell_Cell.lower_plane.y=Shell_Cell.mid_plane.y-Shell_physical_proprties.width*Shell_Cell.normal.y/2;
Shell_Cell.lower_plane.z=Shell_Cell.mid_plane.z-Shell_physical_proprties.width*Shell_Cell.normal.z/2;



%normal of junction only needed in plot
[Shell_Junction.normal.x,Shell_Junction.normal.y,Shell_Junction.normal.z] =...
    surfnorm(Shell_Junction.mid_plane.x,Shell_Junction.mid_plane.y,Shell_Junction.mid_plane.z);

%flip normal if needed needs to point down!
if Shell_Junction.normal.z(1,1)>0
    Shell_Junction.normal.x=Shell_Junction.normal.x*-1;
    Shell_Junction.normal.y=Shell_Junction.normal.y*-1;
    Shell_Junction.normal.z=Shell_Junction.normal.z*-1;   
end

Shell_Junction.upper_plane.x=Shell_Junction.mid_plane.x+Shell_physical_proprties.width*Shell_Junction.normal.x/2;
Shell_Junction.upper_plane.y=Shell_Junction.mid_plane.y+Shell_physical_proprties.width*Shell_Junction.normal.y/2;
Shell_Junction.upper_plane.z=Shell_Junction.mid_plane.z+Shell_physical_proprties.width*Shell_Junction.normal.z/2;


Shell_Junction.lower_plane.x=Shell_Junction.mid_plane.x-Shell_physical_proprties.width*Shell_Junction.normal.x/2;
Shell_Junction.lower_plane.y=Shell_Junction.mid_plane.y-Shell_physical_proprties.width*Shell_Junction.normal.y/2;
Shell_Junction.lower_plane.z=Shell_Junction.mid_plane.z-Shell_physical_proprties.width*Shell_Junction.normal.z/2;


figure(5);
hold on
surf(Shell_Diaphragm.mid_plane.x,Shell_Diaphragm.mid_plane.y,Shell_Diaphragm.mid_plane.z,'FaceColor','none','EdgeColor',Shell_color);
surf(-Shell_Diaphragm.mid_plane.x,Shell_Diaphragm.mid_plane.y,Shell_Diaphragm.mid_plane.z,'FaceColor','none','EdgeColor',Shell_color);
surf(Shell_Diaphragm.mid_plane.x,-Shell_Diaphragm.mid_plane.y,Shell_Diaphragm.mid_plane.z,'FaceColor','none','EdgeColor',Shell_color);
surf(-Shell_Diaphragm.mid_plane.x,-Shell_Diaphragm.mid_plane.y,Shell_Diaphragm.mid_plane.z,'FaceColor','none','EdgeColor',Shell_color);

surf(Shell_Cell.mid_plane.x,Shell_Cell.mid_plane.y,Shell_Cell.mid_plane.z,'FaceColor','none','EdgeColor',Shell_color);
surf(-Shell_Cell.mid_plane.x,Shell_Cell.mid_plane.y,Shell_Cell.mid_plane.z,'FaceColor','none','EdgeColor',Shell_color);
surf(Shell_Cell.mid_plane.x,-Shell_Cell.mid_plane.y,Shell_Cell.mid_plane.z,'FaceColor','none','EdgeColor',Shell_color);
surf(-Shell_Cell.mid_plane.x,-Shell_Cell.mid_plane.y,Shell_Cell.mid_plane.z,'FaceColor','none','EdgeColor',Shell_color);

surf(Shell_Junction.mid_plane.x,Shell_Junction.mid_plane.y,Shell_Junction.mid_plane.z,'FaceColor','none','EdgeColor','b');
surf(-Shell_Junction.mid_plane.x,Shell_Junction.mid_plane.y,Shell_Junction.mid_plane.z,'FaceColor','none','EdgeColor','b');
surf(Shell_Junction.mid_plane.x,-Shell_Junction.mid_plane.y,Shell_Junction.mid_plane.z,'FaceColor','none','EdgeColor','b');
surf(-Shell_Junction.mid_plane.x,-Shell_Junction.mid_plane.y,Shell_Junction.mid_plane.z,'FaceColor','none','EdgeColor','b');


figure(7) %plot at phi=0
hold on
plot(Shell_Diaphragm.mid_plane.x(1,:),Shell_Diaphragm.mid_plane.z(1,:),'--','Color',Shell_color)
plot(Shell_Cell.mid_plane.x(1,:),Shell_Cell.mid_plane.z(1,:),'--','Color',Shell_color)
plot(Shell_Junction.mid_plane.x(1,:),Shell_Junction.mid_plane.z(1,:),'--','Color','b')


fill([Shell_Diaphragm.upper_plane.x(1,:), flip(Shell_Diaphragm.lower_plane.x(1,:))],...
 [Shell_Diaphragm.upper_plane.z(1,:) ,flip(Shell_Diaphragm.lower_plane.z(1,:))],Shell_color,'LineStyle','none','FaceAlpha',0.1);
fill([Shell_Cell.upper_plane.x(1,:), flip(Shell_Cell.lower_plane.x(1,:))],...
 [Shell_Cell.upper_plane.z(1,:) ,flip(Shell_Cell.lower_plane.z(1,:))],Shell_color,'LineStyle','none','FaceAlpha',0.1);
fill([Shell_Junction.upper_plane.x(1,:), flip(Shell_Junction.lower_plane.x(1,:))],...
 [Shell_Junction.upper_plane.z(1,:) ,flip(Shell_Junction.lower_plane.z(1,:))],'b','LineStyle','none','FaceAlpha',0.1);


figure(8) %plot at phi=pi/2
hold on
plot(Shell_Diaphragm.mid_plane.y(res_struc.phi_res,:),Shell_Diaphragm.mid_plane.z(res_struc.phi_res,:),'--','Color',Shell_color)
plot(Shell_Cell.mid_plane.y(res_struc.phi_res,:),Shell_Cell.mid_plane.z(res_struc.phi_res,:),'--','Color',Shell_color)
plot(Shell_Junction.mid_plane.y(res_struc.phi_res,:),Shell_Junction.mid_plane.z(res_struc.phi_res,:),'--','Color','b')

fill([Shell_Diaphragm.upper_plane.y(res_struc.phi_res,:), flip(Shell_Diaphragm.lower_plane.y(res_struc.phi_res,:))],...
 [Shell_Diaphragm.upper_plane.z(res_struc.phi_res,:) ,flip(Shell_Diaphragm.lower_plane.z(res_struc.phi_res,:))],Shell_color,'LineStyle','none','FaceAlpha',0.1);
fill([Shell_Cell.upper_plane.y(res_struc.phi_res,:), flip(Shell_Cell.lower_plane.y(res_struc.phi_res,:))],...
 [Shell_Cell.upper_plane.z(res_struc.phi_res,:) ,flip(Shell_Cell.lower_plane.z(res_struc.phi_res,:))],Shell_color,'LineStyle','none','FaceAlpha',0.1);
fill([Shell_Junction.upper_plane.y(res_struc.phi_res,:), flip(Shell_Junction.lower_plane.y(res_struc.phi_res,:))],...
 [Shell_Junction.upper_plane.z(res_struc.phi_res,:) ,flip(Shell_Junction.lower_plane.z(res_struc.phi_res,:))],'r','LineStyle','none','FaceAlpha',0.1);


title('phi=0');
axis equal;
set(gca,'FontSize',18); 
set(gca,'FontWeight','bold'); 
set(gca,'LineWidth',2); 
set(gcf,'color','white')
set(gcf,'color','white')


    

end

