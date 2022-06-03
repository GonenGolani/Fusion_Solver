function  plot_membrane_shape(mid_plane,lipid_director,lipid_length,dividing_plane,tilt,splay,saddle_splay,res_struct,up_color,down_color,decimation_factor)



tilt_sqaure_up=tilt.up.x.^2+tilt.up.y.^2+tilt.up.z.^2;
tilt_sqaure_down=tilt.down.x.^2+tilt.down.y.^2+tilt.down.z.^2;

figure(1);
hold on
surf(dividing_plane.up.x,dividing_plane.up.y,dividing_plane.up.z,tilt_sqaure_up,'EdgeColor','none')
surf(-dividing_plane.up.x,dividing_plane.up.y,dividing_plane.up.z,tilt_sqaure_up,'EdgeColor','none')
surf(dividing_plane.up.x,-dividing_plane.up.y,dividing_plane.up.z,tilt_sqaure_up,'EdgeColor','none')
surf(-dividing_plane.up.x,-dividing_plane.up.y,dividing_plane.up.z,tilt_sqaure_up,'EdgeColor','none')

surf(dividing_plane.down.x,dividing_plane.down.y,dividing_plane.down.z,tilt_sqaure_down,'EdgeColor','none')
surf(-dividing_plane.down.x,dividing_plane.down.y,dividing_plane.down.z,tilt_sqaure_down,'EdgeColor','none')
surf(dividing_plane.down.x,-dividing_plane.down.y,dividing_plane.down.z,tilt_sqaure_down,'EdgeColor','none')
surf(-dividing_plane.down.x,-dividing_plane.down.y,dividing_plane.down.z,tilt_sqaure_down,'EdgeColor','none')
title('tilt^2');
axis equal
colorbar

figure(2);
hold on
surf(dividing_plane.up.x,dividing_plane.up.y,dividing_plane.up.z,splay.up,'EdgeColor','none')
surf(-dividing_plane.up.x,dividing_plane.up.y,dividing_plane.up.z,splay.up,'EdgeColor','none')
surf(dividing_plane.up.x,-dividing_plane.up.y,dividing_plane.up.z,splay.up,'EdgeColor','none')
surf(-dividing_plane.up.x,-dividing_plane.up.y,dividing_plane.up.z,splay.up,'EdgeColor','none')

surf(dividing_plane.down.x,dividing_plane.down.y,dividing_plane.down.z,splay.down,'EdgeColor','none')
surf(-dividing_plane.down.x,dividing_plane.down.y,dividing_plane.down.z,splay.down,'EdgeColor','none')
surf(dividing_plane.down.x,-dividing_plane.down.y,dividing_plane.down.z,splay.down,'EdgeColor','none')
surf(-dividing_plane.down.x,-dividing_plane.down.y,dividing_plane.down.z,splay.down,'EdgeColor','none')

title('splay');
axis equal
colorbar

figure(3);
hold on
surf(dividing_plane.up.x,dividing_plane.up.y,dividing_plane.up.z,saddle_splay.up,'EdgeColor','none')
surf(-dividing_plane.up.x,dividing_plane.up.y,dividing_plane.up.z,saddle_splay.up,'EdgeColor','none')
surf(dividing_plane.up.x,-dividing_plane.up.y,dividing_plane.up.z,saddle_splay.up,'EdgeColor','none')
surf(-dividing_plane.up.x,-dividing_plane.up.y,dividing_plane.up.z,saddle_splay.up,'EdgeColor','none')

surf(dividing_plane.down.x,dividing_plane.down.y,dividing_plane.down.z,saddle_splay.down,'EdgeColor','none')
surf(-dividing_plane.down.x,dividing_plane.down.y,dividing_plane.down.z,saddle_splay.down,'EdgeColor','none')
surf(dividing_plane.down.x,-dividing_plane.down.y,dividing_plane.down.z,saddle_splay.down,'EdgeColor','none')
surf(-dividing_plane.down.x,-dividing_plane.down.y,dividing_plane.down.z,saddle_splay.down,'EdgeColor','none')

title('saddle splay');
axis equal
colorbar


mid_plane.x=my_downsample(mid_plane.x',decimation_factor)';
mid_plane.y=my_downsample(mid_plane.y',decimation_factor)';
mid_plane.z=my_downsample(mid_plane.z',decimation_factor)';
mid_plane.normal.x=my_downsample(mid_plane.normal.x',decimation_factor)';
mid_plane.normal.y=my_downsample(mid_plane.normal.y',decimation_factor)';
mid_plane.normal.z=my_downsample(mid_plane.normal.z',decimation_factor)';
lipid_director.up.x=my_downsample(lipid_director.up.x',decimation_factor)';
lipid_director.up.y=my_downsample(lipid_director.up.y',decimation_factor)';
lipid_director.up.z=my_downsample(lipid_director.up.z',decimation_factor)';
lipid_director.down.x=my_downsample(lipid_director.down.x',decimation_factor)';
lipid_director.down.y=my_downsample(lipid_director.down.y',decimation_factor)';
lipid_director.down.z=my_downsample(lipid_director.down.z',decimation_factor)';
lipid_length.up=my_downsample(lipid_length.up',decimation_factor)';
lipid_length.down=my_downsample(lipid_length.down',decimation_factor)';


figure(5)
hold on
%mid plane surface
surf(mid_plane.x,mid_plane.y,mid_plane.z,'FaceColor','none');
surf(-mid_plane.x,mid_plane.y,mid_plane.z,'FaceColor','none');
surf(mid_plane.x,-mid_plane.y,mid_plane.z,'FaceColor','none');
surf(-mid_plane.x,-mid_plane.y,mid_plane.z,'FaceColor','none');
axis equal;


%plot lipid director up
quiver3(mid_plane.x,mid_plane.y,mid_plane.z,lipid_director.up.x.*lipid_length.up,lipid_director.up.y.*lipid_length.up,lipid_director.up.z.*lipid_length.up,0,'color',up_color,'ShowArrowHead','off');
quiver3(-mid_plane.x,mid_plane.y,mid_plane.z,-lipid_director.up.x.*lipid_length.up,lipid_director.up.y.*lipid_length.up,lipid_director.up.z.*lipid_length.up,0,'color',up_color,'ShowArrowHead','off');
quiver3(mid_plane.x,-mid_plane.y,mid_plane.z,lipid_director.up.x.*lipid_length.up,-lipid_director.up.y.*lipid_length.up,lipid_director.up.z.*lipid_length.up,0,'color',up_color,'ShowArrowHead','off');
quiver3(-mid_plane.x,-mid_plane.y,mid_plane.z,-lipid_director.up.x.*lipid_length.up,-lipid_director.up.y.*lipid_length.up,lipid_director.up.z.*lipid_length.up,0,'color',up_color,'ShowArrowHead','off');

%plot lipid director down
quiver3(mid_plane.x,mid_plane.y,mid_plane.z,lipid_director.down.x.*lipid_length.down,lipid_director.down.y.*lipid_length.down,lipid_director.down.z.*lipid_length.down,0,'color',down_color,'ShowArrowHead','off');
quiver3(-mid_plane.x,mid_plane.y,mid_plane.z,-lipid_director.down.x.*lipid_length.down,lipid_director.down.y.*lipid_length.down,lipid_director.down.z.*lipid_length.down,0,'color',down_color,'ShowArrowHead','off');
quiver3(mid_plane.x,-mid_plane.y,mid_plane.z,lipid_director.down.x.*lipid_length.down,-lipid_director.down.y.*lipid_length.down,lipid_director.down.z.*lipid_length.down,0,'color',down_color,'ShowArrowHead','off');
quiver3(-mid_plane.x,-mid_plane.y,mid_plane.z,-lipid_director.down.x.*lipid_length.down,-lipid_director.down.y.*lipid_length.down,lipid_director.down.z.*lipid_length.down,0,'color',down_color,'ShowArrowHead','off');

axis equal;


figure(6) %plot only diaphragm
hold on
%mid plane surface
surf(mid_plane.x,mid_plane.y,mid_plane.z,'FaceColor','none');
surf(-mid_plane.x,mid_plane.y,mid_plane.z,'FaceColor','none');
surf(mid_plane.x,-mid_plane.y,mid_plane.z,'FaceColor','none');
surf(-mid_plane.x,-mid_plane.y,mid_plane.z,'FaceColor','none');


%dividing plane surface
%up
surf(dividing_plane.up.x,dividing_plane.up.y,dividing_plane.up.z,'FaceColor','none');
surf(-dividing_plane.up.x,dividing_plane.up.y,dividing_plane.up.z,'FaceColor','none');
surf(dividing_plane.up.x,-dividing_plane.up.y,dividing_plane.up.z,'FaceColor','none');
surf(-dividing_plane.up.x,-dividing_plane.up.y,dividing_plane.up.z,'FaceColor','none');

%down
surf(dividing_plane.down.x,dividing_plane.down.y,dividing_plane.down.z,'FaceColor','none');
surf(-dividing_plane.down.x,dividing_plane.down.y,dividing_plane.down.z,'FaceColor','none');
surf(dividing_plane.down.x,-dividing_plane.down.y,dividing_plane.down.z,'FaceColor','none');
surf(-dividing_plane.down.x,-dividing_plane.down.y,dividing_plane.down.z,'FaceColor','none');

%plot tilt
%up
quiver3(dividing_plane.up.x,dividing_plane.up.y,dividing_plane.up.z,tilt.up.x,tilt.up.y,tilt.up.z,'color','r');
quiver3(-dividing_plane.up.x,dividing_plane.up.y,dividing_plane.up.z,-tilt.up.x,tilt.up.y,tilt.up.z,'color','r');
quiver3(dividing_plane.up.x,-dividing_plane.up.y,dividing_plane.up.z,tilt.up.x,-tilt.up.y,tilt.up.z,'color','r');
quiver3(-dividing_plane.up.x,-dividing_plane.up.y,dividing_plane.up.z,-tilt.up.x,-tilt.up.y,tilt.up.z,'color','r');

%down
quiver3(dividing_plane.down.x,dividing_plane.down.y,dividing_plane.down.z,tilt.up.x,tilt.down.y,tilt.down.z,'color','r');
quiver3(-dividing_plane.down.x,dividing_plane.down.y,dividing_plane.down.z,-tilt.up.x,tilt.down.y,tilt.down.z,'color','r');
quiver3(dividing_plane.down.x,-dividing_plane.down.y,dividing_plane.down.z,tilt.up.x,-tilt.down.y,tilt.down.z,'color','r');
quiver3(-dividing_plane.down.x,-dividing_plane.down.y,dividing_plane.down.z,-tilt.up.x,-tilt.down.y,tilt.down.z,'color','r');
axis equal;




figure(7) %plot at phi=0
hold on
plot(mid_plane.x(1,:),mid_plane.z(1,:),'--k')
plot(dividing_plane.up.x(1,:),dividing_plane.up.z(1,:))
plot(dividing_plane.down.x(1,:),dividing_plane.down.z(1,:))
quiver(mid_plane.x(1,:),mid_plane.z(1,:),lipid_director.up.x(1,:).*lipid_length.up(1,:),lipid_director.up.z(1,:).*lipid_length.up(1,:),0,'color',up_color,'ShowArrowHead','off');
quiver(mid_plane.x(1,:),mid_plane.z(1,:),lipid_director.down.x(1,:).*lipid_length.down(1,:),lipid_director.down.z(1,:).*lipid_length.down(1,:),0,'color',down_color,'ShowArrowHead','off');
title('phi=0');
axis equal;
set(gca,'FontSize',18); 
set(gca,'FontWeight','bold'); 
set(gca,'LineWidth',2); 
set(gcf,'color','white')
set(gcf,'color','white')
xlim([0 5]);
ylim([0 5]);

figure(8) %plot at phi=pi/2
hold on
plot(mid_plane.y(res_struct.phi_res,:),mid_plane.z(res_struct.phi_res,:))
plot(dividing_plane.up.y(res_struct.phi_res,:),dividing_plane.up.z(res_struct.phi_res,:))
plot(dividing_plane.down.y(res_struct.phi_res,:),dividing_plane.down.z(res_struct.phi_res,:))
quiver(mid_plane.y(res_struct.phi_res,:),mid_plane.z(res_struct.phi_res,:),lipid_director.up.y(res_struct.phi_res,:).*lipid_length.up(res_struct.phi_res,:)...
    ,lipid_director.up.z(res_struct.phi_res,:).*lipid_length.up(res_struct.phi_res,:),0,'color',up_color,'ShowArrowHead','off');
quiver(mid_plane.y(res_struct.phi_res,:),mid_plane.z(res_struct.phi_res,:),lipid_director.down.y(res_struct.phi_res,:).*lipid_length.down(res_struct.phi_res,:)...
    ,lipid_director.down.z(res_struct.phi_res,:).*lipid_length.down(res_struct.phi_res,:),0,'color',down_color,'ShowArrowHead','off');
title('phi=pi/2');
axis equal;
set(gca,'FontSize',18); 
set(gca,'FontWeight','bold'); 
set(gca,'LineWidth',2); 
set(gcf,'color','white')
set(gcf,'color','white')




end







