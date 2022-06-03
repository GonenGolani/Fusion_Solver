function  plot_profile(Diaphragm,Cell,Virus,Shell,value,angle_indix,flip_video,up_down_symmetry)
%phi index is  angle that will be plotted 1 is x axis y is the
%max(res_struct.phi_res)
Shell_color=[0.7 0.6 0.3];

if flip_video==1
    flip_sign=-1;
else
    flip_sign=1;
end

Factor=4;


%z shift
Diaphragm.mid_plane.z=Diaphragm.mid_plane.z-Diaphragm.z0;
Diaphragm.mid_plane.normal.z=Diaphragm.mid_plane.normal.z-Diaphragm.z0;
Diaphragm.dividing_plane.up.z=Diaphragm.dividing_plane.up.z-Diaphragm.z0;
Diaphragm.dividing_plane.down.z=Diaphragm.dividing_plane.down.z-Diaphragm.z0;
Cell.mid_plane.z=Cell.mid_plane.z-Diaphragm.z0;
Cell.mid_plane.normal.z=Cell.mid_plane.normal.z-Diaphragm.z0;
Cell.dividing_plane.up.z=Cell.dividing_plane.up.z-Diaphragm.z0;
Cell.dividing_plane.down.z=Cell.dividing_plane.down.z-Diaphragm.z0;

Virus.mid_plane.z=Virus.mid_plane.z-Diaphragm.z0;
Virus.mid_plane.normal.z=Virus.mid_plane.normal.z-Diaphragm.z0;
Virus.dividing_plane.up.z=Virus.dividing_plane.up.z-Diaphragm.z0;
Virus.dividing_plane.down.z=Virus.dividing_plane.down.z-Diaphragm.z0;

%flip picture
Diaphragm.lipid_director.up.z=Diaphragm.lipid_director.up.z*flip_sign;
Diaphragm.lipid_director.down.z=Diaphragm.lipid_director.down.z*flip_sign;
Cell.lipid_director.up.z=Cell.lipid_director.up.z*flip_sign;
Cell.lipid_director.down.z=Cell.lipid_director.down.z*flip_sign;
Virus.lipid_director.up.z=Virus.lipid_director.up.z*flip_sign;
Virus.lipid_director.down.z=Virus.lipid_director.down.z*flip_sign;

Diaphragm.mid_plane.z=Diaphragm.mid_plane.z*flip_sign;
Diaphragm.mid_plane.normal.z=Diaphragm.mid_plane.normal.z*flip_sign;
Diaphragm.dividing_plane.up.z=Diaphragm.dividing_plane.up.z*flip_sign;
Diaphragm.dividing_plane.down.z=Diaphragm.dividing_plane.down.z*flip_sign;
Cell.mid_plane.z=Cell.mid_plane.z*flip_sign;
Cell.mid_plane.normal.z=Cell.mid_plane.normal.z*flip_sign;
Cell.dividing_plane.up.z=Cell.dividing_plane.up.z*flip_sign;
Cell.dividing_plane.down.z=Cell.dividing_plane.down.z*flip_sign;
Virus.mid_plane.z=Virus.mid_plane.z*flip_sign;
Virus.mid_plane.normal.z=Virus.mid_plane.normal.z*flip_sign;
Virus.dividing_plane.up.z=Virus.dividing_plane.up.z*flip_sign;
Virus.dividing_plane.down.z=Virus.dividing_plane.down.z*flip_sign;

%upper bound of line
if Shell.Shell_exist==1
Shell.Diaphragm.upper_plane.x=Shell.Diaphragm.mid_plane.x+Shell.Shell_physical_proprties.width*Shell.Diaphragm.normal.x/2;
Shell.Diaphragm.upper_plane.z=(Shell.Diaphragm.mid_plane.z+Shell.Shell_physical_proprties.width*Shell.Diaphragm.normal.z/2-Diaphragm.z0)*flip_sign;

Shell.Diaphragm.lower_plane.x=Shell.Diaphragm.mid_plane.x-Shell.Shell_physical_proprties.width*Shell.Diaphragm.normal.x/2;
Shell.Diaphragm.lower_plane.z=(Shell.Diaphragm.mid_plane.z-Shell.Shell_physical_proprties.width*Shell.Diaphragm.normal.z/2-Diaphragm.z0)*flip_sign;

Shell.Cell.upper_plane.x=Shell.Cell.mid_plane.x+Shell.Shell_physical_proprties.width*Shell.Cell.normal.x/2;
Shell.Cell.upper_plane.z=(Shell.Cell.mid_plane.z+Shell.Shell_physical_proprties.width*Shell.Cell.normal.z/2-Diaphragm.z0)*flip_sign;

Shell.Cell.lower_plane.x=Shell.Cell.mid_plane.x-Shell.Shell_physical_proprties.width*Shell.Cell.normal.x/2;
Shell.Cell.lower_plane.z=(Shell.Cell.mid_plane.z-Shell.Shell_physical_proprties.width*Shell.Cell.normal.z/2-Diaphragm.z0)*flip_sign;


Shell.Junction.upper_plane.x=Shell.Junction.mid_plane.x+Shell.Shell_physical_proprties.width*Shell.Junction.normal.x/2;
Shell.Junction.upper_plane.z=(Shell.Junction.mid_plane.z+Shell.Shell_physical_proprties.width*Shell.Junction.normal.z/2-Diaphragm.z0)*flip_sign;


Shell.Junction.lower_plane.x=Shell.Junction.mid_plane.x-Shell.Shell_physical_proprties.width*Shell.Junction.normal.x/2;
Shell.Junction.lower_plane.z=(Shell.Junction.mid_plane.z-Shell.Shell_physical_proprties.width*Shell.Junction.normal.z/2-Diaphragm.z0)*flip_sign;
end

if angle_indix==1



    Diaphragm_mid_plane_plot=Diaphragm.mid_plane.x;
    Diaphragm_dividing_plane_up_plot=Diaphragm.dividing_plane.up.x;
    Diaphragm_dividing_plane_down_plot=Diaphragm.dividing_plane.down.x;
    
    Cell_mid_plane_plot=Cell.mid_plane.x;
    Cell_dividing_plane_up_plot=Cell.dividing_plane.up.x;
    Cell_dividing_plane_down_plot=Cell.dividing_plane.down.x;
 
    Virus_mid_plane_plot=Virus.mid_plane.x;
    Virus_dividing_plane_up_plot=Virus.dividing_plane.up.x;
    Virus_dividing_plane_down_plot=Virus.dividing_plane.down.x;    
    
    Diaphragm_lipid_director_up_plot=Diaphragm.lipid_director.up.x;
    Diaphragm_lipid_director_down_plot=Diaphragm.lipid_director.down.x;
    
    Cell_lipid_director_up_plot=Cell.lipid_director.up.x;
    Cell_lipid_director_down_plot=Cell.lipid_director.down.x;

    Virus_lipid_director_up_plot=Virus.lipid_director.up.x;
    Virus_lipid_director_down_plot=Virus.lipid_director.down.x;
    
    % shape of virus mid plane out side of fusion site
    max_x=min(Virus.r_bv,2*Virus.Rv);
    vec_virus_MP_down=linspace(Virus.Rv,max_x,1000);
    vec_virus_MP_up=linspace(0,Virus.r_bv,1000);

    z_vec_virus_MP_down=(-((Virus.r_bv+Virus.r_sv)^2-vec_virus_MP_down.^2).^0.5+Virus.r_sv+Virus.r_bv+Virus.hight_above_plane-Diaphragm.z0);
    z_vec_virus_MP_up=(((Virus.r_bv+Virus.r_sv)^2-vec_virus_MP_up.^2).^0.5+Virus.r_sv+Virus.r_bv+Virus.hight_above_plane-Diaphragm.z0);
    
else % just assume is is y
    Diaphragm_mid_plane_plot=Diaphragm.mid_plane.y;
    Diaphragm_dividing_plane_up_plot=Diaphragm.dividing_plane.up.y;
    Diaphragm_dividing_plane_down_plot=Diaphragm.dividing_plane.down.y;
    
    Cell_mid_plane_plot=Cell.mid_plane.y;
    Cell_dividing_plane_up_plot=Cell.dividing_plane.up.y;
    Cell_dividing_plane_down_plot=Cell.dividing_plane.down.y;
 
    Virus_mid_plane_plot=Virus.mid_plane.y;
    Virus_dividing_plane_up_plot=Virus.dividing_plane.up.y;
    Virus_dividing_plane_down_plot=Virus.dividing_plane.down.y;
    
    Diaphragm_lipid_director_up_plot=Diaphragm.lipid_director.up.y;
    Diaphragm_lipid_director_down_plot=Diaphragm.lipid_director.down.y;
    
    Cell_lipid_director_up_plot=Cell.lipid_director.up.y;
    Cell_lipid_director_down_plot=Cell.lipid_director.down.y;

    Virus_lipid_director_up_plot=Virus.lipid_director.up.y;
    Virus_lipid_director_down_plot=Virus.lipid_director.down.y;

    % shape of virus mid plane out side of fusion site
    max_y=min(Virus.r_sv,2*Virus.Rv);
    vec_virus_MP_down=linspace(Virus.Rv,max_y,1000);
    vec_virus_MP_up=linspace(0,Virus.r_sv,1000);

    z_vec_virus_MP_down=(-(Virus.r_sv^2-vec_virus_MP_down.^2).^0.5-Diaphragm.z0+Virus.r_sv+Virus.hight_above_plane);
    z_vec_virus_MP_up=((Virus.r_sv^2-vec_virus_MP_up.^2).^0.5-Diaphragm.z0+Virus.r_sv+Virus.hight_above_plane);
    
end


% shape of cell mid plane out side of fusion site
x_vec_cell_MP_up=linspace(Cell.R_curv,Cell.R_rim,1000);
x_vec_cell_MP_down=linspace(0,Cell.R_curv,1000);

if Cell.R_curv>0
    z_vec_cell_MP_up=((Cell.R_curv^2-x_vec_cell_MP_up.^2).^0.5-(Cell.R_curv^2-Cell.R_rim^2)^0.5-Diaphragm.z0);
    z_vec_cell_MP_down=(-(Cell.R_curv^2-x_vec_cell_MP_down.^2).^0.5-(Cell.R_curv^2-Cell.R_rim^2)^0.5-Diaphragm.z0);
else
    z_vec_cell_MP=(-(Cell.R_curv^2-x_vec_cell_MP.^2).^0.5+(Cell.R_curv^2-Cell.R_rim^2)^0.5-Diaphragm.z0);
end


%down sample

length_cell=sum((diff(Cell_mid_plane_plot(angle_indix,:)).^2+diff(Cell.mid_plane.z(angle_indix,:)).^2).^0.5);
length_virus=sum((diff(Virus_mid_plane_plot(angle_indix,:)).^2+diff(Virus.mid_plane.z(angle_indix,:)).^2).^0.5);


max_rho_dia=max(abs(Diaphragm.rim.x0),abs(Diaphragm.rim.x0));
[max_rho,~]=max([max_rho_dia ; Cell.R_rim-max_rho_dia ; Virus.Rv-max_rho_dia]);
factor_diaphragm=Factor*ceil(max_rho/max_rho_dia); 
factor_cell=ceil(Factor*max_rho/length_cell);
factor_virus=ceil(Factor*max_rho/length_virus);
if 1==1
Diaphragm.mid_plane.z=my_downsample(Diaphragm.mid_plane.z',factor_diaphragm)';
Diaphragm.mid_plane.normal.z=my_downsample(Diaphragm.mid_plane.normal.z',factor_diaphragm)';
Diaphragm.lipid_director.up.z=my_downsample(Diaphragm.lipid_director.up.z',factor_diaphragm)';
Diaphragm.lipid_director.down.z=my_downsample(Diaphragm.lipid_director.down.z',factor_diaphragm)';
Diaphragm.lipid_length.up=my_downsample(Diaphragm.lipid_length.up',factor_diaphragm)';
Diaphragm.lipid_length.down=my_downsample(Diaphragm.lipid_length.down',factor_diaphragm)';
Diaphragm_mid_plane_plot=my_downsample(Diaphragm_mid_plane_plot',factor_diaphragm)';
Diaphragm_dividing_plane_up_plot=my_downsample(Diaphragm_dividing_plane_up_plot',factor_diaphragm)';
Diaphragm_dividing_plane_down_plot=my_downsample(Diaphragm_dividing_plane_down_plot',factor_diaphragm)';
Diaphragm.dividing_plane.up.z=my_downsample(Diaphragm.dividing_plane.up.z',factor_diaphragm)';
Diaphragm.dividing_plane.down.z=my_downsample(Diaphragm.dividing_plane.down.z',factor_diaphragm)';
Diaphragm_lipid_director_up_plot=my_downsample(Diaphragm_lipid_director_up_plot',factor_diaphragm)';
Diaphragm_lipid_director_down_plot=my_downsample(Diaphragm_lipid_director_down_plot',factor_diaphragm)';

Cell.mid_plane.z=my_downsample(Cell.mid_plane.z',factor_cell)';
Cell.mid_plane.normal.z=my_downsample(Cell.mid_plane.normal.z',factor_cell)';
Cell.lipid_director.up.z=my_downsample(Cell.lipid_director.up.z',factor_cell)';
Cell.lipid_director.down.z=my_downsample(Cell.lipid_director.down.z',factor_cell)';
Cell.lipid_length.up=my_downsample(Cell.lipid_length.up',factor_cell)';
Cell.lipid_length.down=my_downsample(Cell.lipid_length.down',factor_cell)';
Cell_mid_plane_plot=my_downsample(Cell_mid_plane_plot',factor_cell)';
Cell_dividing_plane_up_plot=my_downsample(Cell_dividing_plane_up_plot',factor_cell)';
Cell_dividing_plane_down_plot=my_downsample(Cell_dividing_plane_down_plot',factor_cell)';
Cell.dividing_plane.up.z=my_downsample(Cell.dividing_plane.up.z',factor_cell)';
Cell.dividing_plane.down.z=my_downsample(Cell.dividing_plane.down.z',factor_cell)';
Cell_lipid_director_up_plot=my_downsample(Cell_lipid_director_up_plot',factor_cell)';
Cell_lipid_director_down_plot=my_downsample(Cell_lipid_director_down_plot',factor_cell)';

Virus.mid_plane.z=my_downsample(Virus.mid_plane.z',factor_virus)';
Virus.mid_plane.normal.z=my_downsample(Virus.mid_plane.normal.z',factor_virus)';
Virus.lipid_director.up.z=my_downsample(Virus.lipid_director.up.z',factor_virus)';
Virus.lipid_director.down.z=my_downsample(Virus.lipid_director.down.z',factor_virus)';
Virus.lipid_length.up=my_downsample(Virus.lipid_length.up',factor_virus)';
Virus.lipid_length.down=my_downsample(Virus.lipid_length.down',factor_virus)';
Virus_mid_plane_plot=my_downsample(Virus_mid_plane_plot',factor_virus)';
Virus_dividing_plane_up_plot=my_downsample(Virus_dividing_plane_up_plot',factor_virus)';
Virus_dividing_plane_down_plot=my_downsample(Virus_dividing_plane_down_plot',factor_virus)';
Virus.dividing_plane.up.z=my_downsample(Virus.dividing_plane.up.z',factor_virus)';
Virus.dividing_plane.down.z=my_downsample(Virus.dividing_plane.down.z',factor_virus)';
Virus_lipid_director_up_plot=my_downsample(Virus_lipid_director_up_plot',factor_virus)';
Virus_lipid_director_down_plot=my_downsample(Virus_lipid_director_down_plot',factor_virus)';


end

figure(1) 
hold on
%diaphragm

if (Diaphragm.rim.x0^2+Diaphragm.rim.y0^2)>2*0.15^2
plot(Diaphragm_mid_plane_plot(angle_indix,:),Diaphragm.mid_plane.z(angle_indix,:),'--k')
plot(Diaphragm_dividing_plane_up_plot(angle_indix,:),Diaphragm.dividing_plane.up.z(angle_indix,:),'--k')
plot(Diaphragm_dividing_plane_down_plot(angle_indix,:),Diaphragm.dividing_plane.down.z(angle_indix,:),'--k')

quiver(Diaphragm_mid_plane_plot(angle_indix,:),Diaphragm.mid_plane.z(angle_indix,:),Diaphragm_lipid_director_up_plot(angle_indix,:).*Diaphragm.lipid_length.up(angle_indix,:)...
    ,Diaphragm.lipid_director.up.z(angle_indix,:).*Diaphragm.lipid_length.up(angle_indix,:),0,'color','b','ShowArrowHead','off');
quiver(Diaphragm_mid_plane_plot(angle_indix,:),Diaphragm.mid_plane.z(angle_indix,:),Diaphragm_lipid_director_down_plot(angle_indix,:).*Diaphragm.lipid_length.down(angle_indix,:)...
    ,Diaphragm.lipid_director.down.z(angle_indix,:).*Diaphragm.lipid_length.down(angle_indix,:),0,'color','b','ShowArrowHead','off');
end
%Cell
plot(Cell_mid_plane_plot(angle_indix,:),Cell.mid_plane.z(angle_indix,:),'--k')

plot(Cell_dividing_plane_up_plot(angle_indix,:),Cell.dividing_plane.up.z(angle_indix,:),'--k')

plot(Cell_dividing_plane_down_plot(angle_indix,:),Cell.dividing_plane.down.z(angle_indix,:),'--k')
quiver(Cell_mid_plane_plot(angle_indix,:),Cell.mid_plane.z(angle_indix,:),Cell_lipid_director_up_plot(angle_indix,:).*Cell.lipid_length.up(angle_indix,:)...
    ,Cell.lipid_director.up.z(angle_indix,:).*Cell.lipid_length.up(angle_indix,:),0,'color','r','ShowArrowHead','off');
quiver(Cell_mid_plane_plot(angle_indix,:),Cell.mid_plane.z(angle_indix,:),Cell_lipid_director_down_plot(angle_indix,:).*Cell.lipid_length.down(angle_indix,:)...
    ,Cell.lipid_director.down.z(angle_indix,:).*Cell.lipid_length.down(angle_indix,:),0,'color','b','ShowArrowHead','off');
%virus
if up_down_symmetry~=1
plot(Virus_mid_plane_plot(angle_indix,:),Virus.mid_plane.z(angle_indix,:),'--k')
plot(Virus_dividing_plane_up_plot(angle_indix,:),Virus.dividing_plane.up.z(angle_indix,:),'--k')
plot(Virus_dividing_plane_down_plot(angle_indix,:),Virus.dividing_plane.down.z(angle_indix,:),'--k')
quiver(Virus_mid_plane_plot(angle_indix,:),Virus.mid_plane.z(angle_indix,:),Virus_lipid_director_up_plot(angle_indix,:).*Virus.lipid_length.up(angle_indix,:)...
    ,Virus.lipid_director.up.z(angle_indix,:).*Virus.lipid_length.up(angle_indix,:),0,'color','b','ShowArrowHead','off');
quiver(Virus_mid_plane_plot(angle_indix,:),Virus.mid_plane.z(angle_indix,:),Virus_lipid_director_down_plot(angle_indix,:).*Virus.lipid_length.down(angle_indix,:)...
    ,Virus.lipid_director.down.z(angle_indix,:).*Virus.lipid_length.down(angle_indix,:),0,'color','r','ShowArrowHead','off');
else
    plot(Cell_mid_plane_plot(angle_indix,:),-Cell.mid_plane.z(angle_indix,:),'--k')
    
    plot(Cell_dividing_plane_up_plot(angle_indix,:),-Cell.dividing_plane.up.z(angle_indix,:),'--k')
    
    plot(Cell_dividing_plane_down_plot(angle_indix,:),-Cell.dividing_plane.down.z(angle_indix,:),'--k')
    quiver(Cell_mid_plane_plot(angle_indix,:),-Cell.mid_plane.z(angle_indix,:),Cell_lipid_director_up_plot(angle_indix,:).*Cell.lipid_length.up(angle_indix,:)...
        ,-Cell.lipid_director.up.z(angle_indix,:).*Cell.lipid_length.up(angle_indix,:),0,'color','r','ShowArrowHead','off');
    quiver(Cell_mid_plane_plot(angle_indix,:),-Cell.mid_plane.z(angle_indix,:),Cell_lipid_director_down_plot(angle_indix,:).*Cell.lipid_length.down(angle_indix,:)...
        ,-Cell.lipid_director.down.z(angle_indix,:).*Cell.lipid_length.down(angle_indix,:),0,'color','b','ShowArrowHead','off');

end



if Shell.Shell_exist==1

plot(Shell.Diaphragm.mid_plane.x(angle_indix,:),(Shell.Diaphragm.mid_plane.z(angle_indix,:)-Diaphragm.z0)*flip_sign,'--','Color',Shell_color)
plot(Shell.Cell.mid_plane.x(angle_indix,:),(Shell.Cell.mid_plane.z(angle_indix,:)-Diaphragm.z0)*flip_sign,'--','Color',Shell_color)
plot(Shell.Junction.mid_plane.x(angle_indix,:),(Shell.Junction.mid_plane.z(angle_indix,:)-Diaphragm.z0)*flip_sign,'--','Color','b')

fill([Shell.Diaphragm.upper_plane.x(angle_indix,:), flip(Shell.Diaphragm.lower_plane.x(angle_indix,:))],...
[Shell.Diaphragm.upper_plane.z(angle_indix,:) ,flip(Shell.Diaphragm.lower_plane.z(angle_indix,:))],Shell_color,'LineStyle','none','FaceAlpha',0.4);
fill([Shell.Cell.upper_plane.x(angle_indix,:), flip(Shell.Cell.lower_plane.x(angle_indix,:))],...
[Shell.Cell.upper_plane.z(angle_indix,:) ,flip(Shell.Cell.lower_plane.z(angle_indix,:))],Shell_color,'LineStyle','none','FaceAlpha',0.4);
fill([Shell.Junction.upper_plane.x(angle_indix,:), flip(Shell.Junction.lower_plane.x(angle_indix,:))],...
[Shell.Junction.upper_plane.z(angle_indix,:) ,flip(Shell.Junction.lower_plane.z(angle_indix,:))],'r','LineStyle','none','FaceAlpha',0.4);
end


title(num2str(value,'%.4f'));
axis equal;
%shift x 
xl=xlim;
yl=ylim;

if xl(1)~=0
    xl(2)=xl(2)-xl(1)+3;
    xl(1)=0;
    yl(1)=yl(1)-1.5;
    yl(2)=yl(2)+1.5;

    xlim(xl);
    ylim(yl);

end

% plot MD out of fusion site
plot(x_vec_cell_MP_up,z_vec_cell_MP_up*flip_sign,'--k');
plot(x_vec_cell_MP_down,z_vec_cell_MP_down*flip_sign,'--k');

if up_down_symmetry~=1

    plot(vec_virus_MP_down,z_vec_virus_MP_down*flip_sign,'--k');
    plot(vec_virus_MP_up,z_vec_virus_MP_up*flip_sign,'--k');

else
    plot(x_vec_cell_MP_up,-z_vec_cell_MP_up*flip_sign,'--k');
    plot(x_vec_cell_MP_down,-z_vec_cell_MP_down*flip_sign,'--k');
end

set(gca,'FontSize',18); 
set(gca,'FontWeight','bold'); 
set(gca,'LineWidth',2); 
set(gcf,'color','white')




end

