function [Cell_Radius] = force_cell_volume_constraint(Cell_old,Diaphragm,Minimazation,res_struc)
    % this function change the cell radius so that the intial volume is
    % conserved
    Volume_error=zeros(1,Minimazation.MAX_ITR);
    Cell=Cell_old;
    res_factor=1;
    Cell_Radius=Cell.R_curv;
    int=1;
    while int<=Minimazation.MAX_ITR

        [Cell,Diaphragm] = recalc_only_cell_mid_plane(Cell,Diaphragm,Minimazation,res_struc);
        volume_target=Minimazation.fix_Cell_volume;
        volume_now=find_Cell_volume(Cell,Diaphragm,res_struc);
        delta_volume=volume_target-volume_now;
        Volume_error(int)=abs(delta_volume/volume_target);
        if Volume_error(int)<min(Minimazation.Tolerance*100,0.002) % break if target is achived
            return;
        end

        if delta_volume>0
            change_Rc_vector=linspace(0,abs(delta_volume)^(1/3),1000*res_factor);
        else
            change_Rc_vector=linspace(-(-abs(delta_volume))^(1/3),0,1000*res_factor);
        end

        %find cell area out side of Fusion site assume that area is unchanged
        %in the fusion site
        new_volume_bulk=2*pi/3.*(1+(1-(Cell.R_rim./(Cell.R_curv+change_Rc_vector)).^2).^0.5).*(Cell.R_curv+change_Rc_vector).^3;
        old_volume_bulk=2*pi/3*(1+(1-(Cell.R_rim./(Cell.R_curv)).^2).^0.5).*(Cell.R_curv)^3;

        [~,index]=min(abs(new_volume_bulk-old_volume_bulk-delta_volume));

        %update Cell radius to new optimum 
        Cell.R_curv=abs(Cell.R_curv+change_Rc_vector(index));
        Cell_Radius=Cell.R_curv;
        int=int+1;
        if int>=3
            if abs(Volume_error(int-1)-Volume_error(int-2))/volume_target<Minimazation.Tolerance*100
                res_factor=res_factor*4;
            end
            if Volume_error(int-1)==Volume_error(int-2)
                break;
            end
            if res_factor>100
                break;
            end
        end
        
    end
    to_print=['Volume not fixed, error: ' num2str(Volume_error(int-1),4)];
    disp(to_print);
end

