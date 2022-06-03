function [fix_DOF,DOF_vector] = Force_DOF_constreints_old(DOF_vector,Diaphragm,Cell,Virus,Minimazation,res_struc)
% this function fixes some of the DOF given the constrints specified by the user
% the output is the DOF_vector after constreints fixing and fix_DOF: a
% vector conteining 1 in each DOF fixed and zero where it is not fixed

fix_DOF=DOF_vector*0;
int=1;
max_DOF=length(DOF_vector);
while int<=max_DOF

    if Minimazation.axial_symmetry==1 && int==1  % if there is axial symmetry x0=y0
        DOF_vector(1)=Diaphragm.rim.y0;
        fix_DOF(1)=1;
        int=2;
    end

    if Minimazation.up_down_symmetry==1 && int==3  % if there is up_down fix the dipharagm in the middle
        if Minimazation.fix_inter_membrane_distance~=0
            DOF_vector(4)=Minimazation.fix_inter_membrane_distance; 
            fix_DOF(4)=1;
            Virus.hight_above_plane=Minimazation.fix_inter_membrane_distance;
        end
        DOF_vector(3)=Virus.hight_above_plane/2;
        fix_DOF(3)=1; 
        if Minimazation.HD_rim_fixed==0
            DOF_vector(6)=Virus.Rv;      
            fix_DOF(6)=1;
            int=7;            
        else
        int=5; 
        end
    end
    
    if Minimazation.fix_inter_membrane_distance~=0 && int==4
       DOF_vector(4)=Minimazation.fix_inter_membrane_distance; 
       fix_DOF(4)=1;
       Virus.hight_above_plane=Minimazation.fix_inter_membrane_distance;
       int=5;
    end
    
    if Minimazation.HD_rim_fixed>0 && int==5 
        DOF_vector(5)=Minimazation.HD_rim_fixed;
        DOF_vector(6)=Minimazation.HD_rim_fixed;
        
        fix_DOF(5)=1;
        fix_DOF(6)=1;
       
        int=7;
    end
    
    if Minimazation.HD_rim_fixed_virus>0 && int==5
        DOF_vector(5)=Minimazation.HD_rim_fixed_virus;
        fix_DOF(5)=1; 
        int=6;
    end

    if Minimazation.HD_rim_fixed_cell>0 && int==6
        DOF_vector(6)=Minimazation.HD_rim_fixed_cell;
        fix_DOF(6)=1; 
        int=7;
    end
    
    %set the junction tangent angle of the (Virus.rim.drho_dz_vector=Cell.rim.drho_dz_vector)
    if int>=7+res_struc.poly_degree_rim*2 && int<=7+res_struc.poly_degree_rim*3-1 && Minimazation.up_down_symmetry==1
        DOF_vector(7+res_struc.poly_degree_rim*2:7+res_struc.poly_degree_rim*3-1)=-Cell.rim.drho_dz_vector;
        fix_DOF(7+res_struc.poly_degree_rim*2:7+res_struc.poly_degree_rim*3-1)=1;
        int=7+res_struc.poly_degree_rim*3-1;
    end
    
    %tilt decay length
    if Minimazation.up_down_symmetry==1 && int==7+res_struc.poly_degree_rim*4 && (strcmp(Minimazation.step,'first') || Minimazation.axial_symmetry==1) 
        
       DOF_vector(7+res_struc.poly_degree_rim*4)=Diaphragm.Gaussian_decay_vector_down.k_vector(1); % Diaphragm.Gaussian_decay_vector_up.k_vector
       DOF_vector(7+res_struc.poly_degree_rim*4+2)=Cell.Gaussian_decay_vector_up.k_vector(1); % Virus.Gaussian_decay_vector_up.k_vector
       DOF_vector(7+res_struc.poly_degree_rim*4+3)=Cell.Gaussian_decay_vector_down.k_vector(1); % Virus.Gaussian_decay_vector_down.k_vector
       
       fix_DOF(7+res_struc.poly_degree_rim*4)=1;
       fix_DOF(7+res_struc.poly_degree_rim*4+2)=1;
       fix_DOF(7+res_struc.poly_degree_rim*4+3)=1;
       
    end

    
    diaphragm_radius=(Diaphragm.rim.x0^2+Diaphragm.rim.y0^2)^0.5/2^0.5;
    % if the diaphragm is fixed to flat or very small
    if Minimazation.flat_diaphragm==1 || diaphragm_radius<Minimazation.diaphragm_flat_below  
        if int==7 % diaphragm rim 
            DOF_vector(7)=Diaphragm.z0; %the rim hight: hl_vector(1)
            fix_DOF(7)=1; 

            if res_struc.poly_degree_rim==1
                DOF_vector(8)=0; % Diaphragm.rim.drho_dz_vector(1)=0
                fix_DOF(8)=1; 
                int=9;
            else % no axial symmetric case
                DOF_vector(8:8+res_struc.poly_degree_rim-2)=0; %the hight variance: hl_vector(2:poly_degree_rim)
                DOF_vector(7+res_struc.poly_degree_rim : 7+2*res_struc.poly_degree_rim-1)=0; % Diaphragm.rim.drho_dz_vector(1:poly_degree_rim)=0
                
                fix_DOF(8:8+res_struc.poly_degree_rim-2)=1; 
                fix_DOF(7+res_struc.poly_degree_rim : 7+2*res_struc.poly_degree_rim-1)=1; 
                
                int=7+2*res_struc.poly_degree_rim-1;
            end

        end

        %the diaphragm z matrix coefficiants
        if int>=7+4*res_struc.poly_degree_rim+6 && int<=7+4*res_struc.poly_degree_rim+res_struc.poly_degree_z-4+6 ...
                && (strcmp(Minimazation.step,'first') || Minimazation.axial_symmetry==1)
            
            DOF_vector(7+4*res_struc.poly_degree_rim+6:7+4*res_struc.poly_degree_rim+res_struc.poly_degree_z-5+6)=0;
            fix_DOF(7+4*res_struc.poly_degree_rim+6:7+4*res_struc.poly_degree_rim+res_struc.poly_degree_z-5+6)=1; 
            
            int=int+res_struc.poly_degree_z-4;
        end
        


    end
    

    
    int=int+1;
end  


    

    if Minimazation.do_stalk==1
        fix_DOF(1)=1; %x0
        fix_DOF(2)=1; %y0
        DOF_vector(1)=2*res_struc.pore_open_res;
        DOF_vector(1)=2*res_struc.pore_open_res;
        
        % tanget angle at diaphragm rim
        start_int=7+res_struc.poly_degree_rim;
        end_int=start_int+res_struc.poly_degree_rim*3;
        % diaphragm  rim hight derivative 
        fix_DOF(start_int:end_int)=1;     
        DOF_vector(7+res_struc.poly_degree_rim:7+res_struc.poly_degree_rim-1)=
        
       %fix angle in diaphragm if do_stalk==1
       if strcmp(Minimazation.step,'first') || Minimazation.axial_symmetry==1
          start_int=7+res_struc.poly_degree_rim*4+6+(res_struc.poly_degree_z-4)*3;
          end_int=start_int+2*res_struc.poly_degree_LD_rho-8-1;
          fix_DOF(start_int:end_int)=1;
       end
       if strcmp(Minimazation.step,'second') &&  Minimazation.axial_symmetry==0
           fix_DOF(Minimazation.Kill_DOF_vector_if_flat)=1;
       end
    end
    

end

