function [new_matrix] = angle_coefficiants_polynom(old_matrix,cofficiant_polynum_matrix,first_columme,last_columme,res_struc,polynom_or_hamonic)
% cofficiant_polynum_matrix: contains the coeffciants of the polynum for each line
% first_columme: the first colummes are not used, specify the first collume to be modified 
% last_columme: the MAX_poly of this matrix

    if first_columme<=last_columme
        new_matrix=old_matrix(1:res_struc.phi_res,1:(first_columme-1)); % first first_columme-1 are not tuoched
        phi_vec=linspace(0,pi/2,res_struc.phi_res); %at this range the polynoms will wwe estimated
        int=1;
        while int<=last_columme-first_columme+1
             
             if strcmp(polynom_or_hamonic,'polynom')
                new_matrix(:,first_columme+int-1)=mypolyval(cofficiant_polynum_matrix(int,:),phi_vec);
             end
             if strcmp(polynom_or_hamonic,'hamonic')
                new_matrix(:,first_columme+int-1)=my_harmonic_expanstion(cofficiant_polynum_matrix(int,:),phi_vec);
             end
            int=int+1;
        end
    else
    new_matrix=old_matrix;
    end



end

