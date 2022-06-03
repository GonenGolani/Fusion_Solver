function [filled_matrix] = insert_value_to_matrix_col(vector,number_of_terms_to_insert,phi_res,first_zeros)
    int=1;
    filled_matrix=zeros(phi_res,first_zeros+number_of_terms_to_insert);
    while  int<=number_of_terms_to_insert
        filled_matrix(:,first_zeros+int)=vector(int);
        int=int+1;
    end
end

