function [Array] = vector_to_array(vector,x_length,y_length)
    Array=NaN(y_length,x_length);
    int=1;
    while int<=y_length
        Array(int,:)=vector((int-1)*x_length+1:int*x_length);
        int=int+1;
    end
end

