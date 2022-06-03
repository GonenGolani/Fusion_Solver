function [vector] = array_to_vector(Array,x_length,y_length)
    int=1;
    while int<=y_length
        vector((int-1)*x_length+1:int*x_length)=Array(int,:);
        int=int+1;
    end

end

