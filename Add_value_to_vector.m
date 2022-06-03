function [new_vector] = Add_value_to_vector(old_vector,zeros_num,after_index,value)
    %this function add cells 'zeros_num' to vector filled with 'value'
    % the cells are added 'after_index'
    
    vector_length=length(old_vector);
    if after_index<vector_length && zeros_num>0 % add at the middle
    new_vector(1:after_index)=old_vector(1:after_index);
    new_vector(after_index+1:after_index+zeros_num)=value;
    new_vector(after_index+zeros_num+1:vector_length+zeros_num)=old_vector(after_index+1:vector_length);
    end
    
    if vector_length<=after_index && zeros_num>0 %add at the end
        new_vector(1:vector_length)=old_vector(1:vector_length);
        new_vector(vector_length+1:vector_length+zeros_num)=value;
    end
    
    if zeros_num==0
        new_vector=old_vector;
    end
end

