function [new_array] = double_array(old_array)
% this function takes the old array and plug in new lines, each new line
% equales to the last one. each line is doubled
    old_array_line=length(old_array(1,:));
    old_array_cole=length(old_array(:,1));

    new_array=zeros(2*old_array_cole,old_array_line);
    int=1;
    while int<=old_array_cole
        new_array(2*int-1,:)=old_array(int,:);
        new_array(2*int,:)=old_array(int,:);
        int=int+1;
    end
end

