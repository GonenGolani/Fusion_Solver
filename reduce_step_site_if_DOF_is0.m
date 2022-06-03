function [step_size_vector] = reduce_step_site_if_DOF_is0(step_size_vector,DOF_vector)
    % find where DOF=0
    leng=length(DOF_vector);
    int=1;
    while int<=leng
        if DOF_vector(int)==0
            step_size_vector(int)=step_size_vector(int)/10;
        end
        int=int+1;
    end
end

