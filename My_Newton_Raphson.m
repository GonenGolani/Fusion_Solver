function [Diaphragm,Cell,Virus,DOF_vector,final_Energy,final_step_size] = My_Newton_Raphson...
    (Diaphragm,Cell,Virus,DOF_vector,index,res_struc,initial_steps_size,Minimazation,initial_energy,General_physical_properties)
% this function minimize the energy using the Newton Raphson method X(n+1)=X(n)-f'(Xn)/f''(Xn) 
% it takes on paramter in DOF_vector(index) and minimize the energy up tp Tolerance
% the function recives the inital step size but not its direction, which is
% chosen randomly.
% the output is the new stracture, final energy and the new DOF_vector
% vector
%initilize vectors
energy_vec=NaN(1,Minimazation.MAX_ITR);
steps_size=NaN(1,Minimazation.MAX_ITR);
f_tag=NaN(1,Minimazation.MAX_ITR);
f_tag_tag=NaN(1,Minimazation.MAX_ITR);
pos_vec=NaN(1,Minimazation.MAX_ITR);

% the energy vector record the energy at each step of the iteration
energy_vec(1)=initial_energy;
% the position vector record the DOF value at each step
steps_size(1:2)=initial_steps_size;
%first dirivitive of energy
f_tag(1)=0;
%second dirivitive of energy
f_tag_tag(1:2)=0;
%value of DOF at each iteration
pos_vec(1)=DOF_vector(index);

    int=2; % start from second step
    while int<=Minimazation.MAX_ITR %run until hit MAX_ITR which is global paramter
        DOF_vector(index)=DOF_vector(index)+steps_size(int-1); % change paramter
        pos_vec(int)=DOF_vector(index); 
        % re calc shape and energy
        [Diaphragm,Cell,Virus] = DOF_vector_to_parameters_no_phi(Diaphragm,Cell,Virus,DOF_vector,res_struc);
        [Diaphragm,Cell,Virus] = build_Hemifusion_structure(Diaphragm,Cell,Virus,res_struc,General_physical_properties);
        [Diaphragm,Cell,Virus,~,energy_vec(int)] = find_HD_energy(Diaphragm,Cell,Virus,res_struc);
        %find the energy gradient
         f_tag(int)=(energy_vec(int)-energy_vec(int-1))/steps_size(int);

        
        % if it is the first step size is the same
        if (int<=4)
            if (energy_vec(int)<energy_vec(int-1)) % if the energy decreases keep direction
                steps_size(int+1)=steps_size(int);
            else % if the energy increases flip direction and do half step next time (to avoid same position)
                steps_size(int+1)=-steps_size(int)/2;
            end
        else % if it is second step or later
            f_tag_tag(int)=(f_tag(int)-f_tag(int-1))/steps_size(int);
            % find next step with Newton-Raphson
            steps_size(int+1)=-f_tag_tag(int)/f_tag(int);
        end
        
        %stop condition
        if abs(energy_vec(int)-energy_vec(int-1))<Minimazation.Tolerance
            final_Energy=energy_vec(int);
            final_step_size=steps_size(int);
            break;
        end
        
        int=int+1;
    end


end

