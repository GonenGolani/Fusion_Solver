function [Diaphragm,Virus,Cell,Shell,General_physical_properties,Minimazation,res_struc,intital_tilt_decay_length] = read_conditions(~)

cd('Parameters')
parmeters_path=pwd;

cd('Physical properties');
%general proprties
General_physical_properties.lipid_length=importdata('Lipid length.txt'); %[nm]
General_physical_properties.streatching_modulus=importdata('stretching modulus.txt'); %kappa/nm^2
General_physical_properties.tilt_modulus=importdata('Tilt modulus.txt'); %kappa/nm^2
General_physical_properties.chi=importdata('chi.txt'); %kappa_bar/kappa
General_physical_properties.Virus_distal_J0=importdata('Virus distal J0.txt'); %kappa_bar/kappa
General_physical_properties.Cell_distal_J0=importdata('Cell distal J0.txt'); %kappa_bar/kappa
General_physical_properties.proximal_J0=importdata('Proximal J0.txt'); %kappa_bar/kappa
General_physical_properties.cell_distal_kappa=importdata('cell distal kappa.txt'); %cell distal monolayer bending rigidity
General_physical_properties.virus_distal_kappa=importdata('virus distal kappa.txt'); %virus distal monolayer bending rigidity
General_physical_properties.proximal_kappa=importdata('proximal kappa.txt'); %merged proximal monolayer bending rigidity
General_physical_properties.Splay_grad_square_modulus=importdata('Splay grad square modulus.txt');
General_physical_properties.Splay_grad_tilt_modulus=importdata('Splay grad tilt modulus.txt');
General_physical_properties.Surface_tension=importdata('Surface_tension.txt');




cd(parmeters_path);


%cell proprties
cd('System dimensions');
Cell.R_curv=importdata('Cell radius.txt');  % curvature of cell

Virus.r_sv=importdata('Virus small radius.txt'); % small radius of torus
Virus.r_bv=importdata('Virus big radius.txt'); % large radius of torus

cd(parmeters_path);


%% minimazation parameters
cd('Minimization');
res_struc.phi_res=importdata('angular resolution.txt');  %how many points sampelled in phi between 0 to pi/2
res_struc.rho_res=importdata('radial resolution.txt');  % number of points in rho direction
res_struc.poly_degree_rim=importdata('poly degree rim.txt');   % number of coefficiants in expanstion of the hight z ###min 1
res_struc.poly_degree_z=importdata('poly degree midplane.txt');   % number of coefficiants in expanstion of the hight z ###min 4
res_struc.poly_degree_LD_rho=importdata('poly degree Lipid director rho.txt');   % number of coefficiants in expanstion of lipid director rho direction ###min 2
res_struc.poly_degree_LD_phi=importdata('poly degree Lipid director phi.txt');   % number of coefficiants in expanstion of lipid director phi direction ###min 4
res_struc.pore_open_res=importdata('Pore size resolution.txt');       % pore size resolution in nm
res_struc.polar_angle_polynom=importdata('polar angle polynom.txt');  % only used if non_axial_symmetric_polynom==1, the order of the polynom


Minimazation.Tolerance=importdata('Energy tolerance.txt'); 
Minimazation.MAX_ITR=importdata('Max number of iterations.txt'); 
Minimazation.MAX_DOF_MP=importdata('Max mid plane degree of freedom.txt'); 
Minimazation.MAX_DOF_rim=importdata('Max rim degree of freedom.txt'); 
Minimazation.MAX_DOF_angle_RHO=importdata('Max angle rho degree of freedom.txt'); 
Minimazation.MAX_DOF_angle_PHI=importdata('Max angle phi degree of freedom.txt'); 
Minimazation.max_step=importdata('max step size.txt');
Minimazation.temperature=importdata('start annealing temperature.txt');
Minimazation.number_of_repeats_each_step=importdata('repeats at each step.txt');


cd(parmeters_path);


cd('Handls');
Minimazation.up_down_symmetry=importdata('up down symmetry.txt'); %0 is no 1 is yes 
Minimazation.axial_symmetry=importdata('axial symmetry.txt'); %0 is no 1 is yes 
Minimazation.flat_diaphragm=importdata('flat diaphragm.txt'); %0 is no 1 is yes 
Minimazation.diaphragm_flat_below=importdata('flat diaphragm if below.txt'); % diaphragm size if below is flat
Minimazation.HD_rim_fixed=importdata('HD rim fixed.txt'); % 0 if not fixed. if other number than fixed to that
Minimazation.HD_rim_fixed_cell=importdata('HD rim fixed cell.txt'); % 0 if not fixed. if other number than fixed to that. only cell side is fixed
Minimazation.HD_rim_fixed_virus=importdata('HD rim fixed virus.txt'); % 0 if not fixed. if other number than fixed to that. only virus side is fixed
Minimazation.do_iterative=importdata('iterative calculation.txt'); %0 is no 1 is yes 
Minimazation.relaxed_spherical_VLP=importdata('relaxed spherical VLP.txt'); %0 is no 1 is yes 
Minimazation.relaxed_cylindrical_VLP=importdata('relaxed cylindrical VLP.txt'); %0 is no 1 is yes 
Minimazation.symmetric_J0=importdata('symmetric J0.txt'); %0 is no 1 is yes 
Minimazation.symmetric_distal_J0=importdata('symmetric distal J0.txt'); %0 is no anyother number is the J0 of both distal monolayrs

Minimazation.non_axial_symmetric_polynom=importdata('smooth non_axial symmetric polynom.txt'); % if 1 use polynom approximation to describe 
Minimazation.exponential_decay_tilt=importdata('exponential decay tilt.txt'); % 0 if unsing only polynom to describe if 1 multiply the polynom with exponent. the exponent is a new DOF
Minimazation.fix_Cell_volume=importdata('fix Cell volume.txt'); % 0 if not fixed, any other value will fix the volume by adjusting the Cell radius
Minimazation.fix_inter_membrane_distance=importdata('fix inter membrane distance.txt'); % 0 if not fixed, any other will fix the value

Minimazation.Cell_trans_membrane_added_proprties_exist=importdata('Cell trans membrane added proprties.txt'); %1 if there are added proprties to the cell membrane, 0 if not used
Minimazation.sorting_protein_in_cell_membrane=char(importdata('Add sorting inducing in cell membrane.txt'));   %1 if allow lipid sorting betweeen distal monolayer of cell and diaphragm, 0 if not used
Minimazation.do_stalk=importdata('do stalk.txt');   % if 1 the diaphragm radius size effectivaly zero and it is not included in energy 

Minimazation.full_lipid_flip_flop=importdata('lipid full flip flop.txt'); 
Minimazation.Volume_fixed_based_on_Rc=importdata('Volume_fixed_based_on_Rc.txt');
  
Minimazation.tension_distance_coupleing.exist=importdata('tension_distance_coupleing_exist.txt'); 
Minimazation.tension_distance_coupleing.kappa_b=importdata('tension_distance_coupleing_kappa_b.txt'); 
Minimazation.tension_distance_coupleing.Pext=importdata('tension_distance_coupleing_Pext.txt'); 



cd(parmeters_path);


cd('Shell')
Shell.Shell_exist=importdata('exist shell.txt'); % if 1 use polynom approximation to describe
Minimazation.Shell_exist=Shell.Shell_exist;
res_struc.rho_res_Junction=importdata('radial resolution junction.txt');
res_struc.poly_degree_shell=importdata('poly degree shell.txt');

Shell.Shell_physical_proprties.width=importdata('Shell width.txt');
Shell.Shell_physical_proprties.interaction_distance=importdata('membrane-shell interaction distance.txt');  %in nm
Shell.Shell_physical_proprties.intercation_potential_energy_density=importdata('membrane-shell interaction energy density.txt'); % in kappa/nm^2
Shell.Shell_physical_proprties.Youngs_modulus=importdata('Shell Youngs modulus.txt'); %
Shell.Shell_physical_proprties.J0=importdata('Shell J0.txt'); 
Shell.Shell_physical_proprties.max_shell_curvature=importdata('max shell curvature.txt'); 
Shell.Shell_physical_proprties.P_ratio=importdata('shell poisson ratio.txt'); 
Shell.Shell_physical_proprties.Shell_membrane_interaction=char(importdata('Shell membrane interaction type.txt')); 



res_struc.poly_degree_shell=importdata('poly degree shell.txt');  % polynom dgree of the shell hight function

cd('Shell initial conditions')
Shell.Diaphragm.z_i=importdata('Shell Diphargm zi.txt');
Shell.Diaphragm.z_f=importdata('Shell Diphargm zf.txt');
Shell.Diaphragm.dzdrho_edge=importdata('Shell Diphargm dzdrho.txt');
Shell.Diaphragm.length_from_edge=importdata('Shell Diphargm length from edge.txt');

Shell.Cell.z_i=importdata('Shell Cell zi.txt');
Shell.Cell.z_f=importdata('Shell Cell zf.txt');
Shell.Cell.dzdrho_edge=importdata('Shell Cell dzdrho.txt');
Shell.Cell.length_from_edge=importdata('Shell Cell length from edge.txt');

cd(parmeters_path);


if strcmp(Minimazation.sorting_protein_in_cell_membrane,'affinity_based')...
        || strcmp(Minimazation.sorting_protein_in_cell_membrane,'splay_based')...
        || Minimazation.Cell_trans_membrane_added_proprties_exist==1

    cd('Handls');
    Minimazation.no_sorting=importdata('no sorting.txt'); 
    Minimazation.no_TM_direct_intrinsic_curvature=importdata('no TM direct intrinsic curvature.txt'); 

    cd(parmeters_path);

    cd('TM protein Cell')
    General_physical_properties.Cell_mid_plane_physical_proprties.kappa=importdata('Cell trans memb kappa.txt');
    General_physical_properties.Cell_mid_plane_physical_proprties.kappa_bar=importdata('Cell trans memb kappa_bar.txt');
    General_physical_properties.Cell_mid_plane_physical_proprties.J0=importdata('Cell trans memb J0.txt');

    General_physical_properties.Cell_mid_plane_physical_proprties.sorting.affinity=importdata('affinity.txt');
    General_physical_properties.Cell_mid_plane_physical_proprties.sorting.number_of_intercating_lipids=importdata('number_of_intercating_lipids.txt');


    General_physical_properties.Cell_mid_plane_physical_proprties.sorting.kappa_m_kbT_a=importdata('kappa_m_kbT_a.txt');
    General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_s=importdata('j_s.txt');
    General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_0=importdata('j_0.txt');
    General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_TM=importdata('j_TM.txt');
    General_physical_properties.Cell_mid_plane_physical_proprties.sorting.j_TM_luminal=importdata('j_TM_luminal.txt');
    General_physical_properties.Cell_mid_plane_physical_proprties.sorting.ratio_s_to_0=importdata('ratio_s_to_0.txt'); % ratio between the mole ratio of the two lipids
    General_physical_properties.Cell_mid_plane_physical_proprties.sorting.phi_TM=importdata('phi_TM.txt');  %mole ratio of TM protein in cell membrane

    
cd(parmeters_path);
end


if Minimazation.symmetric_J0==1
    General_physical_properties.Virus_distal_J0=General_physical_properties.proximal_J0;
    General_physical_properties.Cell_distal_J0=General_physical_properties.proximal_J0;
end

if Minimazation.symmetric_distal_J0~=0
    General_physical_properties.Virus_distal_J0=Minimazation.symmetric_distal_J0;
    General_physical_properties.Cell_distal_J0=Minimazation.symmetric_distal_J0;  
end

cd(parmeters_path);

if Minimazation.do_iterative==1
    cd('Scan variable');
    Minimazation.scan_variable.target_value=importdata('target value.txt');
    Minimazation.scan_variable.delta_value=importdata('delta value.txt');
    Minimazation.scan_variable.variable_name=char(importdata('varible name.txt'));
    if exist('value_n.txt','file')==2
        Minimazation.scan_variable.value_n=importdata('value_n.txt');
    else
        Minimazation.scan_variable.value_n=1; % just linear dependence
    end
    
    if exist('target value 2.txt','file')==2 % all the rest added epoch data in same version, no need to check
        Minimazation.scan_variable.target_value_2=importdata('target value 2.txt');
        Minimazation.scan_variable.target_value_3=importdata('target value 3.txt');
        Minimazation.scan_variable.delta_value_2=importdata('delta value 2.txt');
        Minimazation.scan_variable.delta_value_3=importdata('delta value 3.txt');

    else
        Minimazation.scan_variable.target_value_2=0;
        Minimazation.scan_variable.target_value_3=0;
        Minimazation.scan_variable.delta_value_2=0;
        Minimazation.scan_variable.delta_value_3=0;    
    end    
    
else
    Minimazation.scan_variable=0;   
end
cd(parmeters_path);


%% inital DOF guess
Cell.rim.drho_dz_vector=zeros(1, res_struc.poly_degree_rim); %create vector containing the coeffciants of the Diaphragm angle sum(h_l*cos(2*l*phi) 
Virus.rim.drho_dz_vector=zeros(1, res_struc.poly_degree_rim); %create vector containing the coeffciants of the Diaphragm angle 
Diaphragm.rim.hl_vector=zeros(1,res_struc.poly_degree_rim); %create vector containing the coeffciants of the Diaphragm z
Diaphragm.rim.drho_dz_vector=zeros(1,res_struc.poly_degree_rim); %create vector containing the coeffciants of the Diaphragm angle 


cd('Initial configuration')
Diaphragm.r_pore=importdata('initial pore radius.txt');  % pore radius
Diaphragm.rim.x0=importdata('Diaphragm x0.txt');  %ellipse size in x direction 
Diaphragm.rim.y0=importdata('Diaphragm y0.txt'); %ellipse size in y direction 
Diaphragm.rim.hl_vector(1)=importdata('Diaphragm z0.txt'); % set Diaphragm hight
Diaphragm.z0=importdata('Diaphragm z0.txt'); % set Diaphragm hight rho=0 
Virus.hight_above_plane=importdata('virus hight above cell.txt'); %hight of torus edge above z=0
Virus.Rv=importdata('virus fusion site rim radius.txt'); % fusion site size in torus rim
Cell.R_rim=importdata('cell fusion site rim radius.txt');  %fusion site size at cell interface rim
intital_tilt_decay_length=importdata('intital tilt decay length.txt'); %exponenesial tilt decay rate
Cell.rim.drho_dz_vector(1)=-importdata('bifuraction tanagent.txt');  %1 is 45 deg; 
Virus.rim.drho_dz_vector(1)=importdata('bifuraction tanagent.txt'); 
cd(parmeters_path);

if Minimazation.axial_symmetry==1
    Diaphragm.rim.x0=Diaphragm.rim.y0;
end

if Minimazation.do_stalk==1
    Diaphragm.rim.x0=2*res_struc.pore_open_res;
    Diaphragm.rim.y0=2*res_struc.pore_open_res;
    Cell.rim.drho_dz_vector(1)=-1;
    Virus.rim.drho_dz_vector(1)=1;

end

end

