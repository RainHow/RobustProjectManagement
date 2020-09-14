clear% Deterministic Case
clear;
clc;

%******************************
instance_type = 2;
%******************************

for i = 1:240
    
    inputfile = ['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance' num2str(i) '.mat'];
    outputfile = ['DATA/DC' num2str(instance_type) '/Instance_Solution/result_um' num2str(i) '.mat'];

    load(inputfile);
   
    rome_begin;
    h = rome_model('UseMean');
    x = newvar(M,1,'binary');
    y = newvar(N,1,'binary');
    Wf = newvar(N_pr,1,'nonneg');
    Wc = newvar(bar_n,1,'nonneg');
    V = newvar(Q,1,'nonneg');
    C = newvar(N,1,'nonneg');
    
    rome_minimize( Ux'*x + Uy'*y + Gf'*Wf + Gc'*Wc + d'*V );
    rome_constraint( AF*Wf + ACC*Wc + AV*V + AC*C >= AX*x + AY*y + AZ*(Z0+Z*Mu) + B0 );
    rome_constraint( A_xy * [x;y] <= b_xy);
    
    whole_solution_t_start = tic;
    h.solve('CPLEX');
    T_um = toc(whole_solution_t_start);
    
    obj_um = h.objective;
    x_um = h.eval(x);
    y_um = h.eval(y); 
    
    Iteration_um = 1;
    save(outputfile, 'obj_um', 'x_um', 'y_um', 'T_um', 'Iteration_um');
    
end