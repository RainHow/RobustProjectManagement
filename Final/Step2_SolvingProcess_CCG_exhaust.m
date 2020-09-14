% Solving problem by CCG algorithm
clear;
clc;

%*******************************************
instance_type = 2; 
%*******************************************

tolerant = 1e-4;
for i = 1:240
      
    inputfile = ['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance' num2str(i) '.mat'];  
    outputfile = ['DATA/DC' num2str(instance_type) '/Instance_Solution/result_ccg_new' num2str(i) '.mat'];

    load(inputfile);
           
    whole_solution_t_start = tic;
    [obj_ccg_new, x_ccg_new, y_ccg_new, DELTA_ccg_new, UB_ccg_new, LB_ccg_new, Iteration_ccg_new] = CCG_exhaust(tolerant, A_xy,b_xy,N,N_pr,M,K,Q,bar_n, Ux,Uy,Gf,Gc,d,AF,ACC,AV,AC,AX,AY,AZ,B0,Z,Z0,L, REALIZATION);

    T_ccg_new = toc(whole_solution_t_start);
    
    save(outputfile, 'obj_ccg_new', 'x_ccg_new', 'y_ccg_new', 'DELTA_ccg_new', 'T_ccg_new',  'Iteration_ccg_new');
       
end