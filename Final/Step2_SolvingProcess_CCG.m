% Solving problem by CCG algorithm
clear;
clc;

%*******************************************
instance_type = 2; 
%*******************************************

tolerant = 1e-4;
for i = 152
      
    inputfile = ['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance' num2str(i) '.mat'];  
    outputfile = ['DATA/DC' num2str(instance_type) '/Instance_Solution/result_ccg' num2str(i) '.mat'];

    load(inputfile);
           
    whole_solution_t_start = tic;
     [obj_ccg, x_ccg, y_ccg, DELTA_ccg, UB_ccg, LB_ccg, Iteration_ccg, T_master_ccg, T_subproblem_ccg] = CCG(tolerant, A_xy,b_xy,N,N_pr,M,K,Q,bar_n, Ux,Uy,Gf,Gc,d,AF,ACC,AV,AC,AX,AY,AZ,B0,Z,Z0,L, REALIZATION);
%     [obj_ccg, x_ccg, y_ccg, DELTA_ccg, UB_ccg, LB_ccg, Iteration_ccg] = CCG_exhaust(tolerant, A_xy,b_xy,N,N_pr,M,K,Q,bar_n, Ux,Uy,Gf,Gc,d,AF,ACC,AV,AC,AX,AY,AZ,B0,Z,Z0,L, REALIZATION);

    T_ccg = toc(whole_solution_t_start);
    
   save(outputfile, 'obj_ccg', 'x_ccg', 'y_ccg', 'DELTA_ccg', 'T_ccg', 'T_master_ccg', 'T_subproblem_ccg', 'Iteration_ccg');
       
end