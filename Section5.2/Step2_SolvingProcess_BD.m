% Solving problem by BD_approx algorithm
clear;
clc
%*************************************
instance_type = 2;
%*************************************

tolerant = 1.0e-4;     
for i = 1:240
    
 	inputfile = ['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance' num2str(i) '.mat'];
    outputfile = ['DATA/DC' num2str(instance_type) '/Instance_Solution/result_bd' num2str(i) '.mat'];

    load(inputfile);
          
    whole_solution_t_start = tic;
    [obj_bd, x_bd, y_bd, DELTA_bd, LAMBDA_bd, UB, LB, Iteration_bd, T_master_bd] = BD(tolerant, A_xy,b_xy,N,N_pr,M,K,Q,bar_n, Ux,Uy,Gf,Gc,d,AF,ACC,AV,AC,AX,AY,AZ,B0,Z,Z0,L, REALIZATION);
    T_bd = toc(whole_solution_t_start);
    
    save(outputfile, 'obj_bd', 'x_bd', 'y_bd', 'DELTA_bd', 'LAMBDA_bd', 'T_bd', 'T_master_bd', 'Iteration_bd');

        
end