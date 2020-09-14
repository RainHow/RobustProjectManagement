% Solving problem by ARO algorithm
clear;
clc;

%*********************************
instance_type = 2;
%*********************************

for i = 1:240
    
 	inputfile = ['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance' num2str(i) '.mat']; 
    outputfile = ['DATA/DC' num2str(instance_type) '/Instance_Solution/result_eo' num2str(i) '.mat'];
  
    load(inputfile);
    
    whole_solution_t_start = tic;
    [ obj_eo, x_eo, y_eo, DELTA_eo, Iteration_eo] = EO( A_xy, b_xy, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Mu, Sigma, Cor, L, REALIZATION );
    T_eo = toc(whole_solution_t_start);
    

    save(outputfile, 'obj_eo', 'x_eo', 'y_eo', 'DELTA_eo', 'Iteration_eo', 'T_eo' );    
   
end