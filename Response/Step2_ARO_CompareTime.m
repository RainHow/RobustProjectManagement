% Solving problem by ARO algorithm
clear;
clc;

ratio = 1.2;
instance_type = 2;

for i = 70
    
    outputfile1 = ['DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro' num2str(i) '.mat'];
    outputfile2 = ['DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro_rev' num2str(i) '.mat'];
    errorfile = ['DATA/DC' num2str(instance_type) '/Error/DC1Error1.mat'];
    
    load(['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance' num2str(i) '.mat']);
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_eo' num2str(i) '.mat']);
    target = obj_eo * ratio;
       
    whole_solution_t_start = tic;
    [ obj_aro, x_aro, y_aro, DELTA_aro, r_aro, r0_aro, vartheta_aro, rc_aro, Iteration_aro, T_master_aro, T_bilinear_aro ] = ARO(A_xy, b_xy, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Mu, Sigma,Cor, L, REALIZATION,Am,bm,target, errorfile );
    T_aro = toc(whole_solution_t_start);
    
    save(outputfile1, 'obj_aro', 'x_aro', 'y_aro', 'T_aro');
    
    whole_solution_t_start_rev = tic;
    [ obj_aro_rev, x_aro_rev, y_aro_rev, r_aro, r0_aro, vartheta_aro, rc_aro ] = ARO_rev(A_xy, b_xy, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Mu, Sigma,Cor, L, REALIZATION,Am,bm,target, errorfile );
    T_aro_rev = toc(whole_solution_t_start_rev);
    
    save(outputfile2,  'T_aro_rev');
    
    
    
end

