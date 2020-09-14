% Solving problem by ARO algorithm
clear;
clc;

%*********************************
instance_type = 2;
%*********************************
ratio = 1.2;

for i = 1:240
    
    outputfile = ['DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro' num2str(i) '.mat'];
    errorfile = ['DATA/DC' num2str(instance_type) '/Error/DC' num2str(instance_type) 'Error' num2str(i) '.mat'];
    
    load(['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance' num2str(i) '.mat']);
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_eo' num2str(i) '.mat']);
    target = obj_eo * ratio;
    
    whole_solution_t_start = tic;
    [ obj_aro, x_aro, y_aro, DELTA_aro, r_aro, r0_aro, vartheta_aro, rc_aro, Iteration_aro, T_master_aro, T_bilinear_aro ] = ARO(A_xy, b_xy, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Mu, Sigma,Cor, L, REALIZATION,Am,bm,target, errorfile );
    T_aro = toc(whole_solution_t_start);
    
    save(outputfile, 'obj_aro', 'x_aro', 'y_aro', 'DELTA_aro', 'r_aro', 'r0_aro', 'vartheta_aro', 'rc_aro', 'Iteration_aro', 'T_aro', 'T_master_aro', 'T_bilinear_aro', 'ratio', 'target');
    clear obj_eo target
end
