% Solving problem by ARO algorithm
clear;
clc;

%*********************************
instance_type = 1;
%*********************************
ratio = 2

for i =  [1721:1800]
  
    outputfile = ['DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro' num2str(i) '.mat'];
    errorfile = ['DATA/DC' num2str(instance_type) '/Error/DC' num2str(instance_type) 'Error' num2str(i) '.mat'];
    
    load(['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance' num2str(i) '.mat']);
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_eo' num2str(i) '.mat']);
    target = obj_eo * ratio;
    
    whole_solution_t_start = tic;
    [ obj_aro, x_aro, y_aro, DELTA_aro, r_aro, r0_aro, vartheta_aro, rc_aro, Iteration_aro, T_master_aro, T_bilinear_aro ] = ARO(A_xy, b_xy, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Mu, Sigma,Cor, L, REALIZATION,Am,bm,target, errorfile );
    T_aro = toc(whole_solution_t_start);
    
    save(outputfile, 'obj_aro', 'x_aro', 'y_aro', 'DELTA_aro', 'r_aro', 'r0_aro', 'vartheta_aro', 'rc_aro', 'Iteration_aro', 'T_aro', 'T_master_aro', 'T_bilinear_aro', 'ratio', 'target', 'N');
    clear obj_eo target
end

% [  ]

% [1:300 361:660 721:1020 1081:1380 1441:1740] Target

% [1472 1473 1474 1475 1476 1479 1492]
% [1494:1527 1529 1531 1532 1533 1607 1629]

% [615 636 637 641 966 979 983 988 1003 1016 1125 1130 1235 1283 1293 1306 1315 1326 1372]
% [ 217 319 321 606 652 867 1607 1620]