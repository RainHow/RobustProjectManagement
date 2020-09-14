clear;
clc;

ratio = 1.2;
instance_type = 2;

outputfile = ['DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro1_all.mat'];
errorfile = ['DATA/DC' num2str(instance_type) '/Error/DC1Error1.mat'];

load(['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance1.mat']);
load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_eo1.mat']);
target = obj_eo * ratio;


whole_solution_t_start = tic;
[ obj_aro, x_aro, y_aro, r_aro, r0_aro, vartheta_aro, rc_aro ] = ARO_all(A_xy, b_xy, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Mu, Sigma,Cor, L, REALIZATION,Am,bm,target, errorfile );
T_aro = toc(whole_solution_t_start);

save(outputfile, 'obj_aro', 'x_aro', 'y_aro', 'DELTA_aro', 'r0_aro', 'vartheta_aro', 'rc_aro', 'T_aro');

