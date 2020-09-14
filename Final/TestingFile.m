clear;
clc;

%%

instance_type = 2;
Prece = [];

for i = 1:240
      
    A = load(['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance' num2str(i) '.mat']);
    Prece = [ Prece ; A.N_pr];

    
end


%%
% compare solution 'ccg' and 'bd' based on objective function
instance_type = 2;

for i = 1:240
    %load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_ccg'  num2str(i) '.mat']);
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_bd'  num2str(i) '.mat']);
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_ccg_new'  num2str(i) '.mat']);

    if abs(obj_ccg_new-obj_bd) <= 1.0e-7
        disp(['The instance' num2str(i) ' is ok']);
    else
        disp(['The instance' num2str(i) ' is not ok']);        
    end
end

%%
% compare solution 'ccg' and 'bd' based on solution 'x'
instance_type = 2;

for i = 4
    %load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_ccg'  num2str(i) '.mat']);
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_bd'  num2str(i) '.mat']);
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_ccg_new'  num2str(i) '.mat']);

    if any(round(y_bd-y_ccg_new)) == 0
        disp(['The instance' num2str(i) ' is ok']);
    else
        disp(['The instance' num2str(i) ' x is not ok']);        
    end
end
%%
% check sample

% K = 3;
% load([ 'DATA/DC3/Sample_Training.mat' ])

index = zeros(SampleSize,1);

for i = 1:SampleSize
   A = DeltaSample(:,i);
   B = find(A == 6);      % max = 6
   index(i) = length(B);  
end

max(index)
find(index==3)

%%
solution_approach = 'bd';
eval([ 'quantile_index_' solution_approach ' = 0' ])
eval([ 'quantile_' solution_approach ' = cost_' solution_approach '_prob(quantile_index_' solution_approach '+1:end, [1 2] );' ]);
eval([ 'quantile_' solution_approach ' =[ quantile_' solution_approach ', quantile_' solution_approach '(:,2)/ sum(quantile_' solution_approach '(:,2))];' ]);

eval([ 'cost_' solution_approach '_prob = [0,0,0; cost_' solution_approach '_prob];' ]) 
eval([ 'quantile_index_' solution_approach ' = 1;' ])


%%
% compare solution 'ccg' and 'bd'
instance_type = 5;

for i = 1:100
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_ccg'  num2str(i) '.mat']);
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_bd'  num2str(i) '.mat']);
    if abs(obj_ccg-obj_bd) <= 1.0e-7
        disp(['The instance' num2str(i) ' is ok']);
    else
        disp(['The instance' num2str(i) ' is not ok']);        
    end
end


%%
% ??old?new?performance????  
clear all;
clc;
i = 1;
oldfile = ['DATA_Ori/DC1/Instance_Performance/aro_performance' num2str(i) '_beta0.6.mat'];
newfile = ['DATA/DC1/Instance_Performance/aro_performance' num2str(i) '_beta0.6.mat'];

old = load(oldfile);
new = load(newfile);

% ???sample?????

%%
clear;
clc;
Beta = [0.1, 0.6, 1.0];  number_beta = length(Beta);
i = 1; j=1;  beta = Beta(j);
inputfile = ['DATA/DC1/MAT_for_Optimization/instance' num2str(i) '.mat'];
load(inputfile);




%%
clear all;
clc;

i = 1;
tolerant = 0.001;

bench = [];
inputfile = ['DATA/DC1/MAT_for_Optimization/instance' num2str(i) '.mat'];
for j = [0.1 0.6 1]
    inputfile_bench = ['DATA/DC1/Instance_solution/solution' num2str(i) '_beta' num2str(j) '.mat'];
    bench = load(inputfile_bench);
    bench = [bench, bench.obj];
end

load(inputfile);
[obj,obj_sub, x, y, DELTA , UB, LB, Iteration] = CCG( tolerant, A_xy,b_xy,N,N_pr,M,K,Q,bar_n, Ux,Uy,Gf,Gc,d,AF,ACC,AV,AC,AX,AY,AZ,B0,Z,Z0,L, REALIZATION )

obj
bench



%%
%  test whether the optimal point could be in extreme point not bound
x = sdpvar(2,1);
Obj = [1,-2]*x;
Constraints = [ [1;1] <= x <= [2;2]  ];
options = sdpsettings('solver','cplex');
optimize(Constraints, Obj, options);

value(x)
