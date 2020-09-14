
function [obj, x, y, Delta, Lambda, UB, LB, Iteration, master_time] = BD(tolerant, A_xy,b_xy,N,N_pr,M,K,Q,bar_n, Ux,Uy,Gf,Gc,d,AF,ACC,AV,AC,AX,AY,AZ,B0,Z,Z0,L, REALIZATION )
% Benders Dual cutting method
% Algorithm BD

stop_criterion = tolerant;       % tolerant 

limit_extreme_points = 1000;
limit_time = 200;
limit_iteration = 100;

%**** Initialization *********
x = zeros(M,1);  y=zeros(N,1);
LB = -inf;  UB = inf; Iteration = 0;    % bound & record
master_time = []; subproblem_time = [];

Acomb = [AF,ACC,AV,AC];          % combined matrix
Gcomb = [Gf;Gc;d;zeros(N,1)];    % combined matrix

[ ~, Delta, Lambda, UB_oa, LB_oa ] = SubProblem_BD(x,y, Gcomb,Acomb,AX,AY,AZ,B0,Z,Z0, K,L,REALIZATION,limit_extreme_points, limit_time);  % initail Lambda


while UB-LB>tolerant && Iteration<=limit_iteration
    
    Iteration = Iteration +1;    
    % Outer problem solving
    [obj, x, y, eta, time] = MasterProblem_BD(Iteration,Delta,Lambda, A_xy,b_xy,N,M, Ux,Uy,AX,AY,AZ,B0,Z,Z0,limit_time); 
    master_time = [master_time, time];
    LB = obj;   % update lower bound
    [obj_temp, delta_sol , lambda_sol, UB_oa, LB_oa] = SubProblem_BD( x,y, Gcomb,Acomb,AX,AY,AZ,B0,Z,Z0, K,L,REALIZATION, limit_extreme_points, limit_time );
    UB = min(UB , obj_temp + Ux'*x + Uy'*y);
    Delta = [Delta, delta_sol];    Lambda = [Lambda, lambda_sol];  %update
        
end

obj = LB;

end


% ****************** MasterProblem_BD Function ************************
function [ obj, x_sol, y_sol, eta_sol,relaxation_solution_time ] = MasterProblem_BD(Iteration,Delta,Lambda, A_xy,b_xy,N,M, Ux,Uy,AX,AY,AZ,B0,Z,Z0,limit_time )
%MASTERPROBLEM_BD Summary of this function goes here
%   Detailed explanation goes here
  
rome_begin;
h = rome_model('MasterProblem');
x = newvar(M,1,'binary');
y = newvar(N,1,'binary');
eta = newvar(1,1);

rome_minimize( Ux'*x + Uy'*y + eta );

rome_constraint( A_xy * [x;y] <= b_xy );     % original constraint
for i = 1:Iteration
    rome_constraint( eta >= (AX*x+AY*y+AZ*(Z0+Z*Delta(:,i))+B0)'*Lambda(:,i) );
end

relaxation_t_start = tic;
h.solve('CPLEX',limit_time);
relaxation_solution_time = toc(relaxation_t_start);

x_sol = h.eval(x);
y_sol = h.eval(y);
eta_sol = h.eval(eta);
obj = h.objective;

end


% ************************* SubProblem_BD Function ************************
function [ obj, delta_sol, lambda_sol, UB, LB ] = SubProblem_BD( x,y, Gcomb,Acomb,AX,AY,AZ,B0,Z,Z0, K,L,REALIZATION, limit_extreme_points, limit_time )
%SUBPROBLEM_BD Summary of this function goes here

tolerant = 1.0e-4;

Delta = [REALIZATION(:,1)];   % »Ý­n¤@ŸÄinitial delta1
Lambda = [];
LB = -Inf; UB = Inf;    % bound
J = 1;
while J <= limit_extreme_points
    
    [obj_OA1, lambda_sol] = OA1( Delta(:,J),x,y, Gcomb,Acomb,AX,AY,AZ,B0,Z,Z0, limit_time);
    LB = obj_OA1; 
    Lambda = [Lambda, lambda_sol];
    % ****** check whether break out **********
    if UB-LB < tolerant
        break;
    end
    
    [obj_OA2, delta_sol, lambda_sol, eta] = OA2(J,Delta,Lambda,L,K, x,y, Gcomb,Acomb,AX,AY,AZ,B0,Z,Z0, REALIZATION, limit_time);
    UB = obj_OA2;
    Delta = [Delta, delta_sol];
    J = J+1;
    
end

delta_sol = Delta(:,J);
lambda_sol = Lambda(:,J);
obj = LB;



end

function [obj_sol, lambda_sol] = OA1( delta,x,y, Gcomb,Acomb,AX,AY,AZ,B0,Z,Z0, limit_time)

rome_begin;
h = rome_model('OA1');
lambda = newvar(size(Acomb,1),1,'nonneg');
rome_maximize( (AX*x+AY*y+AZ*(Z0+Z*delta)+B0)'*lambda );
rome_constraint(Acomb'*lambda <= Gcomb);

h.solve('CPLEX',limit_time);

lambda_sol = h.eval(lambda);
obj_sol = h.objective;

end

function [obj_sol, delta_sol, lambda_sol, eta_sol ] = OA2(J,Delta,Lambda,L,K, x,y, Gcomb,Acomb,AX,AY,AZ,B0,Z,Z0, REALIZATION, limit_time)

Realization = zeros(sum(L),1);
for k = 1 : K
    Realization(sum(L(1:k-1))+1:sum(L(1:k))) = REALIZATION(k,1:L(k))';
end

rome_begin;
h = rome_model('OA2');

delta = newvar(size(Delta,1),1);
lambda = newvar(size(Lambda,1),1,'nonneg');
eta = newvar(1,1);
aux = newvar(sum(L),1,'binary');    % for discrete case

rome_maximize((AX*x+AY*y+AZ*Z0+B0)'*lambda+eta);

for k = 1:K
    Index_k = sum(L(1:k-1))+1 : sum(L(1:k)); % the index of represent k_th component   
    rome_constraint( sum ( aux(Index_k) ) == 1 );
    delta(k) = aux(Index_k)'*Realization(Index_k);
end

rome_constraint( Acomb'*lambda <= Gcomb );
for i = 1:J
    rome_constraint( eta <= (AZ*Z*Delta(:,i))'*lambda+(AZ*Z*delta)'*Lambda(:,i)-(AZ*Z*Delta(:,i))'*Lambda(:,i) );
end

h.solve('CPLEX',limit_time);

eta_sol = h.eval(eta);
delta_sol = h.eval(delta);
lambda_sol = h.eval(lambda);
obj_sol = h.objective;

end


