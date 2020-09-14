function [obj, x, y, DELTA, UB, LB, Iteration, master_time, subproblem_time ] = CCG( tolerant, A_xy,b_xy,N,N_pr,M,K,Q,bar_n, Ux,Uy,Gf,Gc,d,AF,ACC,AV,AC,AX,AY,AZ,B0,Z,Z0,L, REALIZATION )
% Column-and-constraint generation algorithm
% Algorithm CCG

stop_criterion = tolerant;       % tolerant 
LB = -inf;  UB = inf; 
J = 1;  Index = [1]; 
DELTA = REALIZATION(:,1)';    %start form first realization
Iteration = 0;
master_time = []; subproblem_time = [];

while UB-LB>tolerant
    % Outer problem solving
    [obj, x, y, eta, time] = MasterProblem_CCG(J,Index,DELTA, A_xy,b_xy,N,N_pr,M,K,Q,bar_n, Ux,Uy,Gf,Gc,d,AF,ACC,AV,AC,AX,AY,AZ,B0,Z,Z0);    
    master_time = [master_time, time];
    LB = obj   % update lower bound
    [obj_sub, delta_sol , p_sol, time_sub] = SubProblem_CCG(x,y, N,N_pr,Q,bar_n,K ,Gf,Gc,d,AF,ACC,AV,AC,AX,AY,AZ,B0,Z,Z0, L,REALIZATION);
    subproblem_time = [subproblem_time, time_sub];
    UB = min(UB , obj_sub + Ux'*x + Uy'*y);
    if obj_sub<Inf
        J = J+1;
        DELTA = [DELTA; delta_sol' ]
        Index = [Index J];
    else 
        J = J+1;
        DELTA = [DELTA; delta_sol'];
    end
    Iteration = Iteration+1;  
end

end


% ************** MasterProblem Function ************************
function [ obj, x_sol, y_sol, eta_sol, time_master ] = MasterProblem_CCG(J,Index,DELTA, A_xy,b_xy,N,N_pr,M,K,Q,bar_n, Ux,Uy,Gf,Gc,d,AF,ACC,AV,AC,AX,AY,AZ,B0,Z,Z0)
%MASTERPROBLEM Summary of this function goes here
%   Detailed explanation goes here

E = size(DELTA,1); 

rome_begin;

h = rome_model('MasterProblem');
  
x = newvar(M,1,'binary');
y = newvar(N,1,'binary');
Wf = newvar(N_pr*E,1,'nonneg');
Wc = newvar(bar_n*E,1,'nonneg');
V = newvar(Q*E,1,'nonneg');
C = newvar(N*E,1,'nonneg');
eta = newvar(1,1,'nonneg');

rome_minimize(Ux'*x + Uy'*y + eta);

rome_constraint(A_xy * [x;y] <= b_xy);
for i = Index
    rome_constraint(eta>=Gf'*Wf((i-1)*N_pr+1:i*N_pr) + Gc'*Wc((i-1)*bar_n+1:i*bar_n) + d'*V((i-1)*Q+1:i*Q) );
end

for j = 1:J
    rome_constraint( AF*Wf((j-1)*N_pr+1:j*N_pr)+ACC*Wc((j-1)*bar_n+1:j*bar_n)+AV*V((j-1)*Q+1:j*Q)+AC*C((j-1)*N+1:j*N) >= AX*x+AY*y+AZ*(Z0+Z*DELTA(j,:)')+B0);
end

start = tic;
h.solve('CPLEX');
time_master = toc(start);

x_sol = h.eval(x);
y_sol = h.eval(y);
eta_sol = h.eval(eta);
obj = h.objective;

end


% ******************* SubProblem Function ************************
function [ obj_sub, delta_sol, p_sol, time_sub ] = SubProblem_CCG(x,y, N,N_pr,Q,bar_n,K, Gf,Gc,d,AF,ACC,AV,AC,AX,AY,AZ,B0,Z,Z0, L,REALIZATION )
%SUBPROBLEM Summary of this function goes here
%   variable is (z)

Mbig = 1.0e+15;  % set big M

B = (AZ*Z)';
c = AX*x + AY*y +AZ*Z0 + B0; 

Realization = zeros(sum(L),1);
for k = 1 : K
    Realization(sum(L(1:k-1))+1:sum(L(1:k))) = REALIZATION(k,1:L(k))';
end

rome_begin;
h = rome_model('SubProblem_CGG');

p = newvar(size(AF,1),1,'nonneg');
Phi = newvar(sum(L),1,'binary');
Zeta = newvar(sum(L),1);  

start = tic;
rome_maximize( c'*p + sum(Zeta) );
time_sub = toc(start);

rome_constraint( [AF ACC AV AC]'*p <= [Gf; Gc; d; zeros(N,1)] );
for k = 1 : K
    Index_k = sum(L(1:k-1))+1 : sum(L(1:k));
    rome_constraint( sum ( Phi(Index_k) ) == 1 );       % only 1 realization
    rome_constraint( Zeta(Index_k) <= (Realization(Index_k)*B(k,:))*p + Mbig*(ones(L(k),1)-Phi(Index_k)) );
    rome_constraint( Zeta(Index_k) <= Mbig*Phi(Index_k) );
end

h.solve('CPLEX');

Phi_sol = round(h.eval(Phi));
Chosen_realization = find(Phi_sol);
delta_sol = Realization(Chosen_realization)
p_sol = h.eval(p);

obj_sub = c'*p_sol + delta_sol'*B*p_sol;

% rome_begin;
% h1 = rome_model('check');
% z = newvar(size(b,1),1,'nonneg');
% rome_minimize(b'*z);
% rome_constraint(G*z >= AX*x+AY*y+AZ*(Z0+Z*delta_sol)+B0 );
% 
% h1.solve('CPLEX');
% 
% obj = h1.objective;
% 
% if obj~=obj_sub
%     error('check');
% end

end





