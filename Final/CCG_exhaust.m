function [obj_ms, x, y, DELTA, UB, LB, Iteration ] = CCG_exhasut( tolerant, A_xy,b_xy,N,N_pr,M,K,Q,bar_n, Ux,Uy,Gf,Gc,d,AF,ACC,AV,AC,AX,AY,AZ,B0,Z,Z0,L, REALIZATION )
% Column-and-constraint generation algorithm
% Algorithm CCG

limit_time = 200;
limit_iteration = 1000;

LB = -inf;  UB = inf;         % bound 
J = 1;  Index = [1]; 
DELTA = REALIZATION(:,3)';    %start form first realization
Iteration = 0;

while UB-LB>tolerant && Iteration<=limit_iteration
    % Outer problem solving
    [obj_ms, x, y, eta] = MasterProblem_CCG_exhaust(J,Index,DELTA, A_xy,b_xy,N,N_pr,M,K,Q,bar_n, Ux,Uy,Gf,Gc,d,AF,ACC,AV,AC,AX,AY,AZ,B0,Z,Z0, limit_time);    
    LB = obj_ms;   % update lower bound
    [obj_sub, delta_sol] = SubProblem_CCG_exhaust(x,y, N,N_pr,Q,bar_n,K ,Gf,Gc,d,AF,ACC,AV,AC,AX,AY,AZ,B0,Z,Z0, L,REALIZATION,DELTA, limit_time);
    UB = min(UB , obj_sub + Ux'*x + Uy'*y);
    if obj_sub<Inf
        J = J+1;
        DELTA = [DELTA; delta_sol' ];
        Index = [Index J];
    else 
        J = J+1;
        DELTA = [DELTA; delta_sol'];
    end
    Iteration = Iteration+1;  
    DELTA
end

end


function [ obj, x_sol, y_sol, eta_sol ] = MasterProblem_CCG_exhaust(J,Index,DELTA, A_xy,b_xy,N,N_pr,M,K,Q,bar_n, Ux,Uy,Gf,Gc,d,AF,ACC,AV,AC,AX,AY,AZ,B0,Z,Z0,limit_time)
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
eta = newvar(1,1);

rome_minimize(Ux'*x + Uy'*y + eta);

rome_constraint(A_xy * [x;y] <= b_xy);
for i = 1:J
    rome_constraint( AF*Wf((i-1)*N_pr+1:i*N_pr)+ACC*Wc((i-1)*bar_n+1:i*bar_n)+AV*V((i-1)*Q+1:i*Q)+AC*C((i-1)*N+1:i*N) >= AX*x+AY*y+AZ*(Z0+Z*DELTA(i,:)')+B0);
    if any(i==Index)
        rome_constraint(eta >= Gf'*Wf((i-1)*N_pr+1:i*N_pr) + Gc'*Wc((i-1)*bar_n+1:i*bar_n) + d'*V((i-1)*Q+1:i*Q) );
    end
end

h.solve('CPLEX', limit_time);

x_sol = h.eval(x);
y_sol = h.eval(y);
eta_sol = h.eval(eta);
obj = h.objective;

end


function [ obj_sub, delta_sol ] = SubProblem_CCG_exhaust(x,y, N,N_pr,Q,bar_n,K, Gf,Gc,d,AF,ACC,AV,AC,AX,AY,AZ,B0,Z,Z0, L,REALIZATION,DELTA, limit_time)
%SUBPROBLEM Summary of this function goes here

Mbig = 10000;       % set big M
G = [AF,ACC,AV,AC];          % combine matrix  Gz >= 
b = [Gf;Gc;d;zeros(N,1)];    % simplify b'z 

%%%%%%%%%
%Realization = GeneRealization(REALIZATION);
Realization = REALIZATION(:,end)';
%%%%%%%%%

i=1; obj_sub = -Inf; delta_sol = [];
while i <= size(Realization,1)
    delta = Realization(i,:);
    rome_begin;
    h = rome_model('Identify Problem');
    
    z = newvar(size(b,1),1,'nonneg');   % combine all variable
    
    rome_minimize(b'*z);    
    rome_constraint( G*z >= AX*x+AY*y+AZ*(Z0+Z*delta')+B0 );
     
    h.solve('CPLEX',limit_time);
    
    sol = h.objective;
    if sol>=obj_sub
        delta_sol = delta';
        obj_sub = sol;        
    end
    i=i+1;    
end


end
