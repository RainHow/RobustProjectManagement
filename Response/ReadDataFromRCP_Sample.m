function ReadDataFromRCP_Sample( read_filename , save_filename, Sample )
% This is the REVISED function: Q=2

OriginalData = importdata(read_filename);
RevisedData = OriginalData;

% ************************delete the resoueces-related data************************
RevisedData(:,2:(OriginalData(1,2)+1)) = [];
RevisedData(1:2,:) = [];      

% ************************basic number************************
N = OriginalData(1);             % number of tasks
N_pr = sum(RevisedData(:,2));    % number of relations
M = 5;                           % price levels of general resource
Q = 2;                           % assume 2 piecewise linear function 
U = ones(N,1);                   % assume linear cost, no piecewise
bar_n = sum(U);  
b = RevisedData(:,1);

% ********************** resouce-related number ***************
% Crashing requires one general capacity/per unit
% ith item: resources=floor(i*N*5/M),  price=floor(i*N*2/M)*(1-0.05*(i-1))
Ux = zeros(M,1);  
a = zeros(M,1);
for i = 1 : M
    %a(i) = floor( i*N*5/M );
    a(i) = i*sum(b);
    Ux(i) = 0.01 * a(i) * (1-0.05*(i-1));
end
Uy = 10*ones(N,1);    % price of each unit of specific resource is 10

% no capacity constraint
A_xy = [Ux' Uy'];  b_xy = sum(Ux)+sum(Uy);  


% ************************crashing relevant parameters************************
Gc = zeros(bar_n,1); % the cost of using crashing resources is 0, such that we should use it after reservation
Ec = [ 0; 100*ones(bar_n-2,1); 0];
    % first and last task cannot be crashed since they are auxiliary task
    % for other tasks, set a high upper limit
Eta = ones(N,1);    % special capacity is needed for all tasks
Theta = ones(N,1);  % normalized to one

Precedence = zeros(N_pr,2);
tem_row_number = 0;
for i = 1 : N-1
    Precedence(tem_row_number+1:tem_row_number+RevisedData(i,2),1) = i*ones(RevisedData(i,2),1);
    Precedence(tem_row_number+1:tem_row_number+RevisedData(i,2),2) = (RevisedData(i,3:2+RevisedData(i,2)))';
    tem_row_number = tem_row_number+RevisedData(i,2);
end

%************************** fast tracking parameters*******************
Gf = zeros(N_pr,1); % the cost of using fask tracking is 0
Gamma = zeros(N_pr,1);   % no fast tracking

d = [1;3]; % cost for makespan is higher than the sum of making reservation and use capacity to crash
l = [0; 1000*sum(b)];

rome_begin;
h = rome_model('makespan');
T_completion = newvar(N,1,'nonneg');
rome_minimize(T_completion(N));
for i = 1:N_pr
    source_node = Precedence(i,1);
    destination_node = Precedence(i,2);
    rome_constraint(T_completion(destination_node)>=T_completion(source_node)+b(destination_node));
end
h.solve('CPLEX');

l(1) = 4*(h.objective);

% Information generate
% beta-binomial distribution, n=4, alpha=1, beta=3, mu=1,var=1.2

load(Sample);
K = size(DeltaSample,1);  % 3 factor based model
Z = b*ones(1,K);
Z0 = b;
Mu = mean(DeltaSample,2);
Sigma = mean(DeltaSample.^2,2);
Cor = mean(abs(sum((DeltaSample - Mu*ones(1,SampleSize))./(sqrt(Sigma)*ones(1,SampleSize)),1)));


REALIZATION = [];
for i = 1:K
    REALIZATION = [REALIZATION; unique(DeltaSample(i,:))];
end
L = size(REALIZATION,2)*ones(K,1);

% utility function paramaters
Am = [1,0]; bm = [0,-1];

%target = 2*(K+1) * h.objective;

[ AF, ACC, AV, AC, AX, AY, AZ, B0 ] = ExlicitingParameterForConstraints( N, N_pr, K, U, Q, M, a, Ec, Eta, Theta, Precedence, Gamma, l, b );

save( save_filename, 'N', 'N_pr', 'K', 'Q', 'U', 'bar_n', 'M', 'Ux', 'Uy', 'a', 'A_xy', 'b_xy', 'Gc', 'Ec', 'Eta', 'Theta', 'Precedence', 'Gf', 'Gamma', 'd', 'l', 'b', 'Z', 'Z0', 'Mu', 'Sigma', 'Cor', 'L', 'REALIZATION', 'AF', 'ACC', 'AV', 'AC', 'AX', 'AY', 'AZ', 'B0', 'Am', 'bm');

end


%*****************  Sub function **********************

function [ AF, ACC, AV, AC, AX, AY, AZ, B0 ] = ExlicitingParameterForConstraints( N, N_pr, K, U, Q, M, a, Ec, Eta, Theta, Precedence, Gamma, l, b )
%  Generate the Parameters for S(z) constraints

bar_n = sum(U);

number_general_constraints = N_pr + N + 1 + nnz(Eta) + N_pr + bar_n + Q ;

AF = zeros(number_general_constraints,N_pr);
AF(1:N_pr,:) = eye(N_pr);
AF(N_pr+N+1+nnz(Eta)+1:N_pr+N+1+nnz(Eta)+N_pr,:) = -eye(N_pr);

ACC = zeros(number_general_constraints,bar_n);
for k = 1 : N_pr
    j = Precedence(k,2);
    ACC(k, sum(U(1:j-1))+1:sum(U(1:j)) ) = ones(1,U(j));
end
for j = 1 : N
    ACC(N_pr+N+1,sum(U(1:j-1))+1:sum(U(1:j))) = -Theta(j)*ones(1,U(j));
end
Eta_Nonzeros = find(Eta);
for k = 1 : length(Eta_Nonzeros)
    j = Eta_Nonzeros(k);
    ACC(N_pr+N+1+k,sum(U(1:j-1))+1:sum(U(1:j))) = -ones(1,U(j));
end
ACC(N_pr+N+1+nnz(Eta)+N_pr+1:N_pr+N+1+nnz(Eta)+N_pr+bar_n,:) = -eye(bar_n);


AV = zeros(number_general_constraints,Q);
AV(N_pr+1:N_pr+N,:) = ones(N,Q);
AV(N_pr+N+1+nnz(Eta)+N_pr+bar_n+1:N_pr+N+1+nnz(Eta)+N_pr+bar_n+Q,:) = -eye(Q);


AC = zeros(number_general_constraints,N);
for k = 1 : N_pr
    AC(k,Precedence(k,1)) = -1;
    AC(k,Precedence(k,2)) = 1;
end
AC(N_pr+1:N_pr+N,:) = -eye(N);


AX = zeros(number_general_constraints,M);
AX(N_pr+N+1,:) = -a';

AY = zeros(number_general_constraints,N);
for k = 1 : length(Eta_Nonzeros)
    j = Eta_Nonzeros(k);
    AY(N_pr+N+1+k,j) = -Ec(sum(U(1:j)));
end

AZ = zeros(number_general_constraints,K);
for k = 1 : N_pr
    j = Precedence(k,2);
    AZ(k,j) = 1;
end

B0 = zeros(number_general_constraints,1);
for k = 1 : N_pr
    j = Precedence(k,2);
    B0(k) = b(j);
    B0(N_pr+N+1+nnz(Eta)+k) = -Gamma(k);
end
B0(N_pr+N+1+nnz(Eta)+N_pr+1:N_pr+N+1+nnz(Eta)+N_pr+bar_n) = -Ec;
B0(N_pr+N+1+nnz(Eta)+N_pr+bar_n+1:N_pr+N+1+nnz(Eta)+N_pr+bar_n+Q) = -l;

end
