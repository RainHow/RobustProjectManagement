clc;
clear;

load([ 'DATA/DC1/InstanceAnalysis.mat' ]);
Gamma_bench = 0.4

for i = 1:5
    Gamma_temp = Gamma(360*(i-1)+1:360*i);
    eval([ 'Index' num2str(i) '_small = 360*(i-1) + find(Gamma_temp<Gamma_bench);' ])
    eval([ 'Index' num2str(i) '_large = 360*(i-1) + find(Gamma_temp>=Gamma_bench);' ])   
end

