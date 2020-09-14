
clear;
clc;

%******************************
instance_type = 1;
%******************************


for i = 217
    
    read_filename = ['DATA/DC' num2str(instance_type) '/MAT/mv' num2str(i) '.mat'];
    save_filename = ['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance' num2str(i) '.mat'];
    
    Sample = ['DATA/DC' num2str(instance_type) '/Sample_Training.mat'];
    
    ReadDataFromRCP_Sample( read_filename , save_filename, Sample );    % Process original Data to problem used

end


