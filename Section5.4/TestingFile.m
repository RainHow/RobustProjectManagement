clc;
clear;

%*********************************
instance_type = 3;
%*********************************
ratio = 1.2;

for i = 1:240
    outputfile = ['DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro' num2str(i) '.mat'];
    
    load(outputfile);
    
    save(outputfile, 'obj_aro', 'x_aro', 'y_aro', 'r_aro', 'r0_aro', 'vartheta_aro', 'rc_aro', 'T_aro',   'target');

end
