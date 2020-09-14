clear;
clc

instance_type = 1;
number_instance = 1470;

% instance_type = 2;
% number_instance = 240;


Reservation_x = cell(number_instance,1);
Reservation_y = cell(number_instance,1);

for i = 1 : number_instance

    eval(['load(''DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro' num2str(i) '.mat'');']);
    Reservation_x{i} = x_aro;
    Reservation_y{i} = y_aro;
    
end

eval(['save(''data/DC' num2str(instance_type) '/CapacityAnalysis.mat'', ''Reservation_*'');']);