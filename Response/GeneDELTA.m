
function [ Temp ] = GeneDELTA( REALIZATION )

K = size(REALIZATION,1);

switch K
    case 3
        [X1,X2,X3] = ndgrid( REALIZATION(1,:),REALIZATION(2,:),REALIZATION(3,:));
        Temp = [X1(:) X2(:) X3(:)];        
    case 5
        [X1,X2,X3,X4,X5] = ndgrid( REALIZATION(1,:),REALIZATION(2,:),REALIZATION(3,:), REALIZATION(4,:),REALIZATION(5,:));
        Temp = [X1(:) X2(:) X3(:) X4(:) X5(:)];
    case 7
        [X1,X2,X3,X4,X5,X6,X7] = ndgrid( REALIZATION(1,:),REALIZATION(2,:),REALIZATION(3,:), REALIZATION(4,:),REALIZATION(5,:), REALIZATION(6,:), REALIZATION(7,:) );
        Temp = [X1(:) X2(:) X3(:) X4(:) X5(:) X6(:) X7(:)];        
end

end

