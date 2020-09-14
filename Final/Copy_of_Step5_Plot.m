clear;
clc;

%*****************************
instance_type = 2;
%*****************************

for plot_type_index = 1:4
     
    FontSize = 16;
    BigFontSize = 20;
    
    load(['DATA/DC' num2str(instance_type) '/CostAnalysis_withTarget.mat']);
    
    plot_type_collection = cell(1,4);
    plot_type_collection{1} = 'CCG';
    plot_type_collection{2} = 'BD';
    plot_type_collection{3} = 'EO';
    plot_type_collection{4} = 'UM';

    plot_type = plot_type_collection{plot_type_index};
   
    subplot(2,4,1)
    eval(['histogram(P_' plot_type '_1best)'])
    title('best');
    set(gca,'FontSize',FontSize)
    
    subplot(2,4,2)
    eval(['histogram(P_' plot_type '_2worst)'])
    title('worst')
    set(gca,'FontSize',FontSize)
    
    subplot(2,4,3)
    eval(['histogram(P_' plot_type '_3mean)'])
    title('mean')
    set(gca,'FontSize',FontSize)
    
    subplot(2,4,4)
    eval(['histogram(P_' plot_type '_4variance)'])
    title('variance')
    set(gca,'FontSize',FontSize)
   
    subplot(2,4,5)
    eval(['histogram(P_' plot_type '_7var5)'])
    title('var@5%')
    set(gca,'FontSize',FontSize)
    
    subplot(2,4,6)
    eval(['histogram(P_' plot_type '_8cvar5)'])
    title('cvar@5%')
    set(gca,'FontSize',FontSize)
    
    subplot(2,4,7)
    eval(['histogram(P_' plot_type '_9vp)'])
    title('EP')
    set(gca,'FontSize',FontSize)
    
    subplot(2,4,8)
    eval(['histogram(P_' plot_type '_10el)'])
    title('EL')
    set(gca,'FontSize',FontSize)
   
    plot_title_collection = cell(1,4);
    plot_title_collection{1} = 'Robust-Z';
    plot_title_collection{2} = 'Robust-B';
    plot_title_collection{3} = 'Risk Neutral';
    plot_title_collection{4} = 'Deterministic';

    
    suptitle_info = suptitle(plot_title_collection{plot_type_index});
    
    set(suptitle_info,'FontSize',BigFontSize);
    eval([ 'saveas(gcf,''DATA/DC' num2str(instance_type) '/Results_NoTarget/' plot_type_collection{plot_type_index} '.fig'')' ])
    
end