function [] = plot_sR()
    
    % load variables from corresponding data files
    load('../Output_Data/IDO_Full_12','sR_hf','sR_mf','sR_rnd','Nrepeat','mff','Cost_mf','history_mf');
    load('../Output_Data/IDO_EI_12.mat','objMin','results');


    Budgetmf = [];
    for i = 1:Nrepeat
        Budgetmf = [Budgetmf, Cost_mf{i}'];   
    end    
    Budgetmf = unique(Budgetmf);

    Costmf = [];
    for i = 1:Nrepeat
        Costmf = [Costmf, interp1(Cost_mf{i}, sR_mf{i},Budgetmf)'];
    end   
    
    Costmf = Costmf';
    Costmf_i = cell(10,1);
    
    for i = 1:10
        Costmf_i{i} = interp1([Cost_mf{i}],sR_mf{i},[1:12]);
    end
    
    Costmf_i = cell2mat(Costmf_i);
    
    % Preprocess data before plotting
    Costhf = reshape(cell2mat(sR_hf),size(sR_hf{1},1),Nrepeat)';
    Costrnd = reshape(cell2mat(sR_rnd),size(sR_rnd{1},1),Nrepeat)';
    CostEI = objMin';
    
    Costmfmin = min(min(Costmf));
    Costhfmin = min(min(Costhf));
    Costrndmin = min(min(Costrnd));
    CostEImin = min(min(CostEI));
    
    fstar = min([Costmfmin Costhfmin Costrndmin CostEImin]);
    fmin = -(fstar-mean(Costmf)) - std(Costmf)/sqrt(Nrepeat);
    fmax = -(fstar-mean(Costmf)) + std(Costmf)/sqrt(Nrepeat);
    
    % Plot
    figure; hold on;
    colors = {[0.3010 0.7450 0.9330], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560]};
    p1 = plot(Budgetmf, -(fstar-mean(Costmf)), 'color',colors{1},'LineWidth', 5);
    p2 = errorbar(1:12, -(fstar-mean(Costmf_i)), -(fstar-mean(Costmf_i)) - std(Costmf_i)/sqrt(Nrepeat), -(fstar-mean(Costmf_i)) + std(Costmf_i)/sqrt(Nrepeat), 'color', colors{1},'CapSize', 8, 'LineWidth', 1.25,'LineStyle','none');
    p3 = plot(1:12, -(fstar-mean(Costhf)),'color', colors{2},'LineWidth', 5);
    p4 = errorbar([1:12]+0.025, -(fstar-mean(Costhf)), -(fstar-mean(Costhf)) - std(Costhf)/sqrt(Nrepeat), -(fstar-mean(Costhf)) + std(Costhf)/sqrt(Nrepeat), 'color', colors{2},'CapSize', 8, 'LineWidth', 1.25,'LineStyle','none');
    p5 = plot(1:12, -(fstar-mean(Costrnd)),'color', colors{3},'LineWidth', 5);
    p6 = errorbar([1:12]-0.025, -(fstar-mean(Costrnd)), -(fstar-mean(Costrnd)) - std(Costrnd)/sqrt(Nrepeat), -(fstar-mean(Costrnd)) + std(Costrnd)/sqrt(Nrepeat), 'color', colors{3},'CapSize', 8, 'LineWidth', 1.25,'LineStyle','none');
    p7 = plot(1:12, -(fstar-mean(CostEI)), 'color', colors{4},'LineWidth', 5);
    p8 = errorbar([1:12]+0.05, -(fstar-mean(CostEI)), -(fstar-mean(CostEI)) - std(CostEI)/sqrt(Nrepeat), -(fstar-mean(CostEI)) + std(CostEI)/sqrt(Nrepeat), 'color', colors{4},'CapSize', 8, 'LineWidth', 1.25,'LineStyle','none');
    ylabel('simple regret');
    xlabel('cost');
    xlim([0 13])
    ylim([0.03 40]);
%     titleStr = sprintf('Simple Regret: Costs = %s', mat2str(mff.costs));
%     title(titleStr);
    set(gcf,'color','w');
    set(gca,'yscale','log')
    box on;
    legend([p1 p3 p5 p7],{'MF-GP-UCB','GP-UCB','RAND','EI'});
    set(gca,'FontSize',20)
    
end