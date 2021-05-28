function [] = plot_sR3()
    
    load('../Output_Data/IDC_CaseStudy3_1', 'Cost_mf','history_mf','maxObj_mf')
    C1 = Cost_mf;
    H1 = history_mf;
    obj1 = maxObj_mf;

    load('../Output_Data/IDC_CaseStudy3_2', 'Cost_mf','history_mf','maxObj_mf')
    C2 = Cost_mf;
    H2 = history_mf;
    obj2 = maxObj_mf;

    Nrepeat = size(C1,1);

    %capitol spent at each itteration
    Budget1 = [];
    Budget2 = [];
    for i = 1:Nrepeat
        Budget1 = [Budget1, C1{i}(find(H1{i}.evalFidels==3))'];   
        Budget2 = [Budget2, C2{i}(find(H2{i}.evalFidels==3))'];   
    end    
    Budget1 = unique(Budget1);
    Budget2 = unique(Budget2);
    
    %objective value
    Cost1 = [];
    Cost2 = [];
    for i = 1:Nrepeat
        Cost1 = [Cost1, interp1(C1{i}(find(H1{i}.evalFidels==3)),...
            H1{i}.evalVals(find(H1{i}.evalFidels==3)),Budget1)'];
        Cost2 = [Cost2, interp1(C2{i}(find(H2{i}.evalFidels==3)),...
            H2{i}.evalVals(find(H2{i}.evalFidels==3)),Budget2)'];
    end  
    
    %take best observed
    for i = 2:size(Cost1,1)
       Cost1(i,:) = max(Cost1(1:i,:)); 
    end
    for i = 2:size(Cost2,1)
       Cost2(i,:) = max(Cost2(1:i,:)); 
    end
   
    Cost1_i = cell(10,1);
    Cost2_i = cell(10,1);
    
    for i = 1:Nrepeat
        Cost1_i{i} = interp1(Budget1,Cost1(:,i),[1:30]);
        Cost2_i{i} = interp1(Budget2,Cost2(:,i),[1:30]);
    end
    
    Cost1_i = cell2mat(Cost1_i);
    Cost2_i = cell2mat(Cost2_i);
    
    Cost1max = max(max(Cost1));
    Cost2max = max(max(Cost2));

    fstar = max([Cost1max Cost2max]);
    f1min = - nanstd(Cost1_i);
    f1max = + nanstd(Cost1_i);
    f2min = - nanstd(Cost2_i);
    f2max = + nanstd(Cost2_i);
    
    % Plot
    figure; hold on;
    col = {[0.3010 0.7450 0.9330], [0.8500 0.3250 0.0980],...
                [0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560]};
    p1 = plot([1:30]/30, (fstar-nanmean(Cost1_i)), 'color',col{1},'LineWidth', 5);
    p2 = errorbar([1:30]/30-.0005, (fstar-nanmean(Cost1_i)), f1min,f1max, 'color', ...
            col{1},'CapSize', 8, 'LineWidth', 1.25,'LineStyle','none');
    p3 = plot([1:30]/30, (fstar-nanmean(Cost2_i)),'color', col{2},'LineWidth', 5);
    p4 = errorbar([1:30]/30+.0005, (fstar-nanmean(Cost2_i)), f2min,f2max, ...
            'color', col{2},'CapSize', 8, 'LineWidth', 1.25,'LineStyle','none');
    for i = 1:Nrepeat
        h1 = H1{i}.evalVals(find(H1{i}.evalFidels==3));
        h2 = H2{i}.evalVals(find(H2{i}.evalFidels==3));
        for j = 1:length(h1)
            h1(j) = max(h1(1:j));
        end
        for j = 1:length(h2)
            h2(j) = max(h2(1:j));
        end
%        p5= plot(C1{i}(find(H1{i}.evalFidels==3))'/30,fstar-h1','color',col{1});
%        p6 =plot(C2{i}(find(H2{i}.evalFidels==3))'/30,fstar-h2','color',col{2});
%        p5.Color(4) = .5;
%        p6.Color(4) = .5;
    end
    ylabel('simple regret');
    xlabel('fraction of spent budget');
    xlim([0 1.05])
    ylim([0.04 1]);
    set(gcf,'color','w');
    set(gca,'yscale','log')
    box on;
    legend([p1 p3 ],{'Case 1','Case2'});
    set(gca,'FontSize',20)
    

    
end