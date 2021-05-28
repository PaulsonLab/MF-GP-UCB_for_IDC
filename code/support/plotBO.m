% function [] = plotBO(['legend name 1',...,'legend name n'], fx1,...fxn=[])
% Will plot the cost vs itteration for Baysian Optimization simulations.
% fx data can be an array or table, and any number of datasets can be compared.
% Each row of fx should have the objective evaluations of one BO run. The
% legend should should have as many entries 
%
% If there is only one run, the plot will have blue dots for each objective
% evaluation, and a red line showing the minimum value found at the given
% itteration
%
% if there are multiple runs, the mean will be plotted, and the range of
% the highest and lowset evaluation is shaded around the mean. Legend
% assumes bayesOpt is fx, and fx2 is Turbo
%
function [] = plotBO(legend_list, overlap, varargin)
if nargin<3
    comp =0;
elseif nargin >= 3
    comp = 1;
end

sz= size(varargin{1});   L = sz(1);  W = sz(2);
x = [1:W];

if L == 1
    %plot one simulation
    fx_best = zeros(W,1);
    fx = varargin{1};
    for i = 1:W
        fx_best(i) = min(fx(1:i));
    end
    fx(isnan(fx)) = 0;
    
    figure(); hold on;
    plot(x,fx, '-b','LineWidth',2)
    plot(x,fx_best, 'or','LineWidth',2)
    ylabel('Cost ')
    xlabel('Iteration')
    ylim([min(fx)-10, -min(fx)*2 ])
    legend(legend_list,'Best Point')

elseif (L>=1 && comp == 0)
    
    fx = varargin{1};
    if isa(fx,'table')
        fx = table2array(fx);
    end
    fx(isnan(fx)) = 0;
     
    fbar = mean(fx);
    fmax = fbar + std(fx);
    fmin = fbar - std(fx);
        

	figure(); hold on;
    plot1 = plot(x,fbar, 'b','LineWidth',5);
    for i = 1:L
        plot2 = plot(x,fx(i,:), 'k-');
        plot2.Color(4) = 0.1;
    end
    patch1 = patch([x fliplr(x)], [fmin fliplr(fmax)], 'k');
    alpha(.2)
    ylabel('Cost')
    xlabel('Itteration')
    xlim([1 W])
    
elseif (L>=1 && comp ~= 0)
    color = {'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'};
    aspA  = 1;
    aspB  = 3;
    
    figure(); hold on;
    
    for i2 = 1:nargin-1

        fx = varargin{i2};

        fx(isnan(fx)) = 0;

        fbar = mean(fx);
        fmax = fbar + std(fx);
        fmin = fbar - std(fx);
        
        if ~overlap

            subplot(aspA,aspB,i2); hold on
            % for legend
            plot(NaN,NaN, color{i2},'LineWidth',5);
            plot02 = plot(NaN,NaN, color{i2},'LineWidth',5);
            plot02.Color(4) = .1;
            % make std patch
            patch([x fliplr(x)], [fmin fliplr(fmax)], color{i2});
            alpha(.2)
            % plot actual convergance of all runs 
            for i = 1:L
                plot2 = plot(x,fx(i,:), color{i2});
                plot2.Color(4) = 0.1;
            end
            %plot mean
            plot(x,fbar, color{i2},'LineWidth',5);


            ylabel(legend_list{i2})
            xlabel('Itteration')
            xlim([1 W])
            ylim([-8 0])
            legend('Mean', 'STD')
            box on
            grid on
        else
            % for legend
            for i1 = 1:length(varargin)
                plot(NaN,NaN, color{i1},'LineWidth',5);
            end
            plot02 = plot(NaN,NaN, color{i2},'LineWidth',5);
            plot02.Color(4) = .1;
            % make std patch
            patch([x fliplr(x)], [fmin fliplr(fmax)], color{i2});
            alpha(.2)
            % plot actual convergance of all runs 
            for i = 1:L
                plot2 = plot(x,fx(i,:), color{i2});
                plot2.Color(4) = 0.2;
            end
            %plot mean
            plot(x,fbar, color{i2},'LineWidth',5);


            ylabel('Cost')
            xlabel('Iteration')
            xlim([1 W])
            ylim([-10 0])
            legend(legend_list)
            box on
            grid on
            set(gcf,'color','w');
            set(gca,'FontSize',14)

        end

    end

    
        
else
   'Nothing to plot; verify the dimensions of fx are correct '
    
end



