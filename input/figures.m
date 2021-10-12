%% Foreclosure Dynamics in Complementary Markets
%  Functions for figures
%  Matteo Courthoud
%  25/06/2020

% Run figures.main to produce all figures



classdef figures
    
    % Properties
    properties (Constant)
        
        %policies = ["baseline", "nopredentrypricing", "nopredexitpricing", "nopredentrybundling", "nopredexitbundling"]
        %policies = ["baseline", "nolearning", "nobundling", "datasharing", "limitedbundling"];
        policies = ["baseline", "nolearning", "nomergers", "nobundling", ...
            "datasharing", "limitedmergers",  "limitedbundling", ...
            "nopredentrypricing", "nopredexitpricing", "nopredentrybundling", "nopredexitbundling"];
        statlimits = [5,1,1,1,1,3,3,5,1,1,1,3,3,5]; % Stats min and max
        difflimits = [5,1,1,1,1,1,1,.3,1,1,1,1,1,.3]; % Stats min and max for differences
        varnames = ["Price - Cost (short run)", "Below Cost Pricing (short run)", ...
            "Entry Probability (short run)", "Exit Probability (short run)", "Merger Probability (short run)", ...
            "Total Profits (NPV)", "Consumer Surplus (NPV)", "Total Welfare (NPV)", ...
            "Monopoly Probability (long run)", "Monopoly Probability A (long run)", "Monopoly Probability B (long run)", ...
            "Total Profits (long run)", "Consumer Surplus (long run)", "Total Welfare (long run)"];
        matketnames = ["Monopoly in A&B","B Monopoly in A&B",...
                    "mixed M/Dpoly","B mixed M/Dpoly",...
                    "Duopoly in A&B","B/2 Duopoly in A&B","B Duopoly in A&B"];
        filenames = ["a70g0s30", "a30g0s70"];
        beta = 0.95;
        c = 1;
        gamma = 0;
        p0 = 1.5;
    end
    
    % Methods
    methods (Static)
        
        % Plot one timeline
        function plot_timeline(data, S, T, j, colors, file, policy)
            
            % Init
            colnames = data.Properties.VariableNames(3:end);
                
            % Make figure
            figure();
            set(gcf,'position',[100,1000,900,500]);
            for i=1:size(S,1)
                x = data(data.s==S{i,:} & data.t>0, colnames{j});
                plot(1:max(T), x{:,:}, 'Color', colors(i, :), 'LineWidth',6)
                hold on;
            end   

            % Make pretty
            set(gca, 'FontName', 'SansSerif', 'FontSize', 12);
            xlabel("Time Periods", 'FontSize', 18);
            legend(figures.matketnames, 'Location','best', 'FontSize', 18, 'AutoUpdate','off');
            set(gca, 'Xgrid', 'on', 'Ygrid', 'on', 'Box', 'off', 'TickDir', 'out', 'Gridalpha', 0.05)
            set(gca, 'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3])
            set(gca, 'Linewidth', 2)
            set(gca,'TickLength',[0 0]);
            set(gca,'XColor',[.4 .4 .4],'YColor',[.4 .4 .4]);
            title(colnames{j}, 'FontSize', 28);  

            % Add lines
            v = data(data.t==0, colnames{j});
            yline(v{1,:}, 'r-.', 'LineWidth', 4);
            yline(v{2,:}, 'b-.', 'LineWidth', 4);

            % Print
            rely = (v{1,:} - min(ylim())) / (max(ylim()) - min(ylim()));
            annotation('textbox',[0.02 rely*0.81+0.14 .1 0],'String',['Long-Run', newline, 'B Mpoly'],'FitBoxToText','on','Color','r','EdgeColor','none','FontSize', 14)
            rely = (v{2,:} - min(ylim())) / (max(ylim()) - min(ylim()));
            annotation('textbox',[0.02 rely*0.81+0.14 .1 0],'String',['Long-Run', newline, 'B Dpoly'],'FitBoxToText','on','Color','b','EdgeColor','none','FontSize', 14)

            % Save
            saveas(gca, sprintf('../output/timelines/%s_%s_%s.png', policy, file, colnames{j}));
            close;
        end

        % Plot all the timelines
        function plot_all_timelines()
            
            % Loop over policies
            for file=figures.filenames
                for policy=figures.policies

                    % Import data
                    filename = sprintf('../output/timeseries/%s_%s', policy, file);
                    data = readtable(sprintf("%s.csv", filename));

                    % Get colors
                    T = unique(data.t);
                    S = unique(data(data.t>0, "s"));
                    palette = palettes.("viridis");
                    colors = interp1(linspace(0,1,size(palette,1)),palette,linspace(0,1,size(S,1)));

                    % Make figures for each stat
                    for j=1:size(data, 2)-2
                        figures.plot_timeline(data, S, T, j, colors, file, policy)
                    end
                end
            end
                        
        end
        
        
        
        

        % Plot flows
        function plot_flows(Q, pr_s1, ybars_bottom, colors, I, J)            

            % Loop over rows and columns
            for j=1:J-1
                bottom_rights = ybars_bottom(:,j+1)';
                for i=1:I

                    % Get corners
                    top_lefts = (cumsum(Q(i,:,j)) - Q(i,:,j)) * pr_s1(i,j)*0.9 + ybars_bottom(i,j);
                    bottom_lefts = cumsum(Q(i,:,j)) * pr_s1(i,j)*0.9 + ybars_bottom(i,j);
                    top_rights = bottom_rights;
                    bottom_rights = bottom_rights + bottom_lefts - top_lefts;

                    % Get coordinates
                    [X, Y] = figures.get_coordinates(j, J, top_lefts, bottom_lefts, top_rights, bottom_rights);

                    % Plot flow
                    patch('XData', X, 'YData', Y, 'FaceAlpha', .3, 'Facecolor', colors(i,:), 'EdgeColor', 'none');
                end
            end

        end
        
        % Makes curve between two points
        function [x, y] = get_curves(x1, y1, x2, y2)
            t = linspace(0, pi, 15);
            c = (1-cos(t))./2; 
            Ncurves = numel(y1);
            y = repmat(y1, 15, 1) + repmat(y2 - y1, 15,1) .* repmat(c', 1, Ncurves);
            x = repmat(linspace(x1, x2, 15)', 1, Ncurves);
        end 

        % Get shape coordinates
        function [X, Y] = get_coordinates(j, J, top_lefts, bottom_lefts, top_rights, bottom_rights)

            % Get curves
            w = J/40;
            [bottom_x, bottom_y] = figures.get_curves(j+w, bottom_lefts, j+1-w, bottom_rights);
            [top_x, top_y] = figures.get_curves(j+1-w, top_rights, j+w, top_lefts);

            % Get all oordinates
            X = [bottom_x; top_x];
            Y = [bottom_y; top_y];
        end
       
        % Compute vertical coordinates of bars
        function [ybars_bottom, ybars_top] = get_ybars(pr_s1, I, J)
            ybars_bottom = zeros(size(pr_s1));
            ybars_top = zeros(size(pr_s1));
            for j=1:J
                ybars_bottom(:,j) = cumsum(pr_s1(:,j)*0.9 + 0.1/(I-1)) - 0.1/(I-1) - pr_s1(:,j)*0.9;
                ybars_top(:,j) = cumsum(pr_s1(:,j)*0.9 + 0.1/(I-1)) - 0.1/(I-1);
            end
        end

        % Plot vertical bars
        function plot_bars(ybars_bottom, ybars_top, colors, I, J)

            % Bar width (half)
            w = size(ybars_bottom,2)/40;            

            % Plot bars
            for j=1:J
                set(gca,'ColorOrderIndex',1)
                for i=1:I
                    y_corners = [ybars_bottom(i,j), ybars_bottom(i,j), ybars_top(i,j), ybars_top(i,j)];
                    x_corners = [j-w, j+w, j+w, j-w];
                    patch('X', x_corners, 'Y', y_corners, 'FaceAlpha', .8, 'Facecolor', colors(i,:), 'EdgeColor', 'none');
                end
                hold on
            end
        end

        % Make graph pretty
        function prettify(title, pr_s1, ylabels, xlabels, ymeans, I, J)

            % Y labels
            for i=1:I
                if pr_s1(i,1)>0.05
                    text(1-0.2, ymeans(i,1)-0.02, ylabels(i), 'HorizontalAlignment', 'right','Fontweight', 'Bold', 'Color', [.4 .4 .4])
                    text(1-0.2, ymeans(i,1)+0.02, sprintf("%.2f",pr_s1(i,1)), 'HorizontalAlignment', 'right', 'Color', [.4 .4 .4])
                end
                if pr_s1(i,end)>0.05
                    text(J+0.2, ymeans(i,end)-0.02, ylabels(i),'Fontweight', 'Bold', 'Color', [.4 .4 .4])
                    text(J+0.2, ymeans(i,end)+0.02, sprintf("%.2f",pr_s1(i,end)), 'Color', [.4 .4 .4])
                end
            end

            % X labels
            for j=1:J
                text(j, 1.03, xlabels(j), 'HorizontalAlignment', 'center', 'Color', [.4 .4 .4])
            end
            text((J+1)/2, 1.07, "Periods", 'HorizontalAlignment', 'center', 'Fontsize', 12,'Fontweight', 'Bold', 'Color', [.4 .4 .4])

            % Title
            text((J+1)/2, -0.1, title, 'HorizontalAlignment', 'center', 'Fontsize', 20)

            % Font
            set(gca, 'FontName', 'SansSerif')

        end
        
        % Make alluvial plot
        function plot_alluvial(data, policy, file)
            
            % Init
            I = length(unique(data.s));
            T = unique(data.t);
            J = length(T);
            pr_s1 = ones(I,J+1)/I;
            Q = zeros(I, I, J);
            for t=1:length(T)
               Q(:,:,t) = table2array(data(data.t==T(t), 4:end));
            end
            for j=1:J
                pr_s1(:,j+1) = pr_s1(:,j)' * Q(:,:,j);
            end
            
            % Setup
            xlabels = ["0", "1", "3", "5", "10", "30", "100"];
            title = "Transition Flows";
            
            % Get vertical bars
            I = size(pr_s1,1);
            J = size(pr_s1,2);
            [ybars_bottom, ybars_top] = figures.get_ybars(pr_s1, I, J);
            
            % Get colors
            palette = palettes.("viridis");
            k = 10;
            colors = interp1(linspace(0,1,size(palette,1)),palette,linspace(0,1,3+(size(pr_s1,1)-1)*k));
            colors = colors(2:k:end-1,:);

            % Plot
            figure();
            set(gcf,'position',[100,1000,900,500]);
            ylim([-0.05,1]);
            axis off
            axis ij
            hold on

            % Plot flows
            figures.plot_flows(Q, pr_s1, ybars_bottom, colors, I, J)

            % Plot bars
            figures.plot_bars(ybars_bottom, ybars_top, colors, I, J)

            % Prettify
            ymeans = (ybars_bottom + ybars_top) / 2;
            figures.prettify(title, pr_s1, figures.matketnames, xlabels, ymeans, I, J)
            
            % Save and close
            saveas(gca, sprintf('../output/alluvial/game_%s_%s.png', file, policy));
            close;

        end
        
        % Alluvial graph
        function make_all_alluvialplots()
            
            % Loop over policies
            for file=figures.filenames
                for policy=figures.policies

                    % Import data
                    filename = sprintf('../output/transitions/%s_%s', policy, file);
                    data = readtable(sprintf("%s.csv", filename));

                    %set(gca, 'OuterPosition', [-0.15,0,1.27,1.1])
                    figures.plot_alluvial(data, policy, file);
                end
            end
        end
               
        
        
        
        
        
        
        
        
        % Get compstats
        function compstats = get_compstats(data, label)
            sigmas = sort(unique(data.sigma));
            alphas = sort(unique(data.alpha), 'descend');
            compstats = zeros(length(sigmas), length(alphas));
            for i=1:length(sigmas)
                for j=1:length(alphas)
                    row = (data.sigma==sigmas(i)) .* (data.alpha==alphas(j)) .* ...
                          (data.beta==figures.beta) .* (data.c==figures.c) .* ...
                          (data.gamma==figures.gamma) .* (data.p0==figures.p0) .* ...
                          (data.market==40);
                    compstats(i,j) = table2array(data(logical(row), label));
                end
            end
        end
        
        % Compare policies
        function compstats = compare_games(games, game0, labels, folder)
            
            % Get results
            compstats = zeros(length(games), figures.S);
            for k=1:length(games)
                [~, s0] = ismember({'S2200'}, games(k).labels);
                compstats(k,:) = games(k).sumstats(s0, :);
            end
            
            % Return if graphs are off
            if game0.graphs == false
                return
            end
                           
            % Compute difference
            rel_results = compstats(2:end,:) - compstats(1,:);
            
            % Plot
            figure('visible', game0.showfigures);
            set(gcf,'position',[100,1000,900,400]);
            axes('Position',[0.1 0.15 0.7 0.75]);
            bar(rel_results', 'EdgeColor','none' );

            % Prettify
            set(gca, 'box','off');
            alpha(0.7);
            ylim([-1.5 +1.5]);
            set(gca, 'FontName', 'Helvetica', 'FontSize', 14);
            set(gca,'TickLength',[0 0])
            set(gca,'Xticklabel',figures.statcells)

            % Legend
            hLegend = legend(labels); 
            legend('boxoff');
            set(hLegend, 'position',[0.83 0.5 0.15 0.3]);
            set(hLegend, 'FontName', 'AvantGarde', 'FontSize', 12);

            % Axes
            ax = gca;
            ax.YGrid = 'on';
            ax2 = axes('Position',ax.Position,'XColor',[1 1 1],'YColor',[1 1 1],... 
              'Color','none','XTick',[],'YTick',[]);
            
            % Save figure
            saveas(gcf, sprintf('../output/figures/%s/%s.png', folder, game0.filename))
            close;
            
        end
        
        % Plot comparative statics
        function plot_compstats()
            
            % Loop over policies
            for policy=figures.policies
                disp(policy);
                pause(2);
                
                % Import data
                data = readtable(sprintf("../output/sumstats/%s.csv", policy));
                v1 = unique(data.sigma);
                v2 = unique(data.alpha);

                % Labels 
                labels = data.Properties.VariableNames(8:end);

                % Build subfigures
                for l=1:length(labels)

                    % Get compstats
                    compstats = figures.get_compstats(data, labels{l});

                    % Plot
                    figure();
                    contourf(compstats,20,'linestyle','none'); 
                    figures.prettify_compstats(compstats, l, v1, v2, 0);

                    % Save
                    saveas(gcf, sprintf("../output/compstats/%s/%s.png", policy, labels{l}));
                    close;
                end
            end
            
        end
                
        % Compare policies
        function compare_compstats()
            
            % Loop over policies
            for policy=figures.policies(2:end)
                disp(policy);
                pause(2);
                
                % Import data
                data = readtable(sprintf("../output/sumstats/%s.csv", policy));
                v1 = unique(data.sigma);
                v2 = unique(data.alpha);
                baseline_data = readtable("../output/sumstats/baseline.csv");

                % Labels 
                labels = data.Properties.VariableNames(8:end);

                % Build subfigures
                for l=1:length(labels)

                    % Save data
                    baseline_compstats = figures.get_compstats(baseline_data, labels{l});
                    compstats = figures.get_compstats(data, labels{l});
                    relative_compstats = compstats - baseline_compstats;

                    % Plot
                    figure();
                    contourf(relative_compstats,20,'linestyle','none'); 
                    figures.prettify_compstats(relative_compstats, l, v1, v2, 1);

                    % Save
                    saveas(gcf, sprintf("../output/compstats/%s/diff_%s.png", policy, labels{l}));
                    close;
                end
            end
            
        end
        
        % Prettify compstats
        function prettify_compstats(values, i, v1, v2, diff)
            
            % Set plot
            set(gca, 'OuterPosition', [0,0,1.05,0.95])
            set(gcf,'position',[100,1000,380,360]);
            set(gca,'YDir','normal');
            if diff
                caxis([-1, 1]*figures.difflimits(i));
            else
                negative_min = min(min(values)) < 0;
                positive_max = max(max(values)) > 0;
                caxis([-negative_min, positive_max]*figures.statlimits(i));
            end

            % Prettify
            colormap(centered(20)); 
            shading interp;
            set(gca, 'FontName', 'SansSerif', 'FontSize', 12);
            set(gca,'TickLength',[0 0]);
            set(gca,'XColor',[.4 .4 .4],'YColor',[.4 .4 .4]);

            % Axes
            yticks = [2, (numel(v1)+1)/2, numel(v1)-1];
            yticklabels = round(v1(max(1,floor(yticks))),1);
            xticks = [2, (numel(v2)+1)/2, numel(v2)-1];
            xticklabels = round(v2(max(1,floor(xticks))),1);
            set(gca,'XTick',xticks,'XTickLabel',xticklabels,'YTick',yticks,'YTickLabel',yticklabels);
            ylabel('Product Differentiation: \sigma','FontWeight','bold', 'Color', [.4 .4 .4])
            xlabel('Economies of scale: \alpha','FontWeight','bold', 'Color', [.4 .4 .4])
            
            % Title
            figtitle = figures.varnames(i);
            if diff
                figtitle = strcat("\Delta ", figtitle);
            end
            t = title(figtitle, 'FontSize', 18, 'FontWeight', 'normal');
            titlePos = get(t,'position');
            titlePos(2) = titlePos(2) * 1.05;
            set(t, 'position', titlePos);

            % Colorbar
            c = colorbar('TickLength',0,'YColor',[.4 .4 .4]);
            drawnow
            cdata = c.Face.Texture.CData;
            cdata(end,:) = uint8(0.7 * cdata(end,:));
            c.Face.Texture.ColorType = 'truecoloralpha';
            c.Face.Texture.CData = cdata;
        end
        
        
        
        % Make all figures
        function main()
           
            % Plot timelines
            figures.plot_all_timelines()
            
            % Plot alluvial
            figures.make_all_alluvialplots()
            
            % Plot compstats
            figures.plot_compstats()
            
            % Compare policies
            figures.compare_compstats()
        end
        
        
    end
end