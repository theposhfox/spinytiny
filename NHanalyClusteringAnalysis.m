function a = NHanalyClusteringAnalysis(varargin)

source = varargin{end};

varargin = varargin(1:end-1);

%%%%%%%%%%%%%%%%%%%%%%
%%Color Information%%%
%%%%%%%%%%%%%%%%%%%%%%

    lgray = [0.50 0.51 0.52];   brown = [0.28 0.22 0.14];
    gray = [0.50 0.51 0.52];    lbrown = [0.59 0.45 0.28];
    yellow = [1.00 0.76 0.05];  orange = [0.95 0.40 0.13];
    lgreen = [0.55 0.78 0.25];  green = [0.00 0.43 0.23];
    lblue = [0.00 0.68 0.94];   blue = [0.00 0.33 0.65];
    magenta = [0.93 0.22 0.55]; purple = [0.57 0.15 0.56];
    pink = [0.9 0.6 0.6];       lpurple  = [0.7 0.15 1];
    red = [0.85 0.11 0.14];     black = [0 0 0];
    dred = [0.6 0 0];          dorange = [0.8 0.3 0.03];
    bgreen = [0 0.6 0.7];
    colorj = {red,lblue,green,lgreen,gray,brown,yellow,blue,purple,lpurple,magenta,pink,orange,brown,lbrown};
    rnbo = {dred, red, dorange, orange, yellow, lgreen, green, bgreen, blue, lblue, purple, magenta, lpurple, pink}; 


if strcmpi(source, 'single')

    %%% Set variable lengths %%%
    Correlations = cell(1,14);
    PValues = cell(1,14);
    FarCorrelations = cell(1,14);
    withAPCorrelations = cell(1,14);
    withAP_PValues = cell(1,14);
    CausalCorrelations = cell(1,14);
    CausalPValues = cell(1,14);
    Distances = cell(1,14);
    FarDistances = cell(1,14);
    % FarPValues = cell(1,14);
    NumClusters = cell(1,14);
    ClustSpines = cell(1,14);
    NumSpInClust = cell(1,14);
    ClustDist = cell(1,14);
    CausalClustSpines = cell(1,14);
    NumSpinesInCausCluster = cell(1,14);
    CausalClustDist = cell(1,14);
    NumCausClusters = cell(1,14);
    FilesforSession = zeros(14,1);
    ReferenceTable = cell(1,14);


    for i = 1:length(varargin)
        session = varargin{i}.Session;
        Correlations{session} = [Correlations{session}; varargin{i}.SpineToSpineCorrelation];
        Correlations{session}(isnan(Correlations{session})) = 0;
        PValues{session} = [PValues{session}; varargin{i}.SpineToSpine_PValues];
        withAPCorrelations{session} = [withAPCorrelations{session}; varargin{i}.SpinewithAP_Correlation];
        withAP_PValues{session} = [withAP_PValues{session}; varargin{i}.SpinewithAP_PValues];
        CausalCorrelations{session} = [CausalCorrelations{session}; varargin{i}.CausalCorrelations];
        CausalCorrelations{session}(isnan(CausalCorrelations{session})) = 0;
        CausalPValues{session} = [CausalPValues{session}; varargin{i}.CausalPValues];
        Distances{session} = [Distances{session}; varargin{i}.SpineToSpineDistance];
        if varargin{i}.NumberofDendrites > 1
            FarCorrelations{session} = [FarCorrelations{session}; varargin{i}.FarSpineToSpineCorrelation];
            FarDistances{session} = [FarDistances{session}; varargin{i}.FarSpineToSpineDistance];
    %         FarPValues{session} = [FarPValues{session}; varargin{i}.
        else
            FarCorrelations{session} = [FarCorrelations{session}; nan]; %%% Since you can't compare between dendrites if you only analyze one, several files must have nan!!!!!!!!!!!!
            FarDistances{session} = [FarDistances{session}; nan];
        end
        sesh(i,1) = session;
        if ~isempty(varargin{i}.Clustered_Spines{1})
            NumClusters{session} = [NumClusters{session}; length(varargin{i}.Clustered_Spines)/varargin{i}.NumberofSpines];
        else
            NumClusters{session} = [NumClusters{session}; 0];
        end
        for j = 1:length(varargin{i}.Clustered_Spines)
            if ~isempty(varargin{i}.Clustered_Spines{j})        
                ClustSpines{session} = [ClustSpines{session}; varargin{i}.Clustered_Spines{j}];
                NumSpInClust{session} = [NumSpInClust{session}; length(varargin{i}.Clustered_Spines{j})];
                ClustDist{session} = [ClustDist{session}; varargin{i}.Cluster_Length{j}];
            end
        end
        if ~isempty(varargin{i}.CausalClustered_Spines{1})
            NumCausClusters{session} = [NumCausClusters{session}; length(varargin{i}.CausalClustered_Spines)/varargin{i}.NumberofSpines];
        else
            NumCausClusters{session} = [NumCausClusters{session}; 0];
        end
        for j = 1:length(varargin{i}.CausalClustered_Spines)
            if ~isempty(varargin{i}.CausalClustered_Spines{j})
                CausalClustSpines{session} = [CausalClustSpines{session}; varargin{i}.CausalClustered_Spines{j}];
                NumSpinesInCausCluster{session} = [NumSpinesInCausCluster{session}; length(varargin{i}.CausalClustered_Spines{j})];
                CausalClustDist{session} = [CausalClustDist{session}; varargin{i}.CausalCluster_Length{j}];
            end
        end
        FilesforSession(session,1) = FilesforSession(session, 1) + 1;
        ReferenceTable{session}{FilesforSession(session,1)} = varargin{i}.Filename;
    end

    session = unique(sesh);

    %%%%%%%%%%%%
    binstep = 5;
    %%%%%%%%%%%%

    shufflenum = 1000;
    options = statset('MaxIter', 10000); 

    for i = 1:14
        if isempty(find(session == i))
            continue
        end
        %%% Shuffle Data
        for j = 1:shufflenum
            shuffled{i}(1:length(Correlations{i}),j) = Correlations{i}(randperm(length(Correlations{i})));
            Farshuffled{i}(1:length(FarCorrelations{i}),j) = FarCorrelations{i}(randperm(length(FarCorrelations{i})));
            Causalshuffled{i}(1:length(CausalCorrelations{i}),j) = CausalCorrelations{i}(randperm(length(CausalCorrelations{i})));
        end

        shuffledmedian{i} = nanmedian(shuffled{i},2);
    %     pd = fitdist(shuffledmedian{i}, 'Normal');
    %     ci{i} = paramci(pd);
        [sortedDistances{i} sortedDistIndices{i}] = sort(Distances{i});
        [FarsortDist{i} FarsortDistInd{i}] = sort(FarDistances{i});

        if ~isempty(Distances{i})
            max_distance{i} = max(Distances{i});
            Filt_Corr{i} = Correlations{i}(PValues{i}<0.05 & Correlations{i}>=0.4);
            Filt_Dist{i} = Distances{i}(PValues{i}<0.05 & Correlations{i}>=0.4);
    %         Filt_Corr{i} = Correlations{i}(Correlations{i}>0.5);
    %         Filt_Dist{i} = Distances{i}(Correlations{i}>0.5);
            Filt_wAPCorr{i} = withAPCorrelations{i}(withAP_PValues{i}<0.05);
            Filt_Dist2{i} = Distances{i}(withAP_PValues{i}<0.05); %%% Need to select the distance values differently when pairing when spike-timed events
    %         Filt_Causal{i} = CausalCorrelations{i}(CausalPValues{i}<0.05);
    %         Filt_Dist_Causal{i} = Distances{i}(CausalPValues{i}<0.05);
            Filt_Causal{i} = CausalCorrelations{i}(CausalPValues{i}<0.05 & CausalCorrelations{i}>=0.4);
            Filt_Dist_Causal{i} = Distances{i}(CausalPValues{i}<0.05 & CausalCorrelations{i}>=0.4);
            Correlations_at_bin{i} = cell(1,length([1:binstep:max_distance{i}]));
            X = [];
            X(:,1) = Distances{i}; X(:,2) = Correlations{i};
            try
                obj = []; obj = fitgmdist(X,3, 'Options', options);
                idx = []; idx = cluster(obj,X);
                cluster1{i} = X(idx ==1,:); cluster2{i} = X(idx ==2,:); cluster3{i} = X(idx ==3,:);
                group = [median(cluster1{i}(:,2)) median(cluster2{i}(:,2)) median(cluster3{i}(:,2))];
                poss = [1 2 3];
                [val ind] = max(group);
                eval(['TopGroup{', num2str(i), '} = cluster', num2str(ind), '{', num2str(i), '};']); poss = poss(poss~=ind);
                [val ind] = min(group);
                eval(['BottomGroup{', num2str(i), '} = cluster', num2str(ind), '{', num2str(i), '};']); poss = poss(poss~=ind);
                eval(['MiddleGroup{', num2str(i), '} = cluster', num2str(poss), '{', num2str(i), '};']); 
                catch
                    TopGroup{i}(:,1) = Distances;
                    TopGroup{i}(:,2) = Correlations;
                    MiddleGroup{i}(:,1) = Distances;
                    MiddleGroup{i}(:,2) = Correlations;
                    BottomGroup{i}(:,1) = Distances;
                    BottomGroup{i}(:,2) = Correlations;
            end
            bincount = 1;
            for j = 1:binstep:max_distance{i}
                Correlations_at_bin{i}{bincount} = Correlations{i}(find(Distances{i}>= (j-1) &  Distances{i} < j+binstep));
                Binned{i}(bincount) = nanmedian(Correlations_at_bin{i}{bincount});
                Filt_Binned{i}(bincount) = nanmedian(Filt_Corr{i}(find(Filt_Dist{i}>= (j-1) & Filt_Dist{i} < j+binstep)));
                withAP_Binned{i}(bincount) = nanmedian(withAPCorrelations{i}(find(Distances{i}>=(j-1) & Distances{i} <j+binstep)));
                Filt_withAP_Binned{i}(bincount) = nanmedian(Filt_wAPCorr{i}(find(Filt_Dist2{i}>= (j-1) & Filt_Dist2{i} < j+binstep)));
                Causal_Binned{i}(bincount) = nanmedian(CausalCorrelations{i}(find(Distances{i}>=(j-1) & Distances{i} <j+binstep)));
                Filt_Causal_Binned{i}(bincount) = nanmedian(Filt_Causal{i}(find(Filt_Dist_Causal{i}>= (j-1) & Filt_Dist_Causal{i} < j+binstep)));
                Far_Binned{i}(bincount) = nanmedian(FarCorrelations{i}(find(FarDistances{i}>= (j-1) & FarDistances{i} < j+binstep)));
                try
                    HighCorr_Binned{i}(bincount) = nanmedian(TopGroup{i}(find(TopGroup{i}(:,1)>=(j-1) & TopGroup{i}(:,1) < j+binstep),2));
                    LowCorr_Binned{i}(bincount) = nanmedian(BottomGroup{i}(find(BottomGroup{i}(:,1)>=(j-1) & BottomGroup{i}(:,1) < j+binstep),2));
                catch
                    HighCorrBinned{i}(bincount) = Binned{i}(bincount);
                    LowCorr_Binned{i}(bincount) = Binned{i}(bincount);
                end
                Shuffled_Binned{i}(bincount) = nanmedian(shuffled{i}(find(Distances{i}(:,1)>=(j-1) & Distances{i}(:,1) < j+binstep),2));
                FarShuffled_Binned{i}(bincount) = nanmedian(Farshuffled{i}(find(FarDistances{i}(:,1)>=(j-1) & FarDistances{i}(:,1) < j+binstep),2));
                Causalshuffled_Binned{i}(bincount) = nanmedian(Causalshuffled{i}(find(Distances{i}(:,1)>=(j-1) & Distances{i}(:,1) < j+binstep),2));
                Subtracted_Binned{i}(bincount) = Binned{i}(bincount)-FarShuffled_Binned{i}(bincount);
    %             CausalSubtracted_Binned{i}(bincount) = 
                bincount = bincount+1;
            end
    %         Binned{i} = Binned{i}(Binned{i}~=0);
    %         Filt_Binned{i} = Filt_Binned{i}(Filt_Binned{i}~=0);
    %         withAP_Binned{i} = withAP_Binned{i}(withAP_Binned{i}~= 0);
    %         Filt_withAP_Binned{i} = Filt_withAP_Binned{i}(Filt_withAP_Binned{i} ~= 0);
    %         Causal_Binned{i} = Causal_Binned{i}(Causal_Binned{i}~= 0);
    %         Filt_Causal_Binned{i} = Filt_Causal_Binned{i}(Filt_Causal_Binned{i} ~= 0);
    %         Far_Binned{i} = Far_Binned{i}(Far_Binned{i}~=0);
        else
        end

    end


    %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%
    %%%% Figure 1 %%%%%
    %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%

    scrsz = get(0, 'ScreenSize');
    figure('Position', scrsz)
    extent = 40;


        % [bootx booty] = bootstrp(1000,@(tempD,tempC) fit(tempD,tempC,fittype), tempD,tempC)


        for i = 1:14
            subplot(4,7,i);
            if isempty(find(session == i))
                continue
            end
            if ~isempty(Distances{i})
                col1 = mod(i-1, length(colorj))+1;
                plot(FarDistances{i}, FarCorrelations{i}, '.', 'Color', lgray); hold on;
                plot(Distances{i}(PValues{i}>0.05), Correlations{i}(PValues{i}>0.05), '.k');
                plot(Distances{i}(PValues{i}<0.05), Correlations{i}(PValues{i}<0.05), '.', 'Color', rnbo{col1});
                xlim([-10 max(Distances{i})])
                xlabel('Distance (um)')
                ylim([-0.1 1]);
                ylabel('Spine-Spine Correlation')
                xlim([-5 60])
                title(['Session ', num2str(i)])
                subplot(4,7,14+i);
                %%% Uncomment to use top cluster for bar graphs
        %         bar(1:length(HighCorr_Binned{i}), HighCorr_Binned{i}, 'FaceColor', 'k', 'EdgeColor', rnbo{i}, 'LineWidth', 2); hold on;
                %%% Uncomment to use filtered data (p>0.05 correlation) for bars
%                 bar(1:length(Filt_Binned{i}), Filt_Binned{i},  'FaceColor', 'k', 'EdgeColor', rnbo{i}, 'Linewidth', 2); hold on;
        %         plot(1:length(Binned{i}), median(shuffledmedian{i}), 'color', rnbo{i})
                %%% Uncomment to use all data (no filtering) for bars
                bar(1:length(Binned{i}), Binned{i}, 'FaceColor', 'k', 'EdgeColor', rnbo{col1}, 'Linewidth', 2); hold on;
                %%% Uncomment to use subtracted data 
        %         bar(1:length(Subtracted_Binned{i}), Subtracted_Binned{i}, 'FaceColor', 'k', 'EdgeColor', rnbo{col1}, 'LineWidth', 2); hold on;
        %         bar(0, nanmedian(FarCorrelations{i}), 'FaceColor', gray, 'EdgeColor', 'k', 'LineWidth', 2)
        %         bar(1:length(Far_Binned{i}), Far_Binned{i}, 'FaceColor', gray, 'EdgeColor', 'k', 'Linewidth', 1)
                bar(0, nanmedian(FarCorrelations{i}), 'FaceColor', gray, 'EdgeColor', 'k', 'LineWidth', 2)
                ylim([-0.1 0.5])
                ylabel('median Correlation')
                xlim([-1 round(extent/binstep)+1])
                xtickcounter = [0, binstep*(1:length(Filt_Binned{i}))];
                set(gca, 'XTick', [0:(round(extent/binstep)+1)]);
                set(gca, 'XTickLabel', {'All', xtickcounter})
                xlabel(['Distance Bins (of ', num2str(binstep), ' um)'])
                title(['Session ', num2str(i)])
            else
            end
        end

        a.FarDistances = FarDistances;
        a.FarCorrelations = FarCorrelations;
        a.Distances = Distances;
        a.Correlations = Correlations;


        %%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%
        %%%% Figure 2 %%%%%
        %%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%

        %%%% Cummulative distributions of clustering of causal events (should show
        %%%% whether there is true clustering or if the effect is potentially an artifact of
        %%%% the limitations of imaged distances

        figure('Position', scrsz)

            
        for i = 1:14
            if isempty(find(session == i))
                continue
            end
            place = floor((i-1)/7);
            if place == 0
                row = 7;
            elseif place == 1
                row = 21;
            end
            if i == 7 || i == 14
                shift = 7;
            else
                shift = mod(i,7);
            end
            
            
            subplot(6,14, [row+shift]); hold on;
            for j = 1:shufflenum
                Farshuffled{i}(isnan(Farshuffled{i}(:,j)),j) = 0;
                sortedFarshuffled{i}(:,j) = Farshuffled{i}(FarsortDistInd{i},j)./nansum(Farshuffled{i}(FarsortDistInd{i},j));
                plot(FarsortDist{i},cumsum(sortedFarshuffled{i}(:,j)),'color', [0.5 0.5 0.5]);
            end
            Farshuffledmedian{i} = nanmean(sortedFarshuffled{i},2);
            plot(FarsortDist{i}, cumsum(Farshuffledmedian{i}), 'k', 'LineWidth', 2);
            FarCorrelations_fraction{i} = FarCorrelations{i}(FarsortDistInd{i},1)./nansum(FarCorrelations{i}(FarsortDistInd{i},1));
            Dstat = max(abs(cumsum(Farshuffledmedian{i})-cumsum(FarCorrelations_fraction{i})));   %%% Calculate the D-statistic for manual calculation of Kolmogorov-Smirnov test (using different distributions than the automated function uses)
            plot(FarsortDist{i},cumsum(FarCorrelations_fraction{i}),'color',rnbo{i}, 'LineWidth', 2)
        %     [h_Far p_Far stat_Far] = kstest2(cumsum(Farshuffledmedian{i}), cumsum(FarCorrelations_fraction{i}));
            c_alpha = 1.36; %%% constant based on desired p value of 0.05;
            D_alpha = c_alpha * (sqrt((length(Farshuffledmedian{i})+length(FarCorrelations_fraction{i}))/(length(Farshuffledmedian{i})*length(FarCorrelations_fraction{i}))));
            if Dstat > D_alpha
                text(40, 0.2, '*','FontSize', 12)
            else
                text(40, 0.2, 'n.s.')
            end
            ylim([-0.05 1.05])
            try
                xlim([min(FarsortDist{i})-1 max(FarDistances{i})])
            end

            if place == 0
                row = 27;
                row2 = 41;
            elseif place == 1
                row = 55;
                row2 = 69;
            end
            
            
            subplot(6,14,[row+(2*shift) row+((2*shift)+1) row2+(2*shift) row2+((2*shift)+1)]); hold on;
            for j = 1:shufflenum
                sortedshuffled{i}(:,j) = shuffled{i}(sortedDistIndices{i},j)./nansum(shuffled{i}(sortedDistIndices{i},j));
                plot(sortedDistances{i}, cumsum(sortedshuffled{i}(:,j)),'color', [0.5 0.5 0.5]);
            end
            shuffledmedian{i} = nanmean(sortedshuffled{i},2);
            plot(sortedDistances{i},cumsum(shuffledmedian{i}), 'k', 'LineWidth', 2);
            Correlations_fraction{i} = Correlations{i}(sortedDistIndices{i},1)./nansum(Correlations{i}(sortedDistIndices{i},1));
            plot(sortedDistances{i},cumsum(Correlations_fraction{i}),'color',rnbo{i}, 'LineWidth', 2)
        %     [h p stat] = kstest2(cumsum(shuffledmedian{i}), cumsum(Correlations_fraction{i}));
        %     [h p stat] = kstest2(shuffled{i}, Correlations{i});
            D_alpha = c_alpha * (sqrt((length(shuffledmedian{i})+length(Correlations_fraction{i}))/(length(shuffledmedian{i})*length(Correlations_fraction{i}))));
            if Dstat > D_alpha
                text(40, 0.2, '*','FontSize', 12)
            else
                text(40, 0.2, 'n.s.')
            end
            ylim([-0.05 1.05])
            xlim([-0.05 max(Distances{i})])
            if i == 1 || i == 8
                xlabel('Distance (um)')
                ylabel('Cumulative fraction of correlations')
            end
        end



        %%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%
        %%%% Figure 3 %%%%%
        %%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%

        figure('Position', scrsz)

        for i = 1:14
            subplot(4,7,i);
            if isempty(find(session == i))
                continue
            end
            if ~isempty(Distances{i})
                col1 = mod(i-1, length(colorj))+1;
        %         plot(Distances{i}, withAPCorrelations{i}, '.k');  hold on;
        %         plot(Filt_Dist2{i}, Filt_wAPCorr{i}, '.', 'Color', rnbo{col1})
                plot(Distances{i}, CausalCorrelations{i}, '.k');  hold on;
                plot(Filt_Dist_Causal{i}, Filt_Causal{i}, '.', 'Color', rnbo{col1})
        %         yfit = fit(Distances{i}(Distances{i}<=extent), CausalCorrelations{i}(Distances{i}<=extent), 'exp1');
        %         plot(yfit, '-k'); legend off;
                xlim([-10 max(Distances{i})])
                xlabel('Distance (um)')
                ylim([-0.1 1]);
                ylabel('Spine-Spine Correlation')
                xlim([-5 60])
                title(['Session ', num2str(i)])
                subplot(4,7,14+i);
            %     bar(1:length(Filt_Binned{i}), Filt_Binned{i}, 'c'); hold on;
        %         bar(1:length(withAP_Binned{i}), withAP_Binned{i},'FaceColor', 'k', 'EdgeColor', rnbo{col1}, 'Linewidth', 2); hold on;
                bar(1:length(Filt_Causal_Binned{i}), Filt_Causal_Binned{i},'FaceColor', 'k', 'EdgeColor', rnbo{col1}, 'Linewidth', 2); hold on;
        %         bar(1:length(Causal_Binned{i}), Causal_Binned{i},'FaceColor', 'k', 'EdgeColor', rnbo{col1}, 'Linewidth', 2); hold on;
                %     bar(0, Filt_Far{i}, 'g')
        %         bar(0, nanmedian(FarCorrelations{i}), 'FaceColor', gray, 'EdgeColor', 'k', 'LineWidth', 2)
                bar(0, nanmedian(FarCorrelations{i}), 'FaceColor', gray, 'EdgeColor', 'k', 'LineWidth', 2)
                bar(1:length(Far_Binned{i}), Far_Binned{i}, 'FaceColor', gray, 'EdgeColor', 'k', 'Linewidth', 1)
                ylim([-0.1 0.5])
                ylabel('median Correlation')
                xlim([-1 round(extent/binstep)+1])
                xtickcounter = [0, binstep*(1:length(Filt_Binned{i}))];
                set(gca, 'XTick', [0:(round(extent/binstep)+1)]);
                set(gca, 'XTickLabel', {'All', xtickcounter})
                xlabel(['Distance Bins (of ', num2str(binstep), ' um)'])
                title(['Session ', num2str(i)])
            else
            end
        end

        %%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%
        %%%% Figure 4 %%%%%
        %%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%


        %%%% Cummulative distributions of clustering of causal events (should show
        %%%% whether there is true clustering or if the effect is potentially an artifact of
        %%%% the limitations of imaged distances

        % for i = 1:length(session)
        %     if ~isempty(Distances{i})
        %         try
        %             casualclust_tc(1,i) = Causal_Binned{i}(1);    %%% Estimate the timecourse of clustering over learning
        %             causalclust_tc(2,i) = median([Causal_Binned{i}(2);Causal_Binned{i}(3); Causal_Binned{i}(4)]);
        %             causalclust_tc(3,i) = median([Causal_Binned{i}(5); Causal_Binned{i}(6)]);
        %             causalclust_tc(4,i) = median([Causal_Binned{i}(7); Causal_Binned{i}(8)]);
        %         end
        %     end
        % end
        % 
        % % cluster_tc = log(clust_tc.*100);
        %     clusttc_h = figure('Position', scrsz); subplot(6,14,[1 2 3 4 5 6 7 15 16 17 18 19 20 21]); 
        % try
        %     plot(1:length(session), casualclust_tc(1,:), '.-', 'Color', red, 'LineWidth', 2)
        %     hold on; plot(causalclust_tc(2,:), '.-', 'Color', lblue, 'LineWidth', 2);
        %     plot(causalclust_tc(3,:), '.-','Color', green, 'LineWidth', 2)
        %     plot(causalclust_tc(4,:), '.-','Color',lgreen, 'LineWidth', 2)
        %     
        %     xlabel('Session');
        %     ylabel('median Correlation');
        % 
        %     for j = 1:4
        %         for i = 1:7
        %             causalnewb(j,i) = median(causalclust_tc(j,(2*(i-1)+1):(2*(i-1)+2)));  %%% Average in two-day bins (might show trend more clearly...)
        %         end
        %         plot([1.5 3.5 5.5 7.5 9.5 11.5 13.5], causalnewb(j,:), '*--', 'Color', colorj{j})
        %         legend({['Correlation at ', num2str(binstep), ' um'], ['Correlation from ', num2str(2*binstep), ' to ', num2str(4*binstep), ' um'], ['Correlation from ', num2str(5*binstep), ' to ', num2str(6*binstep), ' um'], ['Correlation from ', num2str(7*binstep), ' to ', num2str(8*binstep), ' um']})
        %     end
        % 
        % catch
        %     spy
        % end
        % 
        % for i = 1:14
        %     if isempty(find(session == i))
        %         continue
        %     end
        %     place = floor((i-1)/7);
        %     if place == 0
        %         row = 7;
        %     elseif place == 1
        %         row = 21;
        %     end
        %     if i == 7 || i == 14
        %         shift = 7;
        %     else
        %         shift = mod(i,7);
        %     end
        %     subplot(6,14, [row+shift]); hold on;
        %     for j = 1:shufflenum
        %         Farshuffled{i}(isnan(Farshuffled{i}(:,j)),j) = 0;
        %         sortedFarshuffled{i}(:,j) = Farshuffled{i}(FarsortDistInd{i},j)./nansum(Farshuffled{i}(FarsortDistInd{i},j));
        %         plot(FarsortDist{i},cumsum(sortedFarshuffled{i}(:,j)),'color', [0.5 0.5 0.5]);
        %     end
        %     Farshuffledmedian{i} = nanmean(sortedFarshuffled{i},2);
        %     plot(FarsortDist{i}, cumsum(Farshuffledmedian{i}), 'k', 'LineWidth', 2);
        %     FarCorrelations_fraction{i} = FarCorrelations{i}(FarsortDistInd{i},1)./nansum(FarCorrelations{i}(FarsortDistInd{i},1));
        %     Dstat = max(abs(cumsum(Farshuffledmedian{i})-cumsum(FarCorrelations_fraction{i})));   %%% Calculate the D-statistic for manual calculation of Kolmogorov-Smirnov test (using different distributions than the automated function uses)
        %     plot(FarsortDist{i},cumsum(FarCorrelations_fraction{i}),'color',rnbo{i}, 'LineWidth', 2)
        % %     [h_Far p_Far stat_Far] = kstest2(cumsum(Farshuffledmedian{i}), cumsum(FarCorrelations_fraction{i}));
        %     c_alpha = 1.36; %%% constant based on desired p value of 0.05;
        %     D_alpha = c_alpha * (sqrt((length(Farshuffledmedian{i})+length(FarCorrelations_fraction{i}))/(length(Farshuffledmedian{i})*length(FarCorrelations_fraction{i}))));
        %     if Dstat > D_alpha
        %         text(40, 0.2, '*','FontSize', 12)
        %     else
        %         text(40, 0.2, 'n.s.')
        %     end
        %     ylim([-0.05 1.05])
        %     try
        %         xlim([min(FarsortDist{i})-1 max(FarDistances{i})])
        %     end
        % 
        %     if place == 0
        %         row = 27;
        %         row2 = 41;
        %     elseif place == 1
        %         row = 55;
        %         row2 = 69;
        %     end
        %     subplot(6,14,[row+(2*shift) row+((2*shift)+1) row2+(2*shift) row2+((2*shift)+1)]); hold on;
        %     for j = 1:shufflenum
        %         sortedcausalshuffled{i}(:,j) = Causalshuffled{i}(sortedDistIndices{i},j)./nansum(Causalshuffled{i}(sortedDistIndices{i},j));
        %         plot(sortedDistances{i}, cumsum(sortedcausalshuffled{i}(:,j)),'color', [0.5 0.5 0.5]);
        %     end
        %     Causalshuffledmedian{i} = nanmean(sortedcausalshuffled{i},2);
        %     plot(sortedDistances{i},cumsum(Causalshuffledmedian{i}), 'k', 'LineWidth', 2);
        %     CausalCorrelations_fraction{i} = CausalCorrelations{i}(sortedDistIndices{i},1)./nansum(CausalCorrelations{i}(sortedDistIndices{i},1));
        %     plot(sortedDistances{i},cumsum(CausalCorrelations_fraction{i}),'color',rnbo{i}, 'LineWidth', 2)
        % %     [h p stat] = kstest2(cumsum(shuffledmedian{i}), cumsum(Correlations_fraction{i}));
        % %     [h p stat] = kstest2(shuffled{i}, Correlations{i});
        %     D_alpha = c_alpha * (sqrt((length(Causalshuffledmedian{i})+length(CausalCorrelations_fraction{i}))/(length(Causalshuffledmedian{i})*length(CausalCorrelations_fraction{i}))));
        %     if Dstat > D_alpha
        %         text(40, 0.2, '*','FontSize', 12)
        %     else
        %         text(40, 0.2, 'n.s.')
        %     end
        %     ylim([-0.05 1.05])
        %     xlim([-0.05 max(Distances{i})])
        %     if i == 1 || i == 8
        %         xlabel('Distance (um)')
        %         ylabel('Cumulative fraction of correlations')
        %     end
        % end
        % 
        % disp('wahh')
        % 
        % figure('Position', scrsz/2)
        % 
        % for i = 1:14
        %     meannumclusters(1,i) = nanmean(NumClusters{i});
        %     numclustersSTD(1,i) = nanstd(NumClusters{i})/sqrt(length(NumClusters{i}));
        %     meannumCclusters(1,i) = nanmean(NumCausClusters{i});
        %     numCclustersSTD(1,i) = nanstd(NumCausClusters{i})/sqrt(length(NumCausClusters{i}));
        %     meanNSPIC(1,i) = nanmean(NumSpInClust{i});
        %     NSPIC_STD(1,i) = nanstd(NumSpInClust{i})/sqrt(length(NumSpInClust{i}));
        %     meanCNSPIC(1,i) = nanmean(NumSpinesInCausCluster{i});
        %     CNSPIC_STD(1,i) = nanstd(NumSpinesInCausCluster{i})/sqrt(length(NumSpinesInCausCluster{i}));
        %     meanDist(1,i) = nanmean(ClustDist{i});
        %     Dist_STD(1,i) = nanstd(ClustDist{i})/sqrt(length(ClustDist{i}));
        %     meanCausDist(1,i) = nanmean(CausalClustDist{i});
        %     CausDist_STD(1,i) = nanstd(CausalClustDist{i})/sqrt(length(CausalClustDist{i}));
        % end
        % 
        % 
        % subplot(1,3,1); 
        % plot(meannumclusters, 'k', 'Linewidth', 2); hold on;
        % plot(meannumCclusters, 'r', 'LineWidth', 2)
        % legend({'Synapse Only', 'Causal'})
        % r_errorbar(1:14, meannumclusters, numclustersSTD, 'k')
        % r_errorbar(1:14, meannumCclusters, numCclustersSTD, 'r')
        % title('Number of Clusters')
        % xlim([0 15])
        % 
        % subplot(1,3,2)
        % plot(meanNSPIC, 'k', 'Linewidth', 2); hold on;
        % plot(meanCNSPIC, 'r', 'Linewidth', 2)
        % legend({'Synapse Only', 'Causal'})
        % r_errorbar(1:14, meanNSPIC, NSPIC_STD, 'k')
        % r_errorbar(1:14, meanCNSPIC, CNSPIC_STD, 'r')
        % title('Mean # of Spines in Each Cluster')
        % xlim([0 15])
        % 
        % subplot(1,3,3)
        % plot(meanDist, 'k', 'Linewidth', 2); hold on;
        % plot(meanCausDist, 'r', 'Linewidth', 2);
        % legend({'Synapse Only', 'Causal'})
        % r_errorbar(1:14, meanDist, Dist_STD, 'k')
        % r_errorbar(1:14, meanCausDist, CausDist_STD, 'r')
        % title('Mean Spatial Extent of Clusters')
        % xlim([0 15])

elseif strcmpi(source, 'multiple')
    
    Distances = cell(1,14);
    Correlations = cell(1,14);
    FarDistances = cell(1,14);
    FarCorrelations = cell(1,14);
    
    for i = 1:length(varargin)
        for j = 1:14
            Distances{j} = [Distances{j}; varargin{i}.Distances{j}];
            Correlations{j} = [Correlations{j}; varargin{i}.Correlations{j}];
            FarDistances{j} = [FarDistances{j}; varargin{i}.FarDistances{j}];
            FarCorrelations{j} = [FarCorrelations{j}; varargin{i}.FarCorrelations{j}];
        end
    end
    
    scrsz = get(0, 'ScreenSize');
    figure('Position', scrsz)
    extent = 40;
    
    for i = 1:14
        ax(i) = subplot(2,7,i);
        plot(Distances{i},Correlations{i}, 'o', 'Color', rnbo{i}); hold on;
        decay = fit(Distances{i},Correlations{i}, 'exp1'); 
        fline = plot(decay); 
        set(fline, 'Color', 'k')
        legend off
        plot(-1/decay.b,0.368*decay.a, 'ok', 'MarkerFaceColor', 'k')
        xlabel('Distance', 'Fontsize', 14)
        ylabel('Correlation', 'Fontsize', 14)
        plot(0:100, zeros(1,101), '--k', 'Linewidth', 2)
        title(['Session ', num2str(i)], 'Fontsize', 14)
    end
    linkaxes(ax);
end



