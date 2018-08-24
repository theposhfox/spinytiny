function [SpectralData] = SpectralClustering(inputData, Correlations, StatClass, Choices)

DendNum = inputData.NumberofDendrites;
session = inputData.Session;

Temporal_Laplacian = cell(1,DendNum);
Temporal_Deg = cell(1,DendNum);
    Dend_Temp_Deg = [];
    Temporal_FirstEigenvector = cell(1,DendNum);
Spatial_Deg = cell(1,DendNum);
    Dend_Spat_Deg = [];
    Spatial_FirstEigenvector = cell(1,DendNum);
Spatiotemporal_Deg = cell(1,DendNum);
Spatiotemporal_Adjacency = cell(1,DendNum);
Spatiotemporal_Laplacian = cell(1,DendNum);
    Dend_SpatTemp_Deg = [];
    SpatioTemporal_FirstEigenvector = cell(1,DendNum);
    
SpectralLengthConstant = 10;

for j = 1:DendNum     
    firstspine = inputData.SpineDendriteGrouping{j}(1);
    lastspine = inputData.SpineDendriteGrouping{j}(end);
    if j == DendNum && lastspine ~= length(inputData.dF_over_F)
        lastspine = length(inputData.dF_over_F);
    end
%     if Choices.UseStatDends
%         if isempty(StatClass{session})
%             SpectralData = [];
%             return
%         end
%         if ~StatClass{session}.MovementDendrites(j)
%             Spatial_Deg{j} = nan(length(firstspine:lastspine),1);
%             Spatiotemporal_Deg{j} = nan(length(firstspine:lastspine),1);
%             Temporal_Deg{j} = nan(length(firstspine:lastspine),1);
%             continue
%         end
%     end
    if firstspine ~= lastspine
        fullDist = full(triu(inputData.DistanceHeatMap));
        ad = fullDist(firstspine:lastspine, firstspine:lastspine);
        ad(ad<1) = 1;
        ad = 1./exp(ad./SpectralLengthConstant);          %%% Adjacency matrix --> 1/e^x, where x = distance
        ad(isnan(ad)) = 0;                                %%% Since the diagonal of the laplacian == the degree, set NaNs in A to be zero to maintain this identity;
%         ad = inputData.AdjacencyMatrix{j};
        ad(ad==0) = nan;
        Spatial_Deg{j} = nanmean(ad,2);
        if strcmpi(Choices.LaplacianToUse, 'Normalized')
            L = inputData.NormalizedLaplacian{j};
        elseif strcmpi(Choices.LaplacianToUse, 'Original')
            L = inputData.LaplacianMatrix{j};
        end
        [eigenvector, ~] = eigs(L, 1, 'lm'); %% Find largest real magnitude eigenvector
            Spatial_FirstEigenvector{j} = abs(eigenvector);

            %%% The temporal correlation matrix will serve as
            %%% the temporal version of the adjacency matrix...
            corrmat = Correlations(Choices.Spine1_Address+firstspine:Choices.Spine1_Address+lastspine, Choices.Spine1_Address+firstspine:Choices.Spine1_Address+lastspine);
            corrmat(corrmat==1) = 0;
            corrmat(isnan(corrmat)) = 0;
            corrmat(corrmat<0) = 0;

%                     Spatial_Deg{session}{j}(Spatial_Deg{session}{j} == 0) = nan;

        Temporal_Deg{j} = nanmean(corrmat,2);                                %%% Temporal degree matrix
        Temporal_Deg{j} = Temporal_Deg{j};
        Temporal_Laplacian{j} = diag(Temporal_Deg{j})-corrmat;     %%% L = D-A;
        if strcmpi(Choices.LaplacianToUse, 'Normalized')
            fullmat = diag(Temporal_Deg{j});
            fullmat(fullmat==0) = eps;
            invTempDeg = inv(fullmat+eye(size(fullmat,1))*1e-4);                         %%% Inverse D for calculating normalized laplacian (add small values to identity to prevent from inv of zeros resulting in inf
            Temporal_Laplacian{j} = invTempDeg*Temporal_Laplacian{j};  %%% Normalize the Laplacian by multiplying by the inverse degree matrix
        else
        end
        [tvecs, ~] = eigs(Temporal_Laplacian{j}, 1, 'lm');
            Temporal_FirstEigenvector{j} = abs(tvecs(:,1));
        Spatiotemporal_Deg{j} = Spatial_Deg{j}.*Temporal_Deg{j};
%                     Spatiotemporal_Deg{session}{j} = Spatiotemporal_Deg{session}{j}(find(Spatiotemporal_Deg{session}{j}));
        Spatiotemporal_Adjacency{j} = inputData.AdjacencyMatrix{j}.*corrmat;
        Spatiotemporal_Adjacency{j}(isnan(Spatiotemporal_Adjacency{j})) = 0;
        Spatiotemporal_Laplacian{j} = full(diag(Spatiotemporal_Deg{j}))-Spatiotemporal_Adjacency{j};
        if strcmpi(Choices.LaplacianToUse, 'Normalized')
            fullmat = diag(Spatiotemporal_Deg{j});
            fullmat(fullmat==0) = eps;
            invSTDeg = inv(fullmat+eye(size(fullmat,1))*1e-4);                          %%% adding a small value to the identity matrix prevents inf values
            Spatiotemporal_Laplacian{j} = invSTDeg*Spatiotemporal_Laplacian{j};
        else
        end
            [stvecs, stvals] = eig(Spatiotemporal_Laplacian{j});
            stvals = diag(stvals);
            stvals(stvals == min(stvals)) = nan;
            [val ind] = min(stvals);
                SpatioTemporalFiedler(1,j) = val;
                SpatioTemporalPartition{j} = stvecs(:,ind);
            [val ind] = max(stvals);
                SpatioTemporal_FirstEigenvector{j} = stvecs(:,ind);

        MovementCorrelationsforAllSpinesonDend{j} = Correlations(Choices.MovementAddress,Choices.Spine1_Address+firstspine:Choices.Spine1_Address+lastspine)';

        MovementCorrelationsforAllSpinesonDend{j} = MovementCorrelationsforAllSpinesonDend{j}.*(MovementCorrelationsforAllSpinesonDend{j});

        if ~isempty(Spatiotemporal_Deg{j})
            address = ~isnan(MovementCorrelationsforAllSpinesonDend{j});
%                         Spatiotemporal_Deg{session}{j} = diag(Spatiotemporal_Deg{session}{j});
                SpatioTemporal_FirstEigenvector{j} = diag(SpatioTemporal_FirstEigenvector{j});
%                         Spatiotemporal_Deg{session}{j} = Spatiotemporal_Deg{session}{j}(address);
                SpatioTemporal_FirstEigenvector{j} = SpatioTemporal_FirstEigenvector{j}(address);
        else
            Spatiotemporal_Deg{j} = nan;
            SpatioTemporal_FirstEigenvector{j} = nan;
        end
        useSmat = Spatial_Deg{j};
%                     useSmat = Spatial_FirstEigenvector{session}{j};
        useTmat = Temporal_Deg{j};
%                     useTmat = Temporal_FirstEigenvector{session}{j};
        useSTmat = Spatiotemporal_Deg{j};
%                     useSTmat = SpatioTemporal_FirstEigenvector{session}{j};
        try
            SpatiotemporalDegree_vs_Movement(1,j) = corr(useSTmat, MovementCorrelationsforAllSpinesonDend{j});
            SpatialDegree_vs_Movement(1,j) = corr(useSmat, MovementCorrelationsforAllSpinesonDend{j});
            TemporalDegree_vs_Movement(1,j) = corr(useTmat, MovementCorrelationsforAllSpinesonDend{j});
            Spatiotemporal_Overlap(1,j) = corr(useSmat, useTmat);
        catch
            SpatiotemporalDegree_vs_Movement(1,j) = NaN;
            SpatialDegree_vs_Movement(1,j) = NaN;
            TemporalDegree_vs_Movement(1,j) = NaN;
        end
    else
        Spatial_Deg{j} = nan;
        Temporal_Deg{j} = nan;
        Spatiotemporal_Deg{j} = nan;

        SpatialDegree_vs_Movement(1,j) = NaN;
        TemporalDegree_vs_Movement(1,j) = NaN;
        SpatiotemporalDegree_vs_Movement(1,j) = NaN;
    end
    Dend_Spat_Deg = [Dend_Spat_Deg; nanmean(Spatial_Deg{j})];
    Dend_Temp_Deg = [Dend_Temp_Deg; nanmean(Temporal_Deg{j})];
    Dend_SpatTemp_Deg = [Dend_SpatTemp_Deg; nanmean(Spatiotemporal_Deg{j})];
end

allspatialspines = cat(1,Spatial_Deg{:});
alltemporalspines = cat(1,Temporal_Deg{:});
allSTspines = cat(1,Spatiotemporal_Deg{:});

SpectralData.Spatial_Deg = Spatial_Deg;
SpectralData.Temporal_Deg = Temporal_Deg;
SpectralData.Spatiotemporal_Deg = Spatiotemporal_Deg;

SpectralData.SpatioTemporalFiedler = SpatioTemporalFiedler;

SpectralData.Temporal_Laplacian = Temporal_Laplacian;

SpectralData.Temporal_Laplacian = Temporal_Laplacian ;

try
    SpectralData.MeanSpatialDegreeofCueSpines = nanmean(allspatialspines(StatClass{session}.CueSpines));
    SpectralData.MeanTemporalDegreeofCueSpines = nanmean(alltemporalspines(StatClass{session}.CueSpines));
    SpectralData.MeanSpatioTemporalDegreeofCueSpines = nanmean(allSTspines(StatClass{session}.CueSpines));
    SpectralData.MeanSpatialDegreeofMovementSpines = nanmean(allspatialspines(StatClass{session}.MovementSpines));
    SpectralData.MeanTemporalDegreeofMovementSpines = nanmean(alltemporalspines(StatClass{session}.MovementSpines));
    SpectralData.MeanSpatioTemporalDegreeofMovementSpines = nanmean(allSTspines(StatClass{session}.MovementSpines));
    SpectralData.MeanSpatialDegreeofMovementDuringCueSpines = nanmean(allspatialspines(StatClass{session}.MovementDuringCueSpines));
    SpectralData.MeanTemporalDegreeofMovementDuringCueSpines = nanmean(alltemporalspines(StatClass{session}.MovementDuringCueSpines));
    SpectralData.MeanSpatioTemporalDegreeofMovementDuringCueSpines = nanmean(allSTspines(StatClass{session}.MovementDuringCueSpines));
    SpectralData.MeanSpatialDegreeofPreSuccessSpines = nanmean(allspatialspines(StatClass{session}.PreSuccessSpines));
    SpectralData.MeanTemporalDegreeofPreSuccessSpines = nanmean(alltemporalspines(StatClass{session}.PreSuccessSpines));
    SpectralData.MeanSpatioTemporalDegreeofPreSuccessSpines = nanmean(allSTspines(StatClass{session}.PreSuccessSpines));
    SpectralData.MeanSpatialDegreeofSuccessSpines = nanmean(allspatialspines(StatClass{session}.SuccessSpines));
    SpectralData.MeanTemporalDegreeofSuccessSpines = nanmean(alltemporalspines(StatClass{session}.SuccessSpines));
    SpectralData.MeanSpatioTemporalDegreeofSuccessSpines = nanmean(allSTspines(StatClass{session}.SuccessSpines));
    SpectralData.MeanSpatialDegreeofRewardSpines = nanmean(allspatialspines(StatClass{session}.RewardSpines));
    SpectralData.MeanTemporalDegreeofRewardSpines = nanmean(alltemporalspines(StatClass{session}.RewardSpines));
    SpectralData.MeanSpatioTemporalDegreeofRewardSpines = nanmean(allSTspines(StatClass{session}.RewardSpines));
catch
    SpectralData.MeanSpatialDegreeofCueSpines = nan;
    SpectralData.MeanTemporalDegreeofCueSpines = nan;
    SpectralData.MeanSpatioTemporalDegreeofCueSpines = nan;
    SpectralData.MeanSpatialDegreeofMovementSpines = nan;
    SpectralData.MeanTemporalDegreeofMovementSpines = nan;
    SpectralData.MeanSpatioTemporalDegreeofMovementSpines = nan;
    SpectralData.MeanSpatialDegreeofMovementDuringCueSpines = nan;
    SpectralData.MeanTemporalDegreeofMovementDuringCueSpines = nan;
    SpectralData.MeanSpatioTemporalDegreeofMovementDuringCueSpines = nan;
    SpectralData.MeanSpatialDegreeofPreSuccessSpines = nan;
    SpectralData.MeanTemporalDegreeofPreSuccessSpines = nan;
    SpectralData.MeanSpatioTemporalDegreeofPreSuccessSpines = nan;
    SpectralData.MeanSpatialDegreeofSuccessSpines = nan;
    SpectralData.MeanTemporalDegreeofSuccessSpines = nan;
    SpectralData.MeanSpatioTemporalDegreeofSuccessSpines = nan;
    SpectralData.MeanSpatialDegreeofRewardSpines = nan;
    SpectralData.MeanTemporalDegreeofRewardSpines = nan;
    SpectralData.MeanSpatioTemporalDegreeofRewardSpines = nan;
end

SpectralData.Spatiotemporal_Overlap = Spatiotemporal_Overlap;

SpectralData.SpatiotemporalDegree_vs_Movement = SpatiotemporalDegree_vs_Movement;
SpectralData.SpatialDegree_vs_Movement = SpatialDegree_vs_Movement;
SpectralData.TemporalDegree_vs_Movement = TemporalDegree_vs_Movement;

SpectralData.Dend_Spat_Deg = Dend_Spat_Deg;
SpectralData.Dend_Temp_Deg = Dend_Temp_Deg;
SpectralData.Dend_SpatTemp_Deg = Dend_SpatTemp_Deg;



