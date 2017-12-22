function ClusterBehaviorCorrelations(varargin)

%%% Input argument list:
%%% If running this on a single animal, inputs are the following:
%%% 1: All Correlation matrices from NHanalyAlignBehavior
%%% 2: All statistically classified structures from NHanalyAlignBehavior
%%% Otherwise, if running with multiple animals, the inputs are the
%%% analyzed files generated from THIS program from each animal

global gui_KomiyamaLabHub

dendsubtract = get(gui_KomiyamaLabHub.figure.handles.DendSubtracted_CheckBox, 'value');
dendexclude = get(gui_KomiyamaLabHub.figure.handles.DendExcluded_CheckBox, 'value');

if isempty(strfind(inputname(1), 'SpineCorrelationTimecourse'))
    %%% Define the inputs
    
    Correlations = varargin{1};
    StatClass = varargin{2};
    firstdatainput = 3;
    
    if length(Correlations)<14
        Correlations{14} = [];
    end
    
    %%%
%     MovementSpineDataToUse = StatClass.MovementSpines;
%     MovementSpineDataToUse = LiberalSpines;
%     CausalMovementDataToUse = StatClass.CausalMovementSpines;
%     CausalMovementDataToUse = CausalLiberalSpines;
    %%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Set addresses for different correlations sets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Cue_address = 1;
    Movement_address = 2; %%% row 2 is binarized movement, row 3 has a larger window
    Presuccess_address = 4; 
    Success_address = 5;  %%% row 5 is binarized rewarded movements, row 6 has a larger window
    MovementDuringCue_address = 7;
    Reward_address = 8;
    Punishment_address = 9;
    Spine1_address = 9;   %%% Actually +1, but easier to add other spine addresses to zero...
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%
    %%% Contingencies %%%
    %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    usespineclassification = 0;
    usecorrminimum = 0;
    usecuecorrminimum = 0;
    corrminimum = 0;
    cuecorrminimum = 0;
    correlatecuewithall = 1;
    useonlyspinesfromstatreldends = 0;         %%% Use only spines on MOVEMENT-RELATED dendrites (applies to spine analysis, not that of dendrites)
        useSTATdends = 0;                      %%% Provides a different contingency option for using only movement-related dendrites for just the "dendritic" portions of the analysis (i.e. not for spines)
    mindendsize = 2;                           %%% Minimum number of spines on a dendrite to be considered for analysis (probably some sensitivity to length of dendrite)
    laplaciantouse = 'Normalized';             %%% Choose 'Normalized' or 'Original'
    
    
    Choices.UseStatisticalClassification = usespineclassification;
    Choices.UseCorrMinimum = usecorrminimum;
    Choices.CorrMinimum = corrminimum;
    Choices.CueCorrMinimum = cuecorrminimum; 
    Choices.CorrelateCueWithAll = correlatecuewithall;
    Choices.UseOnlySpinesFromStatRelDends = useonlyspinesfromstatreldends;
    Choices.UseStatDends = useSTATdends;
    Choices.MinDendSize = mindendsize;
    Choices.LaplacianToUse = laplaciantouse;
    Choices.Spine1_Address = Spine1_address;
    Choices.MovementAddress = Movement_address;
    
    
    %%% Spine data being used: 
    
    %%% Options:
    %%% OverallSpineCorrelations
    %%% SpineCorrelations
    %%% DendSubtractedSpineCorrelations
    %%% CausalCorrelations

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Display above filters in the command window
    
    if usespineclassification
        disp('Clusters ARE being classified by function');
    else
        disp('Clusters are NOT being filtered by statistical relationship with movement');
    end
    %disp('Clusters are being filtered by dominance with cue/movement');
    
    if usecorrminimum
        disp(['Using only correlations above ', num2str(corrminimum), ' for analysis']);
    else
        disp('Not using a minimum correlation')
    end
    
    animal = regexp(varargin{firstdatainput}.Filename, ['[ABCDEFGHIJKLMNOPQRSTUVWXYZ]{2}0+\d+'], 'match');
    
    disp(['File ', animal{1}, ' analyzed'])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% Initialize variables %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    CorrelationofClusters = cell(1,14);
    
    clustcorrwithcue = nan(1,14);
        nonclustcorrwithcue = nan(1,14);
    cueclustcorrwithcue = nan(1,14);
        cuenonclustcorrwithcue = nan(1,14);
        
    clustcorrwithmov = nan(1,14);
        nonclustcorrwithmov = nan(1,14);
    movclustcorrwithmov = nan(1,14);
        movnonclustcorrwithmov = nan(1,14);
        
    clustcorrwithMDC = nan(1,14);
        nonclustcorrwithMDC = nan(1,14);
    MDCclustcorrwithMDC = nan(1,14);
        MDCnonclustcorrwithMDC = nan(1,14);
        
    clustcorrwithsuc = nan(1,14);
        nonclustcorrwithsuc = nan(1,14);
    succlustcorrwithsuc = nan(1,14);
        sucnonclustcorrwithsuc = nan(1,14);
        
    clustcorrwithrew = nan(1,14);
        nonclustcorrwithrew = nan(1,14);
    rewclustcorrwithrew = nan(1,14);
        rewnonclustcorrwithrew = nan(1,14);
        
    Caus_clustcorrwithcue = nan(1,14);
        Caus_nonclustcorrwithcue = nan(1,14);
    Caus_cueclustcorrwithcue = nan(1,14);
        Caus_cuenonclustcorrwithcue = nan(1,14);
        
    Caus_clustcorrwithMDC = nan(1,14);
        Caus_nonclustcorrwithMDC = nan(1,14);
    Caus_MDCclustcorrwithMDC = nan(1,14);
        Caus_MDCnonclustcorrwithMDC = nan(1,14);

        
    Caus_clustcorrwithmov = nan(1,14);
        Caus_nonclustcorrwithmov = nan(1,14);
    Caus_movclustcorrwithmov = nan(1,14);
        Caus_movnonclustcorrwithmov = nan(1,14);
        
    Caus_clustcorrwithsuc = nan(1,14);
        Caus_nonclustcorrwithsuc = nan(1,14);
    Caus_succlustcorrwithsuc = nan(1,14);
        Caus_sucnonclustcorrwithsuc = nan(1,14);
        
    Caus_clustcorrwithrew = nan(1,14);
        Caus_nonclustcorrwithrew = nan(1,14);
    Caus_rewclustcorrwithrew = nan(1,14);
        Caus_rewnonclustcorrwithrew = nan(1,14);
        
    FractionofCueSpinesThatAreClustered = nan(1,14);
    FractionofCueSpinesThatAreNonClustered = nan(1,14);
    FractionofMovementSpinesThatAreClustered = nan(1,14);
    FractionofMovementSpinesThatAreNonClustered = nan(1,14);
    FractionofPreSuccessSpinesThatAreClustered = nan(1,14);
    FractionofSuccessSpinesThatAreClustered = nan(1,14);
    FractionofSuccessSpinesThatAreNonClustered = nan(1,14);
    FractionofMovementDuringCueSpinesThatAreClustered = nan(1,14);
    FractionofRewardSpinesThatAreClustered = nan(1,14);
    FractionofRewardSpinesThatAreNonClustered = nan(1,14);
    
    cluster_freq = nan(1,14);
    nonclustered_freq = nan(1,14);
    cue_cluster_freq = nan(1,14);
    mov_cluster_freq = nan(1,14);
    movduringcue_cluster_freq = nan(1,14);
    presuc_cluster_freq = nan(1,14);
    suc_cluster_freq = nan(1,14);
    rew_cluster_freq = nan(1,14);
    Caus_cluster_freq = nan(1,14);
    Caus_nonclustered_freq = nan(1,14);
    Caus_cue_cluster_freq = nan(1,14);
    Caus_mov_cluster_freq = nan(1,14);
    Caus_movduringcue_cluster_freq = nan(1,14);
    Caus_presuc_cluster_freq = nan(1,14);
    Caus_suc_cluster_freq = nan(1,14);
    Caus_rew_cluster_freq = nan(1,14);
    
    ClustDendFreq = nan(1,14);
    NoClustDendFreq = nan(1,14);
    CueClustDendFreq = nan(1,14);
    MovClustDendFreq = nan(1,14);
    MovDuringCueClustDendFreq = nan(1,14);
    PreSucClustDendFreq = nan(1,14);
    SucClustDendFreq = nan(1,14);
    RewClustDendFreq = nan(1,14);
    NoMovClustDendFreq = nan(1,14);
    DendNumClust = cell(1,14);
    
    
    ClustAmp = nan(1,14);
    NonClusteredAmp = nan(1,14);
    Caus_ClustAmp = nan(1,14);
    Caus_NonClusteredAmp = nan(1,14);
    CueAmp = nan(1,14);
    Caus_CueAmp = nan(1,14);
    MovAmp = nan(1,14);
    Caus_MovAmp = nan(1,14);
    MovDuringCueAmp = nan(1,14);
    Caus_MovDuringCueAmp = nan(1,14);
    PreSucAmp = nan(1,14);
    Caus_PreSucAmp = nan(1,14);
    SucAmp = nan(1,14);
    Caus_SucAmp = nan(1,14);
    RewAmp = nan(1,14);
    Caus_RewAmp = nan(1,14);
    
    NumCueRelSpines = nan(1,14);
    NumMovRelSpines = nan(1,14);
    FractionofMovRelSpinesPerDendrite = arrayfun(@(x) x(:), nan(1,14), 'uni', false);
    NumCueORMovRelSpines = nan(1,14);
    NumPreSucRelSpines = nan(1,14);
    NumSucRelSpines = nan(1,14);
    NumMovDuringCueRelSpines = nan(1,14);
    NumRewRelSpines = nan(1,14);
    NumCausalCueSpines = nan(1,14);
    NumCausalMovSpines = nan(1,14);
    NumCausalSucSpines = nan(1,14);
    NumCausalRewSpines = nan(1,14);
    
    NumClustSpines = nan(1,14);
    NumClustCueSpines = nan(1,14);
    NumClustMovSpines = nan(1,14);
    NumClustMixSpines = nan(1,14);
    NumClustPreSucSpines = nan(1,14);
    NumClustSucSpines = nan(1,14);
    NumClustMovDuringCueSpines = nan(1,14);
    NumClustRewSpines = nan(1,14);
    
    NumFarClustSpines = nan(1,14);
    NumFarClustCueSpines = nan(1,14);
    NumFarClustMovSpines = nan(1,14);
    NumFarClustMixSpines = nan(1,14);
    NumFarClustPreSucSpines = nan(1,14);
    NumFarClustSucSpines = nan(1,14);
    NumFarClustMovDuringCueSpines = nan(1,14);
    NumFarClustRewSpines = nan(1,14);

    NumCausClustSpines = nan(1,14);
    NumCausClustCueSpines = nan(1,14);
    NumCausClustMovSpines = nan(1,14);
    NumCausClustMixSpines = nan(1,14);
    NumCausClustPreSucSpines = nan(1,14);
    NumCausClustSucSpines = nan(1,14);
    NumCausClustMovDuringCueSpines = nan(1,14);
    NumCausClustRewSpines = nan(1,14);
    

    MeanCueClustLength = nan(1,14);
    MaxCueClustLength = nan(1,14);
    AllMovClustLengths = cell(1,14);
    MeanMovClustLength = nan(1,14);
    MaxMovClustLength = nan(1,14);
    MeanMixClustLength = nan(1,14);
    MaxMixClustLength = nan(1,14);
    MeanPreSucClustLength = nan(1,14);
    MaxPreSucClustLength = nan(1,14);
    MeanSucClustLength = nan(1,14);
    MaxSucClustLength = nan(1,14);
    MeanMovDuringCueClustLength = nan(1,14);
    MaxMovDuringCueClustLength = nan(1,14);
    MeanRewClustLength = nan(1,14);
    MaxRewClustLength = nan(1,14);
    MeanAllClustLength = nan(1,14);
    MaxAllClustLength = nan(1,14);
    
    MeanFarCueClustLength = nan(1,14);
    MeanFarMovClustLength = nan(1,14);
    MeanFarMixClustLength = nan(1,14);
    MeanFarPreSucClustLength = nan(1,14);
    MeanFarSucClustLength = nan(1,14);
    MeanFarMovDuringCueClustLength = nan(1,14);
    MeanFarRewClustLength = nan(1,14);
    MeanAllFarClustLength = nan(1,14);

    MeanCausalCueClustLength = nan(1,14);
    MaxCausalCueClustLength = nan(1,14);
    MeanCausalMovClustLength = nan(1,14);
    MaxCausalMovClustLength = nan(1,14);
    MeanCausalSucClustLength = nan(1,14);
    MaxCausalSucClustLength = nan(1,14);
    MeanCausalRewClustLength = nan(1,14);
    MaxCausalRewClustLength = nan(1,14);
    MeanAllCausalClustLength = nan(1,14);
    MaxAllCausalClustLength = nan(1,14);
    
    NearestMovSpine = cell(1,14);
    NextClosest = cell(1,14);
    ThirdClosest = cell(1,14);
    FourthClosest = cell(1,14);
    CorrwithNearestMovSpine = cell(1,14);
    NearestHighCorrMovSpine = cell(1,14);
    NextClosestHighCorrMovSpine = cell(1,14);
    ThirdClosestHighCorrMovSpine = cell(1,14);
    
    DistanceBetweenAllSpines = cell(1,14);
        CorrelationBetweenAllSpines = cell(1,14);
        CorrelationBetweenAllSpinesMovePeriods = cell(1,14);
        CorrelationBetweenAllSpinesStillPeriods = cell(1,14);
        MeanCorrelationBetweenAllSpines = cell(1,14);
    DistanceBetweenCueSpines = cell(1,14);
    MeanDistanceBetweenCueSpines = nan(1,14);
    DistanceBetweenMovementSpines = cell(1,14);
    MeanDistanceBetweenMovementSpines = nan(1,14);
        CorrelationBetweenMovementSpines = cell(1,14);
        CorrelationBetweenMovementSpinesMovePeriods = cell(1,14);
        CorrelationBetweenMovementSpinesStillPeriods = cell(1,14);
        MeanCorrelationBetweenMovementSpines = cell(1,14);
    DistanceBetweenPreSuccessSpines = cell(1,14);
    MeanDistanceBetweenPreSuccessSpines = nan(1,14);
    DistanceBetweenSuccessSpines = cell(1,14);
    MeanDistanceBetweenSuccessSpines = nan(1,14);
    DistanceBetweenMovementDuringCueSpines = cell(1,14);
    MeanDistanceBetweenMovementDuringCueSpines = nan(1,14);
    DistanceBetweenRewardSpines = cell(1,14);
    MeanDistanceBetweenRewardSpines = nan(1,14);
    CorrelationBetweenFarSpines = cell(1,14);
    DistanceBetweenFarSpines = cell(1,14);
    CorrelationBetweenFarMovementSpines = cell(1,14);
    DistanceBetweenFarMovementSpines = cell(1,14);
    
    ClusteredSpines_CorrwithDend = nan(1,14);
    FilteredClusteredSpines_CorrwithDend = nan(1,14);
    NonClusteredSpines_CorrwithDend = nan(1,14);
    CausalClusteredSpines_CorrwithDend = nan(1,14);
    FilteredCausalClusteredSpines_CorrwithDend = nan(1,14);
    NonCausalClusteredSpines_CorrwithDend = nan(1,14);
    CueRelClusteredSpines_CorrwithDend = nan(1,14);
    MovRelClusteredSpines_CorrwithDend = nan(1,14);
    PreSucRelClusteredSpines_CorrwithDend = nan(1,14);
    SucRelClusteredSpines_CorrwithDend = nan(1,14);
    MovDuringCueRelClusteredSpines_CorrwithDend = nan(1,14);
    RewRelClusteredSpines_CorrwithDend = nan(1,14);
    CueRelCausalClusteredSpines_CorrwithDend = nan(1,14);
    MovRelCausalClusteredSpines_CorrwithDend = nan(1,14);
    SucRelCausalClusteredSpines_CorrwithDend = nan(1,14);
    RewRelCausalClusteredSpines_CorrwithDend = nan(1,14);

    Spatial_Deg = cell(1,14);
    Temporal_Deg = cell(1,14);
    Spatiotemporal_Deg = cell(1,14);
    Dend_SpatTemp_Deg = cell(1,14);
    Spatiotemporal_Laplacian = cell(1,14);
    Spatiotemporal_Overlap = cell(1,14);
    MovementCorrelationsforAllSpinesonDend = cell(1,14);
    SpatialDegree_vs_Movement = cell(1,14);
    TemporalDegree_vs_Movement = cell(1,14);
    SpatiotemporalDegree_vs_Movement = cell(1,14);
    DendClust_Deg = cell(1,14);
    MeanSpatialDegreeofCueSpines = nan(1,14);
    MeanTemporalDegreeofCueSpines = nan(1,14);
    MeanSpatioTemporalDegreeofCueSpines = nan(1,14);
    MeanSpatialDegreeofMovementSpines = nan(1,14);
    MeanTemporalDegreeofMovementSpines = nan(1,14);
    MeanSpatioTemporalDegreeofMovementSpines = nan(1,14);
    MeanSpatialDegreeofMDCSpines = nan(1,14);
    MeanTemporalDegreeofMDCSpines = nan(1,14);
    MeanSpatioTemporalDegreeofMDCSpines = nan(1,14);
    MeanSpatialDegreeofPreSuccessSpines = nan(1,14);
    MeanTemporalDegreeofPreSuccessSpines = nan(1,14);
    MeanSpatioTemporalDegreeofPreSuccessSpines = nan(1,14);
    MeanSpatialDegreeofSuccessSpines = nan(1,14);
    MeanTemporalDegreeofSuccessSpines = nan(1,14);
    MeanSpatioTemporalDegreeofSuccessSpines = nan(1,14);
    MeanSpatialDegreeofRewardSpines = nan(1,14);
    MeanTemporalDegreeofRewardSpines = nan(1,14);
    MeanSpatioTemporalDegreeofRewardSpines = nan(1,14);
    
    NumClusters = nan(1,14);
    NumCausalClusters = nan(1,14);
    MeanNumberofSpinesinEachCluster = nan(1,14);
    MeanNumberofSpinesinEachCausalCluster = nan(1,14);
    NumMovClusters = nan(1,14);
    MeanNumberofSpinesinEachMovCluster = nan(1,14);
    
    FractionofClusterThatsMR = nan(1,14);
    FractionofCausalClusterThatsMR = nan(1,14);
    SpinesonCueDend = cell(1,14);
    SpinesonMovementDend = cell(1,14);
    SpinesonPreSuccessDend = cell(1,14);
    SpinesonSuccessDend = cell(1,14);
    SpinesonMovementDuringCueDend = cell(1,14);
    SpinesonRewardDend = cell(1,14);
    PercentCueRelDends = nan(1,14);
    PercentMovRelDends = nan(1,14);
    PercentPreSucRelDends = nan(1,14);
    PercentSucRelDends = nan(1,14);
    PercentMovDuringCueRelDends = nan(1,14);
    PercentRewRelDends = nan(1,14);
    Temporal_Laplacian = cell(1,14);
    Spatial_Deg = cell(1,14);
        Dend_Spat_Deg = cell(1,14);
        Dend_Temp_Deg = cell(1,14);
        Spatial_FirstEigenvector = cell(1,14);
    Temporal_Deg = cell(1,14);
        end_Temp_Deg = cell(1,14);
        Temporal_FirstEigenvector = cell(1,14);
    SpatioTemporalFiedler = cell(1,14);
    SpatioTemporalPartition = cell(1,14);
    SpatioTemporal_FirstEigenvector = cell(1,14);
    
    for i = 1:length(StatClass)
        if ~isempty(StatClass{i})
            StatClass{i}.CueORMovementSpines = logical(StatClass{i}.MovementSpines+StatClass{i}.CueSpines);
            StatClass{i}.DendSub_CueORMovementSpines = logical(StatClass{i}.DendSub_MovementSpines+StatClass{i}.DendSub_CueSpines);
        else
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%
    
    for i = firstdatainput:length(varargin) 
        
        session = varargin{i}.Session;
        SpineToSpineDistance = varargin{i}.DistanceHeatMap;
        tempC = SpineToSpineDistance;
        tempD = SpineToSpineDistance';
        tempC(isnan(tempC) & ~isnan(tempD)) = tempD(isnan(tempC) & ~isnan(tempD));
        DistanceMap = tempC;
        NumClusters(1,session) = 0;
        NumCausalClusters(1,session) = 0;
        
        if varargin{i}.NumberofSpines ~= length(varargin{i}.dF_over_F)
            varargin{i}.NumberofSpines = length(varargin{i}.dF_over_F);
        end
       
        if ~isempty(Correlations{session})
           
            if isfield(varargin{i}, 'SpineDendriteGrouping')
                for d = 1:varargin{i}.NumberofDendrites
                    firstspine = varargin{i}.SpineDendriteGrouping{d}(1);
                    lastspine = varargin{i}.SpineDendriteGrouping{d}(end);
                    SpinesonCueDend{session}(firstspine:lastspine) = ones(length(firstspine:lastspine),1)*StatClass{session}.CueDends(d);
                    SpinesonMovementDend{session}(firstspine:lastspine) = ones(length(firstspine:lastspine),1)*StatClass{session}.MovementDends(d);
                    SpinesonPreSuccessDend{session}(firstspine:lastspine) = ones(length(firstspine:lastspine),1)*StatClass{session}.PreSuccessDends(d);
                    SpinesonSuccessDend{session}(firstspine:lastspine) = ones(length(firstspine:lastspine),1)*StatClass{session}.SuccessDends(d);
                    SpinesonMovementDuringCueDend{session}(firstspine:lastspine) = ones(length(firstspine:lastspine),1)*StatClass{session}.MovementDuringCueDends(d);
                    SpinesonRewardDend{session}(firstspine:lastspine) = ones(length(firstspine:lastspine),1)*StatClass{session}.RewardDends(d);
                end
            else
                if varargin{i}.NumberofDendrites ==1
                    varargin{i}.SpineDendriteGrouping = {1:varargin{i}.NumberofSpines};
                    SpinesonCueDend{session} = ones(varargin{i}.NumberofSpines,1).*StatClass{session}.CueDends;
                    SpinesonMovementDend{session} = ones(varargin{i}.NumberofSpines,1).*StatClass{session}.MovementDends;
                    SpinesonPreSuccessDend{session} = ones(varargin{i}.NumberofSpines,1).*StatClass{session}.PreSuccessDends;
                    SpinesonSuccessDend{session} = ones(varargin{i}.NumberofSpines,1).*StatClass{session}.SuccessDends;
                    SpinesonMovementDuringCueDend{session} = ones(varargin{i}.NumberofSpines,1).*StatClass{session}.MovementDuringCueDends;
                    SpinesonRewardDend{session} = ones(varargin{i}.NumberofSpines,1).*StatClass{session}.RewardDends;
                    disp(['File ', varargin{i}.Filename, ' did not have spine-dendrite grouping information, but can be correctly filled in']);
                else
                    disp(['File ', varargin{i}.Filename, ' does NOT have necessary spine-dendrite grouping info!']);
                end
            end
            
            PercentCueRelDends(1,session) = sum(StatClass{session}.CueDends)/varargin{i}.NumberofDendrites;
            PercentMovRelDends(1,session) = sum(StatClass{session}.MovementDends)/varargin{i}.NumberofDendrites;
            PercentPreSucRelDends(1,session) = sum(StatClass{session}.PreSuccessDends)/varargin{i}.NumberofDendrites;
            PercentSucRelDends(1,session) = sum(StatClass{session}.SuccessDends)/varargin{i}.NumberofDendrites;
            PercentMovDuringCueRelDends(1,session) = sum(StatClass{session}.MovementDuringCueDends)/varargin{i}.NumberofDendrites;
            PercentRewRelDends(1,session) = sum(StatClass{session}.RewardDends)/varargin{i}.NumberofDendrites;
            
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%% Data choice (how to handle bAPs) %%%%%%%%%%%%%%%%
                       
            if dendexclude
                SessionCorrData = Correlations{session}.SpineCorrelations;
                SessionMovCorrData = Correlations{session}.SpineDuringMovePeriods;
                SessionStillCorrData = Correlations{session}.SpineDuringStillPeriods;
                CueSpines = StatClass{session}.CueSpines;
                MovementSpines = StatClass{session}.MovementSpines;
                MovementDuringCueSpines = StatClass{session}.MovementDuringCueSpines;
                PreSuccessSpines = StatClass{session}.PreSuccessSpines;
                SuccessSpines = StatClass{session}.SuccessSpines;
                RewardSpines = StatClass{session}.RewardSpines;
                CueORMovementSpines = StatClass{session}.CueORMovementSpines;
            elseif dendsubtract
                SessionCorrData = Correlations{session}.DendSubtractedSpineCorrelations;
                SessionMovCorrData = Correlations{session}.SpineDuringMovePeriods;
                SessionStillCorrData = Correlations{session}.SpineDuringStillPeriods;
                CueSpines = StatClass{session}.DendSub_CueSpines;
                MovementSpines = StatClass{session}.DendSub_MovementSpines;
                MovementDuringCueSpines = StatClass{session}.DendSub_MovementDuringCueSpines;
                PreSuccessSpines = StatClass{session}.DendSub_PreSuccessSpines;
                SuccessSpines = StatClass{session}.DendSub_SuccessSpines;
                RewardSpines = StatClass{session}.DendSub_RewardSpines;
                CueORMovementSpines = StatClass{session}.DendSub_CueORMovementSpines;
            end
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%% Spectral analysis of clustering %%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [SpectralData] = SpectralClustering(varargin{i}, SessionCorrData, StatClass, Choices);

            

            Spatial_Deg{session} = SpectralData.Spatial_Deg;
            Temporal_Deg{session} = SpectralData.Temporal_Deg;
            Spatiotemporal_Deg{session} = SpectralData.Spatiotemporal_Deg;
            Spatiotemporal_Overlap{session} = SpectralData.Spatiotemporal_Overlap;
            
            SpatioTemporalFiedler{session} = SpectralData.SpatioTemporalFiedler;
            
            Temporal_Laplacian{session} = SpectralData.Temporal_Laplacian;

            MeanSpatialDegreeofCueSpines(1,session) = SpectralData.MeanSpatialDegreeofCueSpines;
            MeanTemporalDegreeofCueSpines(1,session) = SpectralData.MeanTemporalDegreeofCueSpines;
            MeanSpatioTemporalDegreeofCueSpines(1,session) = SpectralData.MeanSpatioTemporalDegreeofCueSpines;
            MeanSpatialDegreeofMovementSpines(1,session) = SpectralData.MeanSpatialDegreeofMovementSpines;
            MeanTemporalDegreeofMovementSpines(1,session) = SpectralData.MeanTemporalDegreeofMovementSpines;
            MeanSpatioTemporalDegreeofMovementSpines(1,session) = SpectralData.MeanSpatioTemporalDegreeofMovementSpines;
            MeanSpatialDegreeofMDCSpines(1,session) = SpectralData.MeanSpatialDegreeofMovementDuringCueSpines;
            MeanTemporalDegreeofMDCSpines(1,session) = SpectralData.MeanTemporalDegreeofMovementDuringCueSpines;
            MeanSpatioTemporalDegreeofMDCSpines(1,session) = SpectralData.MeanSpatioTemporalDegreeofMovementDuringCueSpines;
            MeanSpatialDegreeofPreSuccessSpines(1,session) = SpectralData.MeanSpatialDegreeofPreSuccessSpines;
            MeanTemporalDegreeofPreSuccessSpines(1,session) = SpectralData.MeanTemporalDegreeofPreSuccessSpines;
            MeanSpatioTemporalDegreeofPreSuccessSpines(1,session) = SpectralData.MeanSpatioTemporalDegreeofPreSuccessSpines;
            MeanSpatialDegreeofSuccessSpines(1,session) = SpectralData.MeanSpatialDegreeofSuccessSpines;
            MeanTemporalDegreeofSuccessSpines(1,session) = SpectralData.MeanTemporalDegreeofSuccessSpines;
            MeanSpatioTemporalDegreeofSuccessSpines(1,session) = SpectralData.MeanSpatioTemporalDegreeofSuccessSpines;
            MeanSpatialDegreeofRewardSpines(1,session) = SpectralData.MeanSpatialDegreeofRewardSpines;
            MeanTemporalDegreeofRewardSpines(1,session) = SpectralData.MeanTemporalDegreeofRewardSpines;
            MeanSpatioTemporalDegreeofRewardSpines(1,session) = SpectralData.MeanSpatioTemporalDegreeofRewardSpines;
            
            Dend_Spat_Deg{session} = SpectralData.Dend_Spat_Deg;
            Dend_Temp_Deg{session} = SpectralData.Dend_Temp_Deg;
            Dend_SpatTemp_Deg{session} = SpectralData.Dend_SpatTemp_Deg;      
            
            SpatialDegree_vs_Movement{session} = SpectralData.SpatialDegree_vs_Movement;
            TemporalDegree_vs_Movement{session} = SpectralData.TemporalDegree_vs_Movement;
            SpatiotemporalDegree_vs_Movement{session} = SpectralData.SpatiotemporalDegree_vs_Movement;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%
            
            if useonlyspinesfromstatreldends
                dendcheckC = reshape(SpinesonCueDend{session},1,length(SpinesonCueDend{session}));
                dendcheckM = reshape(SpinesonMovementDend{session},1,length(SpinesonMovementDend{session}));
                dendcheckPS = reshape(SpinesonPreSuccessDend{session},1,length(SpinesonPreSuccessDend{session}));
                dendcheckS = reshape(SpinesonSuccessDend{session},1,length(SpinesonSuccessDend{session}));
                dendcheckMDC = reshape(SpinesonMovementDuringCueDend{session},1,length(SpinesonMovementDuringCueDend{session}));
                dendcheckR = reshape(SpinesonRewardDend{session},1,length(SpinesonRewardDend{session}));
                dendcheck = nansum([dendcheckC; dendcheckM; dendcheckPS; dendcheckS; dendcheckMDC; dendcheckR]);
                dendcheck(dendcheck~=0) = 1;
                dendcheck(dendcheck==0)= nan;
                dendcheckC(dendcheckC==0) = nan;
                dendcheckM(dendcheckM==0) = nan;
                dendcheckPS(dendcheckPS==0) = nan;
                dendcheckS(dendcheckS==0) = nan;
                dendcheckMDC(dendcheckMDC==0) = nan;
                dendcheckR(dendcheckR==0) = nan;
            else
                dendcheck = ones(1,length(varargin{i}.deltaF)); %%% Set this variable to one to allow use of ALL spines on ALL (non just stat-rel) dendrites
                dendcheckC = ones(1,length(varargin{i}.deltaF));
                dendcheckM = ones(1,length(varargin{i}.deltaF));
                dendcheckPS = ones(1,length(varargin{i}.deltaF));
                dendcheckS = ones(1,length(varargin{i}.deltaF));
                dendcheckMDC = ones(1,length(varargin{i}.deltaF));
                dendcheckR = ones(1,length(varargin{i}.deltaF));
            end
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%% Set up correlation matrices %%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                corrmat = SessionCorrData(Spine1_address+1:Spine1_address+length(varargin{i}.dF_over_F),Spine1_address+1:Spine1_address+length(varargin{i}.dF_over_F));
                movcorrmat = SessionMovCorrData;
                stillcorrmat = SessionStillCorrData;
                causalmat = Correlations{session}.CausalCorrelations(Spine1_address+1:Spine1_address+length(varargin{i}.dF_over_F),Spine1_address+1:Spine1_address+length(varargin{i}.dF_over_F));
                
                samedendcorrection = DistanceMap; samedendcorrection(~isnan(samedendcorrection)) = 1;
                separatedendcorrection = DistanceMap; separatedendcorrection(~isnan(separatedendcorrection))=0; separatedendcorrection(isnan(separatedendcorrection)) = 1; 
                    separatedendcorrection(separatedendcorrection==0) = nan; separatedendcorrection = separatedendcorrection.*~eye(size(separatedendcorrection));
                    separatedendcorrection = triu(separatedendcorrection); %%% The distance matrix of spines on different dendrites was converted to a list earlier in analysis, so just taking the upper half makes it easier to find these values
                    separatedendcorrection(separatedendcorrection==0) = nan;
                samedendcorrmat = corrmat.*samedendcorrection; %%% Make all values on the same dendrites == 1, while all others are 0;
                tempcorrmat = corrmat; tempcorrmat(isnan(tempcorrmat)) = 0 ; tempcorrmat = tempcorrmat.*separatedendcorrection;
                separatedendcorrmat = tempcorrmat.*separatedendcorrection; %%% Make vall values on SEPARATE dendrites ==1, while all on the same ==0;
                    separatedendcorrlist = separatedendcorrmat(find(~isnan(separatedendcorrmat)));
                samedendmovcorrmat = movcorrmat.*samedendcorrection;
                samedendstillcorrmat = stillcorrmat.*samedendcorrection;

            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%% Binary Classification of Clusters %%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            

                
                [ClusteredSpines, CorrelationofClusts, FarClusteredSpines, CausalClusteredSpines] = binaryclusterclass(varargin{i}, corrmat, causalmat);
                
            CorrelationofClusters{session} = CorrelationofClusts;
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Find indices for movement- and cue-related clusters
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%% Correlation indices for synapse-only events
            
               [ClusterInfo]...
                = findclusterind(ClusteredSpines, StatClass{session}, Spine1_address);
            
            NumSpinesinCluster = ClusterInfo.NumSpinesinCluster;
            cluster_ind = ClusterInfo.cluster_ind(~isnan(dendcheck(ClusterInfo.cluster_ind-Spine1_address)));
            cue_cluster_ind = ClusterInfo.cue_cluster_ind(~isnan(dendcheckC(ClusterInfo.cue_cluster_ind-Spine1_address)));
            mov_cluster_ind = ClusterInfo.mov_cluster_ind(~isnan(dendcheckM(ClusterInfo.mov_cluster_ind-Spine1_address)));
            mix_cluster_ind = ClusterInfo.mix_cluster_ind(~isnan(dendcheck(ClusterInfo.mix_cluster_ind-Spine1_address)));
            presuc_cluster_ind = ClusterInfo.presuc_cluster_ind(~isnan(dendcheckPS(ClusterInfo.presuc_cluster_ind-Spine1_address)));
            suc_cluster_ind = ClusterInfo.suc_cluster_ind(~isnan(dendcheckS(ClusterInfo.suc_cluster_ind-Spine1_address)));
            movduringcue_cluster_ind = ClusterInfo.movduringcue_cluster_ind(~isnan(dendcheckMDC(ClusterInfo.movduringcue_cluster_ind-Spine1_address)));
            rew_cluster_ind = ClusterInfo.rew_cluster_ind(~isnan(dendcheckR(ClusterInfo.rew_cluster_ind-Spine1_address)));
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Repeat the above for "clusters" on separate dends.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%% Correlation indices for synapse-only events
            
               [ClusterInfo]...
                = findclusterind(FarClusteredSpines, StatClass{session}, Spine1_address);
            
            NumSpinesinFarCluster = ClusterInfo.NumSpinesinCluster;
            Farcluster_ind = ClusterInfo.cluster_ind(~isnan(dendcheck(ClusterInfo.cluster_ind-Spine1_address)));
            Farcue_cluster_ind = ClusterInfo.cue_cluster_ind(~isnan(dendcheckC(ClusterInfo.cue_cluster_ind-Spine1_address)));
            Farmov_cluster_ind = ClusterInfo.mov_cluster_ind(~isnan(dendcheckM(ClusterInfo.mov_cluster_ind-Spine1_address)));
            Farmix_cluster_ind = ClusterInfo.mix_cluster_ind(~isnan(dendcheck(ClusterInfo.mix_cluster_ind-Spine1_address)));
            Farpresuc_cluster_ind = ClusterInfo.presuc_cluster_ind(~isnan(dendcheckPS(ClusterInfo.presuc_cluster_ind-Spine1_address)));
            Farsuc_cluster_ind = ClusterInfo.suc_cluster_ind(~isnan(dendcheckS(ClusterInfo.suc_cluster_ind-Spine1_address)));
            Farmovduringcue_cluster_ind = ClusterInfo.movduringcue_cluster_ind(~isnan(dendcheckMDC(ClusterInfo.movduringcue_cluster_ind-Spine1_address)));
            Farrew_cluster_ind = ClusterInfo.rew_cluster_ind(~isnan(dendcheckR(ClusterInfo.rew_cluster_ind-Spine1_address)));
            
            %%   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            %%%% Repeat the above for causal correlations
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%% Uncomment the following section to separate cue- and
            %%%% movement-correlated clusters
            
               [ClusterInfo]...
                = findclusterind(CausalClusteredSpines, StatClass{session}, Spine1_address);
            
            NumSpinesinCausalCluster = ClusterInfo.NumSpinesinCluster;
            Caus_cluster_ind = ClusterInfo.cluster_ind(~isnan(dendcheck(ClusterInfo.cluster_ind-Spine1_address)));
            Caus_cue_cluster_ind = ClusterInfo.cue_cluster_ind(~isnan(dendcheckC(ClusterInfo.cue_cluster_ind-Spine1_address)));
            Caus_mov_cluster_ind = ClusterInfo.mov_cluster_ind(~isnan(dendcheckM(ClusterInfo.mov_cluster_ind-Spine1_address)));
            Caus_mix_cluster_ind = ClusterInfo.mix_cluster_ind(~isnan(dendcheck(ClusterInfo.mix_cluster_ind-Spine1_address)));
            Caus_presuc_cluster_ind = ClusterInfo.presuc_cluster_ind(~isnan(dendcheckPS(ClusterInfo.presuc_cluster_ind-Spine1_address)));
            Caus_suc_cluster_ind = ClusterInfo.suc_cluster_ind(~isnan(dendcheckS(ClusterInfo.suc_cluster_ind-Spine1_address)));
            Caus_movduringcue_cluster_ind = ClusterInfo.movduringcue_cluster_ind(~isnan(dendcheckMDC(ClusterInfo.movduringcue_cluster_ind-Spine1_address)));
            Caus_rew_cluster_ind = ClusterInfo.rew_cluster_ind(~isnan(dendcheckR(ClusterInfo.rew_cluster_ind-Spine1_address)));

                     
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%% End index-finding section %%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Main variables
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ClusteredSpines = ClusteredSpines(~cellfun(@isempty,ClusteredSpines));
            CausalClusteredSpines = CausalClusteredSpines(~cellfun(@isempty,CausalClusteredSpines));

            NumClusters(1,session) = length(ClusteredSpines);
            NumCausalClusters(1,session) = length(CausalClusteredSpines);
           
            NumSpinesinCluster(NumSpinesinCluster==0) = NaN;
            MeanNumberofSpinesinEachCluster(1,session) = nanmean(NumSpinesinCluster);
            NumSpinesinCausalCluster(NumSpinesinCausalCluster==0) = NaN;
            MeanNumberofSpinesinEachCausalCluster(1,session) = nanmean(NumSpinesinCausalCluster);
            ClusterswithMovRelSpines = ClusteredSpines(cellfun(@(x) x(:), cellfun(@(x) sum(x)>1, cellfun(@(x) ismember(x,mov_cluster_ind-Spine1_address), ClusteredSpines, 'UniformOutput', false), 'Uniformoutput', false)));
                        
            if ~isempty(ClusteredSpines)
                MeanNumberofSpinesinEachMovCluster(1,session) = mean(cell2mat(cellfun(@length, cellfun(@(x,y) x(y), ClusterswithMovRelSpines, cellfun(@(x) ismember(x,mov_cluster_ind-Spine1_address), ClusterswithMovRelSpines, 'Uni', false), 'Uni', false), 'Uni', false)));
                NumMovClusters(1,session) = sum(cell2mat(cellfun(@(x) sum(x)>1, cellfun(@(x) ismember(x,mov_cluster_ind-Spine1_address), ClusteredSpines, 'UniformOutput', false), 'Uniformoutput', false)));  %%% Check if each reported cluster is movement-related
            else
                MeanNumberofSpinesinEachMovCluster(1,session) = NaN;
                NumMovClusters(1,session) = NaN;
            end            
            
            FractionofCueSpinesThatAreClustered(1,session) = length(cue_cluster_ind)/length(find(CueSpines));
            FractionofCueSpinesThatAreNonClustered(1,session) = length((setdiff(find(CueSpines),cue_cluster_ind-Spine1_address)))/length(find(CueSpines));
            FractionofMovementSpinesThatAreClustered(1,session) = length(mov_cluster_ind)/length(find(MovementSpines));
            FractionofMovementSpinesThatAreNonClustered(1,session) = length((setdiff(find(MovementSpines),mov_cluster_ind-Spine1_address)))/length(find(MovementSpines));
            FractionofPreSuccessSpinesThatAreClustered(1,session) = length(presuc_cluster_ind)/length(find(PreSuccessSpines));
            FractionofSuccessSpinesThatAreClustered(1,session) = length(suc_cluster_ind)/length(find(SuccessSpines));
            FractionofSuccessSpinesThatAreNonClustered(1,session) = length((setdiff(find(SuccessSpines),suc_cluster_ind-Spine1_address)))/length(find(SuccessSpines));
            FractionofMovementDuringCueSpinesThatAreClustered(1,session) = length(movduringcue_cluster_ind)/length(find(MovementDuringCueSpines));
            FractionofRewardSpinesThatAreClustered(1,session) = length(rew_cluster_ind)/length(find(RewardSpines));
            FractionofRewardSpinesThatAreNonClustered(1,session) = length((setdiff(find(RewardSpines),rew_cluster_ind-Spine1_address)))/length(find(RewardSpines));
            
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Find number of clustered spines in each category
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            
            NumCueRelSpines(1,session) = sum(CueSpines)/varargin{i}.NumberofSpines;
            NumMovRelSpines(1,session) = sum(MovementSpines)/varargin{i}.NumberofSpines;
                for d = 1:length(varargin{i}.SpineDendriteGrouping)
                    FractionofMovRelSpinesPerDendrite{session}(1,d) = sum(MovementSpines(varargin{i}.SpineDendriteGrouping{d}))/length(varargin{i}.SpineDendriteGrouping{d});
                end
            NumCueORMovRelSpines(1,session) = sum(CueORMovementSpines)/varargin{i}.NumberofSpines;
            NumPreSucRelSpines(1,session) = sum(PreSuccessSpines)/varargin{i}.NumberofSpines;
            NumSucRelSpines(1,session) = sum(SuccessSpines)/varargin{i}.NumberofSpines;
            NumMovDuringCueRelSpines(1,session) = sum(MovementDuringCueSpines)/varargin{i}.NumberofSpines;
            NumRewRelSpines(1,session) = sum(RewardSpines)/varargin{i}.NumberofSpines;
            NumCausalMovSpines(1,session) = sum(StatClass{session}.CausalMovementSpines)/varargin{i}.NumberofSpines;
            
            NumClustSpines(1,session) = length(cluster_ind)/varargin{i}.NumberofSpines;
            NumClustCueSpines(1,session) = length(cue_cluster_ind)/varargin{i}.NumberofSpines;
            NumClustMovSpines(1,session) = length(mov_cluster_ind)/varargin{i}.NumberofSpines;
            NumClustMixSpines(1,session) = length(mix_cluster_ind)/varargin{i}.NumberofSpines;
            NumClustPreSucSpines(1,session) = length(presuc_cluster_ind)/varargin{i}.NumberofSpines;
            NumClustSucSpines(1,session) = length(suc_cluster_ind)/varargin{i}.NumberofSpines;
            NumClustMovDuringCueSpines(1,session) = length(movduringcue_cluster_ind)/varargin{i}.NumberofSpines;
            NumClustRewSpines(1,session) = length(rew_cluster_ind)/varargin{i}.NumberofSpines;
            
            NumFarClustSpines(1,session) = length(Farcluster_ind)/varargin{i}.NumberofSpines;
            NumFarClustCueSpines(1,session) = length(Farcue_cluster_ind)/varargin{i}.NumberofSpines;
            NumFarClustMovSpines(1,session) = length(Farmov_cluster_ind)/varargin{i}.NumberofSpines;
            NumFarClustMixSpines(1,session) = length(Farmix_cluster_ind)/varargin{i}.NumberofSpines;
            NumFarClustPreSucSpines(1,session) = length(Farpresuc_cluster_ind)/varargin{i}.NumberofSpines;
            NumFarClustSucSpines(1,session) = length(Farsuc_cluster_ind)/varargin{i}.NumberofSpines;
            NumFarClustMovDuringCueSpines(1,session) = length(Farmovduringcue_cluster_ind)/varargin{i}.NumberofSpines;
            NumFarClustRewSpines(1,session) = length(Farrew_cluster_ind)/varargin{i}.NumberofSpines;

            NumCausClustSpines(1,session) = length(Caus_cluster_ind)/varargin{i}.NumberofSpines;
            NumCausClustCueSpines(1,session) = length(Caus_cue_cluster_ind)/varargin{i}.NumberofSpines;
            NumCausClustMovSpines(1,session) = length(Caus_mov_cluster_ind)/varargin{i}.NumberofSpines;
            NumCausClustMixSpines(1,session) = length(Caus_mix_cluster_ind)/varargin{i}.NumberofSpines;
            NumCausClustPreSucSpines(1,session) = length(Caus_presuc_cluster_ind)/varargin{i}.NumberofSpines;
            NumCausClustSucSpines(1,session) = length(Caus_suc_cluster_ind)/varargin{i}.NumberofSpines;
            NumCausClustMovDuringCueSpines(1,session) = length(Caus_movduringcue_cluster_ind)/varargin{i}.NumberofSpines;
            NumCausClustRewSpines(1,session) = length(Caus_rew_cluster_ind)/varargin{i}.NumberofSpines;
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Distances between clustered spines
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            val = 1;
            AllClusterDistances = cell(1,length(ClusteredSpines));
            CueClusterDistances = [];
            MovClusterDistances = [];
            MixClusterDistances = [];
            PreSucClusterDistances = [];
            SucClusterDistances = [];
            MovDuringCueClusterDistances = [];
            RewClusterDistances = [];
            
            %%% regular clusters
            
            for CS = 1:length(ClusteredSpines)
                    
                    spinecombos = nchoosek(ClusteredSpines{CS}, 2); %%% Find all possible combinations of spines in a cluster
                    all_clust_distances = zeros(1,size(spinecombos,1));
                    cue_clust_distances = [];
                    mov_clust_distances = [];
                    mix_clust_distances = [];
                    presuc_clust_distances = [];
                    suc_clust_distances = [];
                    movduringcue_clust_distances = [];
                    rew_clust_distances = [];
                    
                    for n = 1:size(spinecombos,1)
                        all_clust_distances(1,n) = DistanceMap(spinecombos(n,1), spinecombos(n,2));
                        if sum(CueSpines(spinecombos(n,1:2))) == 2;
                            cue_clust_distances(1,n) = DistanceMap(spinecombos(n,1), spinecombos(n,2));
                        else
                            cue_clust_distances(1,n) = nan;
                        end
                        if sum(MovementSpines(spinecombos(n,1:2))) == 2;
                            mov_clust_distances(1,n) = DistanceMap(spinecombos(n,1), spinecombos(n,2));
                        else
                            mov_clust_distances(1,n) = nan;
                        end
                        if sum([sum(CueSpines(spinecombos(n,1:2))), sum(MovementSpines(spinecombos(n,1:2)))]) >= 2
                            mix_clust_distances(1,n) = DistanceMap(spinecombos(n,1), spinecombos(n,2));
                        else
                            mix_clust_distances(1,n) = nan;
                        end
                        if sum(PreSuccessSpines(spinecombos(n,1:2))) == 2
                            presuc_clust_distances(1,n) = DistanceMap(spinecombos(n,1),spinecombos(n,2));
                        else
                            presuc_clust_distances(1,n) = nan;
                        end
                        if sum(SuccessSpines(spinecombos(n,1:2))) == 2;
                            suc_clust_distances(1,n) = DistanceMap(spinecombos(n,1), spinecombos(n,2));
                        else
                            suc_clust_distances(1,n) = nan;
                        end
                        if sum(MovementDuringCueSpines(spinecombos(n,1:2))) == 2;
                            movduringcue_clust_distances(1,n) = DistanceMap(spinecombos(n,1), spinecombos(n,2));
                        else
                            movduringcue_clust_distances(1,n) = nan;
                        end
                        if sum(RewardSpines(spinecombos(n,1:2))) == 2;
                            rew_clust_distances(1,n) = DistanceMap(spinecombos(n,1), spinecombos(n,2));
                        else
                            rew_clust_distances(1,n) = nan;
                        end
                    end
                AllClusterDistances{val} = all_clust_distances;
                CueClusterDistances{val} = cue_clust_distances;
                MovClusterDistances{val} = mov_clust_distances;
                MixClusterDistances{val} = mix_clust_distances;
                PreSucClusterDistances{val} = presuc_clust_distances;
                SucClusterDistances{val} = suc_clust_distances;
                MovDuringCueClusterDistances{val} = movduringcue_clust_distances;
                RewClusterDistances{val} = rew_clust_distances;
                val = val+1;
            end
            
            
            MeanCueClustLength(1,session) = nanmean(cell2mat(CueClusterDistances));
            if ~isempty(cell2mat(CueClusterDistances))
                MaxCueClustLength(1,session) = nanmax(cell2mat(CueClusterDistances));
            else
                MaxCueClustLength(1,session) = NaN;
            end
            AllMovClustLengths{session} = cell2mat(MovClusterDistances);
            MeanMovClustLength(1,session) = nanmean(cell2mat(MovClusterDistances));
            if ~isempty(cell2mat(MovClusterDistances))
                MaxMovClustLength(1,session) = nanmax(cell2mat(MovClusterDistances));
            else
                MaxMovClustLength(1,session) = NaN;
            end
            MeanMixClustLength(1,session) = nanmean(cell2mat(MixClusterDistances));
            if ~isempty(cell2mat(MixClusterDistances))
                MaxMixClustLength(1,session) = nanmean(cell2mat(MixClusterDistances));
            else
                MaxMixClustLength(1,session) = NaN;
            end
            MeanPreSucClustLength(1,session) = nanmean(cell2mat(PreSucClusterDistances));
            if ~isempty(cell2mat(PreSucClusterDistances))
                MaxPreSucClustLength(1,session) = nanmax(cell2mat(PreSucClusterDistances));
            else
                MaxPreSucClustLength(1,session) = NaN;
            end
            MeanSucClustLength(1,session) = nanmean(cell2mat(SucClusterDistances));
            if ~isempty(cell2mat(SucClusterDistances))
                MaxSucClustLength(1,session) = nanmax(cell2mat(SucClusterDistances));
            else
                MaxSucClustLength(1,session) = NaN;
            end
            MeanMovDuringCueClustLength(1,session) = nanmean(cell2mat(MovDuringCueClusterDistances));
            if ~isempty(cell2mat(MovDuringCueClusterDistances))
                MaxMovDuringCueClustLength(1,session) = nanmax(cell2mat(MovDuringCueClusterDistances));
            else
                MaxMovDuringCueClustLength(1,session) = NaN;
            end
            MeanRewClustLength(1,session) = nanmean(cell2mat(RewClusterDistances));
            if ~isempty(cell2mat(RewClusterDistances))
                MaxRewClustLength(1,session) = nanmax(cell2mat(RewClusterDistances));
            else
                MaxRewClustLength(1,session) = NaN;
            end
            MeanAllClustLength(1,session) = nanmean(cell2mat(AllClusterDistances));
            if ~isempty(cell2mat(AllClusterDistances))
                MaxAllClustLength(1,session) = nanmax(cell2mat(AllClusterDistances));
            else
                MaxAllClustLength(1,session) = NaN;
            end
            
            %%% Far clusters
            
            val = 1;
            AllFarClusterDistances = [];
            FarCueClusterDistances = [];
            FarMovClusterDistances = [];
            FarMixClusterDistances = [];
            FarPreSucClusterDistances = [];
            FarSucClusterDistances = [];
            FarMovDuringCueClusterDistances = [];
            FarRewClusterDistances = [];
            
            pixpermicron = 4.65;
            if isfield(varargin{i}, 'ZoomValue')
                if varargin{i}.ZoomValue ~= 0
                    pixpermicron = (pixpermicron*varargin{i}.ZoomValue)/12.1;
                end
            end
            
            FarDistanceMap = DistanceMap;
            FarDistanceMap(logical(eye(length(DistanceMap),length(DistanceMap)))) = 0;
            [r c] = find(isnan(triu(FarDistanceMap)));
            
            for m = 1:length(r)
                spine_pos1 = [varargin{i}.ROIPosition{r(m)}(1)+varargin{i}.ROIPosition{r(m)}(3)/2, varargin{i}.ROIPosition{r(m)}(2)+varargin{i}.ROIPosition{r(m)}(4)/2];
                spine_pos2 = [varargin{i}.ROIPosition{c(m)}(1)+varargin{i}.ROIPosition{c(m)}(3)/2, varargin{i}.ROIPosition{c(m)}(2)+varargin{i}.ROIPosition{c(m)}(4)/2];
                FarDistanceMap(r(m),c(m)) = (sqrt((spine_pos1(1)-spine_pos2(1)).^2 +(spine_pos1(2)-spine_pos2(2)).^2))/pixpermicron;
            end
            
            for CS = 1:length(FarClusteredSpines)
                if ~isempty(FarClusteredSpines{CS})
                    
                    spinecombos = nchoosek(FarClusteredSpines{CS}, 2); %%% Find all possible combinations of spines in a cluster
                    all_far_clust_distances = zeros(1,size(spinecombos,1));
                    far_cue_clust_distances = [];
                    far_mov_clust_distances = [];
                    far_mix_clust_distances = [];
                    far_presuc_clust_distances = [];
                    far_suc_clust_distances = [];
                    far_movduringcue_clust_distances = [];
                    far_rew_clust_distances = [];
                    
                    for n = 1:size(spinecombos,1)
                        all_far_clust_distances(1,n) = FarDistanceMap(spinecombos(n,1), spinecombos(n,2));
                        if sum(CueSpines(spinecombos(n,1:2))) == 2;
                            far_cue_clust_distances(1,n) = FarDistanceMap(spinecombos(n,1), spinecombos(n,2));
                        else
                            far_cue_clust_distances(1,n) = nan;
                        end
                        if sum(MovementSpines(spinecombos(n,1:2))) == 2;
                            far_mov_clust_distances(1,n) = FarDistanceMap(spinecombos(n,1), spinecombos(n,2));
                        else
                            far_mov_clust_distances(1,n) = nan;
                        end
                        if sum([sum(CueSpines(spinecombos(n,1:2))), sum(MovementSpines(spinecombos(n,1:2)))]) >= 2
                            far_mix_clust_distances(1,n) = FarDistanceMap(spinecombos(n,1), spinecombos(n,2));
                        else
                            far_mix_clust_distances(1,n) = nan;
                        end
                        if sum(PreSuccessSpines(spinecombos(n,1:2))) == 2
                            far_presuc_clust_distances(1,n) = FarDistanceMap(spinecombos(n,1),spinecombos(n,2));
                        else
                            far_presuc_clust_distances(1,n) = nan;
                        end
                        if sum(SuccessSpines(spinecombos(n,1:2))) == 2;
                            far_suc_clust_distances(1,n) = FarDistanceMap(spinecombos(n,1), spinecombos(n,2));
                        else
                            far_suc_clust_distances(1,n) = nan;
                        end
                        if sum(MovementDuringCueSpines(spinecombos(n,1:2))) == 2;
                            far_movduringcue_clust_distances(1,n) = FarDistanceMap(spinecombos(n,1), spinecombos(n,2));
                        else
                            far_movduringcue_clust_distances(1,n) = nan;
                        end
                        if sum(RewardSpines(spinecombos(n,1:2))) == 2;
                            far_rew_clust_distances(1,n) = FarDistanceMap(spinecombos(n,1), spinecombos(n,2));
                        else
                            far_rew_clust_distances(1,n) = nan;
                        end
                    end
                else
                    far_cue_clust_distances = nan;
                    far_mov_clust_distances = nan;
                    far_mix_clust_distances = nan;
                    far_presuc_clust_distances = nan;
                    far_suc_clust_distances = nan;
                    far_movduringcue_clust_distances = nan;
                    far_rew_clust_distances = nan;
                    all_far_clust_distances = nan;
                end
                AllFarClusterDistances{val} = all_far_clust_distances;
                FarCueClusterDistances{val} = far_cue_clust_distances;
                FarMovClusterDistances{val} = far_mov_clust_distances;
                FarMixClusterDistances{val} = far_mix_clust_distances;
                FarPreSucClusterDistances{val} = far_presuc_clust_distances;
                FarSucClusterDistances{val} = far_suc_clust_distances;
                FarMovDuringCueClusterDistances{val} = far_movduringcue_clust_distances;
                FarRewClusterDistances{val} = far_rew_clust_distances;
                val = val+1;
            end

            MeanFarCueClustLength(1,session) = nanmean(cell2mat(FarCueClusterDistances));
            MeanFarMovClustLength(1,session) = nanmean(cell2mat(FarMovClusterDistances));
            MeanFarMixClustLength(1,session) = nanmean(cell2mat(FarMixClusterDistances));
            MeanFarPreSucClustLength(1,session) = nanmean(cell2mat(FarPreSucClusterDistances));
            MeanFarSucClustLength(1,session) = nanmean(cell2mat(FarSucClusterDistances));
            MeanFarMovDuringCueClustLength(1,session) = nanmean(cell2mat(FarMovDuringCueClusterDistances));
            MeanFarRewClustLength(1,session) = nanmean(cell2mat(FarRewClusterDistances));
            MeanAllFarClustLength(1,session) = nanmean(cell2mat(AllFarClusterDistances));
            
            
            %%% causal clusters
            
            val = 1;
            AllCausClusterDistances = [];
            CausCueClusterDistances = [];
            CausMovClusterDistances = [];
            CausMixClusterDistances = [];
            CausPreSucClusterDistances = [];
            CausSucClusterDistances = [];
            CausMovDuringCueClusterDistances = [];
            CausRewClustDistances = [];
            

            for CS = 1:length(CausalClusteredSpines) 
                spinecombos = nchoosek(CausalClusteredSpines{CS}, 2); %%% Find all possible combinations of spines in a cluster
                all_clust_distances = zeros(1,size(spinecombos,1));
                cue_clust_distances = [];
                mov_clust_distances = [];
                suc_clust_distances = [];
                rew_clust_distances = [];

                for n = 1:size(spinecombos,1)
                    all_clust_distances(1,n) = DistanceMap(spinecombos(n,1), spinecombos(n,2));
                    if sum(MovementSpines(spinecombos(n,1:2))) == 2;
                        mov_clust_distances(1,n) = DistanceMap(spinecombos(n,1), spinecombos(n,2));
                    else
                        mov_clust_distances(1,n) = nan;
                    end
                end
                AllCausClusterDistances{val} = nanmean(all_clust_distances);
                CausCueClusterDistances{val} = nanmean(cue_clust_distances);
                CausMovClusterDistances{val} = nanmean(mov_clust_distances);
                CausSucClusterDistances{val} = nanmean(suc_clust_distances);
                CausRewClustDistances{val} = nanmean(rew_clust_distances); 
                val=val+1;
            end

            
            MeanCausalCueClustLength(1,session) = nanmean(cell2mat(CausCueClusterDistances));
            if ~isempty(cell2mat(CausCueClusterDistances))
                MaxCausalCueClustLength(1,session) = nanmax(cell2mat(CausCueClusterDistances));
            else
                MaxCausalCueClustLength(1,session) = NaN;
            end
            MeanCausalMovClustLength(1,session) = nanmean(cell2mat(CausMovClusterDistances));
            if ~isempty(cell2mat(CausMovClusterDistances))
                MaxCausalMovClustLength(1,session) = nanmax(cell2mat(CausMovClusterDistances));
            else
                MaxCausalMovClustLength(1,session) = NaN;
            end
            MeanCausalSucClustLength(1,session) = nanmean(cell2mat(CausSucClusterDistances));
            if ~isempty(cell2mat(CausSucClusterDistances))
                MaxCausalSucClustLength(1,session) = nanmax(cell2mat(CausSucClusterDistances));
            else
                MaxCausalSucClustLength(1,session) = NaN;
            end
            MeanCausalRewClustLength(1,session) = nanmean(cell2mat(CausRewClustDistances));
            if ~isempty(cell2mat(CausRewClustDistances))
                MaxCausalRewClustLength(1,session) = nanmax(cell2mat(CausRewClustDistances));
            else
                MaxCausalRewClustLength(1,session) = NaN;
            end
            MeanAllCausalClustLength(1,session) = nanmean(cell2mat(AllCausClusterDistances));
            if ~isempty(cell2mat(AllCausClusterDistances))
                MaxAllCausalClustLength(1,session) = nanmax(cell2mat(AllCausClusterDistances));
            else
                MaxAllCausalClustLength(1,session) = NaN;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Distance between functionally relevant spines
            
            allspine_combos = nchoosek(1:length(varargin{i}.dF_over_F),2);
            
            if length(find(CueSpines)) > 1
                cuespine_combos = nchoosek(find(CueSpines),2);
            else
                cuespine_combos = [];
            end
            
            if length(find(MovementSpines)) > 1
                movspine_combos = nchoosek(find(MovementSpines),2);
            else
                movspine_combos = [];
            end
            
            if length(find(PreSuccessSpines)) > 1
                presucspine_combos = nchoosek(find(PreSuccessSpines),2);
            else
                presucspine_combos = [];
            end
            
            if length(find(SuccessSpines)) > 1
                sucspine_combos = nchoosek(find(SuccessSpines),2);
            else
                sucspine_combos = [];
            end
            
            if length(find(MovementDuringCueSpines)) > 1
                movduringcuespine_combos = nchoosek(find(MovementDuringCueSpines),2);
            else
                movduringcuespine_combos = [];
            end
            
            if length(find(RewardSpines)) > 1
                rewspine_combos = nchoosek(find(RewardSpines),2);
            else
                rewspine_combos = [];
            end
            
            dist_between_all_spines = nan(1,size(allspine_combos,1));
                corr_between_all_spines = nan(1,size(allspine_combos,1));
                corr_between_all_spines_movperiods = nan(1,size(allspine_combos,1));
                corr_between_all_spines_stillperiods = nan(1,size(allspine_combos,1));
            dist_between_cue_spines = nan(1,size(cuespine_combos,1));
            dist_between_mov_spines = nan(1,size(movspine_combos,1));
                corr_between_mov_spines = nan(1,size(movspine_combos,1));
                corr_between_mov_spines_movperiods = nan(1,size(movspine_combos,1));
                corr_between_mov_spines_stillperiods = nan(1,size(movspine_combos,1));
            dist_between_presuc_spines = nan(1,size(presucspine_combos,1));
            dist_between_suc_spines = nan(1,size(sucspine_combos,1));
            dist_between_movduringcue_spines = nan(1,size(movduringcuespine_combos,1));
            dist_between_rew_spines = nan(1,size(rewspine_combos,1));
            dist_between_far_mov_spines = nan(1,size(allspine_combos,1));
                corr_between_far_mov_spines = nan(1,size(allspine_combos,1));
            
            tempcorrmat = separatedendcorrmat';
            separatedendcorrmat(isnan(separatedendcorrmat) & ~isnan(tempcorrmat)) = tempcorrmat(isnan(separatedendcorrmat) & ~isnan(tempcorrmat));
            
            
            for n = 1:size(allspine_combos,1)
                dist_between_all_spines(1,n) = DistanceMap(allspine_combos(n,1), allspine_combos(n,2));
                corr_between_all_spines(1,n) = samedendcorrmat(allspine_combos(n,1), allspine_combos(n,2));
                corr_between_all_spines_movperiods(1,n) = samedendmovcorrmat(allspine_combos(n,1), allspine_combos(n,2));
                corr_between_all_spines_stillperiods(1,n) = samedendstillcorrmat(allspine_combos(n,1), allspine_combos(n,2));
            end
            
            for n = 1:size(cuespine_combos,1)
                dist_between_cue_spines(1,n) = DistanceMap(cuespine_combos(n,1), cuespine_combos(n,2));
            end
            
            for n = 1:size(movspine_combos,1)
                dist_between_mov_spines(1,n) = DistanceMap(movspine_combos(n,1), movspine_combos(n,2));
                corr_between_mov_spines(1,n) = samedendcorrmat(movspine_combos(n,1), movspine_combos(n,2));
                corr_between_mov_spines_movperiods(1,n) = samedendmovcorrmat(movspine_combos(n,1), movspine_combos(n,2));
                corr_between_mov_spines_stillperiods(1,n) = samedendstillcorrmat(movspine_combos(n,1), movspine_combos(n,2));
                counter = 1;
                if ~isnan(separatedendcorrmat(movspine_combos(n,1), movspine_combos(n,2)))
                    numpos = find(separatedendcorrlist==separatedendcorrmat(movspine_combos(n,1), movspine_combos(n,2)));
                    if length(numpos)>1
                        separatedendcorrlist(numpos(1)) = NaN;  %%% Set the value in the correlation list to NaN so that it is not found again
                        numpos = numpos(1);
                    else
                        dist_between_far_mov_spines(1,n) = varargin{i}.FarSpineToSpineDistance(numpos);
                        corr_between_far_mov_spines(1,n) = separatedendcorrmat(movspine_combos(n,1), movspine_combos(n,2));
                    end
                else
                    dist_between_far_mov_spines(1,n) = nan;
                    corr_between_far_mov_spines(1,n) = nan;
                end
            end
            for n = 1:size(presucspine_combos,1)
                dist_between_presuc_spines(1,n) = DistanceMap(presucspine_combos(n,1), presucspine_combos(n,2));
            end
            for n = 1:size(sucspine_combos,1)
                dist_between_suc_spines(1,n) = DistanceMap(sucspine_combos(n,1), sucspine_combos(n,2));
            end
            for n = 1:size(movduringcuespine_combos,1)
                dist_between_movduringcue_spines = DistanceMap(movduringcuespine_combos(n,1), movduringcuespine_combos(n,2));
            end
            for n = 1:size(rewspine_combos,1)
                dist_between_rew_spines(1,n) = DistanceMap(rewspine_combos(n,1), rewspine_combos(n,2));
            end
            
           
            
            corr_between_all_spines = corr_between_all_spines(~isnan(dist_between_all_spines));
            corr_between_all_spines_movperiods = corr_between_all_spines_movperiods(~isnan(dist_between_all_spines));
            corr_between_all_spines_stillperiods = corr_between_all_spines_stillperiods(~isnan(dist_between_all_spines));
                CorrelationBetweenAllSpines{session} = corr_between_all_spines;
                CorrelationBetweenAllSpinesMovePeriods{session} = corr_between_all_spines_movperiods;
                CorrelationBetweenAllSpinesStillPeriods{session} = corr_between_all_spines_stillperiods;
                MeanCorrelationBetweenAllSpines{session} = nanmean(corr_between_all_spines);
            dist_between_all_spines = dist_between_all_spines(~isnan(dist_between_all_spines));
                DistanceBetweenAllSpines{session} = dist_between_all_spines;
            dist_between_cue_spines = dist_between_cue_spines(~isnan(dist_between_cue_spines));
                DistanceBetweenCueSpines{session} = dist_between_cue_spines;
                MeanDistanceBetweenCueSpines(1,session) = nanmean(dist_between_cue_spines);
            corr_between_mov_spines = corr_between_mov_spines(~isnan(dist_between_mov_spines)); %%% Needs to be calculated prior to updating (dist_between_mov_spines)!!
            corr_between_mov_spines_movperiods = corr_between_mov_spines_movperiods(~isnan(dist_between_mov_spines));
            corr_between_mov_spines_stillperiods = corr_between_mov_spines_stillperiods(~isnan(dist_between_mov_spines));
                CorrelationBetweenMovementSpines{session} = corr_between_mov_spines;
                CorrelationBetweenMovementSpinesMovePeriods{session} = corr_between_mov_spines_movperiods;
                CorrelationBetweenMovementSpinesStillPeriods{session} = corr_between_mov_spines_stillperiods;
                MeanCorrelationBetweenMovementSpines{1,session} = nanmean(corr_between_mov_spines);
            dist_between_mov_spines = dist_between_mov_spines(~isnan(dist_between_mov_spines));
                DistanceBetweenMovementSpines{session} = dist_between_mov_spines;
                MeanDistanceBetweenMovementSpines(1,session) = nanmean(dist_between_mov_spines);
            dist_between_presuc_spines = dist_between_presuc_spines(~isnan(dist_between_presuc_spines));
                DistanceBetweenPreSuccessSpines{session} = dist_between_presuc_spines;
                MeanDistanceBetweenPreSuccessSpines(1,session) = nanmean(dist_between_presuc_spines);
            dist_between_suc_spines = dist_between_suc_spines(~isnan(dist_between_suc_spines));
                DistanceBetweenSuccessSpines{session} = dist_between_suc_spines;
                MeanDistanceBetweenSuccessSpines(1,session) = nanmean(dist_between_suc_spines);
            dist_between_movduringcue_spines = dist_between_movduringcue_spines(~isnan(dist_between_movduringcue_spines));
                DistanceBetweenMovementDuringCueSpines{session} = dist_between_movduringcue_spines;
                MeanDistanceBetweenMovementDuringCueSpines(1,session) = nanmean(dist_between_movduringcue_spines);
            dist_between_rew_spines = dist_between_rew_spines(~isnan(dist_between_rew_spines));
                DistanceBetweenRewardSpines{session} = dist_between_rew_spines;
                MeanDistanceBetweenRewardSpines(1,session) = nanmean(dist_between_rew_spines);
            
            corr_between_far_mov_spines = corr_between_far_mov_spines(~isnan(dist_between_far_mov_spines));
                CorrelationBetweenFarMovementSpines{session} = corr_between_far_mov_spines;
            dist_between_far_mov_spines = dist_between_far_mov_spines(~isnan(dist_between_far_mov_spines));
                DistanceBetweenFarMovementSpines{session} = dist_between_far_mov_spines;
            CorrelationBetweenFarSpines{session} = separatedendcorrlist;
            if isempty(separatedendcorrlist) && ~isempty(varargin{i}.FarSpineToSpineDistance)
                if sum(varargin{i}.FarSpineToSpineDistance)
                    DistanceBetweenFarSpines{session} = [];
                else
                end
            else
                DistanceBetweenFarSpines{session} = varargin{i}.FarSpineToSpineDistance;
            end
            
            if sum(MovementSpines)>1
                spinelist = find(MovementSpines);
                for s = 1:length(spinelist)
                    currentspine = spinelist(s);            %%% Select current spine on the list
                    otherspines = setdiff(spinelist, currentspine); %%% Make a list of all other spines in the list
                    movspinedistlist = sort(cell2mat(arrayfun(@(x) DistanceMap(currentspine,x), otherspines, 'uni', false)));
                    NearestMovSpine{session}(1,s) = movspinedistlist(1);   %%% Find the smallest distance from the current spine to other spines on the same dendrite
                    if length(movspinedistlist) >1
                        NextClosest{session}(1,s) = movspinedistlist(2);
                    else
                        NextClosest{session}(1,s) = NaN;
                    end
                    if length(movspinedistlist) >2
                        ThirdClosest{session}(1,s) = movspinedistlist(3);
                    else
                        ThirdClosest{session}(1,s) = NaN;
                    end
                    if length(movspinedistlist) > 3
                        FourthClosest{session}(1,s) = movspinedistlist(4);
                    end
                    minloc = find(cell2mat(arrayfun(@(x) DistanceMap(currentspine,x), otherspines, 'uni', false)) == movspinedistlist(1));
                    if ~isempty(minloc)
                        CorrwithNearestMovSpine{session}(1,s) = mean(samedendcorrmat(currentspine, otherspines(minloc)));
                        otherhighcorrspines = mov_cluster_ind-Spine1_address;
                        distancestoHSCs = sort(cell2mat(arrayfun(@(x) DistanceMap(currentspine,x), otherhighcorrspines, 'uni', false)));
                        if ismember(currentspine, otherhighcorrspines)
                            NearestHighCorrMovSpine{session}(1,s) = distancestoHSCs(1);
                        end
                    if length(distancestoHSCs)>1;
                        NextClosestHighCorrMovSpine{session}(1,s) = distancestoHSCs(2);
                    end
                    if length(distancestoHSCs)>2;
                        ThirdClosestHighCorrMovSpine{session}(1,s) = distancestoHSCs(3);
                    end

                    end
                end
            end
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Correlation of Clusters with Different Features %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            allcorrs = [];
            invcheckmat = ones(varargin{i}.NumberofSpines,1);
            
            allspines = Spine1_address+(1:length(varargin{i}.deltaF));
            nonclustered = allspines(~ismember(allspines,cluster_ind));
            causal_nonclustered = allspines(~ismember(allspines,Caus_cluster_ind));
                     
            
            allclusterscorrwithcue = SessionCorrData(Cue_address, cluster_ind);
            nonclusteredcorrwithcue = SessionCorrData(Cue_address, nonclustered);
                cueclusterscorrwithcue = SessionCorrData(Cue_address, cue_cluster_ind);
                cuenonclusteredcorrwithcue = SessionCorrData(Cue_address, setdiff(find(CueSpines),cue_cluster_ind-Spine1_address)+Spine1_address).*dendcheckC(setdiff(find(CueSpines),cue_cluster_ind-Spine1_address));
            allclusterscorrwithmov = SessionCorrData(Movement_address, cluster_ind);
            nonclusteredcorrwithmov = SessionCorrData(Movement_address, nonclustered);
                movclusterscorrwithmov = SessionCorrData(Movement_address, mov_cluster_ind);
                movnonclusteredcorrwithmov = SessionCorrData(Movement_address, setdiff(find(MovementSpines),mov_cluster_ind-Spine1_address)+Spine1_address).*dendcheckM(setdiff(find(MovementSpines),mov_cluster_ind-Spine1_address));
            
            allclusterscorrwithMDC = SessionCorrData(MovementDuringCue_address, cluster_ind);
            nonclusteredcorrwithMDC = SessionCorrData(MovementDuringCue_address, nonclustered);
                MDCclusterscorrwithMDC = SessionCorrData(MovementDuringCue_address, movduringcue_cluster_ind);
                MDCnonclusteredcorrwithMDC = SessionCorrData(MovementDuringCue_address, setdiff(find(MovementDuringCueSpines),movduringcue_cluster_ind-Spine1_address)+Spine1_address).*dendcheckMDC(setdiff(find(MovementDuringCueSpines), movduringcue_cluster_ind-Spine1_address));
                
            allclusterscorrwithsuc = SessionCorrData(Success_address, cluster_ind);
            nonclusteredcorrwithsuc = SessionCorrData(Success_address, nonclustered);
                succlusterscorrwithsuc = SessionCorrData(Success_address, suc_cluster_ind);
                sucnonclusteredcorrwithsuc = SessionCorrData(Success_address, setdiff(find(SuccessSpines),suc_cluster_ind-Spine1_address)+Spine1_address).*dendcheckS(setdiff(find(SuccessSpines),suc_cluster_ind-Spine1_address));
            allclusterscorrwithrew = SessionCorrData(Reward_address, cluster_ind);
            nonclusteredcorrwithrew = SessionCorrData(Reward_address, nonclustered);
                rewclusterscorrwithrew = SessionCorrData(Reward_address, rew_cluster_ind);
                rewnonclusteredcorrwithrew = SessionCorrData(Reward_address, setdiff(find(RewardSpines),rew_cluster_ind-Spine1_address)+Spine1_address).*dendcheckR(setdiff(find(RewardSpines),rew_cluster_ind-Spine1_address));
            
                
                
                
            if usecorrminimum
                allclusterscorrwithcue = allclusterscorrwithcue(find(allclusterscorrwithcue>corrminimum));
                nonclusteredcorrwithcue = nonclusteredcorrwithcue(find(nonclusteredcorrwithcue>corrminimum));
                    cueclusterscorrwithcue = cueclusterscorrwithcue(find(cueclusterscorrwithcue>corrminimum));
                    cuenonclusteredcorrwithcue = cuenonclusteredcorrwithcue(find(cuenonclusteredcorrwithcue>corrminimum));
                allclusterscorrwithmov = allclusterscorrwithmov(find(allclusterscorrwithmov>corrminimum));
                nonclusteredcorrwithmov = nonclusteredcorrwithmov(find(nonclusteredcorrwithmov>corrminimum));
                    movclusterscorrwithmov = movclusterscorrwithmov(find(movclusterscorrwithmov>corrminimum));
                    movnonclusteredcorrwithmov = movnonclusteredcorrwithmov(find(movnonclusteredcorrwithmov>corrminimum));
                allclusterscorrwithMDC = allclusterscorrwithMDC(find(allclusterscorrwithMDC>corrminimum));
                nonclusteredcorrwithMDC = nonclusteredcorrwithMDC(find(nonclusteredcorrwithMDC>corrminimum));
                    MDCclusterscorrwithMDC = MDCclusterscorrwithMDC(find(MDCclusterscorrwithMDC>corrminimum));
                    MDCnonclusteredcorrwithMDC = MDCnonclusteredcorrwithMDC(find(MDCnonclusteredcorrwithMDC>corrminimum));
                allclusterscorrwithsuc = allclusterscorrwithsuc(find(allclusterscorrwithsuc>corrminimum));
                nonclusteredcorrwithsuc = nonclusteredcorrwithsuc(find(nonclusteredcorrwithsuc>corrminimum));
                    succlusterscorrwithsuc = succlusterscorrwithsuc(find(succlusterscorrwithsuc>corrminimum));
                    sucnonclusteredcorrwithsuc = sucnonclusteredcorrwithsuc(find(sucnonclusteredcorrwithsuc));
                allclusterscorrwithrew = allclusterscorrwithrew(find(allclusterscorrwithrew>corrminimum));
                nonclusteredcorrwithrew = nonclusteredcorrwithrew(find(nonclusteredcorrwithrew>corrminimum));
                    rewclusterscorrwithrew = rewclusterscorrwithrew(find(rewclusterscorrwithrew>corrminimum));
                    rewnonclusteredcorrwithrew = rewnonclusteredcorrwithrew(find(rewnonclusteredcorrwithrew));
            end

            clustcorrwithcue(1,session) = nanmean(allclusterscorrwithcue);
            nonclustcorrwithcue(1,session) = nanmean(nonclusteredcorrwithcue);
                cueclustcorrwithcue(1,session) = nanmean(cueclusterscorrwithcue);
                cuenonclustcorrwithcue(1,session) = nanmean(cuenonclusteredcorrwithcue);
            clustcorrwithmov(1,session) = nanmean(allclusterscorrwithmov);
            nonclustcorrwithmov(1,session) = nanmean(nonclusteredcorrwithmov);
                movclustcorrwithmov(1,session) = nanmean(movclusterscorrwithmov);
                movnonclustcorrwithmov(1,session) = nanmean(movnonclusteredcorrwithmov);
                
            clustcorrwithMDC(1,session) = nanmean(allclusterscorrwithMDC);
            nonclustcorrwithMDC(1,session) = nanmean(nonclusteredcorrwithMDC);
                MDCclustcorrwithMDC(1,session) = nanmean(MDCclusterscorrwithMDC);
                MDCnonclustcorrwithMDC(1,session) = nanmean(MDCnonclusteredcorrwithMDC);
            
            clustcorrwithsuc(1,session) = nanmean(allclusterscorrwithsuc);
            nonclustcorrwithsuc(1,session) = nanmean(nonclusteredcorrwithsuc);
                succlustcorrwithsuc(1,session) = nanmean(succlusterscorrwithsuc);
                sucnonclustcorrwithsuc(1,session) = nanmean(sucnonclusteredcorrwithsuc);
            clustcorrwithrew(1,session) = nanmean(allclusterscorrwithrew);
            nonclustcorrwithrew(1,session) = nanmean(nonclusteredcorrwithrew);
                rewclustcorrwithrew(1,session) = nanmean(rewclusterscorrwithrew);
                rewnonclustcorrwithrew(1,session) = nanmean(rewnonclusteredcorrwithrew);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Repeat for Causal Clusters %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
             
            if ~isempty(Caus_cluster_ind)
                Caus_allclusterscorrwithcue = Correlations{session}.CausalCorrelations(Cue_address, Caus_cluster_ind).*dendcheck(Caus_cluster_ind-Spine1_address);
                Caus_nonclusteredcorrwithcue = Correlations{session}.CausalCorrelations(Cue_address, causal_nonclustered).*dendcheck(causal_nonclustered-Spine1_address);
                    Caus_cueclusterscorrwithcue = Correlations{session}.CausalCorrelations(Cue_address, Caus_cue_cluster_ind).*dendcheckC(Caus_cue_cluster_ind-Spine1_address);
                    Caus_cuenonclusteredcorrwithcue = Correlations{session}.CausalCorrelations(Cue_address, setdiff(find(CueSpines),Caus_cue_cluster_ind-Spine1_address)+Spine1_address).*dendcheckC(setdiff(find(CueSpines),Caus_cue_cluster_ind-Spine1_address));
                
                Caus_allclusterscorrwithmov = Correlations{session}.CausalCorrelations(Movement_address, Caus_cluster_ind).*dendcheck(Caus_cluster_ind-Spine1_address);
                Caus_nonclusteredcorrwithmov = Correlations{session}.CausalCorrelations(Movement_address, causal_nonclustered).*dendcheck(causal_nonclustered-Spine1_address);
                    Caus_movclusterscorrwithmov = Correlations{session}.CausalCorrelations(Movement_address, Caus_mov_cluster_ind).*dendcheckM(Caus_mov_cluster_ind-Spine1_address);
                    Caus_movnonclusteredcorrwithmov = Correlations{session}.CausalCorrelations(Movement_address, setdiff(find(MovementSpines),Caus_mov_cluster_ind-Spine1_address)+Spine1_address).*dendcheckM(setdiff(find(MovementSpines),Caus_mov_cluster_ind-Spine1_address));
                
                Caus_allclusterscorrwithMDC = Correlations{session}.CausalCorrelations(MovementDuringCue_address, cluster_ind);
                Caus_nonclusteredcorrwithMDC = Correlations{session}.CausalCorrelations(MovementDuringCue_address, nonclustered).*dendcheck(nonclustered-Spine1_address);
                    Caus_MDCclusterscorrwithMDC = Correlations{session}.CausalCorrelations(MovementDuringCue_address, movduringcue_cluster_ind);
                    Caus_MDCnonclusteredcorrwithMDC = Correlations{session}.CausalCorrelations(MovementDuringCue_address, setdiff(find(MovementDuringCueSpines),movduringcue_cluster_ind-Spine1_address)+Spine1_address).*dendcheckMDC(setdiff(find(MovementDuringCueSpines), movduringcue_cluster_ind-Spine1_address));

                 
                Caus_allclusterscorrwithsuc = Correlations{session}.CausalCorrelations(Success_address, Caus_cluster_ind).*dendcheck(Caus_cluster_ind-Spine1_address);
                Caus_nonclusteredcorrwithsuc = Correlations{session}.CausalCorrelations(Success_address, causal_nonclustered).*dendcheck(causal_nonclustered-Spine1_address);
                    Caus_succlusterscorrwithsuc = Correlations{session}.CausalCorrelations(Success_address, Caus_suc_cluster_ind).*dendcheckS(Caus_suc_cluster_ind-Spine1_address);
                    Caus_sucnonclusteredcorrwithsuc = Correlations{session}.CausalCorrelations(Success_address, setdiff(find(SuccessSpines),Caus_suc_cluster_ind-Spine1_address)+Spine1_address).*dendcheckS(setdiff(find(SuccessSpines),Caus_suc_cluster_ind-Spine1_address));
                Caus_allclusterscorrwithrew = Correlations{session}.CausalCorrelations(Reward_address, Caus_cluster_ind).*dendcheck(Caus_cluster_ind-Spine1_address);
                Caus_nonclusteredcorrwithrew = Correlations{session}.CausalCorrelations(Reward_address, causal_nonclustered).*dendcheck(causal_nonclustered-Spine1_address);
                    Caus_rewclusterscorrwithrew = Correlations{session}.CausalCorrelations(Reward_address, Caus_rew_cluster_ind).*dendcheckR(Caus_rew_cluster_ind-Spine1_address);
                    Caus_rewnonclusteredcorrwithrew = Correlations{session}.CausalCorrelations(Reward_address, setdiff(find(RewardSpines),Caus_rew_cluster_ind-Spine1_address)+Spine1_address).*dendcheckR(setdiff(find(RewardSpines),Caus_rew_cluster_ind-Spine1_address));
            else
                Caus_allclusterscorrwithcue = [];
                Caus_nonclusteredcorrwithcue = [];
                    Caus_cueclusterscorrwithcue = [];
                    Caus_cuenonclusteredcorrwithcue = [];
                Caus_allclusterscorrwithmov = [];
                Caus_nonclusteredcorrwithmov = [];
                    Caus_movclusterscorrwithmov = [];
                    Caus_movnonclusteredcorrwithmov =[];
                Caus_allclusterscorrwithMDC = [];
                Caus_nonclusteredcorrwithMDC = [];
                    Caus_MDCclusterscorrwithMDC = [];
                    Caus_MDCnonclusteredcorrwithMDC = [];
                Caus_allclusterscorrwithsuc = [];
                Caus_nonclusteredcorrwithsuc = [];
                    Caus_succlusterscorrwithsuc = [];
                    Caus_sucnonclusteredcorrwithsuc = [];
                Caus_allclusterscorrwithrew = [];
                Caus_nonclusteredcorrwithrew = [];
                    Caus_rewclusterscorrwithrew = [];
                    Caus_rewnonclusteredcorrwithrew = [];
            end
                
            
            if usecorrminimum
                Caus_allclusterscorrwithcue = Caus_allclusterscorrwithcue(find(Caus_allclusterscorrwithcue>corrminimum));
                Caus_nonclusteredcorrwithcue = Caus_nonclusteredcorrwithcue(find(Caus_nonclusteredcorrwithcue>corrminimum));
                    Caus_cueclusterscorrwithcue = Caus_cueclusterscorrwithcue(find(Caus_cueclusterscorrwithcue>corrminimum));
                    Caus_cuenonclusteredcorrwithcue = Caus_cuenonclusteredcorrwithcue(find(Caus_cuenonclusteredcorrwithcue>corrminimum));
                Caus_allclusterscorrwithmov = Caus_allclusterscorrwithmov(find(Caus_allclusterscorrwithmov>corrminimum));
                Caus_nonclusteredcorrwithmov = Caus_nonclusteredcorrwithmov(find(Caus_nonclusteredcorrwithmov>corrminimum));
                    Caus_movclusterscorrwithmov = Caus_movclusterscorrwithmov(find(Caus_movclusterscorrwithmov>corrminimum));
                    Caus_movnonclusteredcorrwithmov = Caus_movnonclusteredcorrwithmov(find(Caus_movnonclusteredcorrwithmov));
                    
                Caus_allclusterscorrwithMDC = Caus_allclusterscorrwithMDC(find(Caus_allclusterscorrwithMDC>corrminimum));
                Caus_nonclusteredcorrwithMDC = Caus_nonclusteredcorrwithMDC(find(Caus_nonclusteredcorrwithMDC>corrminimum));
                    Caus_MDCclusterscorrwithMDC = Caus_MDCclusterscorrwithMDC(find(Caus_MDCclusterscorrwithMDC>corrminimum));
                    Caus_MDCnonclusteredcorrwithMDC = Caus_MDCnonclusteredcorrwithMDC(find(Caus_MDCnonclusteredcorrwithMDC>corrminimum));
                    
                Caus_allclusterscorrwithsuc = Caus_allclusterscorrwithsuc(find(Caus_allclusterscorrwithsuc>corrminimum));
                Caus_nonclusteredcorrwithsuc = Caus_nonclusteredcorrwithsuc(find(Caus_nonclusteredcorrwithsuc>corrminimum));
                    Caus_succlusterscorrwithsuc = Caus_succlusterscorrwithsuc(find(Caus_succlusterscorrwithsuc>corrminimum));
                    Caus_sucnonclusteredcorrwithsuc = Caus_sucnonclusteredcorrwithsuc(find(Caus_sucnonclusteredcorrwithsuc>corrminimum));
                Caus_allclusterscorrwithrew = Caus_allclusterscorrwithrew(find(Caus_allclusterscorrwithrew>corrminimum));
                Caus_nonclusteredcorrwithrew = Caus_nonclusteredcorrwithrew(find(Caus_nonclusteredcorrwithrew>corrminimum));
                    Caus_rewclusterscorrwithrew = Caus_rewclusterscorrwithrew(find(Caus_rewclusterscorrwithrew>corrminimum));
                    Caus_rewnonclusteredcorrwithrew = Caus_rewnonclusteredcorrwithrew(find(Caus_rewnonclusteredcorrwithrew>corrminimum));
            end

            Caus_clustcorrwithcue(1,session) = nanmean(Caus_allclusterscorrwithcue);
            Caus_nonclustcorrwithcue(1,session) = nanmean(Caus_nonclusteredcorrwithcue);
                Caus_cueclustcorrwithcue(1,session) = nanmean(Caus_cueclusterscorrwithcue);
                Caus_cuenonclustcorrwithcue(1,session) = nanmean(Caus_cuenonclusteredcorrwithcue);
            Caus_clustcorrwithmov(1,session) = nanmean(Caus_allclusterscorrwithmov);
            Caus_nonclustcorrwithmov(1,session) = nanmean(Caus_nonclusteredcorrwithmov);
                Caus_movclustcorrwithmov(1,session) = nanmean(Caus_movclusterscorrwithmov);
                Caus_movnonclustcorrwithmov(1,session) = nanmean(Caus_movnonclusteredcorrwithmov);
            Caus_clustcorrwithMDC(1,session) = nanmean(Caus_allclusterscorrwithMDC);
            Caus_nonclustcorrwithMDC(1,session) = nanmean(Caus_nonclusteredcorrwithMDC);
                Caus_MDCclustcorrwithMDC(1,session) = nanmean(Caus_MDCclusterscorrwithMDC);
                Caus_MDCnonclustcorrwithMDC(1,session) = nanmean(Caus_MDCnonclusteredcorrwithMDC);
            Caus_clustcorrwithsuc(1,session) = nanmean(Caus_allclusterscorrwithsuc);
            Caus_nonclustcorrwithsuc(1,session) = nanmean(Caus_nonclusteredcorrwithsuc);
                Caus_succlustcorrwithsuc(1,session) = nanmean(Caus_succlusterscorrwithsuc);
                Caus_sucnonclustcorrwithsuc(1,session) = nanmean(Caus_sucnonclusteredcorrwithsuc);
            Caus_clustcorrwithrew(1,session) = nanmean(Caus_allclusterscorrwithrew);
            Caus_nonclustcorrwithrew(1,session) = nanmean(Caus_nonclusteredcorrwithrew);
                Caus_rewclustcorrwithrew(1,session) = nanmean(Caus_rewclusterscorrwithrew);
                Caus_rewnonclustcorrwithrew(1,session) = nanmean(Caus_rewnonclusteredcorrwithrew);
                
                
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Plasticity Features
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%% Frequency

            cluster_freq(1,session) = nanmean(varargin{i}.SynapseOnlyFreq(cluster_ind-Spine1_address));
            
            nonclustered_freq(1,session) = nanmean(varargin{i}.SynapseOnlyFreq(nonclustered-Spine1_address));
            
            cue_cluster_freq(1,session) = nanmean(varargin{i}.SynapseOnlyFreq(cue_cluster_ind-Spine1_address));
            
            mov_cluster_freq(1,session) = nanmean(varargin{i}.SynapseOnlyFreq(mov_cluster_ind-Spine1_address));
            
            movduringcue_cluster_freq(1,session) = nanmean(varargin{i}.SynapseOnlyFreq(movduringcue_cluster_ind-Spine1_address));
            
            presuc_cluster_freq(1,session) = nanmean(varargin{i}.SynapseOnlyFreq(presuc_cluster_ind-Spine1_address));
            
            suc_cluster_freq(1,session) = nanmean(varargin{i}.SynapseOnlyFreq(suc_cluster_ind-Spine1_address));
            
            rew_cluster_freq(1,session) = nanmean(varargin{i}.SynapseOnlyFreq(rew_cluster_ind-Spine1_address));
            
            Caus_cluster_freq(1,session) = nanmean(varargin{i}.SpikeTimedEvents(Caus_cluster_ind-Spine1_address));
            
            Caus_nonclustered_freq(1,session) = nanmean(varargin{i}.SpikeTimedEvents(causal_nonclustered-Spine1_address));
            
            Caus_cue_cluster_freq(1,session) = nanmean(varargin{i}.SpikeTimedEvents(Caus_cue_cluster_ind-Spine1_address));
            
            Caus_mov_cluster_freq(1,session) = nanmean(varargin{i}.SpikeTimedEvents(Caus_mov_cluster_ind-Spine1_address));
            
            Caus_movduringcue_cluster_freq(1,session) = nanmean(varargin{i}.SpikeTimedEvents(Caus_movduringcue_cluster_ind-Spine1_address));
            
            Caus_presuc_cluster_freq(1,session) = nanmean(varargin{i}.SpikeTimedEvents(Caus_presuc_cluster_ind-Spine1_address));
            
            Caus_suc_cluster_freq(1,session) = nanmean(varargin{i}.SpikeTimedEvents(Caus_suc_cluster_ind-Spine1_address));
            
            Caus_rew_cluster_freq(1,session) = nanmean(varargin{i}.SpikeTimedEvents(Caus_rew_cluster_ind-Spine1_address));
            
            
            %%% Amplitude
            AmpList = cellfun(@(x) nanmean(x(x<100)), varargin{i}.AllSpineAmpData, 'Uni', false);
            clusterAmp = nanmean(cat(1,AmpList{cluster_ind-Spine1_address}));
            nonclusteredAmp = nanmean(cat(1,AmpList{nonclustered-Spine1_address}));
            cueAmp = nanmean(cat(1,AmpList{cue_cluster_ind-Spine1_address}));
            movAmp = nanmean(cat(1,AmpList{mov_cluster_ind-Spine1_address}));
            movduringcueAmp = nanmean(cat(1,AmpList{movduringcue_cluster_ind-Spine1_address}));
            presucAmp = nanmean(cat(1,AmpList{presuc_cluster_ind-Spine1_address}));
            sucAmp = nanmean(cat(1,AmpList{suc_cluster_ind-Spine1_address}));
            rewAmp = nanmean(cat(1,AmpList{rew_cluster_ind-Spine1_address}));
            Caus_clusterAmp = nanmean(cat(1,AmpList{Caus_cluster_ind-Spine1_address}));
            Caus_nonclusteredAmp = nanmean(cat(1,AmpList{causal_nonclustered-Spine1_address}));
            Caus_cueAmp = nanmean(cat(1,AmpList{Caus_cue_cluster_ind-Spine1_address}));
            Caus_movAmp = nanmean(cat(1,AmpList{Caus_mov_cluster_ind-Spine1_address}));
            Caus_movduringcueAmp = nanmean(cat(1,AmpList{Caus_movduringcue_cluster_ind-Spine1_address}));
            Caus_presucAmp = nanmean(cat(1,AmpList{Caus_presuc_cluster_ind-Spine1_address}));
            Caus_sucAmp = nanmean(cat(1,AmpList{Caus_suc_cluster_ind-Spine1_address}));
            Caus_rewAmp = nanmean(cat(1,AmpList{Caus_rew_cluster_ind-Spine1_address}));

                ClustAmp(1,session) = nanmean(clusterAmp);
                    NonClusteredAmp(1,session) = nanmean(nonclusteredAmp);
                Caus_ClustAmp(1,session) = nanmean(Caus_clusterAmp);
                    Caus_NonClusteredAmp(1,session) = nanmean(Caus_nonclusteredAmp);
                CueAmp(1,session) = nanmean(cueAmp);
                Caus_CueAmp(1,session) = nanmean(Caus_cueAmp);
                MovAmp(1,session) = nanmean(movAmp);
                Caus_MovAmp(1,session) = nanmean(Caus_movAmp);
                MovDuringCueAmp(1,session) = nanmean(movduringcueAmp);
                Caus_MovDuringCueAmp(1,session) = nanmean(Caus_movduringcueAmp);
                PreSucAmp(1,session) = nanmean(presucAmp);
                Caus_PreSucAmp(1,session) = nanmean(Caus_presucAmp);
                SucAmp(1,session) = nanmean(sucAmp);
                Caus_SucAmp(1,session) = nanmean(Caus_sucAmp);
                RewAmp(1,session) = nanmean(rewAmp);
                Caus_RewAmp(1,session) = nanmean(Caus_rewAmp);
                
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%
                                %%% Dendrites %%%
                                %%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%
            %%% Find which dendrites contain clusters
            %%%%%
            
            DendNum = varargin{i}.NumberofDendrites;
            
            DendswithClusts = [];
            DendswithCausClusts = [];
            DendsnoClusts = [];
            DendswithCueClusts = [];
            DendswithMovClusts = [];
            DendswithMovDuringCueClusts = [];
            DendswithPreSucClusts = [];
            DendswithSucClusts = [];
            DendswithRewClusts = [];
            DendswithCausmovClusts = [];
            DendsnomovClusts = [];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Count the number of clusters on each dendrite
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if isfield(varargin{i}, 'SpineDendriteGrouping')
                ClustersoneachDend = cell(length(varargin{i}.SpineDendriteGrouping),1);
                CausalClustersoneachDend = cell(length(varargin{i}.SpineDendriteGrouping),1);
                dcounter = cell(length(varargin{i}.SpineDendriteGrouping),1);
                Cdcounter = cell(length(varargin{i}.SpineDendriteGrouping),1);
                for k = 1:varargin{i}.NumberofDendrites
                    dcounter{k} = 0;
                    Cdcounter{k} = 0;
                end
                for k = 1:length(ClusteredSpines)
                    if ~isempty(ClusteredSpines{k})
                        for m = 1:length(varargin{i}.SpineDendriteGrouping)
                            if ~isempty(find(varargin{i}.SpineDendriteGrouping{m} == ClusteredSpines{k}(1)))
                                ClustersoneachDend{m}{dcounter{m}+1} = ClusteredSpines{k}; %%% Note that this labels spine by ROI number, and thus their actual count, NOT by the addressing scheme!
                                dcounter{m} = dcounter{m}+1;
                            end
                        end
                    end
                end
                for k = 1:length(CausalClusteredSpines)
                    if ~isempty(CausalClusteredSpines{k})
                        for m = 1:length(varargin{i}.SpineDendriteGrouping)
                            if ~isempty(find(varargin{i}.SpineDendriteGrouping{m} == CausalClusteredSpines{k}(1)))
                                CausalClustersoneachDend{m}{dcounter{m}+1} = CausalClusteredSpines{k}; %%% Note that this labels spine by ROI number, and thus their actual count, NOT by the addressing scheme!
                                Cdcounter{m} = Cdcounter{m}+1;
                            end
                        end
                    end
                end
            else
                ClustersoneachDend = {[]};
                CausalClustersoneachDend = {[]};
            end
            

            %%% Find the correlation of clustered vs. non-clustered spines
            %%% with dendritic activity (using pre-subtracted
            %%% spine/dendrite activity correlation heatmaps
             [r_dendsubtracted, ~] = corrcoef([varargin{i}.SynapseOnlyBinarized_DendriteSubtracted', varargin{i}.Dendrite_Binarized']);
            
            corrmattouse = r_dendsubtracted;
            causcorrmattouse = r_dendsubtracted;

            %%% initialize variables
            
            Clusters_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            CausalClusters_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            Nonclustered_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            Noncausal_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            FiltClusts_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            FiltCausClusts_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            
            CueClusts_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            MovClusts_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            PreSucClusts_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            SucClusts_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            MovDuringCueClusts_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            RewClusts_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            
            CueCausClusts_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            MovCausClusts_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            SucCausClusts_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            RewCausClusts_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);

            
            for k = 1:varargin{i}.NumberofDendrites
                if useSTATdends
                    usethisdend = MovementDends(k); %%% Finds the boolean corresponding to whether this dendrite is movement related or not
                else
                    usethisdend = 1;
                end
                if usethisdend
                    clusteredspinesonthisdend = [];
                    causalclusteredspinesonthisdend = [];
                    if ~isempty(ClustersoneachDend{k})
                        for l = 1:length(ClustersoneachDend{k})
                            clusteredspinesonthisdend = [clusteredspinesonthisdend; ClustersoneachDend{k}{l}];
                        end
                    end
                    if ~isempty(CausalClustersoneachDend{k})
                        for l = 1:length(CausalClustersoneachDend{k})
                            causalclusteredspinesonthisdend = [causalclusteredspinesonthisdend; CausalClustersoneachDend{k}{l}];
                        end
                    end
                    otherspinesonthisdend = varargin{i}.SpineDendriteGrouping{k}; %%% Labels spines by actual count/ROI num
                        nonclusteredcausalspinesonthisdend = otherspinesonthisdend(~ismember(otherspinesonthisdend, causalclusteredspinesonthisdend));
                    otherspinesonthisdend = otherspinesonthisdend(~ismember(otherspinesonthisdend, clusteredspinesonthisdend));
                    Clusters_CorrwithDend(1,k) = nanmean(corrmattouse(clusteredspinesonthisdend, length(varargin{i}.deltaF)+k));
                        CausalClusters_CorrwithDend(1,k) = nanmean(causcorrmattouse(causalclusteredspinesonthisdend, length(varargin{i}.deltaF)+k));
                    Nonclustered_CorrwithDend(1,k) = nanmean(corrmattouse(otherspinesonthisdend, length(varargin{i}.deltaF)+k));
                        Noncausal_CorrwithDend(1,k) = nanmean(causcorrmattouse(nonclusteredcausalspinesonthisdend, length(varargin{i}.deltaF)+k));
                    
                    Filtered_Clusts = clusteredspinesonthisdend(ismember(clusteredspinesonthisdend, cluster_ind-Spine1_address)); %%% Correct for addressing scheme           %%% Apply all the filtering criteria established in the first section of this function
                        Filtered_CausalClusts = causalclusteredspinesonthisdend(ismember(causalclusteredspinesonthisdend, cluster_ind-Spine1_address));
                    FiltClusts_CorrwithDend(1,k) = nanmean(corrmattouse(Filtered_Clusts, length(varargin{i}.deltaF)+k));
                        FiltCausClusts_CorrwithDend(1,k) = nanmean(causcorrmattouse(Filtered_CausalClusts, length(varargin{i}.deltaF)+k));
                    CueClusts_CorrwithDend(1,k) = nanmean(corrmattouse(cue_cluster_ind-Spine1_address, length(varargin{i}.deltaF)+k));
                        CueCausClusts_CorrwithDend(1,k) = nanmean(causcorrmattouse(Caus_cue_cluster_ind-Spine1_address, length(varargin{i}.deltaF)+k));
                    MovClusts_CorrwithDend(1,k) = nanmean(corrmattouse(mov_cluster_ind-Spine1_address, length(varargin{i}.deltaF)+k));
                        MovCausClusts_CorrwithDend(1,k) = nanmean(causcorrmattouse(Caus_mov_cluster_ind-Spine1_address, length(varargin{i}.deltaF)+k));
                    PreSucClusts_CorrwithDend(1,k) = nanmean(corrmattouse(presuc_cluster_ind-Spine1_address, length(varargin{i}.deltaF)+k));
                        PreSucCausClusts_CorrwithDend(1,k) = nanmean(causcorrmattouse(Caus_presuc_cluster_ind-Spine1_address, length(varargin{i}.deltaF)+k));
                    SucClusts_CorrwithDend(1,k) = nanmean(corrmattouse(suc_cluster_ind-Spine1_address, length(varargin{i}.deltaF)+k));
                        SucCausClusts_CorrwithDend(1,k) = nanmean(causcorrmattouse(Caus_suc_cluster_ind-Spine1_address, length(varargin{i}.deltaF)+k));
                    MovDuringCueClusts_CorrwithDend(1,k) = nanmean(corrmattouse(movduringcue_cluster_ind-Spine1_address, length(varargin{i}.deltaF)+k));
                        MovDuringCueCausClusts_CorrwithDend(1,k) = nanmean(causcorrmattouse(movduringcue_cluster_ind-Spine1_address, length(varargin{i}.deltaF)+k));
                    RewClusts_CorrwithDend(1,k) = nanmean(corrmattouse(rew_cluster_ind-Spine1_address, length(varargin{i}.deltaF)+k));
                        RewCausClusts_CorrwithDend(1,k) = nanmean(causcorrmattouse(Caus_rew_cluster_ind-Spine1_address, length(varargin{i}.deltaF)+k));
                else
                    continue
                end
            end
            
            ClusteredSpines_CorrwithDend(1,session) = nanmean(Clusters_CorrwithDend);
            FilteredClusteredSpines_CorrwithDend(1,session) = nanmean(FiltClusts_CorrwithDend);
            NonClusteredSpines_CorrwithDend(1,session) = nanmean(Nonclustered_CorrwithDend);
            
            CausalClusteredSpines_CorrwithDend(1,session) = nanmean(CausalClusters_CorrwithDend);
            FilteredCausalClusteredSpines_CorrwithDend(1,session) = nanmean(FiltCausClusts_CorrwithDend);
            NonCausalClusteredSpines_CorrwithDend(1,session) = nanmean(Noncausal_CorrwithDend);

            CueRelClusteredSpines_CorrwithDend(1,session) = nanmean(CueClusts_CorrwithDend);
            MovRelClusteredSpines_CorrwithDend(1,session) = nanmean(MovClusts_CorrwithDend);
            PreSucRelClusteredSpines_CorrwithDend(1,session) = nanmean(PreSucClusts_CorrwithDend);
            SucRelClusteredSpines_CorrwithDend(1,session) = nanmean(SucClusts_CorrwithDend);
            MovDuringCueRelClusteredSpines_CorrwithDend(1,session) = nanmean(MovDuringCueClusts_CorrwithDend);
            RewRelClusteredSpines_CorrwithDend(1,session) = nanmean(RewClusts_CorrwithDend);
            
            CueRelCausalClusteredSpines_CorrwithDend(1,session) = nanmean(CueCausClusts_CorrwithDend);
            MovRelCausalClusteredSpines_CorrwithDend(1,session) = nanmean(MovCausClusts_CorrwithDend);
            SucRelCausalClusteredSpines_CorrwithDend(1,session) = nanmean(SucCausClusts_CorrwithDend);
            RewRelCausalClusteredSpines_CorrwithDend(1,session) = nanmean(RewCausClusts_CorrwithDend);

            if ~isempty(cluster_ind)
                for k = 1:numel(cluster_ind)
                    for m = 1:length(varargin{i}.SpineDendriteGrouping)
                        if ~isempty(find(varargin{i}.SpineDendriteGrouping{m} == cluster_ind(k)-Spine1_address))
                            DendswithClusts = [DendswithClusts; m];
                        end
                    end
                end
            else
            end
            
            if ~isempty(Caus_cluster_ind)
                for k = 1:numel(Caus_cluster_ind)
                    for m = 1:length(varargin{i}.SpineDendriteGrouping)
                        if ~isempty(find(varargin{i}.SpineDendriteGrouping{m} == Caus_cluster_ind(k)-Spine1_address))
                            DendswithCausClusts = [DendswithCausClusts; m];
                        end
                    end
                end
            else
            end
            
            %%% Count the number of cue-related clustered spines on each
            %%% dendrite
            if ~isempty(cue_cluster_ind)
                for k = 1:numel(cue_cluster_ind)
                    for m = 1:length(varargin{i}.SpineDendriteGrouping)
                        if ~isempty(find(varargin{i}.SpineDendriteGrouping{m} == cue_cluster_ind(k)-Spine1_address))
                            DendswithCueClusts = [DendswithCueClusts; m];
                        end
                    end
                end
            else
            end
            %%% Count the number of mov-related clustered spines on each
            %%% dendrite
            if ~isempty(mov_cluster_ind)
                for k = 1:numel(mov_cluster_ind)
                    for m = 1:length(varargin{i}.SpineDendriteGrouping)
                        if ~isempty(find(varargin{i}.SpineDendriteGrouping{m} == mov_cluster_ind(k)-Spine1_address))
                            DendswithMovClusts = [DendswithMovClusts; m];
                        end
                    end
                end
            else
            end
            if ~isempty(movduringcue_cluster_ind)
                for k = 1:numel(movduringcue_cluster_ind)
                    for m = 1:length(varargin{i}.SpineDendriteGrouping)
                        if ~isempty(find(varargin{i}.SpineDendriteGrouping{m} == movduringcue_cluster_ind(k)-Spine1_address))
                            DendswithMovDuringCueClusts = [DendswithMovDuringCueClusts; m];
                        end
                    end
                end
            else
            end
            if ~isempty(presuc_cluster_ind)
                for k = 1:numel(presuc_cluster_ind)
                    for m = 1:length(varargin{i}.SpineDendriteGrouping)
                        if ~isempty(find(varargin{i}.SpineDendriteGrouping{m} == presuc_cluster_ind(k)-Spine1_address))
                            DendswithPreSucClusts = [DendswithPreSucClusts; m];
                        end
                    end
                end
            else
            end
            if ~isempty(suc_cluster_ind)
                for k = 1:numel(suc_cluster_ind)
                    for m = 1:length(varargin{i}.SpineDendriteGrouping)
                        if ~isempty(find(varargin{i}.SpineDendriteGrouping{m} == suc_cluster_ind(k)-Spine1_address))
                            DendswithSucClusts = [DendswithSucClusts; m];
                        end
                    end
                end
            else
            end
            if ~isempty(rew_cluster_ind)
                for k = 1:numel(rew_cluster_ind)
                    for m = 1:length(varargin{i}.SpineDendriteGrouping)
                        if ~isempty(find(varargin{i}.SpineDendriteGrouping{m} == rew_cluster_ind(k)-Spine1_address))
                            DendswithRewClusts = [DendswithRewClusts; m];
                        end
                    end
                end
            else
            end
            if ~isempty(Caus_mov_cluster_ind)
                for k = 1:numel(Caus_mov_cluster_ind)
                    Caus_movAmp(1,k) = nanmean(varargin{i}.AllSpineAmpData{Caus_mov_cluster_ind(k)-Spine1_address});
                    for m = 1:length(varargin{i}.SpineDendriteGrouping)
                        if ~isempty(find(varargin{i}.SpineDendriteGrouping{m} == Caus_mov_cluster_ind(k)-Spine1_address))
                            DendswithCausmovClusts = [DendswithCausmovClusts; m];
                        end
                    end
                end
            else
            end

                               
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Build an array of information for each dendrite and
            %%% associated clustering information:
            %%% Col 1;       Col 2:       Col 3;       Col 4;       
            %%% Dend #       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            DendritesArray = [1:varargin{i}.NumberofDendrites]';                                      %%% All dendrites from the current experiment
            allindex = ismember(DendritesArray, DendswithClusts);
            cueindex = ismember(DendritesArray, DendswithCueClusts);
            movindex = ismember(DendritesArray, DendswithMovClusts);                                %%% Logical array of which dendrites contain clusters
            movduringcueindex = ismember(DendritesArray, DendswithMovDuringCueClusts);
            presucindex = ismember(DendritesArray, DendswithPreSucClusts);
            sucindex = ismember(DendritesArray, DendswithSucClusts);
            rewindex = ismember(DendritesArray, DendswithRewClusts);
            
            DendsnoClusts = DendritesArray(allindex ==  0);
            DendsnocueClusts = DendritesArray(cueindex == 0);
            DendsnomovClusts = DendritesArray(movindex == 0);
            DendsnomovduringcueClusts = DendritesArray(movduringcueindex == 0);
            DendsnopresucClusts = DendritesArray(presucindex == 0);
            DendsnosucClusts = DendritesArray(sucindex == 0);
            DendsnorewClusts = DendritesArray(rewindex == 0);
            
            DendNumClust{session} = DendritesArray;                                                 %%% Will be a 4-column array; column 1 is the dendrites index
            for n = 1:length(DendNumClust{session})
                DendNumClust{session}(n,2) = sum(DendswithClusts == DendNumClust{session}(n,1));    %%% Column 2: Counts number of clustered spines on a given dendrite
                DendNumClust{session}(n,3) = sum(DendswithMovClusts == DendNumClust{session}(n,1)); %%% Column 3: Counts number of movement-related clustered spines on a given dendrite 
                DendNumClust{session}(n,4) = length(ClustersoneachDend{n});                         %%% Column 4: Number of separate clusters on each dendrite
                DendNumClust{session}(n,5) = varargin{i}.Dendritic_Frequency(n);                    %%% Column 5: Frequency of dendrite
            end
            
            
                DendswithClusts = unique(DendswithClusts);
                DendsnoClusts = unique(DendsnoClusts);
                DendswithCueClusts = unique(DendswithCueClusts);
                DendswithMovClusts = unique(DendswithMovClusts);
                DendswithMovDuringCueClusts = unique(DendswithMovDuringCueClusts);
                DendswithPreSucClusts = unique(DendswithPreSucClusts);
                DendswithSucClusts = unique(DendswithSucClusts);
                DendswithRewClusts = unique(DendswithRewClusts);
                DendsnomovClusts = unique(DendsnomovClusts);
               
                
                if useSTATdends
                    DendswithClusts = DendswithClusts.*MovementDends(DendswithClusts);
                        DendswithClusts = DendswithClusts(find(DendswithClusts));
                    DendsnoClusts = DendsnoClusts.*MovementDends(DendsnoClusts);
                        DendsnoClusts = DendsnoClusts(find(DendsnoClusts));
                    DendswithCueClusts = DendswithCueClusts.*CueDends(DendswithCueClusts);
                        DendswithCueClusts = DendswithCueClusts(find(DendswithCueClusts));
                    DendswithMovClusts = DendswithMovClusts.*MovementDends(DendswithMovClusts);
                        DendswithMovClusts = DendswithMovClusts(find(DendswithMovClusts));
                    DendswithMovDuringCueClusts = DendswithmovDuringCueClusts.*MovementDuringCueDends(DendswithMovDuringCueClusts);
                        DendswithMovDuringCueClusts = DendswithMovDuringCueClusts(find(DendswithMovDuringCueClusts));
                    DendswithPreSucClusts = DendswithPreSucClusts.*PreSuccessDends(DenswithPreSuccessClusts)
                    DendswithSucClusts = DendswithSucClusts.*SuccessDends(DendswithSucClusts);
                        DendswithSucClusts = DendswithSucClusts(find(DendswithSucClusts));
                    DendswithRewClusts = DendswithRewClusts.*RewardDends(DendswithRewClusts);
                        DendswithRewClusts = DendswithRewClusts(find(DendswithRewClusts));
                    DendsnomovClusts = DendsnomovClusts.*MovementDends(DendsnomovClusts);
                        DendsnomovClusts = DendsnomovClusts(find(DendsnomovClusts));
                    
                    ClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendswithClusts));
                    NoClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendsnoClusts));
                    CueClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendswithCueClusts));
                    MovClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendswithMovClusts));
                    MovDuringCueClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendswithMovDuringCueClusts));
                    PreSuccessClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendswithPreSucClusts));
                    SuccessClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendswithSucClusts));
                    RewardClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendswithRewClusts));
                    NoMovClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendsnomovClusts));
                    counter = 1;
                    
                    for j = 1:DendNum
                        firstspine = varargin{i}.SpineDendriteGrouping{j}(1);
                        lastspine = varargin{i}.SpineDendriteGrouping{j}(end);

                        if length(varargin{i}.SpineDendriteGrouping{j})>mindendsize
                            if strcmpi(laplaciantouse, 'Normalized')
                                L = varargin{i}.NormalizedLaplacian{j};
                            elseif strcmpi(laplaciantouse, 'Original')
                                L = varargin{i}.LaplacianMatrix{j};
                            end
                            eigenvalues = eig(L);
                            Fiedlerval = min(eigenvalues(~ismember(eigenvalues,min(eigenvalues)))); %%% Finds the second smallest eigenvalue (the Fiedler value, or measure or algebraic connectivity)
                            
                            temporalvalues = eig(Temporal_Laplacian{session}{j});
                            temporalvalues(temporalvalues==min(temporalvalues)) = nan;
                            [TemporalFiedlerval, ~] = min(temporalvalues);
                            
                            if MovementDends(j)
                                DendClust_Deg{session}(counter,1) = real(Fiedlerval);
                                DendClust_Deg{session}(counter,2) = real(TemporalFiedlerval);
                                DendClust_Deg{session}(counter,3) = real(SpatioTemporalFiedler{session}(j));
                                DendClust_Deg{session}(counter,4) = varargin{i}.Dendritic_Frequency(j);
                                DendClust_Deg{session}(counter,5) = SessionCorrData(Movement_address, Spine1_address+length(varargin{i}.deltaF)+j);
                                DendClust_Deg{session}(counter,6) = length(varargin{i}.SpineDendriteGrouping{j});
                                counter = counter+1;
                            end
                        else
                            DendClust_Deg{session}(counter,1) = nan;
                            DendClust_Deg{session}(counter,2) = nan;
                            DendClust_Deg{session}(counter,3) = nan;
                            DendClust_Deg{session}(counter,4) = nan;
                            DendClust_Deg{session}(counter,5) = nan;
                            DendClust_Deg{session}(counter,6) = nan;
                            counter = counter+1;                        
                        end
                    end
                else
                    ClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendswithClusts));
                    NoClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendsnoClusts));
                    CueClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendswithCueClusts));
                    MovClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendswithMovClusts));
                    MovDuringCueClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendswithMovDuringCueClusts));
                    PreSucClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendswithPreSucClusts));
                    SucClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendswithSucClusts));
                    RewClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendswithRewClusts));
                    NoMovClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendsnomovClusts));
                    counter = 1;
                    
                    for j = 1:DendNum
                        firstspine = varargin{i}.SpineDendriteGrouping{j}(1);
                        lastspine = varargin{i}.SpineDendriteGrouping{j}(end);

                        if length(varargin{i}.SpineDendriteGrouping{j})>mindendsize
                            if strcmpi(laplaciantouse, 'Normalized')
                                L = varargin{i}.NormalizedLaplacian{j};
                            elseif strcmpi(laplaciantouse, 'Original')
                                L = varargin{i}.LaplacianMatrix{j};
                            end
                            eigenvalues = eig(L);
                            Fiedlerval = min(eigenvalues(~ismember(eigenvalues,min(eigenvalues)))); %%% Finds the second smallest eigenvalue (the Fiedler value, or measure or algebraic connectivity)
                            
                            temporalvalues = eig(Temporal_Laplacian{session}{j});
                            temporalvalues(temporalvalues==min(temporalvalues)) = nan;
                            [TemporalFiedlerval, ~] = min(temporalvalues);
                            DendClust_Deg{session}(counter,1) = real(Fiedlerval);
                            DendClust_Deg{session}(counter,2) = real(TemporalFiedlerval);
                            DendClust_Deg{session}(counter,3) = real(SpatioTemporalFiedler{session}(j));
                            DendClust_Deg{session}(counter,4) = varargin{i}.Dendritic_Frequency(j);
                            DendClust_Deg{session}(counter,5) = SessionCorrData(Movement_address, Spine1_address+length(varargin{i}.deltaF)+j);
                            DendClust_Deg{session}(counter,6) = length(varargin{i}.SpineDendriteGrouping{j});
                            counter = counter+1;
                        else
                            DendClust_Deg{session}(counter,1) = nan;
                            DendClust_Deg{session}(counter,2) = nan;
                            DendClust_Deg{session}(counter,3) = nan;
                            DendClust_Deg{session}(counter,4) = nan;
                            DendClust_Deg{session}(counter,5) = nan;
                            DendClust_Deg{session}(counter,6) = nan;
                            counter = counter+1;
                        end
                    end
                end
%             FractionofCluster = cellfun(@(x,y) sum(x(y))/length(y), repmat(MovementSpineDataToUse(session),[1 length(ClusteredSpines)]),ClusteredSpines,'Uni', false);
%             FractionofCausalCluster = cellfun(@(x,y) sum(x(y))/length(y), repmat(MovementSpineDataToUse(session),[1 length(CausalClusteredSpines)]),CausalClusteredSpines,'Uni', false);
%             FractionofClusterThatsMR(1,session) = nanmean(cat(2,FractionofCluster{:}));
%             FractionofCausalClusterThatsMR(1,session) = nanmean(cat(2,FractionofCausalCluster{:})); 


        else                %%% When there is no behavioral data for this session...
            
            
            if isfield(varargin{i}, 'SpineDendriteGrouping')
                for d = 1:varargin{i}.NumberofDendrites
                    firstspine = varargin{i}.SpineDendriteGrouping{d}(1);
                    lastspine = varargin{i}.SpineDendriteGrouping{d}(end);
                end
            else
                if varargin{i}.NumberofDendrites == 1
                    varargin{i}.SpineDendriteGrouping = {1:length(varargin{i}.deltaF)};
                    disp(['File ', varargin{i}.Filename, ' did not have spine-dendrite grouping information, but can be correctly filled in']);
                else
                    disp(['File ', varargin{i}.Filename, ' does NOT have necessary spine-dendrite grouping info!']);
                end
            end

            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%% Spectral analysis of dendrites %%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%% When no behavioral data is available, simply use
            %%% correlations of calcium traces alone
            
            if dendexclude
                [rsp psp] = corrcoef(varargin{i}.SynapseOnlyBinarized');
            elseif dendsubtract
                [rsp psp] = corrcoef(varargin{i}.SynapseOnlyBinarized_DendriteSubtracted');
            end
            
            corrmat = nan(varargin{i}.NumberofSpines+Spine1_address, varargin{i}.NumberofSpines+Spine1_address);
            
            corrmat(Spine1_address+1:Spine1_address+size(rsp,1),Spine1_address+1:Spine1_address+size(rsp,2)) = rsp;
            
            [SpectralData] = SpectralClustering(varargin{i}, corrmat, StatClass, Choices);
            

            Spatial_Deg{session} = SpectralData.Spatial_Deg;
            Temporal_Deg{session} = SpectralData.Temporal_Deg;
            Spatiotemporal_Deg{session} = SpectralData.Spatiotemporal_Deg;
            Spatiotemporal_Overlap{session} = SpectralData.Spatiotemporal_Overlap;
            
            SpatioTemporalFiedler{session} = SpectralData.SpatioTemporalFiedler;
            
            Temporal_Laplacian{session} = SpectralData.Temporal_Laplacian;

            MeanSpatialDegreeofCueSpines(1,session) = SpectralData.MeanSpatialDegreeofCueSpines;
            MeanTemporalDegreeofCueSpines(1,session) = SpectralData.MeanTemporalDegreeofCueSpines;
            MeanSpatioTemporalDegreeofCueSpines(1,session) = SpectralData.MeanSpatioTemporalDegreeofCueSpines;
            MeanSpatialDegreeofMovementSpines(1,session) = SpectralData.MeanSpatialDegreeofMovementSpines;
            MeanTemporalDegreeofMovementSpines(1,session) = SpectralData.MeanTemporalDegreeofMovementSpines;
            MeanSpatioTemporalDegreeofMovementSpines(1,session) = SpectralData.MeanSpatioTemporalDegreeofMovementSpines;
            MeanSpatialDegreeofMDCSpines(1,session) = SpectralData.MeanSpatialDegreeofMovementDuringCueSpines;
            MeanTemporalDegreeofMDCSpines(1,session) = SpectralData.MeanTemporalDegreeofMovementDuringCueSpines;
            MeanSpatioTemporalDegreeofMDCSpines(1,session) = SpectralData.MeanSpatioTemporalDegreeofMovementDuringCueSpines;
            MeanSpatialDegreeofPreSuccessSpines(1,session) = SpectralData.MeanSpatialDegreeofPreSuccessSpines;
            MeanTemporalDegreeofPreSuccessSpines(1,session) = SpectralData.MeanTemporalDegreeofPreSuccessSpines;
            MeanSpatioTemporalDegreeofPreSuccessSpines(1,session) = SpectralData.MeanSpatioTemporalDegreeofPreSuccessSpines;
            MeanSpatialDegreeofSuccessSpines(1,session) = SpectralData.MeanSpatialDegreeofSuccessSpines;
            MeanTemporalDegreeofSuccessSpines(1,session) = SpectralData.MeanTemporalDegreeofSuccessSpines;
            MeanSpatioTemporalDegreeofSuccessSpines(1,session) = SpectralData.MeanSpatioTemporalDegreeofSuccessSpines;
            MeanSpatialDegreeofRewardSpines(1,session) = SpectralData.MeanSpatialDegreeofRewardSpines;
            MeanTemporalDegreeofRewardSpines(1,session) = SpectralData.MeanTemporalDegreeofRewardSpines;
            MeanSpatioTemporalDegreeofRewardSpines(1,session) = SpectralData.MeanSpatioTemporalDegreeofRewardSpines;
            
            Dend_Spat_Deg{session} = SpectralData.Dend_Spat_Deg;
            Dend_Temp_Deg{session} = SpectralData.Dend_Temp_Deg;
            Dend_SpatTemp_Deg{session} = SpectralData.Dend_SpatTemp_Deg;      
            
            SpatialDegree_vs_Movement{session} = SpectralData.SpatialDegree_vs_Movement;
            TemporalDegree_vs_Movement{session} = SpectralData.TemporalDegree_vs_Movement;
            SpatiotemporalDegree_vs_Movement{session} = SpectralData.SpatiotemporalDegree_vs_Movement;

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%% Binary Classification of Clusters %%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%% Make full correlation matrix %%%
            
            if dendexclude
                corrmat = corrcoef([varargin{i}.SynapseOnlyBinarized]');
            elseif dendsubtract
                corrmat = corrcoef([varargin{i}.SynapseOnlyBinarized_DendriteSubtracted]');
            end
%             corrmat = corrcoef([varargin{i}.OverallSpineActivity]');
            
            %%% Causal clusters

            Clustnum = 1;
            tempA = varargin{i}.CausalHeatMap;
            tempB = varargin{i}.CausalHeatMap';
            tempA(isnan(tempA) & ~isnan(tempB)) = tempB(isnan(tempA) & ~isnan(tempB));
            causalmat = tempA;
            
            [ClusteredSpines, CorrelationofClusts, FarClusteredSpines, CausalClusteredSpines] = binaryclusterclass(varargin{i},corrmat, causalmat);
            
            CorrelationofClusters{session} = CorrelationofClusts;
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Find indices for clusters
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Find indices for movement- and cue-related clusters
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%% Correlation indices for synapse-only events
            
            StatClass{session}.MovementSpines = zeros(length(varargin{i}.deltaF),1);
            StatClass{session}.CueSpines = zeros(length(varargin{i}.deltaF),1);
            StatClass{session}.CueORMovementSpines = zeros(length(varargin{i}.deltaF),1);
            StatClass{session}.PreSuccessSpines = zeros(length(varargin{i}.deltaF),1);
            StatClass{session}.SuccessSpines = zeros(length(varargin{i}.deltaF),1);
            StatClass{session}.MovementDuringCueSpines = zeros(length(varargin{i}.deltaF),1);
            StatClass{session}.RewardSpines = zeros(length(varargin{i}.deltaF),1);
            
            StatClass{session}.DendSub_MovementSpines = zeros(length(varargin{i}.deltaF),1);
            StatClass{session}.DendSub_CueSpines = zeros(length(varargin{i}.deltaF),1);
            StatClass{session}.DendSub_PreSuccessSpines = zeros(length(varargin{i}.deltaF),1);
            StatClass{session}.DendSub_SuccessSpines = zeros(length(varargin{i}.deltaF),1);
            StatClass{session}.DendSub_MovementDuringCueSpines = zeros(length(varargin{i}.deltaF),1);
            StatClass{session}.DendSub_RewardSpines = zeros(length(varargin{i}.deltaF),1);
            
               [ClusterInfo]...
                = findclusterind(ClusteredSpines, StatClass{session}, Spine1_address);
            
            NumSpinesinCluster = ClusterInfo.NumSpinesinCluster;
            cluster_ind = ClusterInfo.cluster_ind;
            cue_cluster_ind = ClusterInfo.cue_cluster_ind;
            mov_cluster_ind = ClusterInfo.mov_cluster_ind;
            mix_cluster_ind = ClusterInfo.mix_cluster_ind;
            presuc_cluster_ind = ClusterInfo.suc_cluster_ind;
            suc_cluster_ind = ClusterInfo.presuc_cluster_ind;
            movduringcue_cluster_ind = ClusterInfo.movduringcue_cluster_ind;
            rew_cluster_ind = ClusterInfo.rew_cluster_ind;
            
            %%   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            %%%% Find indices for clusters on separate dendrites
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            [ClusterInfo]...
                = findclusterind(FarClusteredSpines, StatClass{session}, Spine1_address);
            
            NumSpinesinFarCluster = ClusterInfo.NumSpinesinCluster;
            Farcluster_ind = ClusterInfo.cluster_ind;
            Farcue_cluster_ind = ClusterInfo.cue_cluster_ind;
            Farmov_cluster_ind = ClusterInfo.mov_cluster_ind;
            Farmix_cluster_ind = ClusterInfo.mix_cluster_ind;
            Farpresuc_cluster_ind = ClusterInfo.suc_cluster_ind;
            Farsuc_cluster_ind = ClusterInfo.presuc_cluster_ind;
            Farmovduringcue_cluster_ind = ClusterInfo.movduringcue_cluster_ind;
            Farrew_cluster_ind = ClusterInfo.rew_cluster_ind;
            
            %%   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            %%%% Find indices for causal clusters
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%% Uncomment the following section to separate cue- and
            %%%% movement-correlated clusters
            
            StatClass{i}.CausalMovementSpines = zeros(length(varargin{i}.deltaF),1);
            
            [ClusterInfo]...
                = findclusterind(CausalClusteredSpines, StatClass{session}, Spine1_address);
            
            NumSpinesinCausalCluster = ClusterInfo.NumSpinesinCluster;
            Caus_cluster_ind = ClusterInfo.cluster_ind;
            Caus_cue_cluster_ind = ClusterInfo.cue_cluster_ind;
            Caus_mov_cluster_ind = ClusterInfo.mov_cluster_ind;
            Caus_mix_cluster_ind = ClusterInfo.mix_cluster_ind;
            Caus_presuc_cluster_ind = ClusterInfo.suc_cluster_ind;
            Caus_suc_cluster_ind = ClusterInfo.presuc_cluster_ind;
            Caus_movduringcue_cluster_ind = ClusterInfo.movduringcue_cluster_ind;
            Caus_rew_cluster_ind = ClusterInfo.rew_cluster_ind;


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%% End index-finding section %%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            

            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%% Main variables
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ClusteredSpines = ClusteredSpines(~cellfun(@isempty,ClusteredSpines));
            CausalClusteredSpines = CausalClusteredSpines(~cellfun(@isempty,CausalClusteredSpines));

            NumClusters(1,session) = length(ClusteredSpines);
            NumCausalClusters(1,session) = length(CausalClusteredSpines);
           
            NumSpinesinCluster(NumSpinesinCluster==0) = NaN;
            MeanNumberofSpinesinEachCluster(1,session) = nanmean(NumSpinesinCluster);
            NumSpinesinCausalCluster(NumSpinesinCausalCluster==0) = NaN;
            MeanNumberofSpinesinEachCausalCluster(1,session) = nanmean(NumSpinesinCausalCluster);
            
            if ~isempty(ClusteredSpines)
                MeanNumberofSpinesinEachMovCluster(1,session) = nanmean(cellfun(@(x) x(:), cellfun(@(x) sum(x)==length(x), cellfun(@(x) ismember(x,mov_cluster_ind-Spine1_address), ClusteredSpines, 'UniformOutput', false), 'Uniformoutput', false)));
                NumMovClusters(1,session) = sum(cellfun(@(x) x(:), cellfun(@(x) sum(x)>1, cellfun(@(x) ismember(x,mov_cluster_ind-Spine1_address), ClusteredSpines, 'UniformOutput', false), 'Uniformoutput', false)));  %%% Check if each reported cluster is movement-related
            else
                MeanNumberofSpinesinEachMovCluster(1,session) = NaN;
                NumMovClusters(1,session) = NaN;
            end
            
            %%% Find number of clustered spines in each category
            
            %%%%%%%%%%%%%%%%%%%            
            %%% All Clustered spines
            %%%%%%%%%%%%%%%%%%%
            
            NumCueRelSpines(1,session) = nan;
            NumMovRelSpines(1,session) = nan;
            FractionofMovRelSpinesPerDendrite{session} = nan;
            NumCueORMovRelSpines(1,session) = nan;
            NumPreSucRelSpines(1,session) = nan;
            NumSucRelSpines(1,session) = nan;
            NumMovDuringCueRelSpines(1,session) = nan;
            NumRewRelSpines(1,session) = nan;
            NumCausalMovSpines(1,session) = nan;
            
            NumClustSpines(1,session) = length(cluster_ind)/varargin{i}.NumberofSpines;
            NumClustCueSpines(1,session) = nan;
            NumClustMovSpines(1,session) = nan;
            NumClustMixSpines(1,session) = nan;
            NumClustPreSucSpines(1,session) = nan;
            NumClustSucSpines(1,session) = nan;
            NumClustMovDuringCueSpines(1,session) = nan;
            NumClustRewSpines(1,session) = nan;
            
            NumFarClustSpines(1,session) = length(Farcluster_ind)/varargin{i}.NumberofSpines;
            NumFarClustCueSpines(1,session) = nan;
            NumFarClustMovSpines(1,session) = nan;
            NumFarClustMixSpines(1,session) = nan;
            NumFarClustPreSucSpines(1,session) = nan;
            NumFarClustSucSpines(1,session) = nan;
            NumFarClustMovDuringCueSpines(1,session) = nan;
            NumFarClustRewSpines(1,session) = nan;

            NumCausClustSpines(1,session) = length(Caus_cluster_ind)/varargin{i}.NumberofSpines;
            NumCausClustCueSpines(1,session) = nan;
            NumCausClustMovSpines(1,session) = nan;
            NumCausClustMixSpines(1,session) = nan;
            NumCausClustPreSucSpines(1,session) = nan;
            NumCausClustSucSpines(1,session) = nan;
            NumCausClustMovDuringCueSpines(1,session) = nan;
            NumCausClustRewSpines(1,session) = nan;
            
            
            %%% Distance between spines (standard clusters)
            if ~isempty(ClusteredSpines)
                allclustlength = nan(1,length(ClusteredSpines));
            else
                allclustlength = nan;
            end
            for CS = 1:length(ClusteredSpines)
                spinecombos = nchoosek(ClusteredSpines{CS}, 2); %%% Find all possible combinations of spines in a cluster
                all_clust_distances = zeros(1,size(spinecombos,1));

                for n = 1:size(spinecombos,1)
                    all_clust_distances(1,n) = DistanceMap(spinecombos(n,1), spinecombos(n,2));
                end
                allclustlength(1,CS) = nanmean(all_clust_distances);
            end

            MeanAllClustLength(1,session) = nanmean(allclustlength);
            MaxAllClustLength(1,session) = nanmax(allclustlength);
            
            %%% Distance between spines ("far" clusters/separate dendrites)
            
            pixpermicron = 4.65;
            if isfield(varargin{i}, 'ZoomValue')
                if varargin{i}.ZoomValue ~= 0
                    pixpermicron = (pixpermicron*varargin{i}.ZoomValue)/12.1;
                end
            end
            
            FarDistanceMap = DistanceMap;
            FarDistanceMap(logical(eye(length(DistanceMap),length(DistanceMap)))) = 0;
            [r c] = find(isnan(triu(FarDistanceMap)));
            
            for m = 1:length(r)
                spine_pos1 = [varargin{i}.ROIPosition{r(m)}(1)+varargin{i}.ROIPosition{r(m)}(3)/2, varargin{i}.ROIPosition{r(m)}(2)+varargin{i}.ROIPosition{r(m)}(4)/2];
                spine_pos2 = [varargin{i}.ROIPosition{c(m)}(1)+varargin{i}.ROIPosition{c(m)}(3)/2, varargin{i}.ROIPosition{c(m)}(2)+varargin{i}.ROIPosition{c(m)}(4)/2];
                FarDistanceMap(r(m),c(m)) = (sqrt((spine_pos1(1)-spine_pos2(1)).^2 +(spine_pos1(2)-spine_pos2(2)).^2))/pixpermicron;
            end

            if ~isempty(FarClusteredSpines)
                for CS = 1:length(FarClusteredSpines)
                    if ~isempty(FarClusteredSpines{CS})

                        spinecombos = nchoosek(FarClusteredSpines{CS}, 2); %%% Find all possible combinations of spines in a cluster
                        all_clust_distances = zeros(1,size(spinecombos,1));

                        for n = 1:size(spinecombos,1)
                            all_clust_distances(1,n) = FarDistanceMap(spinecombos(n,1), spinecombos(n,2));
                        end
                        Faralllength = nanmean(all_clust_distances);
                    else
                        Faralllength = nan;
                    end
                end
            else
                Faralllength = [];
            end
            
            MeanAllFarClustLength(1,session) = nanmean(Faralllength);
            
            
            %%% Distance between spines (causal clusters)
            if ~isempty(CausalClusteredSpines)
                for CS = 1:length(CausalClusteredSpines)
                    if ~isempty(CausalClusteredSpines{CS})

                        spinecombos = nchoosek(CausalClusteredSpines{CS}, 2); %%% Find all possible combinations of spines in a cluster
                        all_clust_distances = zeros(1,size(spinecombos,1));

                        for n = 1:size(spinecombos,1)
                            all_clust_distances(1,n) = DistanceMap(spinecombos(n,1), spinecombos(n,2));
                        end
                        Causalllength = nanmean(all_clust_distances);
                    else
                        Causalllength = nan;
                    end
                end
            else
                Causalllength = [];
            end
            
            MeanAllCausalClustLength(1,session) = nanmean(Causalllength);
            if ~isempty(CausalClusteredSpines)
                MaxAllCausalClustLength(1,session) = nanmax(Causalllength);
            end

            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Plasticity Features
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            allspines = Spine1_address+(1:length(varargin{i}.deltaF));
            nonclustered = allspines(~ismember(allspines,cluster_ind));
            causal_nonclustered = allspines(~ismember(allspines,Caus_cluster_ind));
            
            %%% Frequency
            cluster_freq(1,session) = nanmean(varargin{i}.SynapseOnlyFreq(cluster_ind-Spine1_address));
            
            nonclustered_freq(1,session) = nanmean(varargin{i}.SynapseOnlyFreq(nonclustered-Spine1_address));
                         
            Caus_cluster_freq(1,session) = nanmean(varargin{i}.SpikeTimedEvents(Caus_cluster_ind-Spine1_address));
            
            Caus_nonclustered_freq(1,session) = nanmean(varargin{i}.SpikeTimedEvents(causal_nonclustered-Spine1_address));
                        
            
            %%% Amplitude
            AmpList = cellfun(@(x) nanmean(x(x<100)), varargin{i}.AllSpineAmpData, 'Uni', false);
            clusterAmp = nanmean(cat(1,AmpList{cluster_ind-Spine1_address}));
            nonclusteredAmp = nanmean(cat(1,AmpList{nonclustered-Spine1_address}));
            Caus_clusterAmp = nanmean(cat(1,AmpList{Caus_cluster_ind-Spine1_address}));
            Caus_nonclusteredAmp = nanmean(cat(1,AmpList{causal_nonclustered-Spine1_address}));

                ClustAmp(1,session) = nanmean(clusterAmp);
                    NonClusteredAmp(1,session) = nanmean(nonclusteredAmp);
                Caus_ClustAmp(1,session) = nanmean(Caus_clusterAmp);
                    Caus_NonClusteredAmp(1,session) = nanmean(Caus_nonclusteredAmp);
            
            %%
            %%%%%%%%%%%%%%%%%
            %%% Dendrites %%%
            %%%%%%%%%%%%%%%%%
            
            DendNum = varargin{i}.NumberofDendrites; 
            
            %%%%%
            %%% Find which dendrites contain clusters
            %%%%%
            
            DendswithClusts = [];
            DendswithCausClusts = [];
            DendsnoClusts = [];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Count the number of clusters on each dendrite
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if isfield(varargin{i}, 'SpineDendriteGrouping')
                ClustersoneachDend = cell(length(varargin{i}.SpineDendriteGrouping),1);
                CausalClustersoneachDend = cell(length(varargin{i}.SpineDendriteGrouping),1);
                dcounter = cell(length(varargin{i}.SpineDendriteGrouping),1);
                Cdcounter = cell(length(varargin{i}.SpineDendriteGrouping),1);
                for k = 1:varargin{i}.NumberofDendrites
                    dcounter{k} = 0;
                    Cdcounter{k} = 0;
                end
                for k = 1:length(ClusteredSpines)
                    if ~isempty(ClusteredSpines{k})
                        for m = 1:length(varargin{i}.SpineDendriteGrouping)
                            if ~isempty(find(varargin{i}.SpineDendriteGrouping{m} == ClusteredSpines{k}(1)))
                                ClustersoneachDend{m}{dcounter{m}+1} = ClusteredSpines{k}; %%% Note that this labels spine by ROI number, and thus their actual count, NOT by the addressing scheme!
                                dcounter{m} = dcounter{m}+1;
                            end
                        end
                    end
                end
                for k = 1:length(CausalClusteredSpines)
                    if ~isempty(CausalClusteredSpines{k})
                        for m = 1:length(varargin{i}.SpineDendriteGrouping)
                            if ~isempty(find(varargin{i}.SpineDendriteGrouping{m} == CausalClusteredSpines{k}(1)))
                                CausalClustersoneachDend{m}{dcounter{m}+1} = CausalClusteredSpines{k}; %%% Note that this labels spine by ROI number, and thus their actual count, NOT by the addressing scheme!
                                Cdcounter{m} = Cdcounter{m}+1;
                            end
                        end
                    end
                end
            else
                ClustersoneachDend = {[]};
                CausalClustersoneachDend = {[]};
            end
            

            %%% Find the correlation of clustered vs. non-clustered spines
            %%% with dendritic activity (using pre-subtracted
            %%% spine/dendrite activity correlation heatmaps
            
            spinedatatouse = varargin{i}.SynapseOnlyBinarized_DendriteSubtracted;
            
             [r_dendsubtracted, ~] = corrcoef([spinedatatouse', varargin{i}.Dendrite_Binarized']);
            
            corrmattouse = r_dendsubtracted;
            causcorrmattouse = r_dendsubtracted;

            %%% initialize variables
            Clusters_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            CausalClusters_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            Nonclustered_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            Noncausal_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            FiltClusts_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            FiltCausClusts_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            MovClusts_CorrwithDend = nan(1,varargin{i}.NumberofDendrites);
            MovCausClusts_CorrwithDend = nan(1,varargin{i}.NumberofDendrites); 
            
            for k = 1:varargin{i}.NumberofDendrites
                    clusteredspinesonthisdend = [];
                    causalclusteredspinesonthisdend = [];
                    if ~isempty(ClustersoneachDend{k})
                        for l = 1:length(ClustersoneachDend{k})
                            clusteredspinesonthisdend = [clusteredspinesonthisdend; ClustersoneachDend{k}{l}];
                        end
                    end
                    if ~isempty(CausalClustersoneachDend{k})
                        for l = 1:length(CausalClustersoneachDend{k})
                            causalclusteredspinesonthisdend = [causalclusteredspinesonthisdend; CausalClustersoneachDend{k}{l}];
                        end
                    end
                    otherspinesonthisdend = varargin{i}.SpineDendriteGrouping{k}; %%% Labels spines by actual count/ROI num
                        nonclusteredcausalspinesonthisdend = otherspinesonthisdend(~ismember(otherspinesonthisdend, causalclusteredspinesonthisdend));
                    otherspinesonthisdend = otherspinesonthisdend(~ismember(otherspinesonthisdend, clusteredspinesonthisdend));
                    Clusters_CorrwithDend(1,k) = nanmean(corrmattouse(clusteredspinesonthisdend, length(varargin{i}.deltaF)+k));
                        CausalClusters_CorrwithDend(1,k) = nanmean(causcorrmattouse(causalclusteredspinesonthisdend, length(varargin{i}.deltaF)+k));
                    Nonclustered_CorrwithDend(1,k) = nanmean(corrmattouse(otherspinesonthisdend, length(varargin{i}.deltaF)+k));
                        Noncausal_CorrwithDend(1,k) = nanmean(causcorrmattouse(nonclusteredcausalspinesonthisdend, length(varargin{i}.deltaF)+k));
                    
                    Filtered_Clusts = clusteredspinesonthisdend(ismember(clusteredspinesonthisdend, cluster_ind-Spine1_address)); %%% Correct for addressing scheme           %%% Apply all the filtering criteria established in the first section of this function
                        Filtered_CausalClusts = causalclusteredspinesonthisdend(ismember(causalclusteredspinesonthisdend, cluster_ind-Spine1_address));
                    FiltClusts_CorrwithDend(1,k) = nanmean(corrmattouse(Filtered_Clusts, length(varargin{i}.deltaF)+k));
                        FiltCausClusts_CorrwithDend(1,k) = nanmean(causcorrmattouse(Filtered_CausalClusts, length(varargin{i}.deltaF)+k));
            end
            
            ClusteredSpines_CorrwithDend(1,session) = nanmean(Clusters_CorrwithDend);
            FilteredClusteredSpines_CorrwithDend(1,session) = nanmean(FiltClusts_CorrwithDend);
            NonClusteredSpines_CorrwithDend(1,session) = nanmean(Nonclustered_CorrwithDend);
            
            CausalClusteredSpines_CorrwithDend(1,session) = nanmean(CausalClusters_CorrwithDend);
            FilteredCausalClusteredSpines_CorrwithDend(1,session) = nanmean(FiltCausClusts_CorrwithDend);
            NonCausalClusteredSpines_CorrwithDend(1,session) = nanmean(Noncausal_CorrwithDend);
            
            %%% Find amplitude for each cluster, and associate with its
            %%% parent dendrite
            
            if ~isempty(cluster_ind)
                for k = 1:numel(cluster_ind)
                    clusterAmp(1,k) = nanmean(varargin{i}.AllSpineAmpData{cluster_ind(k)-Spine1_address});
                    for m = 1:length(varargin{i}.SpineDendriteGrouping)
                        if ~isempty(find(varargin{i}.SpineDendriteGrouping{m} == cluster_ind(k)-Spine1_address))
                            DendswithClusts = [DendswithClusts; m];
                        end
                    end
                end
            else
                clusterAmp = nan;
            end
            
            if ~isempty(Caus_cluster_ind)
                for k = 1:numel(Caus_cluster_ind)
                    Caus_clusterAmp(1,k) = nanmean(varargin{i}.AllSpineAmpData{Caus_cluster_ind(k)-Spine1_address});
                    for m = 1:length(varargin{i}.SpineDendriteGrouping)
                        if ~isempty(find(varargin{i}.SpineDendriteGrouping{m} == Caus_cluster_ind(k)-Spine1_address))
                            DendswithCausClusts = [DendswithCausClusts; m];
                        end
                    end
                end
            else
                Caus_clusterAmp = nan;
            end
            
            %%% Count the number of mov-related clustered spines on each
            %%% dendrite
            
                               
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Build an array of information for each dendrite and
            %%% associated clustering information:
            %%% Col 1;       Col 2:       Col 3;       Col 4;       
            %%% Dend #       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            DendritesArray = [1:varargin{i}.NumberofDendrites]';                                    %%% All dendrites from the current experiment
            allindex = ismember(DendritesArray, DendswithClusts);
            DendsnoClusts = DendritesArray(allindex ==  0);
            
            DendNumClust{session} = DendritesArray;                                                 %%% Will be a 4-column array; column 1 is the dendrites index
            for n = 1:length(DendNumClust{session})
                DendNumClust{session}(n,2) = sum(DendswithClusts == DendNumClust{session}(n,1));    %%% Column 2: Counts number of clustered spines on a given dendrite
                DendNumClust{session}(n,3) = nan;                                                   %%% Column 3: Counts number of movement-related clustered spines on a given dendrite 
                DendNumClust{session}(n,4) = length(ClustersoneachDend{n});                         %%% Column 4: Number of separate clusters on each dendrite
                DendNumClust{session}(n,5) = varargin{i}.Dendritic_Frequency(n);                    %%% Column 5: Frequency of dendrite
            end
            
            
                DendswithClusts = unique(DendswithClusts);
                DendsnoClusts = unique(DendsnoClusts);               
                
                ClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendswithClusts));
                NoClustDendFreq(1,session) = nanmean(varargin{i}.Dendritic_Frequency(DendsnoClusts));
                CueClustDendFreq(1,session) = nan;
                MovClustDendFreq(1,session) = nan;
                SucClustDendFreq(1,session) = nan;
                RewClustDendFreq(1,session) = nan;
                NoMovClustDendFreq(1,session) = nan;
                counter = 1;

                for j = 1:DendNum
                    firstspine = varargin{i}.SpineDendriteGrouping{j}(1);
                    lastspine = varargin{i}.SpineDendriteGrouping{j}(end);

                    if length(varargin{i}.SpineDendriteGrouping{j})>mindendsize
                        if strcmpi(laplaciantouse, 'Normalized')
                            L = varargin{i}.NormalizedLaplacian{j};
                        elseif strcmpi(laplaciantouse, 'Original')
                            L = varargin{i}.LaplacianMatrix{j};
                        end
                        eigenvalues = eig(L);
                        Fiedlerval = min(eigenvalues(~ismember(eigenvalues,min(eigenvalues)))); %%% Finds the second smallest eigenvalue (the Fiedler value, or measure or algebraic connectivity)

                        temporalvalues = eig(Temporal_Laplacian{session}{j});
                        temporalvalues(temporalvalues==min(temporalvalues)) = nan;
                        [TemporalFiedlerval, ~] = min(temporalvalues);
                        DendClust_Deg{session}(counter,1) = real(Fiedlerval);
                        DendClust_Deg{session}(counter,2) = real(TemporalFiedlerval);
                        DendClust_Deg{session}(counter,3) = real(SpatioTemporalFiedler{session}(j));
                        DendClust_Deg{session}(counter,4) = varargin{i}.Dendritic_Frequency(j);
                        DendClust_Deg{session}(counter,5) = nan;
                        DendClust_Deg{session}(counter,6) = length(varargin{i}.SpineDendriteGrouping{j});
                        counter = counter+1;
                    else
                        DendClust_Deg{session}(counter,1) = nan;
                        DendClust_Deg{session}(counter,2) = nan;
                        DendClust_Deg{session}(counter,3) = nan;
                        DendClust_Deg{session}(counter,4) = nan;
                        DendClust_Deg{session}(counter,5) = nan;
                        DendClust_Deg{session}(counter,6) = nan;
                        counter = counter+1;
                    end
                end
        end
    end
    

    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% Color Information %%
    %%%%%%%%%%%%%%%%%%%%%%%%

    lgray = [0.50 0.51 0.52];   brown = [0.28 0.22 0.14];
    gray = [0.50 0.51 0.52];    lbrown = [0.59 0.45 0.28];
    yellow = [1.00 0.76 0.05];  orange = [0.95 0.40 0.13];
    lgreen = [0.45 0.8 0.35];  green = [0.00 0.43 0.23];
    lblue = [0.30 0.65 0.94];   blue = [0.00 0.33 0.65];
    magenta = [0.93 0.22 0.55]; purple = [0.57 0.15 0.56];
    pink = [0.9 0.6 0.6];       lpurple  = [0.7 0.15 1];
    red = [0.85 0.11 0.14];     black = [0.1 0.1 0.15];
    dred = [0.6 0 0];          dorange = [0.8 0.3 0.03];
    bgreen = [0 0.6 0.7];
    colorj = {red,lblue,green,lgreen,gray,brown,yellow,blue,purple,lpurple,magenta,pink,orange,brown,lbrown};
    rnbo = {dred, red, dorange, orange, yellow, lgreen, green, bgreen, blue, lblue, purple, magenta, lpurple, pink}; 

    %%%%%%%%%%%%%%%%%%%%%%%


    scrsz = get(0, 'ScreenSize');
    
    %%%%%%%%%%%%%%%%
    %%% Figure 1 %%%
    %%% Correlation of Clustered vs. Non-clustered spines with behavioral
    %%% features %%%
    %%%%%%%%%%%%%%%%
    %%
    figure('Position', scrsz);
    subplot(2,5,1);
    plot(clustcorrwithcue, '-o','Color', black, 'MarkerFaceColor', 'k', 'Linewidth', 2); hold on; % plot(cueanticorr, '-o', 'Color', gray, 'MarkerFaceColor', gray, 'LineWidth', 2); 
    plot(nonclustcorrwithcue, '-o', 'MarkerFaceColor', dred,'Color', dred, 'Linewidth', 2);
    plot(cueclustcorrwithcue, '-o', 'MarkerFaceColor', blue,'Color', blue, 'Linewidth', 2);
    plot(cuenonclustcorrwithcue, '-o', 'MarkerFaceColor', gray, 'Color', gray, 'Linewidth', 2);
    title('Correlation with cue');
    legend({'All Clustered', 'Nonclustered', 'Clust Mov-rel', 'NClust mov-rel'}, 'Location', 'NW')
    
    %%%%%%%%%%%%%%%
    %%
    
    subplot(2,5,2)
    plot(clustcorrwithMDC, '-o', 'Color', black, 'MarkerFaceColor', 'k', 'Linewidth', 2); hold on;
    plot(nonclustcorrwithMDC, '-o', 'Color', dred, 'MarkerFaceColor', dred, 'Linewidth', 2);
    plot(MDCclustcorrwithMDC, '-o', 'Color', blue, 'MarkerFaceColor', blue, 'Linewidth', 2)
    plot(MDCnonclustcorrwithMDC, '-o', 'Color', gray, 'MarkerFaceColor', gray, 'Linewidth', 2)
    title('Correlation with MDC');
    legend({'All Clustered', 'Nonclustered', 'Clust Mov-rel', 'NClust mov-rel'}, 'Location', 'NW')

    %%%%%%%%%%%%%%%
    %%
   
    
    subplot(2,5,3);
    plot(clustcorrwithmov, '-o','Color', black, 'MarkerFaceColor', 'k', 'Linewidth', 2); hold on; % plot(clustanticorr, '-o', 'Color', gray, 'MarkerFaceColor', gray, 'LineWidth', 2);
    plot(nonclustcorrwithmov, '-o', 'MarkerFaceColor',dred, 'Color', dred, 'Linewidth', 2)
    plot(movclustcorrwithmov, '-o', 'MarkerFaceColor',blue, 'Color', blue, 'Linewidth', 2)
    plot(movnonclustcorrwithmov, '-o', 'MarkerFaceColor', gray, 'Color', gray, 'Linewidth', 2)
    title('Correlation with movement')
    legend({'All Clustered', 'Nonclustered', 'Clust Mov-rel', 'NClust mov-rel'}, 'Location', 'NW')
    
    
    %%
   
    
    subplot(2,5,4);
    plot(clustcorrwithsuc, '-o','Color', black, 'MarkerFaceColor', black, 'Linewidth', 2); hold on;
    plot(nonclustcorrwithsuc, '-o','Color', dred, 'MarkerFaceColor',dred, 'Linewidth', 2)
    plot(succlustcorrwithsuc, '-o','Color', blue, 'MarkerFaceColor',blue, 'Linewidth', 2)
    plot(sucnonclustcorrwithsuc, '-o', 'Color', gray, 'MarkerFaceColor', gray, 'Linewidth', 2)
    title('Correlation with rewarded presses');
    legend({'All Clustered', 'Nonclustered', 'Clust Mov-rel', 'NClust mov-rel'}, 'Location', 'NW')
    
    %%
    
    subplot(2,5,5)
    plot(clustcorrwithrew, '-o','Color', black, 'MarkerFaceColor', black, 'Linewidth', 2); hold on;
    plot(nonclustcorrwithrew, '-o','Color', dred, 'MarkerFaceColor', dred, 'Linewidth', 2);
    plot(rewclustcorrwithrew, '-o', 'MarkerFaceColor',blue, 'Color', blue, 'Linewidth', 2);
    plot(rewnonclustcorrwithrew, '-o', 'MarkerFaceColor', gray, 'Color', gray, 'Linewidth', 2);
    title('Correlation with reward');
    legend({'All Clustered', 'Nonclustered', 'Clust Mov-rel', 'NClust mov-rel'}, 'Location', 'NW')
    
    %%
    
    subplot(2,5,6)
    plot(Caus_clustcorrwithcue, '-o','Color', black, 'MarkerFaceColor', black, 'Linewidth', 2); hold on;% plot(causalclustanticorr, '-o', 'Color', gray, 'MarkerFaceColor', gray, 'LineWidth', 2);
    plot(Caus_nonclustcorrwithcue, '-o', 'MarkerFaceColor',dred, 'Color', dred, 'Linewidth', 2)
    plot(Caus_cueclustcorrwithcue, '-o', 'MarkerFaceColor',blue, 'Color', blue, 'Linewidth', 2)
    plot(Caus_cuenonclustcorrwithcue, '-o', 'MarkerFaceColor', gray, 'Color', gray, 'Linewidth', 2);
    title('Causal Correlation with Cue');
    
    
    %%
    
    subplot(2,5,7)
    plot(Caus_clustcorrwithMDC,  '-o','Color', black, 'MarkerFaceColor', black, 'Linewidth', 2); hold on;
    plot(Caus_nonclustcorrwithMDC, '-o', 'MarkerFaceColor', dred, 'Color', dred, 'Linewidth', 2)
    plot(Caus_MDCclustcorrwithMDC, '-o', 'MarkerFaceColor', blue, 'Color', blue, 'Linewidth', 2);
    plot(Caus_MDCnonclustcorrwithMDC, '-o', 'MarkerFaceColor', gray, 'Color', gray', 'Linewidth', 2);
    
    %%
    
    subplot(2,5,8)
    plot(Caus_clustcorrwithmov, '-o','Color', black, 'MarkerFaceColor', black, 'Linewidth', 2); hold on;
    plot(Caus_nonclustcorrwithmov, '-o', 'MarkerFaceColor',dred, 'Color', dred, 'Linewidth', 2)
    plot(Caus_movclustcorrwithmov, '-o', 'MarkerFaceColor',blue, 'Color', blue, 'Linewidth', 2);
    plot(Caus_movnonclustcorrwithmov, '-o', 'MarkerFaceColor', gray, 'Color', gray, 'Linewidth', 2);
    title('Causal Correlation with Mov');
    
    %%
    
    subplot(2,5,9)
    plot(Caus_clustcorrwithsuc, '-o','Color', black, 'MarkerFaceColor', black, 'Linewidth', 2); hold on;
    plot(Caus_nonclustcorrwithsuc, '-o', 'MarkerFaceColor',dred, 'Color', dred, 'Linewidth', 2)
    plot(Caus_succlustcorrwithsuc, '-o', 'MarkerFaceColor',blue, 'Color', blue, 'Linewidth', 2);
    plot(Caus_sucnonclustcorrwithsuc, '-o', 'MarkerFaceColor', gray, 'Color', gray, 'Linewidth', 2);
    title('Causal Correlation with Rewarded Presses');
    
    %%
    
    subplot(2,5,10)
    plot(Caus_clustcorrwithrew, '-o','Color', black, 'MarkerFaceColor', black, 'Linewidth', 2); hold on;
    plot(Caus_nonclustcorrwithrew, '-o', 'MarkerFaceColor',dred, 'Color', dred, 'Linewidth', 2)
    plot(Caus_rewclustcorrwithrew, '-o', 'MarkerFaceColor',blue, 'Color', blue, 'Linewidth', 2);
    plot(Caus_rewnonclustcorrwithrew, '-o', 'MarkerFaceColor', gray, 'Color', gray, 'Linewidth', 2);
    title('Causal Correlation with Reward');
    
    %%
    %%%%%%%%%%%%
    %%% Figure 2
    %%% Plasticity features
    %%%%%%%%%%%%
    
    figure('Position', scrsz); subplot(2,3,1); 
    plot(cluster_freq, '-o', 'Color', black, 'MarkerFaceColor', black, 'Linewidth', 2)
    hold on; plot(Caus_cluster_freq, '-o', 'Color', dred,'MarkerFaceColor', dred, 'Linewidth', 2);
    hold on; plot(nonclustered_freq, '-o', 'Color', lpurple, 'MarkerFaceColor', lpurple, 'Linewidth', 2);
    plot(Caus_nonclustered_freq, '-o', 'Color', lgreen, 'MarkerFaceColor', lgreen, 'Linewidth', 2);
    legend({'Clustered Spines', 'Causal Clusters', 'NonClustered', 'C nonclustered'})
    ylabel('Frequency of Clustered Spines')
    xlabel('Session')
    
    subplot(2,3,2); 
    plot(ClustAmp, '-o', 'Color', black, 'MarkerFaceColor', black, 'Linewidth', 2); hold on; plot(MovAmp, '-o', 'Color', dred, 'MarkerFaceColor', dred, 'Linewidth', 2); 
    plot(Caus_ClustAmp, '--o', 'Color', black, 'MarkerFaceColor', black, 'Linewidth', 2); 
    plot(NonClusteredAmp, '-o', 'Color', dred, 'MarkerFaceColor', dred, 'Linewidth', 2);
    plot(Caus_NonClusteredAmp, '--o', 'Color', dred, 'MarkerFaceColor', dred, 'Linewidth', 2);
    legend({'All', 'Mvmnt Rel.', 'Caus', 'Caus mov'})
    ylabel('Amp. of Clustered Spines');
    xlabel('Session')
    
    subplot(2,3,3); plot(ClustDendFreq, '-o','Color', gray, 'MarkerFaceColor', gray, 'Linewidth', 2); hold on; 
    plot(NoClustDendFreq, '-o', 'Color', dred, 'MarkerFaceColor', dred, 'Linewidth', 2)
    plot(CueClustDendFreq, '-o', 'Color', lgreen, 'MarkerFaceColor', lgreen, 'Linewidth', 2);
    plot(MovClustDendFreq, '-o', 'Color', black, 'MarkerFaceColor', black, 'Linewidth', 2); % plot(NoMovClustDendFreq, '-o', 'Color', lpurple, 'MarkerFaceColor', lpurple, 'Linewidth', 2); 
    plot(SucClustDendFreq, '-o', 'Color', lblue, 'MarkerFaceColor', lblue, 'Linewidth', 2);
    plot(RewClustDendFreq, '-o', 'Color', purple, 'MarkerFaceColor', purple, 'Linewidth', 2);
    legend({'Dends w/ Clusts', 'Dends no Clusts', 'Dends w/ CueClusts', 'Dends w/ MovClusts', 'Dends w/ SucClusts', 'Dends w/ RewClusts'})
    ylabel('Dendrite Frequency')
    xlabel('Session')
    
    subplot(2,3,4)
    for s = 1:length(DendNumClust)
        if ~isempty(DendNumClust{s})
            for p = 1:size(DendNumClust{s},1)
                plot(DendNumClust{s}(p,2), DendNumClust{s}(p,5), 'o', 'Color', black, 'MarkerFaceColor', black, 'Linewidth', 2); hold on;
            end
        end
    end
    xlabel('# of Clustered Spines')
    ylabel('Dendritic Freq.')
    
    subplot(2,3,5)
    for s = 1:length(DendNumClust)
        if ~isempty(DendNumClust{s})
            for p = 1:size(DendNumClust{s},1)
                plot(DendNumClust{s}(p,3), DendNumClust{s}(p,5), 'o', 'Color', black, 'MarkerFaceColor', black, 'Linewidth', 2); hold on;
            end
        end
    end
    xlabel('# of Mov-related Clustered Spines')
    ylabel('Dendritic Freq.')
    
    subplot(2,3,6)
    for s = 1:length(DendNumClust)
        if ~isempty(DendNumClust{s})
            for p = 1:size(DendNumClust{s},1)
                plot(DendNumClust{s}(p,4), DendNumClust{s}(p,5), 'o', 'Color', black, 'MarkerFaceColor', black, 'Linewidth', 2); hold on;
            end
        end
    end
    xlabel('# of Clusters on Dend')
    ylabel('Dendritic Freq.')
    
    %%
    %%%%%%%%%%%%%
    %%% Figure 3
    %%% Number of spine types
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure('Position', scrsz)
   
    sub1 = 3;
    sub2 = 3;
    
        subplot(sub1,sub2,1)
    plot(NumCueRelSpines, '-o', 'Color', lgreen, 'MarkerFaceColor', lgreen,'Linewidth', 2); hold on;
    plot(NumMovRelSpines, '-o', 'Color', black, 'MarkerFaceColor', black, 'Linewidth', 2);
    plot(cell2mat(cellfun(@nanmean, FractionofMovRelSpinesPerDendrite, 'uni', false)), '-o', 'Color', 'r', 'Linewidth', 2);
    plot(NumCueORMovRelSpines, '--o', 'Color', red, 'MarkerFaceColor', red, 'Linewidth', 2); 
    plot(NumPreSucRelSpines, '-o', 'Color', bgreen, 'MarkerFaceColor', bgreen, 'Linewidth', 2);
    plot(NumSucRelSpines, '-o', 'Color', lblue, 'MarkerFaceColor', lblue, 'Linewidth', 2);
    plot(NumMovDuringCueRelSpines, '-o', 'Color', green, 'MarkerFaceColor', green, 'Linewidth', 2);
    plot(NumRewRelSpines, '-o', 'Color', purple, 'MarkerFaceColor', purple, 'Linewidth', 2);
    plot(NumCausalMovSpines, '-o', 'Color', dred, 'MarkerFaceColor', dred, 'Linewidth', 2)
    legend({'Cue Rel', 'Mov Rel', 'MovORCue Rel', 'PreSuc Rel', 'Suc Rel', 'MovDuringCue','Rew Rel','Caus Mov Rel'}, 'Location', 'NorthEast');
    xlabel('Session','Fontsize',14);
    ylabel('Fraction of spines','Fontsize',14);
    title('Classes of Spines', 'Fontsize',14);
    
        subplot(sub1,sub2,2)
    plot(NumClustSpines, '-o','Color', black, 'MarkerFaceColor', black, 'Linewidth', 2); hold on;
    plot(NumCausClustSpines, '-o', 'Color', dred, 'MarkerFaceColor', dred, 'Linewidth', 2);
    plot(NumFarClustSpines, '-o', 'Color', gray, 'MarkerFaceColor', gray, 'Linewidth', 2);
    legend({'Clustered Spines', 'Caus Clustered Spines', 'Far Clustered Spines'}, 'Location', 'NorthEast');
    xlabel('Session','Fontsize',14);
    ylabel('Fraction of spines','Fontsize',14);
    title('Number of clustered spines', 'Fontsize', 14);
    
        subplot(sub1,sub2,3)
    plot(MeanNumberofSpinesinEachCluster, '-o', 'Color', purple, 'Linewidth', 2, 'MarkerFaceColor', purple); hold on
    plot(MeanNumberofSpinesinEachCausalCluster, '-o', 'Color', gray, 'Linewidth', 2, 'MarkerFaceColor', gray)
    plot(NumClusters, '-ok', 'Linewidth', 2, 'MarkerFaceColor', 'k')
    plot(NumCausalClusters, '-o', 'Color', lblue, 'Linewidth', 2, 'MarkerFaceColor', lblue)
    xlabel('Session', 'Fontsize', 14);
    ylabel('Raw Number', 'Fontsize', 14);
    legend({'Number of spines in each cluster', 'Num sp in each caus clust', 'Number of Clusters', 'Num Caus Clust'})
    title('Quantification of clusters', 'Fontsize', 14)
    
        subplot(sub1,sub2,4)
    plot(PercentCueRelDends, '-o', 'Color', lgreen, 'MarkerFaceColor', lgreen, 'Linewidth', 2); hold on;
    plot(PercentMovRelDends, '-o', 'Color', black, 'MarkerFaceColor', black,'Linewidth', 2)
    plot(PercentPreSucRelDends, '-o', 'Color', bgreen, 'MarkerFaceColor', bgreen, 'Linewidth', 2);
    plot(PercentSucRelDends, '-o', 'Color', lblue, 'MarkerFaceColor', lblue, 'Linewidth', 2)
    plot(PercentMovDuringCueRelDends, '-o', 'Color', green, 'MarkerFaceColor', green, 'Linewidth', 2)
    plot(PercentRewRelDends, '-o', 'Color', purple, 'MarkerFaceColor', purple, 'Linewidth', 2)
    xlabel('Session', 'Fontsize', 14)
    ylabel('Fraction of Dendrites', 'Fontsize', 14)
    legend({'Cue rel dends', 'Mov rel dends', 'Suc rel dends', 'Rew rel dends'})
    title([{'Fraction of dendrites that'}, {'are (function)-related'}], 'Fontsize', 14)
    
        subplot(sub1,sub2,5)
    plot(NumClustCueSpines, '-o','Color', lgreen, 'MarkerFaceColor', lgreen, 'Linewidth', 2); hold on;
        plot(NumCausClustCueSpines, '--o','Color', lgreen, 'MarkerFaceColor', lgreen, 'Linewidth', 2); hold on;
    plot(NumClustMovSpines, '-o','Color', black, 'MarkerFaceColor', black, 'LineWidth', 2)
        plot(NumCausClustMovSpines, '--o','Color', black, 'MarkerFaceColor', black, 'LineWidth', 2)
    plot(NumClustMixSpines, '-o','Color', red, 'MarkerFaceColor', red, 'Linewidth', 2);
        plot(NumCausClustMixSpines, '--o','Color', red, 'MarkerFaceColor', red, 'LineWidth', 2)
    plot(NumClustPreSucSpines, '-o', 'Color', bgreen, 'MarkerFaceColor', bgreen, 'Linewidth', 2)
        plot(NumCausClustPreSucSpines, '--o', 'Color', bgreen, 'MarkerFaceColor', bgreen, 'Linewidth', 2)
    plot(NumClustSucSpines, '-o','Color', lblue, 'MarkerFaceColor', lblue, 'Linewidth', 2)
        plot(NumCausClustSucSpines, '--o', 'Color', lblue, 'MarkerFaceColor', lblue, 'Linewidth', 2)
    plot(NumClustMovDuringCueSpines, '-o', 'Color', green, 'MarkerFaceColor', green, 'Linewidth', 2)
        plot(NumCausClustMovDuringCueSpines, '--o', 'Color', green, 'MarkerFaceColor', green, 'Linewidth', 2)
    plot(NumClustRewSpines, '-o','Color', purple, 'MarkerFaceColor', purple, 'Linewidth', 2)
        plot(NumCausClustRewSpines, '--o', 'Color', purple, 'MarkerFaceColor', purple, 'Linewidth', 2)
    xlabel('Session', 'Fontsize', 14)
    ylabel('Fraction of Spines', 'Fontsize',14)
    title('Classes of Clustered Spines', 'Fontsize', 14)
    legend({'Cue Clust', 'Mov Clust', 'Mix Clust','PreSuc Clust', 'Suc Clust', 'MovDuringCue Clust', 'Rew Clust'})
    
    
            subplot(sub1, sub2, 6)
    plot(MeanNumberofSpinesinEachMovCluster, '-o', 'Color', purple, 'Linewidth', 2, 'MarkerFaceColor', purple); hold on
    plot(NumMovClusters, '-ok', 'Linewidth', 2, 'MarkerFaceColor', 'k')
    xlabel('Session', 'Fontsize', 14);
    ylabel('Raw Number', 'Fontsize', 14);
    legend({'Number of spines in each mov cluster', 'Number of Clusters'})
    title('Quantification of movement clusters', 'Fontsize', 14)
    
            subplot(sub1,sub2,7)
    plot(NumFarClustCueSpines, '-o','Color', lgreen, 'MarkerFaceColor', lgreen, 'Linewidth', 2); hold on;
    plot(NumFarClustMovSpines, '-o','Color', black, 'MarkerFaceColor', black, 'LineWidth', 2)
    plot(NumFarClustMixSpines, '-o','Color', red, 'MarkerFaceColor', red, 'Linewidth', 2);
    plot(NumFarClustPreSucSpines, '-o', 'Color', bgreen, 'MarkerFaceColor', bgreen, 'Linewidth', 2)
    plot(NumFarClustSucSpines, '-o','Color', lblue, 'MarkerFaceColor', lblue, 'Linewidth', 2)
    plot(NumFarClustMovDuringCueSpines, '-o', 'Color', green, 'MarkerFaceColor', green, 'Linewidth', 2)
    plot(NumFarClustRewSpines, '-o','Color', purple, 'MarkerFaceColor', purple, 'Linewidth', 2)
    xlabel('Session', 'Fontsize', 14)
    ylabel('Fraction of Spines', 'Fontsize',14)
    title('Classes of Far-Clustered Spines', 'Fontsize', 14)
    legend({'Cue Clust', 'Mov Clust', 'Mix Clust','PreSuc Clust', 'Suc Clust', 'MovDuringCue Clust', 'Rew Clust'})
    
        subplot(sub1,sub2,8)
    FractionofCueSpinesThatAreClustered(FractionofCueSpinesThatAreClustered == Inf) = NaN;
    FractionofMovementSpinesThatAreClustered(FractionofMovementSpinesThatAreClustered == Inf) = NaN;
    FractionofPreSuccessSpinesThatAreClustered(FractionofPreSuccessSpinesThatAreClustered == Inf) = NaN;
    FractionofSuccessSpinesThatAreClustered(FractionofSuccessSpinesThatAreClustered == Inf) = NaN;
    FractionofMovementDuringCueSpinesThatAreClustered(FractionofMovementDuringCueSpinesThatAreClustered == Inf) = NaN;
    FractionofRewardSpinesThatAreClustered(FractionofRewardSpinesThatAreClustered == Inf) = NaN;
    
    plot(FractionofCueSpinesThatAreClustered, 'Color', lgreen, 'Linewidth', 2); hold on;
    plot(FractionofMovementSpinesThatAreClustered, 'Color', black, 'Linewidth', 2); 
    plot(FractionofPreSuccessSpinesThatAreClustered, 'Color', bgreen, 'Linewidth', 2);
    plot(FractionofSuccessSpinesThatAreClustered, 'Color', lblue, 'Linewidth', 2);
    plot(FractionofMovementDuringCueSpinesThatAreClustered, 'Color', green, 'Linewidth', 2);
    plot(FractionofRewardSpinesThatAreClustered, 'Color', purple, 'Linewidth', 2);
    xlabel('Session', 'Fontsize', 14)
    ylabel('Fraction of Spines', 'Fontsize',14)
    title([{'Fraction of (function) spines'},{'that are clustered'}], 'Fontsize', 14)
    legend({'Cue Clust', 'Mov Clust', 'PreSuc Clust', 'Suc Clust', 'MovDuringCue Clust', 'Rew Clust'})

    
    %%
    %%%%%%%%%%%%
    %%% Figure 4:
    %%%
    %%% Correlation of spine types with dendrites
    %%%%%%%%%%%%%
    
    figure('Position', scrsz);hold on;
    
    subplot(2,2,1)
    plot(cell2mat(cellfun(@nanmean, CorrelationofClusters, 'Uni', false)), 'Color', black, 'MarkerFaceColor', black, 'Linewidth', 2)
    
    xlabel('Session', 'Fontsize', 14)
    ylabel('Correlation of Clustered Spines', 'Fontsize', 14)
    title('Mean Correlation of Clustered Spines', 'Fontsize', 14)
    
    subplot(2,2,2)
    plot(ClusteredSpines_CorrwithDend, '-ok', 'LineWidth', 3, 'MarkerFaceColor', 'k'); hold on
    plot(FilteredClusteredSpines_CorrwithDend, '--o', 'Color', black, 'LineWidth', 2, 'MarkerFaceColor', orange);
    plot(CausalClusteredSpines_CorrwithDend, '-o', 'Color', red, 'Linewidth', 2, 'MarkerFaceColor', red);
    plot(FilteredCausalClusteredSpines_CorrwithDend, '--o', 'Color', red, 'Linewidth', 2, 'MarkerFaceColor', red);
    plot(NonClusteredSpines_CorrwithDend, '-o','Color', gray, 'LineWidth', 2, 'LineWidth', 2, 'MarkerFaceColor', gray)

    legend({'All Clustered', 'Filtered Clustered', 'Causal Clustered','Filtered Causal Clustered', 'Non-Clustered'})
    xlabel('Session')
    ylabel('Correlation with Dendrite', 'Fontsize', 14)
    title('Correlation of Different Spine Types with Dendritic Activity', 'Fontsize', 14)
    
    subplot(2,2,3)
    plot(CueRelClusteredSpines_CorrwithDend, '-o', 'Color', lgreen, 'Linewidth', 2, 'MarkerFaceColor', lgreen); hold on;
    plot(MovRelClusteredSpines_CorrwithDend, '-o', 'Color', black, 'LineWidth', 2, 'MarkerFaceColor', black);
    plot(PreSucRelClusteredSpines_CorrwithDend, '-o', 'Color', bgreen, 'Linewidth', 2, 'MarkerFaceColor', bgreen);
    plot(SucRelClusteredSpines_CorrwithDend, '-o', 'Color', lblue, 'Linewidth', 2, 'MarkerFaceColor', lblue)
    plot(MovDuringCueRelClusteredSpines_CorrwithDend, '-o', 'Color', green, 'Linewidth', 2, 'MarkerFaceColor', green);
    plot(RewRelClusteredSpines_CorrwithDend, '-o', 'Color', purple, 'Linewidth', 2, 'MarkerFaceColor', purple)
    legend({'Cue-rel Clusters', 'Mov-rel clusters', 'PreSuc','Suc-rel clusters', 'MovDuringCue', 'Rew-rel clusters'})
    xlabel('Session')
    ylabel('Correlation with Dendrite')
    title('Correlation of Different Spine Types with Dendritic Activity', 'Fontsize', 14)

    %%
    %%% Figure 5;
    %%%
    %%% Spectral analysis of clusters
    
    figure('Position', scrsz)
    MeanSpatialDegree = nan(1,14);
    MeanTemporalDegree = nan(1,14);
    MeanSpatioTemporalDegree = nan(1,14);
    DendriteClusteringDegree = nan(1,14);
    TemporalClusteringDegree = nan(1,14);
    SpatioTemporalClusteringDegree = nan(1,14);
    
    for i = 1:14
        MeanSpatialDegree(1,i) = nanmean(Dend_Spat_Deg{i});
        MeanTemporalDegree(1,i) = nanmean(Dend_Temp_Deg{i});
        MeanSpatioTemporalDegree(1,i) = nanmean(Dend_SpatTemp_Deg{i});
        if ~isempty(DendClust_Deg{i})
            DendriteClusteringDegree(1,i) = nanmean(DendClust_Deg{i}(:,1));
            TemporalClusteringDegree(1,i) = nanmean(DendClust_Deg{i}(:,2));
            SpatioTemporalClusteringDegree(1,i) = nanmean(DendClust_Deg{i}(:,3));
        else
            DendriteClusteringDegree(1,i) = nan;
            TemporalClusteringDegree(1,i) = nan;
            SpatioTemporalClusteringDegree(1,i) = nan;
        end
        MeanSpatialDegree_vs_Movement(1,i) = nanmean(SpatialDegree_vs_Movement{i});
        MeanTemporalDegree_vs_Movement(1,i) = nanmean(TemporalDegree_vs_Movement{i});
        MeanSpatioTemporalDegree_vs_Movement(1,i) = nanmean(SpatiotemporalDegree_vs_Movement{i});
        MeanSpatioTemporalOverlap(1,i) = nanmean(Spatiotemporal_Overlap{i});
    end
    subplot(2,4,1)
    plot(1:14, MeanSpatialDegree, '-ok', 'MarkerFaceColor', 'k'); hold on;
    plot(1:14, MeanTemporalDegree, '-or', 'MarkerFaceColor', 'r');
    plot(1:14, MeanSpatioTemporalDegree, '-og', 'MarkerFaceColor', 'g');
    legend({'Spatial', 'Temporal', 'Spatiotemporal'});
    xlabel('Session', 'Fontsize', 12)
    ylabel('Mean Degree', 'Fontsize', 12)
    title('Mean Graphical Degree of Dendrites', 'Fontsize', 14)
    
    subplot(2,4,2)
    plot(1:14,DendriteClusteringDegree,'-ok', 'MarkerFaceColor', 'k'); hold on;
    plot(1:14,TemporalClusteringDegree,'-or', 'MarkerFaceColor', 'r');
    plot(1:14,SpatioTemporalClusteringDegree, '-og', 'MarkerFaceColor', 'g');
    legend({'Spatial', 'Temporal'});
    
    ylabel('Mean algebraic connectivity of dendrites', 'Fontsize', 12)
    xlabel('Session', 'Fontsize', 12)
    title('Dendritic Fiedler Values', 'Fontsize', 14)
    
    subplot(2,4,3)
    plot(1:14, MeanSpatialDegree_vs_Movement, '-ok', 'MarkerFaceColor', 'k'); hold on;
    plot(1:14, MeanTemporalDegree_vs_Movement, '-or', 'MarkerFaceColor', 'r'); 
    plot(1:14, MeanSpatioTemporalDegree_vs_Movement, '-og', 'MarkerFaceColor', 'g');
    legend({'Spatial v Movment', 'Temporal v Movement', 'ST v Movement'})
    
    ylabel('Correlation of Degree and Movement Correlation')
    xlabel('Session')
    title('Correlation of Degree Measurements with Movement Correlation Profile')
    
    subplot(2,4,4)
    plot(1:14, MeanSpatioTemporalOverlap, '-ok', 'MarkerFaceColor', 'k');
    ylabel('Correlation')
    xlabel('Session')
    title('Correlation between Spatial and Temporal Profiles')
    
    subplot(2,4,5)
    plot(1:14, MeanSpatialDegreeofCueSpines, '-o', 'Color', lgreen, 'MarkerFaceColor', lgreen, 'Linewidth', 2); hold on;
    plot(1:14, MeanSpatialDegreeofMovementSpines, '-o', 'Color', black, 'MarkerFaceColor', black, 'Linewidth', 2)
    plot(1:14, MeanSpatialDegreeofMDCSpines, '-o', 'Color', green, 'MarkerFaceColor', green, 'Linewidth', 2)
    plot(1:14, MeanSpatialDegreeofPreSuccessSpines, '-o', 'Color', bgreen, 'MarkerFaceColor', bgreen , 'Linewidth', 2)
    plot(1:14, MeanSpatialDegreeofSuccessSpines, '-o', 'Color', lblue, 'MarkerFaceColor', lblue, 'Linewidth', 2)
    plot(1:14, MeanSpatialDegreeofRewardSpines, '-o', 'Color', purple, 'MarkerFaceColor', purple, 'Linewidth', 2);
    xlim([0 15])
    xlabel('Session', 'Fontsize', 14);
    title('Mean Spatial Degree of Spine Types')
    
    subplot(2,4,6)
    plot(1:14, MeanTemporalDegreeofCueSpines, '-o', 'Color', lgreen, 'MarkerFaceColor', lgreen, 'Linewidth', 2); hold on;
    plot(1:14, MeanTemporalDegreeofMovementSpines, '-o', 'Color', black, 'MarkerFaceColor', black, 'Linewidth', 2)
    plot(1:14, MeanTemporalDegreeofMDCSpines, '-o', 'Color', green, 'MarkerFaceColor', green, 'Linewidth', 2)
    plot(1:14, MeanTemporalDegreeofPreSuccessSpines, '-o', 'Color', bgreen, 'MarkerFaceColor', bgreen, 'Linewidth', 2)
    plot(1:14, MeanTemporalDegreeofSuccessSpines, '-o', 'Color', lblue, 'MarkerFaceColor', lblue, 'Linewidth', 2)
    plot(1:14, MeanTemporalDegreeofRewardSpines, '-o', 'Color', purple, 'MarkerFaceColor', purple, 'Linewidth', 2);
    xlim([0 15])
    xlabel('Session', 'Fontsize', 14);
    title('Mean Temporal Degree')
    
    subplot(2,4,7)
    plot(1:14, MeanSpatioTemporalDegreeofCueSpines, '-o', 'Color', lgreen, 'MarkerFaceColor', lgreen, 'Linewidth', 2); hold on;
    plot(1:14, MeanSpatioTemporalDegreeofMovementSpines, '-o', 'Color', black, 'MarkerFaceColor', black, 'Linewidth', 2)
    plot(1:14, MeanSpatioTemporalDegreeofMDCSpines, '-o', 'Color', green, 'MarkerFaceColor', green, 'Linewidth', 2)
    plot(1:14, MeanSpatioTemporalDegreeofPreSuccessSpines, '-o', 'Color', bgreen, 'MarkerFaceColor', bgreen, 'Linewidth', 2)
    plot(1:14, MeanSpatioTemporalDegreeofSuccessSpines, '-o', 'Color', lblue, 'MarkerFaceColor', lblue, 'Linewidth', 2)
    plot(1:14, MeanSpatioTemporalDegreeofRewardSpines, '-o', 'Color', purple, 'MarkerFaceColor', purple, 'Linewidth', 2);
    xlim([0 15])
    xlabel('Session', 'Fontsize', 14);
    title('Mean SpatioTemporal Degree');
    
    %%
    %%%%%%%%%%%%
    %%% Figure 6
    %%%
    %%% Length of Clusters
    %%%%%%%%%%%%%
    
    figure('Position', scrsz)
    subplot(1,4,1)
        plot(MeanCueClustLength, '-o', 'Color', lgreen, 'LineWidth', 2, 'MarkerFaceColor', lgreen); hold on;
        plot(MeanMovClustLength, '-o', 'Color', black, 'LineWidth', 2, 'MarkerFaceColor', black);
        plot(MeanMixClustLength, '-o', 'Color', red, 'Linewidth', 2, 'MarkerFaceColor', red)
        plot(MeanPreSucClustLength, '-o', 'Color', bgreen', 'Linewidth', 2, 'MarkerFaceColor', bgreen);
        plot(MeanSucClustLength, '-o', 'Color', lblue, 'Linewidth', 2, 'MarkerFaceColor', lblue);
        plot(MeanMovDuringCueClustLength, '-o', 'Color', green, 'Linewidth', 2, 'MarkerFaceColor', green);
        plot(MeanRewClustLength, '-o', 'Color', purple, 'Linewidth', 2, 'MarkerFaceColor', purple);
        plot(MeanAllClustLength, '-o', 'Color', gray, 'Linewidth', 2, 'MarkerFaceColor', gray);
        ylabel('Mean Distance b/w Spines in Clusters', 'Fontsize', 14);
        xlabel('Session','Fontsize', 14)
        xlim([0 15])
        legend({'Cue clust', 'Mov Clust', 'Mix Clust', 'Suc. Clust', 'Rew Clust', 'All Clust'}, 'Location', 'Northwest')
        title('Mean Cluster Length')
            pos = get(gca,'Position');
            axes('Position', [pos(1)+0.7*pos(3), pos(2)+0.7*pos(4), 0.25*pos(3), 0.25*pos(4)]);
            hist(cell2mat(AllMovClustLengths));
            
        
    subplot(1,4,2)
        plot(MaxCueClustLength, '-o', 'Color', lgreen, 'LineWidth', 2, 'MarkerFaceColor', lgreen); hold on;
        plot(MaxMovClustLength, '-o', 'Color', black, 'Linewidth', 2, 'MarkerFaceColor', black);
        plot(MaxMixClustLength, '-o', 'Color', red, 'Linewidth', 2, 'MarkerFaceColor', red)
        plot(MaxSucClustLength, '-o', 'Color', lblue, 'Linewidth', 2, 'MarkerFaceColor', lblue);
        plot(MaxRewClustLength, '-o', 'Color', purple, 'Linewidth', 2, 'MarkerFaceColor', purple);
        plot(MaxAllClustLength, '-o', 'Color', gray, 'Linewidth', 2, 'MarkerFaceColor', gray);
        ylabel('Max Distance b/w Spines in Clusters', 'Fontsize', 14);
        xlabel('Session', 'Fontsize', 14)
        xlim([0 15])
        legend({'Cue clust', 'Mov Clust', 'Mix Clust', 'Suc. Clust', 'Rew Clust', 'All Clust'})
        title('Max Cluster Length')
        
    subplot(1,4,3)
        plot(MeanDistanceBetweenCueSpines, '-o', 'Color', lgreen, 'Linewidth', 2, 'MarkerFaceColor', lgreen); hold on;
        plot(MeanDistanceBetweenMovementSpines, '-o', 'Color', black, 'Linewidth', 2, 'MarkerFaceColor', black)
        plot(MeanDistanceBetweenPreSuccessSpines, '-o', 'Color', bgreen, 'Linewidth', 2, 'MarkerFaceColor', bgreen)
        plot(MeanDistanceBetweenSuccessSpines, '-o', 'Color', lblue, 'Linewidth', 2, 'MarkerFaceColor', lblue)
        plot(MeanDistanceBetweenMovementDuringCueSpines, '-o', 'Color', green, 'Linewidth', 2, 'MarkerFaceColor', green)
        plot(MeanDistanceBetweenRewardSpines, '-o', 'Color', purple, 'Linewidth', 2, 'MarkerFaceColor', purple)
        ylabel('Distance (um)', 'Fontsize', 14)
        xlabel('Session', 'Fontsize', 14)
        xlim([0 15])
        title('Distance between (function) spines')
        
    subplot(1,4,4)
        plot(MeanFarCueClustLength, '-o', 'Color', lgreen, 'LineWidth', 2, 'MarkerFaceColor', lgreen); hold on;
        plot(MeanFarMovClustLength, '-o', 'Color', black, 'LineWidth', 2, 'MarkerFaceColor', black);
        plot(MeanFarMixClustLength, '-o', 'Color', red, 'Linewidth', 2, 'MarkerFaceColor', red)
        plot(MeanFarPreSucClustLength, '-o', 'Color', bgreen', 'Linewidth', 2, 'MarkerFaceColor', bgreen);
        plot(MeanFarSucClustLength, '-o', 'Color', lblue, 'Linewidth', 2, 'MarkerFaceColor', lblue);
        plot(MeanFarMovDuringCueClustLength, '-o', 'Color', green, 'Linewidth', 2, 'MarkerFaceColor', green);
        plot(MeanFarRewClustLength, '-o', 'Color', purple, 'Linewidth', 2, 'MarkerFaceColor', purple);
        plot(MeanAllFarClustLength, '-o', 'Color', gray, 'Linewidth', 2, 'MarkerFaceColor', gray);
        ylabel('Mean Distance b/w Spines in Clusters', 'Fontsize', 14);
        xlabel('Session','Fontsize', 14)
        xlim([0 15])
        legend({'Cue clust', 'Mov Clust', 'Mix Clust', 'Suc. Clust', 'Rew Clust', 'All Clust'})
        title('Mean Inter-dendrite cluster Length')
        
        
        
        %%% Figure 7
    figure('Position', scrsz); 
    try
        subplot(2,4,1)
        plot(DistanceBetweenMovementSpines{1}, CorrelationBetweenMovementSpines{1}, 'ok', 'MarkerFaceColor', 'k')
            xlim([0 100])
            ylim([-0.05 1])
            xlabel('Distance (\mum)', 'FontSize', 14)
            ylabel('Correlation', 'FontSize', 14)
            title('All Spines, Session 1')
        subplot(2,4,2)
        plot(DistanceBetweenFarMovementSpines{1}, CorrelationBetweenFarMovementSpines{1}, 'o', 'MarkerFaceColor', gray)
            xlim([0 100])
            ylim([-0.05 1])
            xlabel('Distance (\mum)', 'FontSize', 14)
            ylabel('Correlation', 'FontSize', 14)
            title('All Spines, Session 1')
        subplot(2,4,3)
        plot(DistanceBetweenMovementSpines{10}, CorrelationBetweenMovementSpines{10}, 'ok', 'MarkerFaceColor', 'k')
            xlim([0 100])
            ylim([-0.05 1])
            xlabel('Distance (\mum)', 'FontSize', 14)
            ylabel('Correlation', 'FontSize', 14)
            title('All Spines, Session 10')
        subplot(2,4,4)
        plot(DistanceBetweenFarMovementSpines{10}, CorrelationBetweenFarMovementSpines{10}, 'o', 'MarkerFaceColor', gray)
            xlim([0 100])
            ylim([-0.05 1])
            xlabel('Distance (\mum)', 'FontSize', 14)
            ylabel('Correlation', 'FontSize', 14)
            title('All Spines, Session 10')
        subplot(2,4,5)
        plot(DistanceBetweenAllSpines{1}, CorrelationBetweenAllSpines{1}, 'o', 'MarkerEdgeColor', dred, 'MarkerFaceColor', dred)
            xlim([0 100])
            ylim([-0.05 1])
            xlabel('Distance (\mum)', 'FontSize', 14)
            ylabel('Correlation', 'FontSize', 14)
            title('All Spines, Session 1')
        subplot(2,4,6)
        plot(DistanceBetweenFarSpines{1}, CorrelationBetweenFarSpines{1}, 'o', 'MarkerFaceColor', gray)
            xlim([0 100])
            ylim([-0.05 1])
            xlabel('Distance (\mum)', 'FontSize', 14)
            ylabel('Correlation', 'FontSize', 14)
            title('All Spines, Session 1')
        subplot(2,4,7)
        plot(DistanceBetweenAllSpines{10}, CorrelationBetweenAllSpines{10}, 'o', 'MarkerEdgeColor', dred, 'MarkerFaceColor', dred)
            xlim([0 100])
            ylim([-0.05 1])
            xlabel('Distance (\mum)', 'FontSize', 14)
            ylabel('Correlation', 'FontSize', 14)
            title('All Spines, Session 10')
        subplot(2,4,8)
        plot(DistanceBetweenFarSpines{10}, CorrelationBetweenFarSpines{10}, 'o', 'MarkerFaceColor', gray)
            xlim([0 100])
            ylim([-0.05 1])
            xlabel('Distance (\mum)', 'FontSize', 14)
            ylabel('Correlation', 'FontSize', 14)
            title('All Spines, Session 10')        
    catch
    end
        
    a.ClustCorrwithCue = clustcorrwithcue;
    a.NonClusteredCorrwithCue = nonclustcorrwithcue;
    a.ClustCorrwithMovementDuringCue = clustcorrwithMDC;
    a.NonClusteredCorrwithMovementDuringCue = nonclustcorrwithMDC;
    a.ClustCorrwithMovement = clustcorrwithmov;
    a.NonClusteredCorrwithMovement = nonclustcorrwithmov;
    a.ClustCorrwithSuccess = clustcorrwithsuc;
    a.NonClusteredCorrwithSuccess = nonclustcorrwithsuc;
    a.ClustCorrwithReward = clustcorrwithrew;
    a.NonClusteredCorrwithReward = nonclustcorrwithrew;
    a.CorrelationofClusters = CorrelationofClusters;
       
    a.CausalClustCorrwithCue = Caus_clustcorrwithcue;
    a.CausalNonClusteredCorrwithCue = Caus_nonclustcorrwithcue;
    a.CausalClustCorrwithMovementDuringCue = Caus_clustcorrwithMDC;
    a.CausalNonClusteredCorrwithMovementDuringCue = Caus_nonclustcorrwithMDC;
    a.CausalClustCorrwithMovement = Caus_clustcorrwithmov;
    a.CausalNonClusteredCorrwithMovement = Caus_nonclustcorrwithmov;
    a.CausalClustCorrwithSuccess = Caus_clustcorrwithsuc;
    a.CausalNonClusteredCorrwithSuccess = Caus_nonclustcorrwithsuc;
    a.CausalClustCorrwithReward = Caus_clustcorrwithrew;
    a.CausalNonClusteredCorrwithReward = Caus_nonclustcorrwithrew;
    
    a.CueClustersCorrwithCue = cueclustcorrwithcue;
    a.CueNonClustCorrwithCue = cuenonclustcorrwithcue;
    a.MovementDuringCueClustersCorrwithMovementDuringCue = MDCclustcorrwithMDC;
    a.MovementDuringCueNonClustCorrwithMovementDuringCue = MDCnonclustcorrwithMDC;
    a.MovementClustersCorrwithMovement = movclustcorrwithmov;
    a.MovementNonClustCorrwithMovement = movnonclustcorrwithmov;
    a.SuccessClustersCorrwithSuccess = succlustcorrwithsuc;
    a.SuccessNonClustCorrwithSuccess = sucnonclustcorrwithsuc;
    a.RewardClustersCorrwithReward = rewclustcorrwithrew;
    a.RewardNonClustCorrwithReward = rewnonclustcorrwithrew;
    
    a.CausalCueClustersCorrwithCue = Caus_cueclustcorrwithcue;
    a.CausalCueNonClustCorrwithCue = Caus_cuenonclustcorrwithcue;
    a.CausalMovementDuringCueClustersCorrwithMovementDuringCue = Caus_MDCclustcorrwithMDC;
    a.CausalMovementDuringCueNonClustCorrwithMovementDuringCue = Caus_MDCnonclustcorrwithMDC;
    a.CausalMovementClustersCorrwithMovement = Caus_movclustcorrwithmov;
    a.CausalMovementNonClustCorrwithMovement = Caus_movnonclustcorrwithmov;
    a.CausalSuccessClustersCorrwithSuccess = Caus_succlustcorrwithsuc;
    a.CausalSuccessNonClustCorrwithSuccess = Caus_sucnonclustcorrwithsuc;
    a.CausalRewardClustersCorrwithReward = Caus_rewclustcorrwithrew;
    a.CausalRewardNonClustCorrwithReward = Caus_rewnonclustcorrwithrew;
    
    a.FractionofCueSpinesThatAreClustered = FractionofCueSpinesThatAreClustered;
    a.FractionofCueSpinesThatAreNonClustered = FractionofCueSpinesThatAreNonClustered;
    a.FractionofMovementSpinesThatAreClustered = FractionofMovementSpinesThatAreClustered;
    a.FractionofMovementSpinesThatAreNonClustered = FractionofMovementSpinesThatAreNonClustered;
    a.FractionofPreSuccessSpinesThatAreClustered = FractionofPreSuccessSpinesThatAreClustered; 
    a.FractionofSuccessSpinesThatAreClustered = FractionofSuccessSpinesThatAreClustered;
    a.FractionofSuccessSpinesThatAreNonClustered = FractionofSuccessSpinesThatAreNonClustered;
    a.FractionofMovementDuringCueSpinesThatAreClustered = FractionofMovementDuringCueSpinesThatAreClustered;
    a.FractionofRewardSpinesThatAreClustered = FractionofRewardSpinesThatAreClustered;
    a.FractionofRewardSpinesThatAreNonClustered = FractionofRewardSpinesThatAreNonClustered;
    
    a.ClusterFrequency = cluster_freq;
    a.NonClusteredFrequency = nonclustered_freq;
    a.CueClusterFrequency = cue_cluster_freq;
    a.MovementClusterFrequency = mov_cluster_freq;
    a.MovementDuringCueClusterFrequency = movduringcue_cluster_freq;
    a.PreSuccessClusterFrequency = presuc_cluster_freq;
    a.SuccessClusterFrequency = suc_cluster_freq;
    a.RewardClusterFrequency = rew_cluster_freq;
    a.CausalClusterFrequency = Caus_cluster_freq;
    a.NonClusteredCausalFrequency = Caus_nonclustered_freq;
    a.CausalCueClusterFrequency = Caus_cue_cluster_freq;
    a.CausalMovementClusterFrequency = Caus_mov_cluster_freq;
    a.CausalMovementDuringCueClusterFrequency = Caus_movduringcue_cluster_freq;
    a.CausalPreSuccessClusterFrequency = Caus_presuc_cluster_freq;
    a.CausalSuccessClusterFrequency = Caus_suc_cluster_freq;
    a.CausalRewardClusterFrequency = Caus_rew_cluster_freq;
    
    a.DendriteswithMovementClusters = DendswithMovClusts;
    a.DendriteswithoutMovementClusters = DendsnomovClusts;
    a.DendriteswithClustersFrequency = ClustDendFreq;
    a.DendriteswithoutClustersFrequency = NoClustDendFreq;
    a.DendriteswithCueClustersFrequency = CueClustDendFreq;
    a.DendriteswithMovClustersFrequency = MovClustDendFreq;
    a.DendriteswithMovDuringCueClustersFrequency = MovDuringCueClustDendFreq;
    a.DendriteswithPreSucClustersFrequency = PreSucClustDendFreq;
    a.DendriteswithSucClustersFrequency = SucClustDendFreq;
    a.DendriteswithRewClustersFrequency = RewClustDendFreq;
    a.DendriteswithoutMovClustersFrequency = NoMovClustDendFreq;
    a.DendriticClusterNumberInfo = DendNumClust;
    
    a.ClusteredSpineAmp = ClustAmp;
    a.NonClusteredSpineAmp = NonClusteredAmp;
    a.ClusteredCueSpineAmp = CueAmp;
    a.ClusteredMovSpineAmp = MovAmp;
    a.ClusteredMovDuringCueSpineAmp = MovDuringCueAmp;
    a.ClusteredPreSuccessSpineAmp = PreSucAmp;
    a.ClusteredSuccessSpineAmp = SucAmp;
    a.ClusteredRewardSpineAmp = RewAmp;
    a.CausalClusteredSpineAmp = Caus_ClustAmp;
    a.CausalNonClusteredSpineAmp = Caus_NonClusteredAmp;
    a.CausalClusteredCueSpineAmp = Caus_CueAmp;
    a.CausalClusteredMovSpineAmp = Caus_MovAmp;
    a.CausalClusteredMovDuringCueSpineAmp = Caus_MovDuringCueAmp;
    a.CausalClusteredPreSuccessSpineAmp = Caus_PreSucAmp;
    a.CausalClusteredSuccessSpineAmp = Caus_SucAmp;
    a.CausalClusteredRewardSpineAmp = Caus_RewAmp;
    
    
    a.NumberofCueSpines = NumCueRelSpines;
    a.NumberofMovementRelatedSpines = NumMovRelSpines;
    a.FractionofMovementRelatedSpinesPerDendrite = FractionofMovRelSpinesPerDendrite;
    a.NumberofCueORMovementRelatedSpines = NumCueORMovRelSpines;
    a.NumberofPreSuccessSpines = NumPreSucRelSpines;
    a.NumberofSuccessSpines = NumSucRelSpines;
    a.NumberofMovementDuringCueSpines = NumMovDuringCueRelSpines;
    a.NumberofRewardSpines = NumRewRelSpines;
    a.NumberofCausalCueSpines = NumCausalCueSpines;
    a.NumberofCausalMvmntSpines = NumCausalMovSpines;
    a.NumberofCausalSuccessSpines = NumCausalSucSpines;
    a.NumberofCausalRewardSpines = NumCausalRewSpines;
    
    a.NumberofClusteredSpines = NumClustSpines;
    a.NumberofClusteredCueSpines = NumClustCueSpines;
    a.NumberofClusteredMovementSpines = NumClustMovSpines;
    a.NumberofClusteredMixedFunctionSpines = NumClustMixSpines;
    a.NumberofClusteredPreSuccessSpines = NumClustPreSucSpines;
    a.NumberofClusteredSuccessSpines = NumClustSucSpines;
    a.NumberofClusteredMovementDuringCueSpines = NumClustMovDuringCueSpines;
    a.NumberofClusteredRewardSpines = NumClustRewSpines;
    
    a.NumberofFarClusteredSpines = NumFarClustSpines;
    a.NumberofFarClusteredCueSpines = NumFarClustCueSpines;
    a.NumberofFarClusteredMovementSpines = NumFarClustMovSpines;
    a.NumberofFarClusteredMixedFunctionSpines = NumFarClustMixSpines;
    a.NumberofFarClusteredPreSuccessSpines = NumFarClustPreSucSpines;
    a.NumberofFarClusteredSuccessSpines = NumFarClustSucSpines;
    a.NumberofFarClusteredMovementDuringCueSpines = NumFarClustMovDuringCueSpines;
    a.NumberofFarClusteredRewardSpines = NumFarClustRewSpines;
    
    a.NumberofCausalClusteredSpines = NumCausClustSpines;
    a.NumberofCausalClusteredCueSpines = NumCausClustCueSpines;
    a.NumberofCausalClusteredMovementDuringCueSpines = NumCausClustMovDuringCueSpines;
    a.NumberofCausalClusteredMovementSpines = NumCausClustMovSpines;
    a.NumberofCausalClusteredPreSuccessSpines = NumCausClustPreSucSpines;
    a.NumberofCausalClusteredSuccessSpines = NumCausClustSucSpines;
    a.NumberofCausalClusteredRewardSpines = NumCausClustRewSpines;
    
    a.MeanCueClustLength = MeanCueClustLength;
    a.MaxCueClustLength = MaxCueClustLength;
    a.MeanMovClustLength = MeanMovClustLength;
    a.AllMovClustLengths = AllMovClustLengths;
    a.MaxMovClustLength = MaxMovClustLength;
    a.MeanMixClustLength = MeanMixClustLength;
    a.MaxMixClustLength = MaxMixClustLength;
    a.MeanPreSucClustLength = MeanPreSucClustLength;
    a.MaxPreSucClustLength = MaxPreSucClustLength;
    a.MeanSucClustLength = MeanSucClustLength;
    a.MaxSucClustLength = MaxSucClustLength;
    a.MeanMovDuringCueClustLength = MeanMovDuringCueClustLength;
    a.MaxMovDuringCueClustLength = MaxMovDuringCueClustLength;
    a.MeanRewClustLength = MeanRewClustLength;
    a.MaxRewClustLength = MaxRewClustLength;
    a.MeanAllClustLength = MeanAllClustLength;
    a.MaxAllClustLength = MaxAllClustLength;
    
    a.MeanFarCueClustLength = MeanFarCueClustLength;
    a.MeanFarMovClustLength = MeanFarMovClustLength;
    a.MeanFarMixClustLength = MeanFarMixClustLength;
    a.MeanFarPreSucClustLength = MeanFarPreSucClustLength;
    a.MeanFarSucClustLength = MeanFarSucClustLength;
    a.MeanFarMovDuringCueClustLength = MeanFarMovDuringCueClustLength;
    a.MeanFarRewClustLength = MeanFarRewClustLength;
    a.MeanAllFarClustLength = MeanAllFarClustLength;
    
    a.MeanCausalCueClustLength = MeanCausalCueClustLength;
    a.MaxCausalCueClustLength = MaxCausalCueClustLength;
    a.MeanCausalMovClustLength = MeanCausalMovClustLength;
    a.MaxCausalMovClustLength = MaxCausalMovClustLength;
    a.MeanCausalSucClustLength = MeanCausalSucClustLength;
    a.MaxCausalSucClustLength = MaxCausalSucClustLength;
    a.MeanCausalRewClustLength = MeanCausalRewClustLength;
    a.MaxCausalRewClustLength = MaxCausalRewClustLength;
    a.MeanCausalAllClustLength = MeanAllCausalClustLength;
    a.MaxCausalAllClustLength = MaxAllCausalClustLength;
    
    a.NearestMovementRelatedSpine = NearestMovSpine;
    a.NextClosestMovementRelatedSpine = NextClosest;
    a.ThirdClosestMovementRelatedSpine = ThirdClosest;
    a.FourthClosestMovementRelatedSpine = FourthClosest;
    a.CorrelationwithNearestMovementRelatedSpine = CorrwithNearestMovSpine; 
    a.NearestHighlyCorrelatedMovementRelatedSpine = NearestHighCorrMovSpine;
    a.NextClosestHighlyCorrelatedMovementRelatedSpine = NextClosestHighCorrMovSpine;
    a.ThirdClosestHighlyCorrelatedMovementRelatedSpine = ThirdClosestHighCorrMovSpine;
    a.DistanceBetweenAllSpines = DistanceBetweenAllSpines;
        a.CorrelationBetweenAllSpines = CorrelationBetweenAllSpines;
        a.CorrelationBetweenAllSpinesMovementPeriods = CorrelationBetweenAllSpinesMovePeriods;
        a.CorrelationBetweenAllSpinesStillPeriods = CorrelationBetweenAllSpinesStillPeriods;
        a.MeanCorrelationBetweenAllSpines = MeanCorrelationBetweenAllSpines;
    a.DistanceBetweenCueSpines = DistanceBetweenCueSpines;
    a.MeanDistanceBetweenCueSpines = MeanDistanceBetweenCueSpines;
    a.DistanceBetweenMovementSpines = DistanceBetweenMovementSpines;
    a.MeanDistanceBetweenMovementSpines = MeanDistanceBetweenMovementSpines;
        a.CorrelationBetweenMovementSpines = CorrelationBetweenMovementSpines;
        a.CorrelationBetweenMovementSpinesMovementPeriods = CorrelationBetweenMovementSpinesMovePeriods;
        a.CorrelationBetweenMovementSpinesStillPeriods = CorrelationBetweenMovementSpinesStillPeriods;
        a.MeanCorrelationBetweenMovementSpines = MeanCorrelationBetweenMovementSpines;
    a.DistanceBetweenPreSuccessSpines = DistanceBetweenPreSuccessSpines;
    a.MeanDistanceBetweenPreSuccessSpines = MeanDistanceBetweenPreSuccessSpines;
    a.DistanceBetweenSuccessSpines = DistanceBetweenSuccessSpines;
    a.MeanDistanceBetweenSuccessSpines = MeanDistanceBetweenSuccessSpines;
    a.DistanceBetweenMovementDuringCueSpines = DistanceBetweenMovementDuringCueSpines;
    a.MeanDistanceBetweenMovementDuringCueSpines = MeanDistanceBetweenMovementDuringCueSpines;
    a.DistanceBetweenRewardSpines = DistanceBetweenRewardSpines;
    a.MeanDistanceBetweenRewardSpines = MeanDistanceBetweenRewardSpines;
    a.DistanceBetweenFarSpines = DistanceBetweenFarSpines;
    a.CorrelationBetweenFarSpines = CorrelationBetweenFarSpines;
    a.DistanceBetweenFarMovementSpines = DistanceBetweenFarMovementSpines;
    a.CorrelationBetweenFarMovementSpines = CorrelationBetweenFarMovementSpines;
    
    a.MeanNumberofSpinesinEachCluster = MeanNumberofSpinesinEachCluster;
    a.MeanNumberofSpinesinEachCausalCluster = MeanNumberofSpinesinEachCausalCluster;
    a.NumberofClusters = NumClusters;
    a.NumberofCausalClusters = NumCausalClusters;
    a.MeanNumberofSpinesinEachMovCluster = MeanNumberofSpinesinEachMovCluster;
    a.NumberofMovClusters = NumMovClusters;
    
    a.ClusteredSpines_CorrwithDend = ClusteredSpines_CorrwithDend;
    a.FilteredClusteredSpines_CorrwithDend = FilteredClusteredSpines_CorrwithDend;
    a.NonClusteredSpines_CorrwithDend = NonClusteredSpines_CorrwithDend;
    a.CausalClusteredSpines_CorrwithDend = CausalClusteredSpines_CorrwithDend;
    a.FilteredCausalClusteredSpines_CorrwithDend = FilteredCausalClusteredSpines_CorrwithDend;
    a.NonCausalClusteredSpines_CorrwithDend = NonCausalClusteredSpines_CorrwithDend;
    a.CueRelClusteredSpines_CorrwithDend = CueRelClusteredSpines_CorrwithDend;
    a.MovRelClusteredSpines_CorrwithDend = MovRelClusteredSpines_CorrwithDend;
    a.PreSucRelClusteredSpines_CorrwithDend = PreSucRelClusteredSpines_CorrwithDend;
    a.SucRelClusteredSpines_CorrwithDend = SucRelClusteredSpines_CorrwithDend;
    a.MovDuringCueRelClusteredSpines_CorrwithDend = MovDuringCueRelClusteredSpines_CorrwithDend;
    a.RewRelClusteredSpines_CorrwithDend = RewRelClusteredSpines_CorrwithDend;
    a.CueRelCausalClusteredSpines_CorrwithDend = CueRelCausalClusteredSpines_CorrwithDend;
    a.MovRelCausalClusteredSpines_CorrwithDend = MovRelCausalClusteredSpines_CorrwithDend;
    a.SucRelCausalClusteredSpines_CorrwithDend = SucRelCausalClusteredSpines_CorrwithDend;
    a.RewRelCausalClusteredSpines_CorrwithDend = RewRelCausalClusteredSpines_CorrwithDend;

    
    a.SpatialDegree = Spatial_Deg;
    a.MeanDendriticSpatialDegree = MeanSpatialDegree;
    a.TemporalDegree = Temporal_Deg;
    a.MeanDendriticTemporalDegree = MeanTemporalDegree;
    a.SpatioTemporalDegree = Spatiotemporal_Deg;
    a.MeanDendriticSpatioTemporalDegree = MeanSpatioTemporalDegree;
    a.SpatioTemporalOverlap = Spatiotemporal_Overlap;
    a.SpectralMovementCorrelation = MovementCorrelationsforAllSpinesonDend;
    a.SpatialDegreevsMovement = SpatialDegree_vs_Movement;
    a.TemporalDegreevsMovement = TemporalDegree_vs_Movement;
    a.SpatioTemporalDegreevsMovement = SpatiotemporalDegree_vs_Movement;
    a.SpatialDegreevsMovement = SpatialDegree_vs_Movement;
    a.TemporalDegreevsMovement = TemporalDegree_vs_Movement;
    a.SpectralDendriteInformation = DendClust_Deg;                  %%% Col 1-3 = Eigenvalue/Fiedler value, col 4 = dend freq., col 5 = dend corr with movment, col 6 = number of spines
    a.SpatioTemporalFiedlerValues = SpatioTemporalFiedler;
    a.SpatioTemporalPartitionVector  = SpatioTemporalPartition;     %%% The Fiedler vector can be used as a partition for the graph (positive in one group, negative in another)
    a.SpatioTemporal_FirstEigenVector = SpatioTemporal_FirstEigenvector;
    a.MeanSpatialDegreeofCueSpines = MeanSpatialDegreeofCueSpines;
    a.MeanTemporalDegreeofCueSpines = MeanTemporalDegreeofCueSpines;
    a.MeanSpatioTemporalDegreeofCueSpines = MeanSpatioTemporalDegreeofCueSpines;
    a.MeanSpatialDegreeofMovementSpines = MeanSpatialDegreeofMovementSpines;
    a.MeanTemporalDegreeofMovementSpines = MeanTemporalDegreeofMovementSpines;
    a.MeanSpatioTemporalDegreeofMovementSpines = MeanSpatioTemporalDegreeofMovementSpines;
    a.MeanSpatialDegreeofMovementDuringCueSpines = MeanSpatialDegreeofMDCSpines;
    a.MeanTemporalDegreeofMovementDuringCueSpines = MeanTemporalDegreeofMDCSpines;
    a.MeanSpatioTemporalDegreeofMovementDuringCueSpines = MeanSpatioTemporalDegreeofMDCSpines;
    a.MeanSpatialDegreeofPreSuccessSpines = MeanSpatialDegreeofPreSuccessSpines;
    a.MeanTemporalDegreeofPreSuccessSpines = MeanTemporalDegreeofPreSuccessSpines;
    a.MeanSpatioTemporalDegreeofPreSuccessSpines = MeanSpatioTemporalDegreeofPreSuccessSpines;
    a.MeanSpatialDegreeofSuccessSpines = MeanSpatialDegreeofSuccessSpines;
    a.MeanTemporalDegreeofSuccessSpines = MeanTemporalDegreeofSuccessSpines;
    a.MeanSpatioTemporalDegreeofSuccessSpines = MeanSpatioTemporalDegreeofSuccessSpines;
    a.MeanSpatialDegreeofRewardSpines = MeanSpatialDegreeofRewardSpines;
    a.MeanTemporalDegreeofRewardSpines = MeanTemporalDegreeofRewardSpines;
    a.MeanSpatioTemporalDegreeofRewardSpines = MeanSpatioTemporalDegreeofRewardSpines;

%     a.AverageFractionofaClusterThatsMovementRelated = FractionofClusterThatsMR;
%     a.AverageFractionofaCausalClusterThatsMovementRelated = FractionofCausalClusterThatsMR;
    a.PercentofCueRelatedDendrites = PercentCueRelDends;
    a.PercentofMovementRelatedDendrites = PercentMovRelDends;
    a.PercentofPreSuccessRelatedDendrites = PercentPreSucRelDends;
    a.PercentofSuccessRelatedDendrites = PercentSucRelDends;
    a.PercentofMovementDuringCueRelatedDendrites = PercentMovDuringCueRelDends;
    a.PercentofRewardRelatedDendrites = PercentRewRelDends;
    
    fname = inputname(1);
    fname = fname(1:5);
    fname = [fname, '_SpineCorrelationTimecourse'];
    eval([fname, '= a;'])
    cd('C:\Users\Komiyama\Desktop\Output Data');
    save(fname, fname);
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Averaging %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%
    %%% Collect Data from input
    %%%
    
    FractionofMovementRelatedSpinesPerDendrite = cell(1,14);
    SpatialDegree = cell(1,14);
    TemporalDegree = cell(1,14);
    SpatioTemporalDegree = cell(1,14);
    SpatialMovementCorrelation = cell(1,14);
    TemporalMovementCorrelation = cell(1,14);
    SpatioTemporalMovementCorrelation = cell(1,14);
    DendClusteringDegree = cell(1,14);
    SpatioTemporalOverlap = cell(1,14);
    CorrelationofClusters = cell(1,14);
    AllDistancesBetweenAllSpines = cell(1,14);
    AllMovClustLengths = cell(1,14);
        AllDistancesBetweenFarSpines = cell(1,14);
    CorrelationBetweenAllSpines = cell(1,14);
    CorrelationBetweenAllSpinesMovePeriods = cell(1,14);
    CorrelationBetweenAllSpinesStillPeriods = cell(1,14);
        CorrelationBetweenFarSpines = cell(1,14);
    AllDistancesBetweenMovementSpines = cell(1,14);
        AllDistancesBetweenFarMovementSpines = cell(1,14);
    CorrelationBetweenMovementSpines = cell(1,14);
    CorrelationBetweenMovementSpinesMovePeriods = cell(1,14);
    CorrelationBetweenMovementSpinesStillPeriods = cell(1,14);
        CorrelationBetweenFarMovementSpines = cell(1,14);
    NearestMovementRelatedSpine = cell(1,14);
    NextClosestMovementRelatedSpine= cell(1,14);
    ThirdClosestMovementRelatedSpine = cell(1,14);
    FourthClosestMovementRelatedSpine = cell(1,14);
    CorrelationwithNearestMovementRelatedSpine = cell(1,14);
    NearestHighlyCorrelatedMovementRelatedSpine = cell(1,14);
    NextClosestHighlyCorrelatedMovementRelatedSpine = cell(1,14);
    ThirdClosestHighlyCorrelatedMovementRelatedSpine = cell(1,14);
    
    for i = 1:length(varargin)
        AllClustersCorrwithCue(i,1:14) = varargin{i}.ClustCorrwithCue;
        NonClusteredCorrwithCue(i,1:14) = varargin{i}.NonClusteredCorrwithCue;
        AllClustersCorrwithMDC(i,1:14) = varargin{i}.ClustCorrwithMovementDuringCue;
        NonClusteredCorrwithMDC(i,1:14) = varargin{i}.NonClusteredCorrwithMovementDuringCue;
        AllClustersCorrwithMovement(i,1:14) = varargin{i}.ClustCorrwithMovement;
        NonClusteredCorrwithMovement(i,1:14) = varargin{i}.NonClusteredCorrwithMovement;
        AllClustersCorrwithSuccess(i,1:14) = varargin{i}.ClustCorrwithSuccess;
        NonClusteredCorrwithSuccess(i,1:14) = varargin{i}.NonClusteredCorrwithSuccess;
        AllClustCorrwithReward(i,1:14) = varargin{i}.ClustCorrwithReward;
        NonClusteredCorrwithReward(i,1:14) = varargin{i}.NonClusteredCorrwithReward;
        

        AllCausalClustersCorrwithCue(i,1:14) = varargin{i}.CausalClustCorrwithCue;
        CausalNonClusteredCorrwithCue(i,1:14) = varargin{i}.CausalNonClusteredCorrwithCue;
        AllCausalClustersCorrwithMovement(i,1:14) = varargin{i}.CausalClustCorrwithMovement;
        CausalNonClusteredCorrwithMovement(i,1:14) = varargin{i}.CausalNonClusteredCorrwithMovement;
        AllCausalClustersCorrwithMDC(i,1:14) = varargin{i}.CausalClustCorrwithMovementDuringCue;
        CausalNonClusteredCorrwithMDC(i,1:14) = varargin{i}.CausalNonClusteredCorrwithMovementDuringCue;
        AllCausalClustersCorrwithSuccess(i,1:14) = varargin{i}.CausalClustCorrwithSuccess;
        CausalNonClusteredCorrwithSuccess(i,1:14) = varargin{i}.CausalNonClusteredCorrwithSuccess;
        AllCausalClustCorrwithReward(i,1:14) = varargin{i}.CausalClustCorrwithReward;
        CausalNonClusteredCorrwithReward(i,1:14) = varargin{i}.CausalNonClusteredCorrwithReward;
        
        CueRelatedClustersCorrwithCue(i,1:14) = varargin{i}.CueClustersCorrwithCue;
            CueRelatedNonClusteredCorrwithCue(i,1:14) = varargin{i}.CueNonClustCorrwithCue;
        MDCRelatedClustersCorrwithMDC(i,1:14) = varargin{i}.MovementDuringCueClustersCorrwithMovementDuringCue;
            MDCRelatedNonClusteredCorrwithMDC(i,1:14) = varargin{i}.MovementDuringCueNonClustCorrwithMovementDuringCue;
        MovementRelatedClustersCorrwithMovement(i,1:14) = varargin{i}.MovementClustersCorrwithMovement;
            MovementRelatedNonClusteredCorrwithMovement(i,1:14) = varargin{i}.MovementNonClustCorrwithMovement;
        SuccessRelatedClustersCorrwithSuccess(i,1:14) = varargin{i}.SuccessClustersCorrwithSuccess;
            SuccessRelatedNonClusteredCorrwithSuccess(i,1:14) = varargin{i}.SuccessNonClustCorrwithSuccess;
        RewardRelatedClustersCorrwithReward(i,1:14) = varargin{i}.RewardClustersCorrwithReward; 
            RewardRelatedNonClusteredCorrwithReward(i,1:14) = varargin{i}.RewardNonClustCorrwithReward; 
            
        CausalCueRelatedClustersCorrwithCue(i,1:14) = varargin{i}.CausalCueClustersCorrwithCue;
            CausalCueRelatedNonClusteredCorrwithCue(i,1:14) = varargin{i}.CausalCueNonClustCorrwithCue;
        CausalMDCRelatedClustersCorrwithMDC(i,1:14) = varargin{i}.CausalMovementDuringCueClustersCorrwithMovementDuringCue;
            CausalMDCRelatedNonClusteredCorrwithMDC(i,1:14) = varargin{i}.CausalMovementDuringCueNonClustCorrwithMovementDuringCue;
        CausalMovementRelatedClustersCorrwithMovement(i,1:14) = varargin{i}.CausalMovementClustersCorrwithMovement;
            CausalMovementRelatedNonClusteredCorrwithMovement(i,1:14) = varargin{i}.CausalMovementNonClustCorrwithMovement;
        CausalSuccessRelatedClustersCorrwithSuccess(i,1:14) = varargin{i}.CausalSuccessClustersCorrwithSuccess;
            CausalSuccessRelatedNonClusteredCorrwithSuccess(i,1:14) = varargin{i}.CausalSuccessNonClustCorrwithSuccess;
        CausalRewardRelatedClustersCorrwithReward(i,1:14) = varargin{i}.CausalRewardClustersCorrwithReward; 
            CausalRewardRelatedNonClusteredCorrwithReward(i,1:14) = varargin{i}.CausalRewardNonClustCorrwithReward; 

        if any(cell2mat(cellfun(@isempty, varargin{i}.MeanCorrelationBetweenMovementSpines, 'Uni', false)))
            varargin{i}.MeanCorrelationBetweenMovementSpines(cell2mat(cellfun(@isempty, varargin{i}.MeanCorrelationBetweenMovementSpines, 'Uni', false))) = {NaN};
        end
        MeanCorrelationBetweenMovementSpines(i,1:14) = cell2mat(varargin{i}.MeanCorrelationBetweenMovementSpines);
                
        FractionofCueSpinesThatAreClustered(i,1:14) = varargin{i}.FractionofCueSpinesThatAreClustered;
        FractionofMovementSpinesThatAreClustered(i,1:14) = varargin{i}.FractionofMovementSpinesThatAreClustered;
        FractionofPreSuccessSpinesThatAreClustered(i,1:14) = varargin{i}.FractionofPreSuccessSpinesThatAreClustered;
        FractionofSuccessSpinesThatAreClustered(i,1:14) = varargin{i}.FractionofSuccessSpinesThatAreClustered;
        FractionofMovementDuringCueSpinesThatAreClustered(i,1:14) = varargin{i}.FractionofMovementDuringCueSpinesThatAreClustered;
        FractionofRewardSpinesThatAreClustered(i,1:14) = varargin{i}.FractionofRewardSpinesThatAreClustered;
        
    
        for j = 1:14    %%% Sessions
            for k = 1:length(varargin{i}.SpatioTemporalDegree{j})
%                 SpatialDegree{j} = [SpatialDegree{j}; varargin{i}.SpectralDendriteInformation{j}
                SpatioTemporalDegree{j} = [SpatioTemporalDegree{j}; varargin{i}.SpatioTemporalDegree{j}{k}];
                SpatialMovementCorrelation{j} = [SpatialMovementCorrelation{j}; varargin{i}.SpatialDegreevsMovement{j}'];
                TemporalMovementCorrelation{j} = [TemporalMovementCorrelation{j}; varargin{i}.TemporalDegreevsMovement{j}'];
                SpatioTemporalMovementCorrelation{j} = [SpatioTemporalMovementCorrelation{j}; varargin{i}.SpatioTemporalDegreevsMovement{j}'];
            end
            if ~isempty(varargin{i}.SpectralDendriteInformation{j})
                DendClusteringDegree{j} = [DendClusteringDegree{j}; varargin{i}.SpectralDendriteInformation{j}];
            end
%             AssociatedDendFreq{i} = [AssociatedDendFreq{i}; varargin{i}.SpectralDendriteInformation{j}];
%             AssociatedDendCorr{i} = [AssociatedDendCorr{i}; varargin{i}.SpectralDendriteInformation{j}];
            SpatioTemporalOverlap{j} = [SpatioTemporalOverlap{j}; nanmean(varargin{i}.SpatioTemporalOverlap{j})];
            for k = 1:length(varargin{i}.SpatialDegree{j})
                SpatialDegree{j} = [SpatialDegree{j}; varargin{i}.SpatialDegree{j}{k}];
                TemporalDegree{j} = [TemporalDegree{j}; varargin{i}.TemporalDegree{j}{k}];
            end
            CorrelationofClusters{j} = [CorrelationofClusters{j}; reshape(varargin{i}.CorrelationofClusters{j},length(varargin{i}.CorrelationofClusters{j}),1)];
        end
        
        SpatialDegreeofCueSpines(i,1:14) = varargin{i}.MeanSpatialDegreeofCueSpines;
        TemporalDegreeofCueSpines(i,1:14) = varargin{i}.MeanTemporalDegreeofCueSpines;
        SpatioTemporalDegreeofCueSpines(i,1:14) = varargin{i}.MeanSpatioTemporalDegreeofCueSpines;
        SpatialDegreeofMovementSpines(i,1:14) = varargin{i}.MeanSpatialDegreeofMovementSpines;
        TemporalDegreeofMovementSpines(i,1:14) = varargin{i}.MeanTemporalDegreeofMovementSpines;
        SpatioTemporalDegreeofMovementSpines(i,1:14) = varargin{i}.MeanSpatioTemporalDegreeofMovementSpines;
        SpatialDegreeofMovementDuringCueSpines(i,1:14) = varargin{i}.MeanSpatialDegreeofMovementDuringCueSpines;
        TemporalDegreeofMovementDuringCueSpines(i,1:14) = varargin{i}.MeanTemporalDegreeofMovementDuringCueSpines;
        SpatioTemporalDegreeofMovementDuringCueSpines(i,1:14) = varargin{i}.MeanSpatioTemporalDegreeofMovementDuringCueSpines;
        SpatialDegreeofPreSuccessSpines(i,1:14) = varargin{i}.MeanSpatialDegreeofPreSuccessSpines;
        TemporalDegreeofPreSuccessSpines(i,1:14) = varargin{i}.MeanTemporalDegreeofPreSuccessSpines;
        SpatioTemporalDegreeofPreSuccessSpines(i,1:14) = varargin{i}.MeanSpatioTemporalDegreeofPreSuccessSpines;
        SpatialDegreeofSuccessSpines(i,1:14) = varargin{i}.MeanSpatialDegreeofSuccessSpines;
        TemporalDegreeofSuccessSpines(i,1:14) = varargin{i}.MeanTemporalDegreeofSuccessSpines;
        SpatioTemporalDegreeofSuccessSpines(i,1:14) = varargin{i}.MeanSpatioTemporalDegreeofSuccessSpines;
        SpatialDegreeofRewardSpines(i,1:14) = varargin{i}.MeanSpatialDegreeofRewardSpines;
        TemporalDegreeofRewardSpines(i,1:14) = varargin{i}.MeanTemporalDegreeofRewardSpines;
        SpatioTemporalDegreeofRewardSpines(i,1:14) = varargin{i}.MeanSpatioTemporalDegreeofRewardSpines;
        
        ClusterFreq(i,1:14) = varargin{i}.ClusterFrequency;
        NonClusteredFreq(i,1:14) = varargin{i}.NonClusteredFrequency;
        CueClusterFrequency(i,1:14) = varargin{i}.CueClusterFrequency;
        MovementClusterFrequency(i,1:14) = varargin{i}.MovementClusterFrequency;
        MovementDuringCueClusterFrequency(i,1:14) = varargin{i}.MovementDuringCueClusterFrequency;
        PreSuccessClusterFrequency(i,1:14) = varargin{i}.PreSuccessClusterFrequency;
        SuccessClusterFrequency(i,1:14) = varargin{i}.SuccessClusterFrequency;
        RewardClusterFrequency(i,1:14) = varargin{i}.RewardClusterFrequency;
        CausalClusterFreq(i,1:14) = varargin{i}.CausalClusterFrequency;
        NonClusteredCausalFreq(i,1:14) = varargin{i}.NonClusteredCausalFrequency;
        CausalCueClusterFrequency(i,1:14) = varargin{i}.CausalCueClusterFrequency;
        CausalMovementClusterFrequency(i,1:14) = varargin{i}.CausalMovementClusterFrequency;
        CausalMovementDuringCueClusterFrequency(i,1:14) = varargin{i}.CausalMovementDuringCueClusterFrequency;
        CausalPreSuccessClusterFrequency(i,1:14) = varargin{i}.CausalPreSuccessClusterFrequency;
        CausalSuccessClusterFrequency(i,1:14) = varargin{i}.CausalSuccessClusterFrequency;
        CausalRewardClusterFrequency(i,1:14) = varargin{i}.CausalRewardClusterFrequency;

        ClusteredSpineAmp(i,1:14) = varargin{i}.ClusteredSpineAmp;
        NonClusteredSpineAmp(i,1:14) = varargin{i}.NonClusteredSpineAmp;
        ClusteredCueSpineAmp(i,1:14) = varargin{i}.ClusteredCueSpineAmp;
        ClusteredMoveSpineAmp(i,1:14) = varargin{i}.ClusteredMovSpineAmp;
        ClusteredMovDuringCueSpineAmp(i,1:14) = varargin{i}.ClusteredMovDuringCueSpineAmp;
        ClusteredPreSuccessSpineAmp(i,1:14) = varargin{i}.ClusteredPreSuccessSpineAmp;
        ClusteredSuccessSpineAmp(i,1:14) = varargin{i}.ClusteredSuccessSpineAmp;
        ClusteredRewardSpineAmp(i,1:14) = varargin{i}.ClusteredRewardSpineAmp;
        CausalClusteredSpineAmp(i,1:14) = varargin{i}.CausalClusteredSpineAmp;
        CausalNonClusteredSpineAmp(i,1:14) = varargin{i}.CausalNonClusteredSpineAmp;
        CausalClusteredCueSpineAmp(i,1:14) = varargin{i}.CausalClusteredCueSpineAmp;
        CausalClusteredMoveSpineAmp(i,1:14) = varargin{i}.CausalClusteredMovSpineAmp;
        CausalClusteredMovDuringCueSpineAmp(i,1:14) = varargin{i}.CausalClusteredMovDuringCueSpineAmp;
        CausalClusteredPreSuccessSpineAmp(i,1:14) = varargin{i}.CausalClusteredPreSuccessSpineAmp;
        CausalClusteredSuccessSpineAmp(i,1:14) = varargin{i}.CausalClusteredSuccessSpineAmp;
        CausalClusteredRewardSpineAmp(i,1:14) = varargin{i}.CausalClusteredRewardSpineAmp;
        
        ClustDendFreq(i,1:14) = varargin{i}.DendriteswithClustersFrequency;
        NonClustDendFreq(i,1:14) = varargin{i}.DendriteswithoutClustersFrequency;
        CueClustDendFreq(i,1:14) = varargin{i}.DendriteswithCueClustersFrequency;
        MovClustDendFreq(i,1:14) = varargin{i}.DendriteswithMovClustersFrequency;
        MovDuringCueClustDendFreq(i,1:14) = varargin{i}.DendriteswithMovDuringCueClustersFrequency;
        PreSucClustDendFreq(i,1:14) = varargin{i}.DendriteswithPreSucClustersFrequency;
        SucClustDendFreq(i,1:14) = varargin{i}.DendriteswithSucClustersFrequency;
        RewClustDendFreq(i,1:14) = varargin{i}.DendriteswithRewClustersFrequency;
        NonMovClustDendFreq(i,1:14) = varargin{i}.DendriteswithoutMovClustersFrequency;        
        
        NumCueRelSpines(i,1:14) = varargin{i}.NumberofCueSpines;
        NumMovRelSpines(i,1:14) = varargin{i}.NumberofMovementRelatedSpines;
        FractionofMovementRelatedSpinesPerDendrite(1:14) = cellfun(@(x,y) [x,y], FractionofMovementRelatedSpinesPerDendrite, varargin{i}.FractionofMovementRelatedSpinesPerDendrite, 'Uni', false);
        NumCueORMovRelSpines(i,1:14) = varargin{i}.NumberofCueORMovementRelatedSpines;
        NumPreSucRelSpines(i,1:14) = varargin{i}.NumberofPreSuccessSpines;
        NumSucRelSpines(i,1:14) = varargin{i}.NumberofSuccessSpines;
        NumMovDuringCueRelSpines(i,1:14) = varargin{i}.NumberofMovementDuringCueSpines;
        NumRewRelSpines(i,1:14) = varargin{i}.NumberofRewardSpines;
        NumCausalMovSpines(i,1:14) = varargin{i}.NumberofCausalMvmntSpines;
        NumCausalSucSpines(i,1:14) = varargin{i}.NumberofCausalSuccessSpines;
        NumCausalCueSpines(i,1:14) = varargin{i}.NumberofCausalCueSpines;
        
        NumClustSpines(i,1:14) = varargin{i}.NumberofClusteredSpines;
        NumClustCueSpines(i,1:14) = varargin{i}.NumberofClusteredCueSpines;
        NumClustMovSpines(i,1:14) = varargin{i}.NumberofClusteredMovementSpines;
        NumClustMixSpines(i,1:14) = varargin{i}.NumberofClusteredMixedFunctionSpines;
        NumClustPreSucSpines(i,1:14) = varargin{i}.NumberofClusteredPreSuccessSpines;
        NumClustSucSpines(i,1:14) = varargin{i}.NumberofClusteredSuccessSpines;
        NumClustMovDuringCueSpines(i,1:14) = varargin{i}.NumberofClusteredMovementDuringCueSpines;
        NumClustRewSpines(i,1:14) = varargin{i}.NumberofClusteredRewardSpines;
        
        NumFarClustSpines(i,1:14) = varargin{i}.NumberofFarClusteredSpines;
        NumFarClustCueSpines(i,1:14) = varargin{i}.NumberofFarClusteredCueSpines;
        NumFarClustMovSpines(i,1:14) = varargin{i}.NumberofFarClusteredMovementSpines;
        NumFarClustMixSpines(i,1:14) = varargin{i}.NumberofClusteredMixedFunctionSpines;
        NumFarClustPreSucSpines(i,1:14) = varargin{i}.NumberofFarClusteredPreSuccessSpines;
        NumFarClustSucSpines(i,1:14) = varargin{i}.NumberofFarClusteredSuccessSpines;
        NumFarClustMovDuringCueSpines(i,1:14) = varargin{i}.NumberofFarClusteredMovementDuringCueSpines;
        NumFarClustRewSpines(i,1:14) = varargin{i}.NumberofFarClusteredRewardSpines;
        
        NumCausClustSpines(i,1:14) = varargin{i}.NumberofCausalClusteredSpines;
        NumCausClustCueSpines(i,1:14) = varargin{i}.NumberofCausalClusteredCueSpines;
        NumCausClustMovSpines(i,1:14) = varargin{i}.NumberofCausalClusteredMovementSpines;
        NumCausClustMovDuringCueSpines(i,1:14) = varargin{i}.NumberofCausalClusteredMovementDuringCueSpines;
        NumCausClustPreSucSpines(i,1:14) = varargin{i}.NumberofCausalClusteredSuccessSpines;
        NumCausClustSucSpines(i,1:14) = varargin{i}.NumberofCausalClusteredSuccessSpines;
        NumCausClustRewSpines(i,1:14) = varargin{i}.NumberofCausalClusteredRewardSpines;

        NumberofClusters(i,1:14) = varargin{i}.NumberofClusters;
        NumberofCausalClusters(i,1:14) = varargin{i}.NumberofCausalClusters;
        NumberofSpinesinEachCluster(i,1:14) = varargin{i}.MeanNumberofSpinesinEachCluster;
        NumberofSpinesinEachCausalCluster(i,1:14) = varargin{i}.MeanNumberofSpinesinEachCausalCluster;
        NumberofMovClusters(i,1:14) = varargin{i}.NumberofMovClusters;
        NumberofSpinesinEachMovCluster(i,1:14) = varargin{i}.MeanNumberofSpinesinEachMovCluster;
        
        CueClusterLength(i,1:14) = varargin{i}.MeanCueClustLength;
        CueClusterMax(i,1:14) = varargin{i}.MaxCueClustLength;
        MovClusterLength(i,1:14) = varargin{i}.MeanMovClustLength;
            AllMovClustLengths(1:14) = cellfun(@(x,y) [x,y], AllMovClustLengths, varargin{i}.AllMovClustLengths, 'Uni', false);
        MovClusterMax(i,1:14) = varargin{i}.MaxMovClustLength;
        MixClusterLength(i,1:14) = varargin{i}.MeanMixClustLength;
        MixClusterMax(i,1:14) = varargin{i}.MaxMixClustLength;
        PreSucClusterLength(i,1:14) = varargin{i}.MeanPreSucClustLength;
        PreSucClusterMax(i,1:14) = varargin{i}.MaxPreSucClustLength;
        SucClusterLength(i,1:14) = varargin{i}.MeanSucClustLength;
        SucClusterMax(i,1:14) = varargin{i}.MaxSucClustLength;
        MovDuringCueClusterLength(i,1:14) = varargin{i}.MeanMovDuringCueClustLength;
        MovDuringCueClusterMax(i,1:14) = varargin{i}.MaxMovDuringCueClustLength;
        RewClusterLength(i,1:14) = varargin{i}.MeanRewClustLength;
        RewClusterMax(i,1:14) = varargin{i}.MaxRewClustLength;
        
        AllClusterLength(i,1:14) = varargin{i}.MeanAllClustLength;
        AllClusterMax(i,1:14) = varargin{i}.MaxAllClustLength;
        CausalCueClusterLength(i,1:14) = varargin{i}.MeanCausalCueClustLength;
        CausalCueClusterMax(i,1:14) = varargin{i}.MaxCausalCueClustLength;
        CausalMovClusterLength(i,1:14) = varargin{i}.MeanCausalMovClustLength;
        CausalMovClusterMax(i,1:14) = varargin{i}.MaxCausalMovClustLength;
        CausalSucClusterLength(i,1:14) = varargin{i}.MeanCausalSucClustLength;
        CausalSucClusterMax(i,1:14) = varargin{i}.MaxCausalSucClustLength;
        CausalRewClusterLength(i,1:14) = varargin{i}.MeanCausalRewClustLength;
        CausalRewClusterMax(i,1:14) = varargin{i}.MaxCausalRewClustLength;
        AllCausalClusterLength(i,1:14) = varargin{i}.MeanCausalAllClustLength;
        AllCausalClusterMax(i,1:14) = varargin{i}.MaxCausalAllClustLength;
        
        AllFarClusterLength(i,1:14) = varargin{i}.MeanAllFarClustLength;
        FarCueClusterLength(i,1:14) = varargin{i}.MeanFarCueClustLength;
        FarMovClusterLength(i,1:14) = varargin{i}.MeanFarMovClustLength;
        FarMixClusterLength(i,1:14) = varargin{i}.MeanFarMixClustLength;
        FarPreSucClusterLength(i,1:14) = varargin{i}.MeanFarPreSucClustLength;
        FarSucClusterLength(i,1:14) = varargin{i}.MeanFarSucClustLength;
        FarMovDuringCueClusterLength(i,1:14) = varargin{i}.MeanFarMovDuringCueClustLength;
        FarRewClusterLength(i,1:14) = varargin{i}.MeanFarRewClustLength;
        
            AllDistancesBetweenAllSpines(1:14) = cellfun(@(x,y) [x,y], AllDistancesBetweenAllSpines, varargin{i}.DistanceBetweenAllSpines, 'Uni', false);
            CorrelationBetweenAllSpines(1:14) = cellfun(@(x,y) [x,y], CorrelationBetweenAllSpines, varargin{i}.CorrelationBetweenAllSpines, 'Uni', false);
            CorrelationBetweenAllSpinesMovePeriods(1:14) = cellfun(@(x,y) [x,y], CorrelationBetweenAllSpinesMovePeriods, varargin{i}.CorrelationBetweenAllSpinesMovementPeriods, 'Uni', false);
            CorrelationBetweenAllSpinesStillPeriods(1:14) = cellfun(@(x,y) [x,y], CorrelationBetweenAllSpinesStillPeriods, varargin{i}.CorrelationBetweenAllSpinesStillPeriods, 'Uni', false);
            if any(cell2mat(cellfun(@(x,y) length(x)~=length(y), varargin{i}.CorrelationBetweenFarSpines, varargin{i}.DistanceBetweenFarSpines, 'Uni', false)))
                problemdays = find(cell2mat(cellfun(@(x,y) length(x)~=length(y), varargin{i}.CorrelationBetweenFarSpines, varargin{i}.DistanceBetweenFarSpines, 'Uni', false)));
                file = inputname(i); file = file(1:5);
                for c = 1:length(problemdays)
                    fprintf('Correlation and distance vectors \n not equal for session %d \n from input %2s \n', problemdays(c), file)
                end
            end
                    AllDistancesBetweenFarSpines(1:14) = cellfun(@(x,y) [x;y], AllDistancesBetweenFarSpines, varargin{i}.DistanceBetweenFarSpines, 'Uni', false);
                    CorrelationBetweenFarSpines(1:14) = cellfun(@(x,y) [x;y], CorrelationBetweenFarSpines, varargin{i}.CorrelationBetweenFarSpines, 'Uni', false);
                    NearestMovementRelatedSpine(1:14) = cellfun(@(x,y) [x,y], NearestMovementRelatedSpine, varargin{i}.NearestMovementRelatedSpine, 'Uni', false);
                    NextClosestMovementRelatedSpine(1:14) = cellfun(@(x,y) [x,y], NextClosestMovementRelatedSpine, varargin{i}.NextClosestMovementRelatedSpine, 'Uni', false);
                    ThirdClosestMovementRelatedSpine(1:14) = cellfun(@(x,y) [x,y], ThirdClosestMovementRelatedSpine, varargin{i}.ThirdClosestMovementRelatedSpine, 'Uni', false);
                    FourthClosestMovementRelatedSpine(1:14) = cellfun(@(x,y) [x,y], FourthClosestMovementRelatedSpine, varargin{i}.FourthClosestMovementRelatedSpine, 'Uni', false);
                    CorrelationwithNearestMovementRelatedSpine(1:14) = cellfun(@(x,y) [x,y], CorrelationwithNearestMovementRelatedSpine, varargin{i}.CorrelationwithNearestMovementRelatedSpine, 'Uni', false);
                    NearestHighlyCorrelatedMovementRelatedSpine(1:14) =  cellfun(@(x,y) [x,y], NearestHighlyCorrelatedMovementRelatedSpine, varargin{i}.NearestHighlyCorrelatedMovementRelatedSpine, 'Uni', false);
                    NextClosestHighlyCorrelatedMovementRelatedSpine(1:14) =  cellfun(@(x,y) [x,y], NextClosestHighlyCorrelatedMovementRelatedSpine, varargin{i}.NextClosestHighlyCorrelatedMovementRelatedSpine, 'Uni', false);
                    ThirdClosestHighlyCorrelatedMovementRelatedSpine(1:14) =  cellfun(@(x,y) [x,y], ThirdClosestHighlyCorrelatedMovementRelatedSpine, varargin{i}.ThirdClosestHighlyCorrelatedMovementRelatedSpine, 'Uni', false);
        DistanceBetweenCueSpines(i,1:14) = varargin{i}.MeanDistanceBetweenCueSpines;
        DistanceBetweenMovementSpines(i,1:14) = varargin{i}.MeanDistanceBetweenMovementSpines;
            AllDistancesBetweenMovementSpines(1:14) = cellfun(@(x,y) [x,y], AllDistancesBetweenMovementSpines, varargin{i}.DistanceBetweenMovementSpines, 'Uni', false);
            CorrelationBetweenMovementSpines(1:14) = cellfun(@(x,y) [x,y], CorrelationBetweenMovementSpines, varargin{i}.CorrelationBetweenMovementSpines, 'Uni', false);
            CorrelationBetweenMovementSpinesMovePeriods(1:14) = cellfun(@(x,y) [x,y], CorrelationBetweenMovementSpinesMovePeriods, varargin{i}.CorrelationBetweenMovementSpinesMovementPeriods, 'Uni', false);
            CorrelationBetweenMovementSpinesStillPeriods(1:14) = cellfun(@(x,y) [x,y], CorrelationBetweenMovementSpinesStillPeriods, varargin{i}.CorrelationBetweenMovementSpinesStillPeriods, 'Uni', false);
                AllDistancesBetweenFarMovementSpines(1:14) = cellfun(@(x,y) [x,y], AllDistancesBetweenFarMovementSpines, varargin{i}.DistanceBetweenFarMovementSpines, 'Uni', false);
                CorrelationBetweenFarMovementSpines(1:14) = cellfun(@(x,y) [x,y], CorrelationBetweenFarMovementSpines, varargin{i}.CorrelationBetweenFarMovementSpines, 'Uni', false);
                MeanCorrelationBetweenFarMovementSpines(i,1:14) = cell2mat(cellfun(@nanmean, varargin{i}.CorrelationBetweenFarMovementSpines, 'Uni', false));
        DistanceBetweenPreSuccessSpines(i,1:14) = varargin{i}.MeanDistanceBetweenPreSuccessSpines;
        DistanceBetweenSuccessSpines(i,1:14) = varargin{i}.MeanDistanceBetweenSuccessSpines;
        DistanceBetweenMovementDuringCueSpines(i,1:14) = varargin{i}.MeanDistanceBetweenMovementDuringCueSpines;
        DistanceBetweenRewardSpines(i,1:14) = varargin{i}.MeanDistanceBetweenRewardSpines;
        
        ClusteredSpines_CorrwithDend(i,1:14) = varargin{i}.ClusteredSpines_CorrwithDend;
        FilteredClusteredSpines_CorrwithDend(i,1:14) = varargin{i}.FilteredClusteredSpines_CorrwithDend;
        NonClusteredSpines_CorrwithDend(i,1:14) = varargin{i}.NonClusteredSpines_CorrwithDend;
        CausalClusteredSpines_CorrwithDend(i,1:14) = varargin{i}.CausalClusteredSpines_CorrwithDend;
        FilteredCausalClusteredSpines_CorrwithDend(i,1:14) = varargin{i}.FilteredCausalClusteredSpines_CorrwithDend;
        NonCausalClusteredSpines_CorrwithDend(i,1:14) = varargin{i}.NonCausalClusteredSpines_CorrwithDend;
        CueRelClusteredSpines_CorrwithDend(i,1:14) = varargin{i}.CueRelClusteredSpines_CorrwithDend;
        MovRelClusteredSpines_CorrwithDend(i,1:14) = varargin{i}.MovRelClusteredSpines_CorrwithDend;
        SucRelClusteredSpines_CorrwithDend(i,1:14) = varargin{i}.SucRelClusteredSpines_CorrwithDend;
        RewRelClusteredSpines_CorrwithDend(i,1:14) = varargin{i}.RewRelClusteredSpines_CorrwithDend;
        CueRelCausalClusteredSpines_CorrwithDend(i,1:14) = varargin{i}.CueRelCausalClusteredSpines_CorrwithDend;
        MovRelCausalClusteredSpines_CorrwithDend(i,1:14) = varargin{i}.MovRelCausalClusteredSpines_CorrwithDend;
        SucRelCausalClusteredSpines_CorrwithDend(i,1:14) = varargin{i}.SucRelCausalClusteredSpines_CorrwithDend;
        RewRelCausalClusteredSpines_CorrwithDend(i,1:14) = varargin{i}.RewRelCausalClusteredSpines_CorrwithDend;
        
%         FractionofClusterThatsMovementRelated(i,1:14) = varargin{i}.AverageFractionofaClusterThatsMovementRelated;
%         FractionofCausalClusterThatsMovementRelated(i,1:14) = varargin{i}.AverageFractionofaCausalClusterThatsMovementRelated;
    
        PercentofCueRelatedDendrites(i,1:14) = varargin{i}.PercentofCueRelatedDendrites;
        PercentofMovementRelatedDendrites(i,1:14) = varargin{i}.PercentofMovementRelatedDendrites;
        PercentofPreSuccessRelatedDendrites(i,1:14) = varargin{i}.PercentofPreSuccessRelatedDendrites;
        PercentofSuccessRelatedDendrites(i,1:14) = varargin{i}.PercentofSuccessRelatedDendrites;
        PercentofMovementDuringCueRelatedDendrites(i,1:14) = varargin{i}.PercentofMovementDuringCueRelatedDendrites;
        PercentofRewardRelatedDendrites(i,1:14) = varargin{i}.PercentofRewardRelatedDendrites;
        
%         Spatial_Degree(i,1:14) = varargin{i}.MeanDendriticSpatialDegree;
%         Temporal_Degree(i,1:14) = varargin{i}.MeanDendriticTemporalDegree;
%         ST_Degree(i,1:14) = varargin{i}.MeanDendriticSpatioTemporalDegree;
    end
    
    
    
    %%% Save Organized data for stats %%%
        
        a.AllClustersCorrwithCue = AllClustersCorrwithCue;
        a.NonClusteredCorrwithCue = NonClusteredCorrwithCue;
        a.AllClustersCorrwithMovement = AllClustersCorrwithMovement;
        a.NonClusteredCorrwithMovement = NonClusteredCorrwithMovement;
        a.AllClustersCorrwithMDC = AllClustersCorrwithMDC;
        a.NonClusteredCorrwithMDC = NonClusteredCorrwithMDC;
        a.AllClustersCorrwithSuccess = AllClustersCorrwithSuccess;
        a.NonClusteredCorrwithSuccess = NonClusteredCorrwithSuccess;
        a.AllClustCorrwithReward = AllClustCorrwithReward;
        a.NonClusteredCorrwithReward = NonClusteredCorrwithReward;
        
        a.CorrelationofClusters = CorrelationofClusters;
        
        a.AllCausalClustersCorrwithMovement = AllCausalClustersCorrwithMovement;
        a.CausalNonClusteredCorrwithMovement = CausalNonClusteredCorrwithMovement;
        a.AllCausalClustersCorrwithCue = AllCausalClustersCorrwithCue;
        a.CausalNonClusteredCorrwithCue = CausalNonClusteredCorrwithCue;
        a.AllCausalClustersCorrwithMDC = AllCausalClustersCorrwithMDC;
        a.CausalNonClusteredCorrwithMDC = CausalNonClusteredCorrwithMDC;
        a.AllCausalClustersCorrwithSuccess = AllCausalClustersCorrwithSuccess;
        a.CausalNonClusteredCorrwithSuccess = CausalNonClusteredCorrwithSuccess;
        a.AllCausalClustCorrwithReward = AllCausalClustCorrwithReward;
        a.CausalNonClusteredCorrwithReward = CausalNonClusteredCorrwithReward;
        
        a.CueRelatedClustersCorrwithCue = CueRelatedClustersCorrwithCue;
        a.CueRelatedNonClusteredCorrwithCue = CueRelatedNonClusteredCorrwithCue;
        a.MovementRelatedClustersCorrwithMovement = MovementRelatedClustersCorrwithMovement;
        a.MovementRelatedNonClusteredCorrwithMovement = MovementRelatedNonClusteredCorrwithMovement;
        a.MDCRelatedClustersCorrwithMDC = MDCRelatedClustersCorrwithMDC;
        a.MDCRelatedNonClusteredCorrwithMDC = MDCRelatedNonClusteredCorrwithMDC; 
        a.SuccessRelatedClustersCorrwithSuccess = SuccessRelatedClustersCorrwithSuccess;
        a.SuccessRelatedNonClsuteredCorrwithSuccess = SuccessRelatedNonClusteredCorrwithSuccess;
        a.RewardRelatedClutersCorrwithReward = RewardRelatedClustersCorrwithReward;
        a.RewardRelatedNonClusteredCorrwithReward = RewardRelatedNonClusteredCorrwithReward;
        
        a.CausalCueRelatedClustersCorrwithCue = CausalCueRelatedClustersCorrwithCue;
        a.CausalCueRelatedNonClusteredCorrwithCue = CausalCueRelatedNonClusteredCorrwithCue;
        a.CausalMovementRelatedClustersCorrwithMovement = CausalMovementRelatedClustersCorrwithMovement;
        a.CausalMovementRelatedNonClusteredCorrwithMovement = CausalMovementRelatedNonClusteredCorrwithMovement;
        a.CausalMDCRelatedClustersCorrwithMDC = CausalMDCRelatedClustersCorrwithMDC;
        a.CausalMDCRelatedNonClusteredCorrwithMDC = CausalMDCRelatedNonClusteredCorrwithMDC; 

        a.CausalSuccessRelatedClustersCorrwithSuccess = CausalSuccessRelatedClustersCorrwithSuccess;
        a.CausalSuccessRelatedNonClsuteredCorrwithSuccess = CausalSuccessRelatedNonClusteredCorrwithSuccess;
        a.CausalRewardRelatedClutersCorrwithReward = CausalRewardRelatedClustersCorrwithReward;
        a.CausalRewardRelatedNonClusteredCorrwithReward = CausalRewardRelatedNonClusteredCorrwithReward;
        
        a.FractionofCueSpinesThatAreClustered = FractionofCueSpinesThatAreClustered;
        a.FractionofMovementSpinesThatAreClustered = FractionofMovementSpinesThatAreClustered;
        a.FractionofPreSuccessSpinesThatAreClustered = FractionofPreSuccessSpinesThatAreClustered;
        a.FractionofSuccessSpinesThatAreClustered = FractionofSuccessSpinesThatAreClustered;
        a.FractionofMDCSpinesThatAreClustered = FractionofMovementDuringCueSpinesThatAreClustered;
        a.FractionofRewardSpinesThatAreClustred = FractionofRewardSpinesThatAreClustered;
        
    
        a.SpatioTemporalDegree = SpatioTemporalDegree;
        a.SpatialMovementCorrelation = SpatialMovementCorrelation;
        a.TemporalMovementCorrelation = TemporalMovementCorrelation;
        a.SpatioTemporalMovementCorrelation = SpatioTemporalMovementCorrelation;
        a.DendClusteringDegree = DendClusteringDegree;
        a.DendClusteringDegree = SpatioTemporalOverlap;
        a.SpatialDegree = SpatialDegree;
        a.TemporalDegree = TemporalDegree;
        a.SpatialDegreeofCueSpines = SpatialDegreeofCueSpines;
        a.TemporalDegreeofCueSpines = TemporalDegreeofCueSpines;
        a.SpatioTemporalDegreeofCueSpines = SpatioTemporalDegreeofCueSpines;
        a.SpatialDegreeofMovementSpines = SpatialDegreeofMovementSpines;
        a.TemporalDegreeofMovementSpines = TemporalDegreeofMovementSpines;
        a.SpatioTemporalDegreeofMovementSpines = SpatioTemporalDegreeofMovementSpines;
        a.SpatialDegreeofMovementDuringCueSpines = SpatialDegreeofMovementDuringCueSpines;
        a.TemporalDegreeofMovementDuringCueSpines = TemporalDegreeofMovementDuringCueSpines;
        a.SpatioTemporalDegreeofMovementDuringCueSpines = SpatioTemporalDegreeofMovementDuringCueSpines;
        a.SpatialDegreeofPreSuccessSpines = SpatialDegreeofPreSuccessSpines;
        a.TemporalDegreeofPreSuccessSpines = TemporalDegreeofPreSuccessSpines;
        a.SpatioTemporalDegreeofPreSuccessSpines = SpatioTemporalDegreeofPreSuccessSpines;
        a.SpatialDegreeofSuccessSpines = SpatialDegreeofSuccessSpines;
        a.TemporalDegreeofSuccessSpines = TemporalDegreeofSuccessSpines;
        a.SpatioTemporalDegreeofSuccessSpines = SpatioTemporalDegreeofSuccessSpines;
        a.SpatialDegreeofRewardSpines = SpatialDegreeofRewardSpines;
        a.TemporalDegreeofRewardSpines = TemporalDegreeofRewardSpines;
        a.SpatioTemporalDegreeofRewardSpines = SpatioTemporalDegreeofRewardSpines;
        
        a.ClusterFreq = ClusterFreq;
        a.NonClusterFreq = NonClusteredFreq;
        a.CueClusterFrequency = CueClusterFrequency;
        a.MovementClusterFrequency = MovementClusterFrequency;
        a.MovementDuringCueClusterFrequency = MovementDuringCueClusterFrequency;
        a.PreSuccessClusterFrequency = PreSuccessClusterFrequency;
        a.SuccessClusterFrequency = SuccessClusterFrequency;
        a.RewardClusterFrequency = RewardClusterFrequency;
        a.CausalClusterFreq = CausalClusterFreq;
        a.NonClusteredCausalFreq = NonClusteredCausalFreq;
        a.CausalCueClusterFrequency = CausalCueClusterFrequency;
        a.CausalMovementClusterFrequency = CausalMovementClusterFrequency;
        a.CausalMovementDuringCueClusterFrequency = CausalMovementDuringCueClusterFrequency;
        a.CausalPreSuccessClusterFrequency = CausalPreSuccessClusterFrequency;
        a.CausalSuccessClusterFrequency = CausalSuccessClusterFrequency;
        a.CausalRewardClusterFrequency = CausalRewardClusterFrequency;

        a.ClusteredSpineAmp = ClusteredSpineAmp;
        a.NonClusteredSpineAmp = NonClusteredSpineAmp;
        a.ClusteredCueSpineAmp = ClusteredCueSpineAmp;
        a.ClusteredMoveSpineAmp = ClusteredMoveSpineAmp;
        a.ClusteredMovDuringCueSpineAmp = ClusteredMovDuringCueSpineAmp;
        a.ClusteredPreSuccessSpineAmp = ClusteredPreSuccessSpineAmp;
        a.ClusteredSuccessSpineAmp = ClusteredSuccessSpineAmp;
        a.ClusteredRewardSpineAmp = ClusteredRewardSpineAmp;
        a.CausalClusteredSpineAmp = CausalClusteredSpineAmp;
        a.CausalNonClusteredSpineAmp = CausalNonClusteredSpineAmp;
        a.CausalClusteredCueSpineAmp = CausalClusteredCueSpineAmp;
        a.CausalClusteredMoveSpineAmp = CausalClusteredMoveSpineAmp;
        a.CausalClusteredMovDuringCueSpineAmp = CausalClusteredMovDuringCueSpineAmp;
        a.CausalClusteredPreSuccessSpineAmp = CausalClusteredPreSuccessSpineAmp;
        a.CausalClusteredSuccessSpineAmp = CausalClusteredSuccessSpineAmp;
        a.CausalClusteredRewardSpineAmp = CausalClusteredRewardSpineAmp;
        
        a.ClustDendFreq = ClustDendFreq;
        a.NonClustDendFreq =NonClustDendFreq;
        a.MovClustDendFreq = MovClustDendFreq;
        a.NonMovClustDendFreq = NonMovClustDendFreq;
        
        a.NumCueRelSpines = NumCueRelSpines;
        a.NumMvmtSpines = NumMovRelSpines;
        a.NumCueORMovRelSpines = NumCueORMovRelSpines;
        a.NumPreSucRelSpines = NumPreSucRelSpines;
        a.NumSucRelSpines = NumSucRelSpines;
        a.NumMovDuringCueRelSpines = NumMovDuringCueRelSpines;
        a.NumRewRelSpines = NumRewRelSpines;
        a.NumCausalMovSpines = NumCausalMovSpines;
        a.NumCausalSucSpines = NumCausalSucSpines;
        a.NumCausalCueSpines = NumCausalCueSpines;
        a.NumClustSpines = NumClustSpines;
        a.NumClustCueSpines = NumClustCueSpines;
        a.NumClustMovSpines = NumClustMovSpines;
        a.NumClustMixSpines = NumClustMixSpines;
        a.NumClustPreSucSpines = NumClustPreSucSpines;
        a.NumClustSucSpines = NumClustSucSpines;
        a.NumClustMovDuringCueSpines = NumClustMovDuringCueSpines;
        a.NumClustRewSpines = NumClustRewSpines;
        a.NumFarClustSpines = NumFarClustSpines;
        a.NumFarClustCueSpines = NumFarClustCueSpines;
        a.NumFarClustMovSpines = NumFarClustMovSpines;
        a.NumFarClustMixSpines = NumFarClustMixSpines;
        a.NumFarClustPreSucSpines = NumFarClustPreSucSpines;
        a.NumFarClustSucSpines = NumFarClustSucSpines;
        a.NumFarClustMovDuringCueSpines = NumFarClustMovDuringCueSpines;
        a.NumFarClustRewSpines = NumFarClustRewSpines;
        a.NumCausClustSpines = NumCausClustSpines;
        a.NumCausClustCueSpines = NumCausClustCueSpines;
        a.NumCausClustMovSpines = NumCausClustMovSpines;
        a.NumCausClustMovDuringCueSpines = NumCausClustMovDuringCueSpines;
        a.NumCausClustPreSucSpines = NumCausClustPreSucSpines;
        a.NumCausClustSucSpines = NumCausClustSucSpines;
        a.NumCausClustRewSpines = NumCausClustRewSpines;

        a.NumberofClusters = NumberofClusters;
        a.NumberofCausalClusters = NumberofCausalClusters;
        a.NumberofSpinesinEachCluster = NumberofSpinesinEachCluster;
        a.NumberofSpinesinEachCausalCluster = NumberofSpinesinEachCausalCluster;
        a.NumberofMovClusters = NumberofMovClusters;
        a.NumberofSpinesinEachMovCluster= NumberofSpinesinEachMovCluster;
        
        a.CueClusterLength = CueClusterLength;
        a.CueClusterMax = CueClusterMax;
        a.MovClusterLength = MovClusterLength;
        a.MovClusterMax = MovClusterMax;
        a.MixClusterLength = MixClusterLength;
        a.MixClusterMax = MixClusterMax;
        a.SucClusterLength = SucClusterLength;
        a.SucClusterMax = SucClusterMax;
        a.RewClusterLength = RewClusterLength;
        a.RewClusterMax = RewClusterMax;
        a.AllClusterLength = AllClusterLength;
        a.AllClusterLength = AllClusterMax;
        a.CausalCueClusterLength = CausalCueClusterLength;
        a.CausalCueClusterMax = CausalCueClusterMax;
        a.CausalMovClusterLength = CausalMovClusterLength;
        a.CausalMovClusterMax = CausalMovClusterMax;
        a.CausalSucClusterLength = CausalSucClusterLength;
        a.CausalSucClusterMax = CausalSucClusterMax;
        a.CausalRewClusterLength = CausalRewClusterLength;
        a.CausalRewClusterMax = CausalRewClusterMax;
        a.AllCausalClusterLength = AllCausalClusterLength;
        a.AllCausalClusterMax = AllCausalClusterMax;
        a.FarCueClusterLength = FarCueClusterLength;
        a.FarMovClusterLength = FarMovClusterLength;
        a.FarMixClusterLength = FarMixClusterLength;
        a.FarPreSucClusterLength = FarPreSucClusterLength;
        a.FarSucClusterLength = FarSucClusterLength;
        a.FarMovDuringCueClusterLength = FarMovDuringCueClusterLength;
        a.FarRewClusterLength = FarRewClusterLength;
        
        a.DistanceBetweenCueSpines = DistanceBetweenCueSpines;
        a.DistanceBetweenMovementSpines = DistanceBetweenMovementSpines;
            a.AllDistancesBetweenMovementSpines = AllDistancesBetweenMovementSpines;
            a.CorrelationBetweenMovementSpines = CorrelationBetweenMovementSpines;
            a.MeanCorrelationBetweenMovementSpines = MeanCorrelationBetweenMovementSpines;
            a.MeanCorrelationBetweenFarMovementSpines = MeanCorrelationBetweenFarMovementSpines;
        a.DistanceBetweenSuccessSpines = DistanceBetweenSuccessSpines;
        a.DistanceBetweenRewardSpines = DistanceBetweenRewardSpines;
        
        a.ClusteredSpines_CorrwithDend = ClusteredSpines_CorrwithDend;
        a.FilteredClusteredSpines_CorrwithDend = FilteredClusteredSpines_CorrwithDend;
        a.NonClusteredSpines_CorrwithDend = NonClusteredSpines_CorrwithDend;
        a.CausalClusteredSpines_CorrwithDend = CausalClusteredSpines_CorrwithDend;
        a.FilteredCausalClusteredSpines_CorrwithDend = FilteredCausalClusteredSpines_CorrwithDend;
        a.NonCausalClusteredSpines_CorrwithDend = NonCausalClusteredSpines_CorrwithDend;
        a.CueRelClusteredSpines_CorrwithDend = CueRelClusteredSpines_CorrwithDend;
        a.MovRelClusteredSpines_CorrwithDend = MovRelClusteredSpines_CorrwithDend;
        a.SucRelClusteredSpines_CorrwithDend = SucRelClusteredSpines_CorrwithDend;
        a.RewRelClusteredSpines_CorrwithDend = RewRelClusteredSpines_CorrwithDend;
        a.CueRelCausalClusteredSpines_CorrwithDend = CueRelCausalClusteredSpines_CorrwithDend;
        a.MovRelCausalClusteredSpines_CorrwithDend = MovRelCausalClusteredSpines_CorrwithDend;
        a.SucRelCausalClusteredSpines_CorrwithDend = SucRelCausalClusteredSpines_CorrwithDend;
        a.RewRelCausalClusteredSpines_CorrwithDend = RewRelCausalClusteredSpines_CorrwithDend; 

%         a.FractionofClusterThatsMovementRelated = FractionofClusterThatsMovementRelated;
%         a.FractionofCausalClusterThatsMovementRelated = FractionofCausalClusterThatsMovementRelated;
        
        a.PercentofCueRelatedDendrites = PercentofCueRelatedDendrites;
        a.PercentofMovementRelatedDendrites = PercentofMovementRelatedDendrites;
        a.PercentofPreSuccessRelatedDendrites = PercentofPreSuccessRelatedDendrites;
        a.PercentofSuccessRelatedDendrites = PercentofSuccessRelatedDendrites;
        a.PercentofMovementDuringCueRelatedDendrites = PercentofMovementDuringCueRelatedDendrites;
        a.PercentofRewardRelatedDendrites = PercentofRewardRelatedDendrites;
        
    eval('ClusteringAllData = a;')
    save('ClusteringAllData', 'ClusteringAllData');
        
    %%%
    %%% Average all input data
    %%%
        MeanAllClustersCorrwithCue = nanmean(AllClustersCorrwithCue,1);
        MeanNonClusteredCorrwithCue = nanmean(NonClusteredCorrwithCue,1);
        MeanAllClustersCorrwithMDC = nanmean(AllClustersCorrwithMDC, 1);
        MeanNonClusteredCorrwithMDC = nanmean(NonClusteredCorrwithMDC,1);
        MeanAllClustersCorrwithMovement = nanmean(AllClustersCorrwithMovement,1);
        MeanNonClusteredCorrwithMovement = nanmean(NonClusteredCorrwithMovement,1);
        MeanAllClustersCorrwithSuccess = nanmean(AllClustersCorrwithSuccess,1);
        MeanNonClusteredCorrwithSuccess = nanmean(NonClusteredCorrwithSuccess,1);
        MeanAllClustCorrwithReward = nanmean(AllClustCorrwithReward,1);
        MeanNonClusteredCorrwithReward = nanmean(NonClusteredCorrwithReward,1);
        
        MeanCorrelationofClusters = cellfun(@nanmean, CorrelationofClusters);

        MeanAllCausalClustersCorrwithMovement = nanmean(AllCausalClustersCorrwithMovement,1);
        MeanCausalNonClusteredCorrwithMovement = nanmean(CausalNonClusteredCorrwithMovement,1);
        MeanAllCausalClustersCorrwithMDC = nanmean(AllCausalClustersCorrwithMDC,1);
        MeanCausalNonClusteredCorrwithMDC = nanmean(CausalNonClusteredCorrwithMDC,1);
        MeanAllCausalClustersCorrwithCue = nanmean(AllCausalClustersCorrwithCue,1);
        MeanCausalNonClusteredCorrwithCue = nanmean(CausalNonClusteredCorrwithCue,1);
        MeanAllCausalClustersCorrwithSuccess = nanmean(AllCausalClustersCorrwithSuccess,1);
        MeanCausalNonClusteredCorrwithSuccess = nanmean(CausalNonClusteredCorrwithSuccess,1);
        MeanAllCausalClustCorrwithReward = nanmean(AllCausalClustCorrwithReward,1);
        MeanCausalNonClusteredCorrwithReward = nanmean(CausalNonClusteredCorrwithReward,1);
        
        MeanCueRelatedClustersCorrwithCue = nanmean(CueRelatedClustersCorrwithCue,1);
        MeanCueRelatedNonClusteredCorrwithCue = nanmean(CueRelatedNonClusteredCorrwithCue, 1);
        MeanMDCRelatedClustersCorrwithMDC = nanmean(MDCRelatedClustersCorrwithMDC,1);
        MeanMDCRelatedNonClusteredCorrwithMDC = nanmean(MDCRelatedNonClusteredCorrwithMDC, 1);
        MeanMovementRelatedClustersCorrwithMovement = nanmean(MovementRelatedClustersCorrwithMovement,1);
        MeanMovementRelatedNonClusteredCorrwithMovement = nanmean(MovementRelatedNonClusteredCorrwithMovement,1);
        MeanSuccessRelatedClustersCorrwithSuccess = nanmean(SuccessRelatedClustersCorrwithSuccess,1);
        MeanSuccessRelatedNonClusteredCorrwithSuccess = nanmean(SuccessRelatedNonClusteredCorrwithSuccess);
        MeanRewardRelatedClustersCorrwithReward = nanmean(RewardRelatedClustersCorrwithReward,1);
        MeanRewardRelatedNonClusteredCorrwithReward = nanmean(RewardRelatedNonClusteredCorrwithReward, 1);
        
        MeanCausalCueRelatedClustersCorrwithCue = nanmean(CausalCueRelatedClustersCorrwithCue,1);
        MeanCausalCueRelatedNonClusteredCorrwithCue = nanmean(CausalCueRelatedNonClusteredCorrwithCue, 1);
        MeanCausalMDCRelatedClustersCorrwithMDC = nanmean(CausalMDCRelatedClustersCorrwithMDC, 1);
        MeanCausalMDCRelatedNonClusteredCorrwithMDC = nanmean(CausalMDCRelatedNonClusteredCorrwithMDC,1);
        MeanCausalMovementRelatedClustersCorrwithMovement = nanmean(CausalMovementRelatedClustersCorrwithMovement,1);
        MeanCausalMovementRelatedNonClusteredCorrwithMovement = nanmean(CausalMovementRelatedNonClusteredCorrwithMovement,1);
        MeanCausalSuccessRelatedClustersCorrwithSuccess = nanmean(CausalSuccessRelatedClustersCorrwithSuccess,1);
        MeanCausalSuccessRelatedNonClusteredCorrwithSuccess = nanmean(CausalSuccessRelatedNonClusteredCorrwithSuccess);
        MeanCausalRewardRelatedClustersCorrwithReward = nanmean(CausalRewardRelatedClustersCorrwithReward,1);
        MeanCausalRewardRelatedNonClusteredCorrwithReward = nanmean(CausalRewardRelatedNonClusteredCorrwithReward, 1);
        
        MeanFractionofCueSpinesThatAreClustered = nanmean(FractionofCueSpinesThatAreClustered,1);
        MeanFractionofMovementSpinesThatAreClustered = nanmean(FractionofMovementSpinesThatAreClustered,1);
        MeanFractionofPreSuccessSpinesThatAreClustered = nanmean(FractionofPreSuccessSpinesThatAreClustered,1);
        MeanFractionofSuccessSpinesThatAreClustered = nanmean(FractionofSuccessSpinesThatAreClustered,1);
        MeanFractionofMovementDuringCueSpinesThatAreClustered = nanmean(FractionofMovementDuringCueSpinesThatAreClustered,1);
        MeanFractionofRewardSpinesThatAreClustered = nanmean(FractionofRewardSpinesThatAreClustered,1);
        
    
        MeanSpatioTemporalDegree = cellfun(@nanmean, SpatioTemporalDegree);
        MeanSpatialMovementCorrelation = cellfun(@nanmean, SpatialMovementCorrelation);
        MeanTemporalMovementCorrelation = cellfun(@nanmean, TemporalMovementCorrelation);
        MeanSpatioTemporalMovementCorrelation = cellfun(@nanmean, SpatioTemporalMovementCorrelation);
        MeanDendSpatialFiedlerValue = cellfun(@(x) nanmean(x(:,1)),DendClusteringDegree);
        MeanDendTemporalFiedlerValue = cellfun(@(x) nanmean(x(:,2)),DendClusteringDegree);
        MeanDendSpatioTemporalFiedlerValue = cellfun(@(x) nanmean(x(:,3)),DendClusteringDegree);
        MeanSpatioTemporalOverlap = cellfun(@nanmean, SpatioTemporalOverlap);
        MeanSpatialDegree = cellfun(@nanmean, SpatialDegree);
        MeanTemporalDegree = cellfun(@nanmean, TemporalDegree);
        MeanSpatialDegreeofCueSpines = nanmean(SpatialDegreeofCueSpines,1);
        MeanTemporalDegreeofCueSpines = nanmean(TemporalDegreeofCueSpines,1);
        MeanSpatioTemporalDegreeofCueSpines = nanmean(SpatioTemporalDegreeofCueSpines,1);
        MeanSpatialDegreeofMovementSpines = nanmean(SpatialDegreeofMovementSpines,1);
        MeanTemporalDegreeofMovementSpines = nanmean(TemporalDegreeofMovementSpines,1);
        MeanSpatioTemporalDegreeofMovementSpines = nanmean(SpatioTemporalDegreeofMovementSpines,1);
        MeanSpatialDegreeofMovementDuringCueSpines = nanmean(SpatialDegreeofMovementDuringCueSpines,1);
        MeanTemporalDegreeofMovementDuringCueSpines = nanmean(TemporalDegreeofMovementDuringCueSpines,1);
        MeanSpatioTemporalDegreeofMovementDuringCueSpines = nanmean(SpatioTemporalDegreeofMovementDuringCueSpines,1);
        MeanSpatialDegreeofPreSuccessSpines = nanmean(SpatialDegreeofPreSuccessSpines,1);
        MeanTemporalDegreeofPreSuccessSpines = nanmean(TemporalDegreeofPreSuccessSpines,1);
        MeanSpatioTemporalDegreeofPreSuccessSpines = nanmean(SpatioTemporalDegreeofPreSuccessSpines,1);
        MeanSpatialDegreeofSuccessSpines = nanmean(SpatialDegreeofSuccessSpines,1);
        MeanTemporalDegreeofSuccessSpines = nanmean(TemporalDegreeofSuccessSpines,1);
        MeanSpatioTemporalDegreeofSuccessSpines = nanmean(SpatioTemporalDegreeofSuccessSpines,1);
        MeanSpatialDegreeofRewardSpines = nanmean(SpatialDegreeofRewardSpines,1);
        MeanTemporalDegreeofRewardSpines = nanmean(TemporalDegreeofRewardSpines,1);
        MeanSpatioTemporalDegreeofRewardSpines = nanmean(SpatioTemporalDegreeofRewardSpines,1);
       
        
        MeanClusterFreq = nanmean(ClusterFreq,1);
        MeanNonClusteredFreq = nanmean(NonClusteredFreq,1);
        MeanCueClusterFrequency = nanmean(CueClusterFrequency,1);
        MeanMovementClusterFrequency = nanmean(MovementClusterFrequency,1);
        MeanMovementDuringCueClusterFrequency = nanmean(MovementDuringCueClusterFrequency, 1);
        MeanPreSuccessClusterFrequency = nanmean(PreSuccessClusterFrequency,1);
        MeanSuccessClusterFrequency = nanmean(SuccessClusterFrequency,1);
        MeanRewardClusterFrequency = nanmean(RewardClusterFrequency,1);
        MeanCausalClusterFreq = nanmean(CausalClusterFreq,1);
        MeanNonClusteredCausalFreq = nanmean(NonClusteredCausalFreq,1);
        MeanCausalCueClusterFrequency = nanmean(CausalCueClusterFrequency,1);
        MeanCausalMovementClusterFrequency = nanmean(CausalMovementClusterFrequency,1);
        MeanCausalMovementDuringCueClusterFrequency = nanmean(CausalMovementDuringCueClusterFrequency,1);
        MeanCausalPreSuccessClusterFrequency = nanmean(CausalPreSuccessClusterFrequency,1);
        MeanCausalSuccessClusterFrequency = nanmean(CausalSuccessClusterFrequency,1);
        MeanCausalRewardClusterFrequency = nanmean(CausalRewardClusterFrequency,1);

        MeanClusteredSpineAmp = nanmean(ClusteredSpineAmp,1);
        MeanNonClusteredSpineAmp = nanmean(NonClusteredSpineAmp,1);
        MeanClusteredCueSpineAmp = nanmean(ClusteredCueSpineAmp,1);
        MeanClusteredMoveSpineAmp = nanmean(ClusteredMoveSpineAmp,1);
        MeanClusteredMovDuringCueSpineAmp = nanmean(ClusteredMovDuringCueSpineAmp,1);
        MeanClusteredPreSuccessSpineAmp = nanmean(ClusteredPreSuccessSpineAmp,1);
        MeanClusteredSuccessSpineAmp = nanmean(ClusteredSuccessSpineAmp,1);
        MeanClusteredRewardSpineAmp = nanmean(ClusteredRewardSpineAmp,1);
        MeanCausalClusteredSpineAmp = nanmean(CausalClusteredSpineAmp,1);
        MeanCausalNonClusteredSpineAmp = nanmean(CausalNonClusteredSpineAmp, 1);
        MeanCausalClusteredCueSpineAmp = nanmean(CausalClusteredCueSpineAmp,1);
        MeanCausalClusteredMoveSpineAmp = nanmean(CausalClusteredMoveSpineAmp,1);
        MeanCausalClusteredMovDuringCueSpineAmp = nanmean(CausalClusteredMovDuringCueSpineAmp,1);
        MeanCausalClusteredPreSuccessSpineAmp = nanmean(CausalClusteredPreSuccessSpineAmp,1);
        MeanCausalClusteredSuccessSpineAmp = nanmean(CausalClusteredSuccessSpineAmp,1);
        MeanCausalClusteredRewardSpineAmp = nanmean(CausalClusteredRewardSpineAmp,1);
        
        MeanClustDendFreq = nanmean(ClustDendFreq,1);
        MeanNonClustDendFreq = nanmean(NonClustDendFreq,1);
        MeanCueClustDendFreq = nanmean(CueClustDendFreq,1);
        MeanMovClustDendFreq = nanmean(MovClustDendFreq,1);
        MeanMovDuringCueClustDendFreq = nanmean(MovDuringCueClustDendFreq,1);
        MeanPreSucClustDendFreq = nanmean(PreSucClustDendFreq,1);
        MeanSucClustDendFreq = nanmean(SucClustDendFreq,1);
        MeanRewClustDendFreq = nanmean(RewClustDendFreq,1);
        MeanNonMovClustDendFreq = nanmean(NonMovClustDendFreq,1);
        
        MeanNumCueRelSpines = nanmean(NumCueRelSpines,1);
        MeanNumMovRelSpines = nanmean(NumMovRelSpines,1);
        MeanNumCueORMovRelSpines = nanmean(NumCueORMovRelSpines,1);
        MeanNumPreSucRelSpines = nanmean(NumPreSucRelSpines,1);
        MeanNumSucRelSpines = nanmean(NumSucRelSpines,1);
        MeanNumMovDuringCueRelSpines = nanmean(NumMovDuringCueRelSpines,1);
        MeanNumRewRelSpines = nanmean(NumRewRelSpines,1);
        MeanNumCausalMovSpines = nanmean(NumCausalMovSpines,1);
        MeanNumCausalSucSpines = nanmean(NumCausalSucSpines,1);
        MeanNumCausalCueSpines = nanmean(NumCausalCueSpines,1);
        
        MeanNumClustSpines = nanmean(NumClustSpines,1);
        MeanNumClustCueSpines = nanmean(NumClustCueSpines,1);
        MeanNumClustMovSpines = nanmean(NumClustMovSpines,1);
        MeanNumClustMixSpines = nanmean(NumClustMixSpines,1);
        MeanNumClustPreSucSpines = nanmean(NumClustPreSucSpines,1);
        MeanNumClustSucSpines = nanmean(NumClustSucSpines,1);
        MeanNumClustMovDuringCueSpines = nanmean(NumClustMovDuringCueSpines,1);
        MeanNumClustRewSpines = nanmean(NumClustRewSpines,1);
        MeanNumFarClustSpines = nanmean(NumFarClustSpines,1);
        MeanNumFarClustCueSpines = nanmean(NumFarClustCueSpines,1);
        MeanNumFarClustMovSpines = nanmean(NumFarClustMovSpines,1);
        MeanNumFarClustMixSpines = nanmean(NumFarClustMixSpines,1);
        MeanNumFarClustPreSucSpines = nanmean(NumFarClustPreSucSpines,1);
        MeanNumFarClustSucSpines = nanmean(NumFarClustSucSpines,1);
        MeanNumFarClustMovDuringCueSpines = nanmean(NumFarClustMovDuringCueSpines,1);
        MeanNumFarClustRewSpines = nanmean(NumFarClustRewSpines,1);
        MeanNumCausClustSpines = nanmean(NumCausClustSpines,1);
        MeanNumCausClustCueSpines = nanmean(NumCausClustCueSpines,1);
        MeanNumCausClustMovSpines = nanmean(NumCausClustMovSpines,1);
        MeanNumCausClustMovDuringCueSpines = nanmean(NumCausClustMovDuringCueSpines,1);
        MeanNumCausClustPreSucSpines = nanmean(NumCausClustPreSucSpines,1);
        MeanNumCausClustSucSpines = nanmean(NumCausClustSucSpines,1);
        MeanNumCausClustRewSpines = nanmean(NumCausClustRewSpines,1);

        
        MeanNumberofClusters = nanmean(NumberofClusters,1);
        MeanNumberofCausalClusters = nanmean(NumberofCausalClusters,1);
        MeanNumberofSpinesinEachCluster = nanmean(NumberofSpinesinEachCluster,1);
        MeanNumberofSpinesinEachCausalCluster = nanmean(NumberofSpinesinEachCausalCluster,1);
        MeanNumberofMovClusters = nanmean(NumberofMovClusters,1);
        MeanNumberofSpinesinEachMovCluster = nanmean(NumberofSpinesinEachMovCluster,1);
        
        MeanCueClusterLength = nanmean(CueClusterLength,1);
        MeanCueClusterMax = nanmean(CueClusterMax,1);
        MeanMovClusterLength = nanmean(MovClusterLength,1);
        MeanMovClusterMax = nanmean(MovClusterMax,1);
        MeanMixClusterLength = nanmean(MixClusterLength,1);
        MeanMixClusterMax = nanmean(MixClusterMax,1);
        MeanPreSucClusterLength = nanmean(PreSucClusterLength,1);
        MeanPreSucClusterMax = nanmean(PreSucClusterMax,1);
        MeanSucClusterLength = nanmean(SucClusterLength,1);
        MeanSucClusterMax = nanmean(SucClusterMax,1);
        MeanMovDuringCueClusterLength = nanmean(MovDuringCueClusterLength,1);
        MeanMovDuringCueClusterMax = nanmean(MovDuringCueClusterMax, 1);
        MeanRewClusterLength = nanmean(RewClusterLength,1);
        MeanRewClusterMax = nanmean(RewClusterMax,1);
        MeanAllClusterLength = nanmean(AllClusterLength,1);
        MeanAllClusterMax = nanmean(AllClusterMax,1);
        MeanCausalCueClusterLength = nanmean(CausalCueClusterLength,1);
        MeanCausalCueClusterMax = nanmean(CausalCueClusterMax,1);
        MeanCausalMovClusterLength = nanmean(CausalMovClusterLength,1);
        MeanCausalMovClusterMax = nanmean(CausalMovClusterMax,1);
        MeanCausalSucClusterLength = nanmean(CausalSucClusterLength,1);
        MeanCausalSucClusterMax = nanmean(CausalSucClusterMax,1);
        MeanCausalRewClusterLength = nanmean(CausalRewClusterLength,1);
        MeanCausalRewClusterMax = nanmean(CausalRewClusterMax,1);
        MeanAllCausalClusterLength = nanmean(AllCausalClusterLength,1);
        MeanAllCausalClusterMax = nanmean(AllCausalClusterMax,1);
        
        MeanFarCueClusterLength = nanmean(FarCueClusterLength,1);
        MeanFarMovClusterLength = nanmean(FarMovClusterLength,1);
        MeanFarMixClusterLength = nanmean(FarMixClusterLength,1);
        MeanFarPreSucClusterLength = nanmean(FarPreSucClusterLength,1);
        MeanFarSucClusterLength = nanmean(FarSucClusterLength,1);
        MeanFarMovDuringCueClusterLength = nanmean(FarMovDuringCueClusterLength,1);
        MeanFarRewClusterLength = nanmean(FarRewClusterLength,1);
        MeanAllFarClusterLength = nanmean(AllFarClusterLength,1);

        
        MeanDistanceBetweenCueSpines = nanmean(DistanceBetweenCueSpines,1);
        MeanDistanceBetweenMovementSpines = nanmean(DistanceBetweenMovementSpines,1);
        MeanDistanceBetweenPreSuccessSpines = nanmean(DistanceBetweenPreSuccessSpines,1);
        MeanDistanceBetweenSuccessSpines = nanmean(DistanceBetweenSuccessSpines,1);
        MeanDistanceBetweenMovementDuringCueSpines = nanmean(DistanceBetweenMovementDuringCueSpines,1);
        MeanDistanceBetweenRewardSpines = nanmean(DistanceBetweenRewardSpines,1);
        
        MeanClusteredSpines_CorrwithDend = nanmean(ClusteredSpines_CorrwithDend,1);
        MeanFilteredClusteredSpines_CorrwithDend = nanmean(FilteredClusteredSpines_CorrwithDend,1);
        MeanNonClusteredSpines_CorrwithDend = nanmean(NonClusteredSpines_CorrwithDend,1);
        MeanCausalClusteredSpines_CorrwithDend = nanmean(CausalClusteredSpines_CorrwithDend,1);
        MeanFilteredCausalClusteredSpines_CorrwithDend = nanmean(FilteredCausalClusteredSpines_CorrwithDend,1);
        MeanNonCausalClusteredSpines_CorrwithDend = nanmean(NonCausalClusteredSpines_CorrwithDend,1);
        MeanCueRelClusteredSpines_CorrwithDend = nanmean(CueRelClusteredSpines_CorrwithDend,1);
        MeanMovRelClusteredSpines_CorrwithDend = nanmean(MovRelClusteredSpines_CorrwithDend,1);
        MeanSucRelClusteredSpines_CorrwithDend = nanmean(SucRelClusteredSpines_CorrwithDend,1);
        MeanRewRelClusteredSpines_CorrwithDend = nanmean(RewRelClusteredSpines_CorrwithDend,1);
        MeanCueRelCausalClusteredSpines_CorrwithDend = nanmean(CueRelCausalClusteredSpines_CorrwithDend,1);
        MeanMovRelCausalClusteredSpines_CorrwithDend = nanmean(MovRelCausalClusteredSpines_CorrwithDend,1);
        MeanSucRelCausalClusteredSpines_CorrwithDend = nanmean(SucRelCausalClusteredSpines_CorrwithDend,1);
        MeanRewRelCausalClusteredSpines_CorrwithDend = nanmean(RewRelCausalClusteredSpines_CorrwithDend,1);
        
        AllMeanCorrelationBetweenFarMovementSpines = nanmean(MeanCorrelationBetweenFarMovementSpines,1);
        
%         MeanFractionofClusterThatsMovementRelated = nanmean(FractionofClusterThatsMovementRelated,1);
%         MeanFractionofCausalClusterThatsMovementRelated = nanmean(FractionofCausalClusterThatsMovementRelated,1);

        MeanPercentofCueRelatedDendrites = nanmean(PercentofCueRelatedDendrites,1);
        MeanPercentofMovementRelatedDendrites = nanmean(PercentofMovementRelatedDendrites,1);
        MeanPercentofPreSuccessRelatedDendrites = nanmean(PercentofPreSuccessRelatedDendrites,1);
        MeanPercentofSuccessRelatedDendrites = nanmean(PercentofSuccessRelatedDendrites,1);
        MeanPercentofMovementDuringCueRelatedDendrites = nanmean(PercentofMovementDuringCueRelatedDendrites, 1);
        MeanPercentofRewardRelatedDendrites = nanmean(PercentofRewardRelatedDendrites,1);

    for i = 1:14
        MeanSpatialDegree(1,i) = nanmean(SpatialDegree{i});
        MeanTemporalDegree(1,i) = nanmean(TemporalDegree{i});
    end
    
    %%%
    %%% Determine error of input data
    %%%
    
    for i = 1:14
        AllClustersCorrwithMovementSEM(1,i) = nanstd(AllClustersCorrwithMovement(:,i),0,1)/sqrt(sum(~isnan(AllClustersCorrwithMovement(:,i))));
        NonClusteredCorrwithMovementSEM(1,i) = nanstd(NonClusteredCorrwithMovement(:,i),0,1)/sqrt(sum(~isnan(NonClusteredCorrwithMovement(:,i))));
        AllClustersCorrwithCueSEM(1,i) = nanstd(AllClustersCorrwithCue(:,i),0,1)/sqrt(sum(~isnan(AllClustersCorrwithCue(:,i))));
        NonClusteredCorrwithCueSEM(1,i) = nanstd(NonClusteredCorrwithCue(:,i),0,1)/sqrt(sum(~isnan(NonClusteredCorrwithCue(:,i))));
        AllClustersCorrwithMDCSEM(1,i) = nanstd(AllClustersCorrwithMDC(:,i),0,1)/sqrt(sum(~isnan(AllClustersCorrwithMDC(:,i))));
        NonClusteredCorrwithMDCSEM(1,i) = nanstd(NonClusteredCorrwithMDC(:,i),0,1)/sqrt(sum(~isnan(NonClusteredCorrwithCue(:,i))));
        AllClustersCorrwithSuccessSEM(1,i) = nanstd(AllClustersCorrwithSuccess(:,i),0,1)/sqrt(sum(~isnan(AllClustersCorrwithSuccess(:,i))));
        NonClusteredCorrwithSuccessSEM(1,i) = nanstd(NonClusteredCorrwithSuccess(:,i),0,1)/sqrt(sum(~isnan(NonClusteredCorrwithSuccess(:,i))));
        AllClustCorrwithRewardSEM(1,i) = nanstd(AllClustCorrwithReward(:,i),0,1)/sqrt(sum(~isnan(AllClustCorrwithReward(:,i))));
        NonClusteredCorrwithRewardSEM(1,i) = nanstd(NonClusteredCorrwithReward(:,i),0,1)/sqrt(sum(~isnan(NonClusteredCorrwithReward(:,i))));
        
        CorrelationofClustersSEM(1,i) = nanstd(CorrelationofClusters{i},0,1)/sqrt(sum(~isnan(CorrelationofClusters{i})));
        
        AllCausalClustersCorrwithMovementSEM(1,i) = nanstd(AllCausalClustersCorrwithMovement(:,i),0,1)/sqrt(sum(~isnan(AllCausalClustersCorrwithMovement(:,i))));
        CausalNonClusteredCorrwithMovementSEM(1,i) = nanstd(CausalNonClusteredCorrwithMovement(:,i),0,1)/sqrt(sum(~isnan(CausalNonClusteredCorrwithMovement(:,i))));
        AllCausalClustersCorrwithCueSEM(1,i) = nanstd(AllCausalClustersCorrwithCue(:,i),0,1)/sqrt(sum(~isnan(AllCausalClustersCorrwithCue(:,i))));
        CausalNonClusteredCorrwithCueSEM(1,i) = nanstd(CausalNonClusteredCorrwithCue(:,i),0,1)/sqrt(sum(~isnan(CausalNonClusteredCorrwithCue(:,i))));
        AllCausalClustersCorrwithMDCSEM(1,i) = nanstd(AllCausalClustersCorrwithMDC(:,i),0,1)/sqrt(sum(~isnan(AllCausalClustersCorrwithMDC(:,i))));
        CausalNonClusteredCorrwithMDCSEM(1,i) = nanstd(CausalNonClusteredCorrwithMDC(:,i),0,1)/sqrt(sum(~isnan(CausalNonClusteredCorrwithMDC(:,i))));
        AllCausalClustersCorrwithSuccessSEM(1,i) = nanstd(AllCausalClustersCorrwithSuccess(:,i),0,1)/sqrt(sum(~isnan(AllCausalClustersCorrwithSuccess(:,i))));
        CausalNonClusteredCorrwithSuccessSEM(1,i) = nanstd(CausalNonClusteredCorrwithSuccess(:,i),0,1)/sqrt(sum(~isnan(CausalNonClusteredCorrwithSuccess(:,i))));
        AllCausalClustCorrwithRewardSEM(1,i) = nanstd(AllCausalClustCorrwithReward(:,i),0,1)/sqrt(sum(~isnan(AllCausalClustCorrwithReward(:,i))));
        CausalNonClusteredCorrwithRewardSEM(1,i) = nanstd(CausalNonClusteredCorrwithReward(:,i),0,1)/sqrt(sum(~isnan(CausalNonClusteredCorrwithReward(:,i))));
        
        CueRelatedClustersCorrwithCueSEM(1,i) = nanstd(CueRelatedClustersCorrwithCue(:,i),0,1)/sqrt(sum(~isnan(CueRelatedClustersCorrwithCue(:,i))));
        CueRelatedNonClusteredCorrwithCueSEM(1,i) = nanstd(CueRelatedNonClusteredCorrwithCue(:,i),0,1)/sqrt(sum(~isnan(CueRelatedNonClusteredCorrwithCue(:,i))));
        MovementRelatedClustersCorrwithMovementSEM(1,i) = nanstd(MovementRelatedClustersCorrwithMovement(:,i),0,1)/sqrt(sum(~isnan(MovementRelatedClustersCorrwithMovement(:,i))));
        MovementRelatedNonClusteredCorrwithMovementSEM(1,i) = nanstd(MovementRelatedNonClusteredCorrwithMovement(:,i),0,1)/sqrt(sum(~isnan(MovementRelatedNonClusteredCorrwithMovement(:,i))));
        MDCRelatedClustersCorrwithMDCSEM(1,i) = nanstd(MDCRelatedClustersCorrwithMDC(:,i),0,1)/sqrt(sum(~isnan(MDCRelatedClustersCorrwithMDC(:,i))));
        MDCRelatedNonClusteredCorrwithMDCSEM(1,i) = nanstd(MDCRelatedNonClusteredCorrwithMDC(:,i),0,1)/sqrt(sum(~isnan(MDCRelatedNonClusteredCorrwithMDC(:,i))));
        SuccessRelatedClustersCorrwithSuccessSEM(1,i) = nanstd(SuccessRelatedClustersCorrwithSuccess(:,i),0,1)/sqrt(sum(~isnan(SuccessRelatedClustersCorrwithSuccess(:,i))));
        SuccessRelatedNonClusteredCorrwithSuccessSEM(1,i) = nanstd(SuccessRelatedNonClusteredCorrwithSuccess(:,i),0,1)/sqrt(sum(~isnan(SuccessRelatedNonClusteredCorrwithSuccess(:,i))));
        RewardRelatedClustersCorrwithRewardSEM(1,i) = nanstd(RewardRelatedClustersCorrwithReward(:,i),0,1)/sqrt(sum(~isnan(RewardRelatedClustersCorrwithReward(:,i))));
        RewardRelatedNonClusteredCorrwithRewardSEM(1,i) = nanstd(RewardRelatedNonClusteredCorrwithReward(:,i),0,1)/sqrt(sum(~isnan(RewardRelatedNonClusteredCorrwithReward(:,i))));
        
        CausalCueRelatedClustersCorrwithCueSEM(1,i) = nanstd(CausalCueRelatedClustersCorrwithCue(:,i),0,1)/sqrt(sum(~isnan(CausalCueRelatedClustersCorrwithCue(:,i))));
        CausalCueRelatedNonClusteredCorrwithCueSEM(1,i) = nanstd(CausalCueRelatedNonClusteredCorrwithCue(:,i),0,1)/sqrt(sum(~isnan(CausalCueRelatedNonClusteredCorrwithCue(:,i))));
        CausalMovementRelatedClustersCorrwithMovementSEM(1,i) = nanstd(CausalMovementRelatedClustersCorrwithMovement(:,i),0,1)/sqrt(sum(~isnan(CausalMovementRelatedClustersCorrwithMovement(:,i))));
        CausalMovementRelatedNonClusteredCorrwithMovementSEM(1,i) = nanstd(CausalMovementRelatedNonClusteredCorrwithMovement(:,i),0,1)/sqrt(sum(~isnan(CausalMovementRelatedNonClusteredCorrwithMovement(:,i))));
        CausalMDCRelatedClustersCorrwithMDCSEM(1,i) = nanstd(CausalMDCRelatedClustersCorrwithMDC(:,i),0,1)/sqrt(sum(~isnan(CausalMDCRelatedClustersCorrwithMDC(:,i))));
        CausalMDCRelatedNonClusteredCorrwithMDCSEM(1,i) = nanstd(CausalMDCRelatedNonClusteredCorrwithMDC(:,i),0,1)/sqrt(sum(~isnan(CausalMDCRelatedNonClusteredCorrwithMDC(:,i))));
        CausalSuccessRelatedClustersCorrwithSuccessSEM(1,i) = nanstd(CausalSuccessRelatedClustersCorrwithSuccess(:,i),0,1)/sqrt(sum(~isnan(CausalSuccessRelatedClustersCorrwithSuccess(:,i))));
        CausalSuccessRelatedNonClusteredCorrwithSuccessSEM(1,i) = nanstd(CausalSuccessRelatedNonClusteredCorrwithSuccess(:,i),0,1)/sqrt(sum(~isnan(CausalSuccessRelatedNonClusteredCorrwithSuccess(:,i))));
        CausalRewardRelatedClustersCorrwithRewardSEM(1,i) = nanstd(CausalRewardRelatedClustersCorrwithReward(:,i),0,1)/sqrt(sum(~isnan(CausalRewardRelatedClustersCorrwithReward(:,i))));
        CausalRewardRelatedNonClusteredCorrwithRewardSEM(1,i) = nanstd(CausalRewardRelatedNonClusteredCorrwithReward(:,i),0,1)/sqrt(sum(~isnan(CausalRewardRelatedNonClusteredCorrwithReward(:,i))));
        
        FractionofCueSpinesThatAreClusteredSEM(1,i) = nanstd(FractionofCueSpinesThatAreClustered(:,i),0,1)/sqrt(sum(~isnan(FractionofCueSpinesThatAreClustered(:,i))));
        FractionofMovementSpinesThatAreClusteredSEM(1,i) = nanstd(FractionofMovementSpinesThatAreClustered(:,i),0,1)/sqrt(sum(~isnan(FractionofMovementSpinesThatAreClustered(:,i))));
        FractionofPreSuccessSpinesThatAreClusteredSEM(1,i) = nanstd(FractionofPreSuccessSpinesThatAreClustered(:,i),0,1)/sqrt(sum(~isnan(FractionofPreSuccessSpinesThatAreClustered(:,i))));
        FractionofSuccessSpinesThatAreClusteredSEM(1,i) = nanstd(FractionofSuccessSpinesThatAreClustered(:,i),0,1)/sqrt(sum(~isnan(FractionofSuccessSpinesThatAreClustered(:,i))));
        FractionofMovementDuringCueSpinesThatAreClusteredSEM(1,i) = nanstd(FractionofMovementDuringCueSpinesThatAreClustered(:,i),0,1)/sqrt(sum(~isnan(FractionofMovementDuringCueSpinesThatAreClustered(:,i))));
        FractionofRewardSpinesThatAreClusteredSEM(1,i) = nanstd(FractionofRewardSpinesThatAreClustered(:,i),0,1)/sqrt(sum(~isnan(FractionofRewardSpinesThatAreClustered(:,i))));
        
        SpatioTemporalDegreeSEM(1,i) = nanstd(SpatioTemporalDegree{i},0,1)/sqrt(sum(~isnan(SpatioTemporalDegree{i})));
        SpatialMovementCorrelationSEM(1,i) = nanstd(SpatialMovementCorrelation{i},0,1)/sqrt(sum(~isnan(SpatialMovementCorrelation{i})));
        TemporalMovementCorrelationSEM(1,i) = nanstd(TemporalMovementCorrelation{i},0,1)/sqrt(sum(~isnan(TemporalMovementCorrelation{i})));
        SpatioTemporalMovementCorrelationSEM(1,i) = nanstd(SpatioTemporalMovementCorrelation{i},0,1)/sqrt(sum(~isnan(SpatioTemporalMovementCorrelation{i})));
        DendSpatialFiedlerValueSEM(1,i) = nanstd(DendClusteringDegree{i}(:,1),0,1)/sqrt(sum(~isnan(DendClusteringDegree{i}(:,1))));
        DendTemporalFiedlerValueSEM(1,i) = nanstd(DendClusteringDegree{i}(:,2),0,1)/sqrt(sum(~isnan(DendClusteringDegree{i}(:,2))));
        DendSpatioTemporalFiedlerValueSEM(1,i) = nanstd(DendClusteringDegree{i}(:,3),0,1)/sqrt(sum(~isnan(DendClusteringDegree{i}(:,3))));
        SpatialDegreeSEM(1,i) = nanstd(SpatialDegree{i},0,1)/sqrt(sum(~isnan(SpatialDegree{i})));
        TemporalDegreeSEM(1,i) = nanstd(TemporalDegree{i},0,1)/sqrt(sum(~isnan(TemporalDegree{i})));
        
        ClusterFreqSEM(1,i) = nanstd(ClusterFreq(:,i),0,1)/sqrt(sum(~isnan(ClusterFreq(:,i))));
        NonClusteredFreqSEM(1,i) = nanstd(NonClusteredFreq(:,i),0,1)/sqrt(sum(~isnan(NonClusteredFreq(:,i))));
        CueClusterFrequencySEM(1,i) = nanstd(CueClusterFrequency(:,i),0,1)/sqrt(sum(~isnan(CueClusterFrequency(:,i))));
        MovementClusterFrequencySEM(1,i) = nanstd(MovementClusterFrequency(:,i),0,1)/sqrt(sum(~isnan(MovementClusterFrequency(:,i))));
        MovementDuringCueClusterFrequencySEM(1,i) = nanstd(MovementDuringCueClusterFrequency(:,i),0,1)/sqrt(sum(~isnan(MovementDuringCueClusterFrequency(:,i))));
        PreSuccessClusterFrequencySEM(1,i) = nanstd(PreSuccessClusterFrequency(:,i),0,1)/sqrt(sum(~isnan(PreSuccessClusterFrequency(:,i))));
        SuccessClusterFrequencySEM(1,i) = nanstd(SuccessClusterFrequency(:,i),0,1)/sqrt(sum(~isnan(SuccessClusterFrequency(:,i))));
        RewardClusterFrequencySEM(1,i) = nanstd(RewardClusterFrequency(:,i),0,1)/sqrt(sum(~isnan(RewardClusterFrequency(:,i))));
        CausalClusterFreqSEM(1,i) = nanstd(CausalClusterFreq(:,i),0,1)/sqrt(sum(~isnan(CausalClusterFreq(:,i))));
        NonClusteredCausalFreqSEM(1,i) = nanstd(NonClusteredCausalFreq(:,i),0,1)/sqrt(sum(~isnan(NonClusteredCausalFreq(:,i))));
        CausalCueClusterFrequencySEM(1,i) = nanstd(CausalCueClusterFrequency(:,i),0,1)/sqrt(sum(~isnan(CausalCueClusterFrequency(:,i))));
        CausalMovementClusterFrequencySEM(1,i) = nanstd(CausalMovementClusterFrequency(:,i),0,1)/sqrt(sum(~isnan(CausalMovementClusterFrequency(:,i))));
        CausalMovementDuringCueClusterFrequencySEM(1,i) = nanstd(CausalMovementDuringCueClusterFrequency(:,i),0,1)/sqrt(sum(~isnan(CausalMovementDuringCueClusterFrequency(:,i))));
        CausalPreSuccessClusterFrequencySEM(1,i) = nanstd(CausalPreSuccessClusterFrequency(:,i),0,1)/sqrt(sum(~isnan(CausalPreSuccessClusterFrequency(:,i))));
        CausalSuccessClusterFrequencySEM(1,i) = nanstd(CausalSuccessClusterFrequency(:,i),0,1)/sqrt(sum(~isnan(CausalSuccessClusterFrequency(:,i))));
        CausalRewardClusterFrequencySEM(1,i) = nanstd(CausalRewardClusterFrequency(:,i),0,1)/sqrt(sum(~isnan(CausalRewardClusterFrequency(:,i))));

        ClusteredSpineAmpSEM(1,i) = nanstd(ClusteredSpineAmp(:,i),0,1)/sqrt(sum(~isnan(ClusteredSpineAmp(:,i))));
        NonClusteredSpineAmpSEM(1,i) = nanstd(NonClusteredSpineAmp(:,i),0,1)/sqrt(sum(~isnan(NonClusteredSpineAmp(:,i))));
        ClusteredCueSpineAmpSEM(1,i) = nanstd(ClusteredCueSpineAmp(:,i),0,1)/sqrt(sum(~isnan(ClusteredCueSpineAmp(:,i))));
        ClusteredMoveSpineAmpSEM(1,i) = nanstd(ClusteredMoveSpineAmp(:,i),0,1)/sqrt(sum(~isnan(ClusteredMoveSpineAmp(:,i))));
        ClusteredSuccessSpineAmpSEM(1,i) = nanstd(ClusteredSuccessSpineAmp(:,i),0,1)/sqrt(sum(~isnan(ClusteredSuccessSpineAmp(:,i))));
        ClusteredRewardSpineAmpSEM(1,i) = nanstd(ClusteredRewardSpineAmp(:,i),0,1)/sqrt(sum(~isnan(ClusteredRewardSpineAmp(:,i))));
        CausalClusteredSpineAmpSEM(1,i) = nanstd(CausalClusteredSpineAmp(:,i),0,1)/sqrt(sum(~isnan(CausalClusteredSpineAmp(:,i))));
        CausalNonClusteredSpineAmpSEM(1,i) = nanstd(CausalNonClusteredSpineAmp(:,i),0,1)/sqrt(sum(~isnan(CausalNonClusteredSpineAmp(:,i))));
        CausalClusteredCueSpineAmpSEM(1,i) = nanstd(CausalClusteredCueSpineAmp(:,i),0,1)/sqrt(sum(~isnan(CausalClusteredCueSpineAmp(:,i))));
        CausalClusteredMoveSpineAmpSEM(1,i) = nanstd(CausalClusteredMoveSpineAmp(:,i),0,1)/sqrt(sum(~isnan(CausalClusteredMoveSpineAmp(:,i))));
        CausalClusteredSuccessSpineAmpSEM(1,i) = nanstd(CausalClusteredSuccessSpineAmp(:,i),0,1)/sqrt(sum(~isnan(CausalClusteredSuccessSpineAmp(:,i))));
        CausalClusteredRewardSpineAmpSEM(1,i) = nanstd(CausalClusteredRewardSpineAmp(:,i),0,1)/sqrt(sum(~isnan(CausalClusteredRewardSpineAmp(:,i))));
        
        ClustDendFreqSEM(1,i) = nanstd(ClustDendFreq(:,i),0,1)/sqrt(sum(~isnan(ClustDendFreq(:,i))));
        NonClustDendFreqSEM(1,i) = nanstd(NonClustDendFreq(:,i),0,1)/sqrt(sum(~isnan(NonClustDendFreq(:,i))));
        CueClustDendFreqSEM(1,i) = nanstd(CueClustDendFreq(:,i),0,1)/sqrt(sum(~isnan(CueClustDendFreq(:,i))));
        MovClustDendFreqSEM(1,i) = nanstd(MovClustDendFreq(:,i),0,1)/sqrt(sum(~isnan(MovClustDendFreq(:,i))));
        MovDuringCueClustDendFreqSEM(1,i) = nanstd(MovDuringCueClustDendFreq(:,i),0,1)/sqrt(sum(~isnan(MovDuringCueClustDendFreq(:,i))));
        PreSucClustDendFreqSEM(1,i) = nanstd(PreSucClustDendFreq(:,i),0,1)/sqrt(sum(~isnan(PreSucClustDendFreq(:,i))));
        SucClustDendFreqSEM(1,i) = nanstd(SucClustDendFreq(:,i),0,1)/sqrt(sum(~isnan(SucClustDendFreq(:,i))));
        RewClustDendFreqSEM(1,i) = nanstd(RewClustDendFreq(:,i),0,1)/sqrt(sum(~isnan(RewClustDendFreq(:,i))));
        NonMovClustDendFreqSEM(1,i) = nanstd(NonMovClustDendFreq(:,i),0,1)/sqrt(sum(~isnan(NonMovClustDendFreq(:,i))));
        
        NumCueRelSpinesSEM(1,i) = nanstd(NumCueRelSpines(:,i),0,1)/sqrt(sum(~isnan(NumCueRelSpines(:,i))));
        NumMovRelSpinesSEM(1,i) = nanstd(NumMovRelSpines(:,i),0,1)/sqrt(sum(~isnan(NumMovRelSpines(:,i))));
        NumCueORMovRelSpinesSEM(1,i) = nanstd(NumCueORMovRelSpines(:,i),0,1)/sqrt(sum(~isnan(NumCueORMovRelSpines(:,i))));
        NumPreSucRelSpinesSEM(1,i) = nanstd(NumPreSucRelSpines(:,i),0,1)/sqrt(sum(~isnan(NumPreSucRelSpines(:,i))));
        NumSucRelSpinesSEM(1,i) = nanstd(NumSucRelSpines(:,i),0,1)/sqrt(sum(~isnan(NumSucRelSpines(:,i))));
        NumMovDuringCueRelSpinesSEM(1,i) = nanstd(NumMovDuringCueRelSpines(:,i),0,1)/sqrt(sum(~isnan(NumMovDuringCueRelSpines(:,i))));
        NumRewRelSpinesSEM(1,i) = nanstd(NumRewRelSpines(:,i),0,1)/sqrt(sum(~isnan(NumRewRelSpines(:,i))));
        NumCausalMovSpinesSEM(1,i) = nanstd(NumCausalMovSpines(:,i),0,1)/sqrt(sum(~isnan(NumCausalMovSpines(:,i))));
        NumCausalSucSpinesSEM(1,i) = nanstd(NumCausalSucSpines(:,i),0,1)/sqrt(sum(~isnan(NumCausalSucSpines(:,i))));
        NumCausalCueSpinesSEM(1,i) = nanstd(NumCausalCueSpines(:,i),0,1)/sqrt(sum(~isnan(NumCausalCueSpines(:,i))));
        
        NumClustSpinesSEM(1,i) = nanstd(NumClustSpines(:,i),0,1)/sqrt(sum(~isnan(NumClustSpines(:,i))));
        NumClustCueSpinesSEM(1,i) = nanstd(NumClustCueSpines(:,i),0,1)/sqrt(sum(~isnan(NumClustCueSpines(:,i))));
        NumClustMovSpinesSEM(1,i) = nanstd(NumClustMovSpines(:,i),0,1)/sqrt(sum(~isnan(NumClustMovSpines(:,i))));
        NumClustMixSpinesSEM(1,i) = nanstd(NumClustMixSpines(:,i),0,1)/sqrt(sum(~isnan(NumClustMixSpines(:,i))));
        NumClustPreSucSpinesSEM(1,i) = nanstd(NumClustPreSucSpines(:,i),0,1)/sqrt(sum(~isnan(NumClustPreSucSpines(:,i))));
        NumClustSucSpinesSEM(1,i) = nanstd(NumClustSucSpines(:,i),0,1)/sqrt(sum(~isnan(NumClustSucSpines(:,i))));
        NumClustMovDuringCueSpinesSEM(1,i) = nanstd(NumClustMovDuringCueSpines(:,i),0,1)/sqrt(sum(~isnan(NumClustMovDuringCueSpines(:,i))));
        NumClustRewSpinesSEM(1,i) = nanstd(NumClustRewSpines(:,i),0,1)/sqrt(sum(~isnan(NumClustRewSpines(:,i))));
        
        NumFarClustSpinesSEM(1,i) = nanstd(NumFarClustSpines(:,i),0,1)/sqrt(sum(~isnan(NumFarClustSpines(:,i))));
        NumFarClustCueSpinesSEM(1,i) = nanstd(NumFarClustCueSpines(:,i),0,1)/sqrt(sum(~isnan(NumFarClustCueSpines(:,i))));
        NumFarClustMovSpinesSEM(1,i) = nanstd(NumFarClustMovSpines(:,i),0,1)/sqrt(sum(~isnan(NumFarClustMovSpines(:,i))));
        NumFarClustMixSpinesSEM(1,i) = nanstd(NumFarClustMixSpines(:,i),0,1)/sqrt(sum(~isnan(NumFarClustMixSpines(:,i))));
        NumFarClustPreSucSpinesSEM(1,i) = nanstd(NumFarClustPreSucSpines(:,i),0,1)/sqrt(sum(~isnan(NumFarClustPreSucSpines(:,i))));
        NumFarClustSucSpinesSEM(1,i) = nanstd(NumFarClustSucSpines(:,i),0,1)/sqrt(sum(~isnan(NumFarClustSucSpines(:,i))));
        NumFarClustMovDuringCueSpinesSEM(1,i) = nanstd(NumFarClustMovDuringCueSpines(:,i),0,1)/sqrt(sum(~isnan(NumFarClustMovDuringCueSpines(:,i))));
        NumFarClustRewSpinesSEM(1,i) = nanstd(NumFarClustRewSpines(:,i),0,1)/sqrt(sum(~isnan(NumFarClustRewSpines(:,i))));

        NumCausClustSpinesSEM(1,i) = nanstd(NumCausClustSpines(:,i),0,1)/sqrt(sum(~isnan(NumCausClustSpines(:,i))));
        NumCausClustCueSpinesSEM(1,i) = nanstd(NumCausClustCueSpines(:,i),0,1)/sqrt(sum(~isnan(NumCausClustCueSpines(:,i))));
        NumCausClustMovSpinesSEM(1,i) = nanstd(NumCausClustMovSpines(:,i),0,1)/sqrt(sum(~isnan(NumCausClustMovSpines(:,i))));
        NumCausClustMovDuringCueSpinesSEM(1,i) = nanstd(NumCausClustMovDuringCueSpines(:,i),0,1)/sqrt(sum(~isnan(NumCausClustMovDuringCueSpines(:,i))));
        NumCausClustPreSucSpinesSEM(1,i) = nanstd(NumCausClustPreSucSpines(:,i),0,1)/sqrt(sum(~isnan(NumCausClustPreSucSpines(:,i))));
        NumCausClustSucSpinesSEM(1,i) = nanstd(NumCausClustSucSpines(:,i),0,1)/sqrt(sum(~isnan(NumCausClustSucSpines(:,i))));
        NumCausClustRewSpinesSEM(1,i) = nanstd(NumCausClustRewSpines(:,i),0,1)/sqrt(sum(~isnan(NumCausClustRewSpines(:,i))));

        NumberofClustersSEM(1,i) = nanstd(NumberofClusters(:,i),0,1)/sqrt(sum(~isnan(NumberofClusters(:,i))));
        NumberofCausalClustersSEM(1,i) = nanstd(NumberofCausalClusters(:,i),0,1)/sqrt(sum(~isnan(NumberofCausalClusters(:,i))));
        NumberofSpinesinEachClusterSEM(1,i) = nanstd(NumberofSpinesinEachCluster(:,i),0,1)/sqrt(sum(~isnan(NumberofSpinesinEachCluster(:,i))));
        NumberofSpinesinEachCausalClusterSEM(1,i) = nanstd(NumberofSpinesinEachCausalCluster(:,i),0,1)/sqrt(sum(~isnan(NumberofSpinesinEachCausalCluster(:,i))));
        NumberofMovClustersSEM(1,i) = nanstd(NumberofMovClusters(:,i),0,1)/sqrt(sum(~isnan(NumberofMovClusters(:,i))));
        NumberofSpinesinEachMovClusterSEM(1,i) = nanstd(NumberofSpinesinEachMovCluster(:,i),0,1)/sqrt(sum(~isnan(NumberofSpinesinEachMovCluster(:,i))));
        
        CueClusterLengthSEM(1,i) = nanstd(CueClusterLength(:,i),0,1)/sqrt(sum(~isnan(CueClusterLength(:,i))));
        CueClusterMaxSEM(1,i) = nanstd(CueClusterMax(:,i),0,1)/sqrt(sum(~isnan(CueClusterMax(:,i))));
        MovClusterLengthSEM(1,i) = nanstd(MovClusterLength(:,i),0,1)/sqrt(sum(~isnan(MovClusterLength(:,i))));
        MovClusterMaxSEM(1,i) = nanstd(MovClusterMax(:,i),0,1)/sqrt(sum(~isnan(MovClusterMax(:,i))));
        MixClusterLengthSEM(1,i) = nanstd(MixClusterLength(:,i),0,1)/sqrt(sum(~isnan(MixClusterLength(:,i))));
        MixClusterMaxSEM(1,i) = nanstd(MixClusterMax(:,i),0,1)/sqrt(sum(~isnan(MixClusterMax(:,i))));
        PreSucClusterLengthSEM(1,i) = nanstd(PreSucClusterLength(:,i),0,1)/sqrt(sum(~isnan(PreSucClusterLength(:,i))));
        PreSucClusterMaxSEM(1,i) = nanstd(PreSucClusterMax(:,i),0,1)/sqrt(sum(~isnan(PreSucClusterMax(:,i))));
        SucClusterLengthSEM(1,i) = nanstd(SucClusterLength(:,i),0,1)/sqrt(sum(~isnan(SucClusterLength(:,i))));
        SucClusterMaxSEM(1,i) = nanstd(SucClusterMax(:,i),0,1)/sqrt(sum(~isnan(SucClusterMax(:,i))));
        MovDuringCueClusterLengthSEM(1,i) = nanstd(MovDuringCueClusterLength(:,i),0,1)/sqrt(sum(~isnan(MovDuringCueClusterLength(:,i))));
        MovDuringCueClusterMaxSEM(1,i) = nanstd(MovDuringCueClusterMax(:,i),0,1)/sqrt(sum(~isnan(MovDuringCueClusterMax(:,i))));
        RewClusterLengthSEM(1,i) = nanstd(RewClusterLength(:,i),0,1)/sqrt(sum(~isnan(RewClusterLength(:,i))));
        RewClusterMaxSEM(1,i) = nanstd(RewClusterMax(:,i),0,1)/sqrt(sum(~isnan(RewClusterMax(:,i))));
        AllClusterLengthSEM(1,i) = nanstd(AllClusterLength(:,i),0,1)/sqrt(sum(~isnan(AllClusterLength(:,i))));
        AllClusterMaxSEM(1,i) = nanstd(AllClusterMax(:,i),0,1)/sqrt(sum(~isnan(AllClusterMax(:,i))));
        CausalCueClusterLengthSEM(1,i) = nanstd(CausalCueClusterLength(:,i),0,1)/sqrt(sum(~isnan(CausalCueClusterLength(:,i))));
        CausalCueClusterMaxSEM(1,i) = nanstd(CausalCueClusterMax(:,i),0,1)/sqrt(sum(~isnan(CausalCueClusterMax(:,i))));
        CausalMovClusterLengthSEM(1,i) = nanstd(CausalMovClusterLength(:,i),0,1)/sqrt(sum(~isnan(CausalMovClusterLength(:,i))));
        CausalMovClusterMaxSEM(1,i) = nanstd(CausalMovClusterMax(:,i),0,1)/sqrt(sum(~isnan(CausalMovClusterMax(:,i))));
        AllCausalClusterLengthSEM(1,i) = nanstd(AllCausalClusterLength(:,i),0,1)/sqrt(sum(~isnan(AllCausalClusterLength(:,i))));
        AllCausalClusterMaxSEM(1,i) = nanstd(AllCausalClusterMax(:,i),0,1)/sqrt(sum(~isnan(AllCausalClusterMax(:,i))));
        
        FarCueClusterLengthSEM(1,i) = nanstd(FarCueClusterLength(:,i),0,1)/sqrt(sum(~isnan(FarCueClusterLength(:,i))));
        FarMovClusterLengthSEM(1,i) = nanstd(FarMovClusterLength(:,i),0,1)/sqrt(sum(~isnan(FarMovClusterLength(:,i))));
        FarMixClusterLengthSEM(1,i) = nanstd(FarMixClusterLength(:,i),0,1)/sqrt(sum(~isnan(FarMixClusterLength(:,i))));
        FarPreSucClusterLengthSEM(1,i) = nanstd(FarPreSucClusterLength(:,i),0,1)/sqrt(sum(~isnan(FarPreSucClusterLength(:,i))));
        FarSucClusterLengthSEM(1,i) = nanstd(FarSucClusterLength(:,i),0,1)/sqrt(sum(~isnan(FarSucClusterLength(:,i))));
        FarMovDuringCueClusterLengthSEM(1,i) = nanstd(FarMovDuringCueClusterLength(:,i),0,1)/sqrt(sum(~isnan(FarMovDuringCueClusterLength(:,i))));
        FarRewClusterLengthSEM(1,i) = nanstd(FarRewClusterLength(:,i),0,1)/sqrt(sum(~isnan(FarRewClusterLength(:,i))));
        AllFarClusterLengthSEM(1,i) = nanstd(AllFarClusterLength(:,i),0,1)/sqrt(sum(~isnan(AllFarClusterLength(:,i))));

        
        DistanceBetweenCueSpinesSEM(1,i) = nanstd(DistanceBetweenCueSpines(:,i),0,1)/sqrt(sum(~isnan(DistanceBetweenCueSpines(:,i))));
        DistanceBetweenMovementSpinesSEM(1,i) = nanstd(DistanceBetweenMovementSpines(:,i),0,1)/sqrt(sum(~isnan(DistanceBetweenMovementSpines(:,i))));
        DistanceBetweenPreSuccessSpinesSEM(1,i) = nanstd(DistanceBetweenPreSuccessSpines(:,i),0,1)/sqrt(sum(~isnan(DistanceBetweenPreSuccessSpines(:,i))));
        DistanceBetweenSuccessSpinesSEM(1,i) = nanstd(DistanceBetweenSuccessSpines(:,i),0,1)/sqrt(sum(~isnan(DistanceBetweenSuccessSpines(:,i))));
        DistanceBetweenMovementDuringCueSpinesSEM(1,i) = nanstd(DistanceBetweenMovementDuringCueSpines(:,i),0,1)/sqrt(sum(~isnan(DistanceBetweenMovementDuringCueSpines(:,i))));
        DistanceBetweenRewardSpinesSEM(1,i) = nanstd(DistanceBetweenRewardSpines(:,i),0,1)/sqrt(sum(~isnan(DistanceBetweenRewardSpines(:,i))));
        
        ClusteredSpines_CorrwithDendSEM(1,i) = nanstd(ClusteredSpines_CorrwithDend(:,i),0,1)/sqrt(sum(~isnan(ClusteredSpines_CorrwithDend(:,i))));
        FilteredClusteredSpines_CorrwithDendSEM(1,i) = nanstd(FilteredClusteredSpines_CorrwithDend(:,i),0,1)/sqrt(sum(~isnan(FilteredClusteredSpines_CorrwithDend(:,i))));
        NonClusteredSpines_CorrwithDendSEM(1,i) = nanstd(NonClusteredSpines_CorrwithDend(:,i),0,1)/sqrt(sum(~isnan(NonClusteredSpines_CorrwithDend(:,i))));
        CausalClusteredSpines_CorrwithDendSEM(1,i) = nanstd(CausalClusteredSpines_CorrwithDend(:,i),0,1)/sqrt(sum(~isnan(CausalClusteredSpines_CorrwithDend(:,i))));
        FilteredCausalClusteredSpines_CorrwithDendSEM(1,i) = nanstd(FilteredCausalClusteredSpines_CorrwithDend(:,i),0,1)/sqrt(sum(~isnan(FilteredCausalClusteredSpines_CorrwithDend(:,i))));
        NonCausalClusteredSpines_CorrwithDendSEM(1,i) = nanstd(NonCausalClusteredSpines_CorrwithDend(:,i),0,1)/sqrt(sum(~isnan(NonCausalClusteredSpines_CorrwithDend(:,i))));
        CueRelClusteredSpines_CorrwithDendSEM(1,i) = nanstd(CueRelClusteredSpines_CorrwithDend(:,i),0,1)/sqrt(sum(~isnan(CueRelClusteredSpines_CorrwithDend(:,i))));
        MovRelClusteredSpines_CorrwithDendSEM(1,i) = nanstd(MovRelClusteredSpines_CorrwithDend(:,i),0,1)/sqrt(sum(~isnan(MovRelClusteredSpines_CorrwithDend(:,i))));
        SucRelClusteredSpines_CorrwithDendSEM(1,i) = nanstd(SucRelClusteredSpines_CorrwithDend(:,i),0,1)/sqrt(sum(~isnan(SucRelClusteredSpines_CorrwithDend(:,i))));
        RewRelClusteredSpines_CorrwithDendSEM(1,i) = nanstd(RewRelClusteredSpines_CorrwithDend(:,i),0,1)/sqrt(sum(~isnan(RewRelClusteredSpines_CorrwithDend(:,i))));
        CueRelCausalClusteredSpines_CorrwithDendSEM(1,i) = nanstd(CueRelCausalClusteredSpines_CorrwithDend(:,i),0,1)/sqrt(sum(~isnan(CueRelCausalClusteredSpines_CorrwithDend(:,i))));
        MovRelCausalClusteredSpines_CorrwithDendSEM(1,i) = nanstd(MovRelCausalClusteredSpines_CorrwithDend(:,i),0,1)/sqrt(sum(~isnan(MovRelCausalClusteredSpines_CorrwithDend(:,i))));
        SucRelCausalClusteredSpines_CorrwithDendSEM(1,i) = nanstd(SucRelCausalClusteredSpines_CorrwithDend(:,i),0,1)/sqrt(sum(~isnan(SucRelCausalClusteredSpines_CorrwithDend(:,i))));
        RewRelCausalClusteredSpines_CorrwithDendSEM(1,i) = nanstd(RewRelCausalClusteredSpines_CorrwithDend(:,i),0,1)/sqrt(sum(~isnan(RewRelCausalClusteredSpines_CorrwithDend(:,i))));
        
        MeanCorrelationBetweenMovementSpinesSEM(1,i) = nanstd(MeanCorrelationBetweenMovementSpines(:,i), 0,1)/sqrt(sum(~isnan(MeanCorrelationBetweenMovementSpines(:,i))));
        MeanCorrelationBetweenFarMovementSpinesSEM(1,i) = nanstd(MeanCorrelationBetweenFarMovementSpines(:,i),0,1)/sqrt(sum(~isnan(MeanCorrelationBetweenFarMovementSpines(:,i))));
        
%         FractionofClusterThatsMovementRelatedSEM(1,i) = nanstd(FractionofClusterThatsMovementRelated(:,i),0,1)/sqrt(sum(~isnan(FractionofClusterThatsMovementRelated(:,i))));
%         FractionofCausalClusterThatsMovementRelatedSEM(1,i) = nanstd(FractionofCausalClusterThatsMovementRelated(:,i),0,1)/sqrt(sum(~isnan(FractionofCausalClusterThatsMovementRelated(:,i))));

        PercentofCueRelatedDendritesSEM(1,i) = nanstd(PercentofCueRelatedDendrites(:,i),0,1)/sqrt(sum(~isnan(PercentofCueRelatedDendrites(:,i))));
        PercentofMovementRelatedDendritesSEM(1,i) = nanstd(PercentofMovementRelatedDendrites(:,i),0,1)/sqrt(sum(~isnan(PercentofMovementRelatedDendrites(:,i))));
        PercentofPreSuccessRelatedDendritesSEM(1,i) = nanstd(PercentofPreSuccessRelatedDendrites(:,i),0,1)/sqrt(sum(~isnan(PercentofPreSuccessRelatedDendrites(:,i))));
        PercentofSuccessRelatedDendritesSEM(1,i) = nanstd(PercentofSuccessRelatedDendrites(:,i),0,1)/sqrt(sum(~isnan(PercentofSuccessRelatedDendrites(:,i))));        
        PercentofMovementDuringCueRelatedDendritesSEM(1,i) = nanstd(PercentofMovementDuringCueRelatedDendrites(:,i),0,1)/sqrt(sum(~isnan(PercentofMovementDuringCueRelatedDendrites(:,i))));
        PercentofRewardRelatedDendritesSEM(1,i) = nanstd(PercentofRewardRelatedDendrites(:,i),0,1)/sqrt(sum(~isnan(PercentofRewardRelatedDendrites(:,i))));

        try
            SpatialDegreeSEM(1,i) = nanstd(SpatialDegree{i},0,1)/sqrt(sum(~isnan(SpatialDegree{i})));
            TemporalDegreeSEM(1,i) = nanstd(TemporalDegree{i},0,1)/sqrt(sum(~isnan(TemporalDegree{i})));
        catch
            SpatialDegreeSEM(1,i) = nan;
            TemporalDegreeSEM(1,i) = nan;
        end
        
        SpatialDegreeofCueSpinesSEM(1,i) = nanstd(SpatialDegreeofCueSpines(:,i),0,1)/sqrt(sum(~isnan(SpatialDegreeofCueSpines(:,i))));
        TemporalDegreeofCueSpinesSEM(1,i) = nanstd(TemporalDegreeofCueSpines(:,i),0,1)/sqrt(sum(~isnan(TemporalDegreeofCueSpines(:,i))));
        SpatioTemporalDegreeofCueSpinesSEM(1,i) = nanstd(SpatioTemporalDegreeofCueSpines(:,i),0,1)/sqrt(sum(~isnan(SpatioTemporalDegreeofCueSpines(:,i))));
        SpatialDegreeofMovementSpinesSEM(1,i) = nanstd(SpatialDegreeofMovementSpines(:,i),0,1)/sqrt(sum(~isnan(SpatialDegreeofMovementSpines(:,i))));
        TemporalDegreeofMovementSpinesSEM(1,i) = nanstd(TemporalDegreeofMovementSpines(:,i),0,1)/sqrt(sum(~isnan(TemporalDegreeofMovementSpines(:,i))));
        SpatioTemporalDegreeofMovementSpinesSEM(1,i) = nanstd(SpatioTemporalDegreeofMovementSpines(:,i),0,1)/sqrt(sum(~isnan(SpatioTemporalDegreeofMovementSpines(:,i))));
        SpatialDegreeofMovementDuringCueSpinesSEM(1,i) = nanstd(SpatialDegreeofMovementDuringCueSpines(:,i),0,1)/sqrt(sum(~isnan(SpatialDegreeofMovementDuringCueSpines(:,i))));
        TemporalDegreeofMovementDuringCueSpinesSEM(1,i) = nanstd(TemporalDegreeofMovementDuringCueSpines(:,i),0,1)/sqrt(sum(~isnan(TemporalDegreeofMovementDuringCueSpines(:,i))));
        SpatioTemporalDegreeofMovementDuringCueSpinesSEM(1,i) = nanstd(SpatioTemporalDegreeofMovementDuringCueSpines(:,i),0,1)/sqrt(sum(~isnan(SpatioTemporalDegreeofMovementDuringCueSpines(:,i))));  
        SpatialDegreeofPreSuccessSpinesSEM(1,i) = nanstd(SpatialDegreeofPreSuccessSpines(:,i),0,1)/sqrt(sum(~isnan(SpatialDegreeofPreSuccessSpines(:,i))));
        TemporalDegreeofPreSuccessSpinesSEM(1,i) = nanstd(TemporalDegreeofPreSuccessSpines(:,i),0,1)/sqrt(sum(~isnan(TemporalDegreeofPreSuccessSpines(:,i))));
        SpatioTemporalDegreeofPreSuccessSpinesSEM(1,i) = nanstd(SpatioTemporalDegreeofPreSuccessSpines(:,i),0,1)/sqrt(sum(~isnan(SpatioTemporalDegreeofPreSuccessSpines(:,i))));
        SpatialDegreeofSuccessSpinesSEM(1,i) = nanstd(SpatialDegreeofSuccessSpines(:,i),0,1)/sqrt(sum(~isnan(SpatialDegreeofSuccessSpines(:,i))));
        TemporalDegreeofSuccessSpinesSEM(1,i) = nanstd(TemporalDegreeofSuccessSpines(:,i),0,1)/sqrt(sum(~isnan(TemporalDegreeofSuccessSpines(:,i))));
        SpatioTemporalDegreeofSuccessSpinesSEM(1,i) = nanstd(SpatioTemporalDegreeofSuccessSpines(:,i),0,1)/sqrt(sum(~isnan(SpatioTemporalDegreeofSuccessSpines(:,i))));
        SpatialDegreeofRewardSpinesSEM(1,i) = nanstd(SpatialDegreeofRewardSpines(:,i),0,1)/sqrt(sum(~isnan(SpatialDegreeofRewardSpines(:,i))));
        TemporalDegreeofRewardSpinesSEM(1,i) = nanstd(TemporalDegreeofRewardSpines(:,i),0,1)/sqrt(sum(~isnan(TemporalDegreeofRewardSpines(:,i))));
        SpatioTemporalDegreeofRewardSpinesSEM(1,i) = nanstd(SpatioTemporalDegreeofRewardSpines(:,i),0,1)/sqrt(sum(~isnan(SpatioTemporalDegreeofRewardSpines(:,i))));
        
    end
    
    SpatioTemporalOverlapSEM = cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))), SpatioTemporalOverlap);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Plots %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    red = [0.85 0.11 0.14];     black = [0.1 0.1 0.15];
    dred = [0.6 0 0];          dorange = [0.8 0.3 0.03];
    bgreen = [0 0.6 0.7];
    colorj = {red,lblue,green,lgreen,gray,brown,yellow,blue,purple,lpurple,magenta,pink,orange,brown,lbrown};
    rnbo = {dred, red, dorange, orange, yellow, lgreen, green, bgreen, blue, lblue, purple, magenta, lpurple, pink}; 

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Figure 1: Spine-Movement correlation using different separation
    %%%           criteria
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    scrsz = get(0, 'ScreenSize');
    figure('Position', scrsz); 
    subplot(2,5,1)
        plot(MeanAllClustersCorrwithCue, 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanNonClusteredCorrwithCue, 'Color', dred, 'Linewidth', 2);
        plot(MeanCueRelatedClustersCorrwithCue, 'Color', blue, 'Linewidth', 2)
        plot(MeanCueRelatedNonClusteredCorrwithCue, 'Color', gray, 'Linewidth', 2)
        title('Synaptic Events with Cue', 'Fontsize', 14)
        xlim([0 15])
        r_errorbar(1:14, MeanAllClustersCorrwithCue, AllClustersCorrwithCueSEM, black)
        r_errorbar(1:14, MeanNonClusteredCorrwithCue, NonClusteredCorrwithCueSEM, dred)
        r_errorbar(1:14, MeanCueRelatedClustersCorrwithCue, CueRelatedClustersCorrwithCueSEM, blue)
        r_errorbar(1:14, MeanCueRelatedNonClusteredCorrwithCue, CueRelatedNonClusteredCorrwithCueSEM, gray)
        xlabel('Session', 'Fontsize', 14)
        set(gca, 'XTick', 0:15); set(gca, 'XTickLabel', 0:15)
        ylabel('Correlation', 'Fontsize', 14)
    subplot(2,5,2)
        plot(MeanAllClustersCorrwithMDC, 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanNonClusteredCorrwithMDC, 'Color', dred, 'Linewidth', 2);
        plot(MeanMDCRelatedClustersCorrwithMDC, 'Color', blue, 'Linewidth', 2);
        plot(MeanMDCRelatedNonClusteredCorrwithMDC, 'Color', gray, 'Linewidth', 2);
        title('Synaptic Events with Cue', 'Fontsize', 14)
        xlim([0 15])
        r_errorbar(1:14, MeanAllClustersCorrwithMDC, AllClustersCorrwithMDCSEM, black)
        r_errorbar(1:14, MeanNonClusteredCorrwithMDC, NonClusteredCorrwithMDCSEM, dred)
        r_errorbar(1:14, MeanMDCRelatedClustersCorrwithMDC, MDCRelatedClustersCorrwithMDCSEM, blue)
        r_errorbar(1:14, MeanMDCRelatedNonClusteredCorrwithMDC, MDCRelatedNonClusteredCorrwithMDCSEM, gray)
        xlabel('Session', 'Fontsize', 14)
        set(gca, 'XTick', 0:15); set(gca, 'XTickLabel', 0:15)
        ylabel('Correlation', 'Fontsize', 14)
    subplot(2,5,3)
        plot(MeanAllClustersCorrwithMovement, 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanNonClusteredCorrwithMovement, 'Color', dred, 'Linewidth', 2);
        plot(MeanMovementRelatedClustersCorrwithMovement, 'Color', blue, 'Linewidth', 2);
        plot(MeanMovementRelatedNonClusteredCorrwithMovement, 'Color', gray, 'Linewidth', 2);
        title('Synaptic Events with Movement', 'Fontsize', 14)
        xlim([0 15])
        r_errorbar(1:14, MeanAllClustersCorrwithMovement, AllClustersCorrwithMovementSEM, black)
        r_errorbar(1:14, MeanNonClusteredCorrwithMovement, NonClusteredCorrwithMovementSEM, dred)
        r_errorbar(1:14, MeanMovementRelatedClustersCorrwithMovement, MovementRelatedClustersCorrwithMovementSEM, blue)
        r_errorbar(1:14, MeanMovementRelatedNonClusteredCorrwithMovement, MovementRelatedNonClusteredCorrwithMovementSEM, gray)
        xlabel('Session', 'Fontsize', 14)
        set(gca, 'XTick', 0:15); set(gca, 'XTickLabel', 0:15)
        ylabel('Correlation', 'Fontsize', 14)
    subplot(2,5,4)
        plot(MeanAllClustersCorrwithSuccess, 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanNonClusteredCorrwithSuccess, 'Color', dred', 'Linewidth', 2);
        plot(MeanSuccessRelatedClustersCorrwithSuccess, 'Color', blue, 'Linewidth', 2);
        plot(MeanSuccessRelatedNonClusteredCorrwithSuccess, 'Color', gray, 'Linewidth', 2);
        title([{'Synaptic Events with'}, {'Successful Movements'}], 'Fontsize', 14)
        xlim([0 15])
        r_errorbar(1:14, MeanAllClustersCorrwithSuccess, AllClustersCorrwithSuccessSEM, black)
        r_errorbar(1:14, MeanNonClusteredCorrwithSuccess, NonClusteredCorrwithSuccessSEM, dred)
        r_errorbar(1:14, MeanSuccessRelatedClustersCorrwithSuccess, SuccessRelatedClustersCorrwithSuccessSEM, blue)
        r_errorbar(1:14, MeanSuccessRelatedNonClusteredCorrwithSuccess, SuccessRelatedNonClusteredCorrwithSuccessSEM, gray);
        xlabel('Session', 'Fontsize', 14)
        set(gca, 'XTick', 0:15); set(gca, 'XTickLabel', 0:15)
        ylabel('Correlation', 'Fontsize', 14)
    subplot(2,5,5)
        plot(MeanAllClustCorrwithReward, 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanNonClusteredCorrwithReward, 'Color', dred', 'Linewidth', 2);
        plot(MeanRewardRelatedClustersCorrwithReward, 'Color', blue, 'Linewidth', 2);
        plot(MeanRewardRelatedNonClusteredCorrwithReward, 'Color', gray, 'Linewidth', 2)
        title('Synaptic Events with Reward', 'Fontsize', 14)
        xlim([0 15])
        legend({'Clustered Spines', 'Nonclustered spines', '(Function)-related clustered spines', '(Function)-related nonclustered spines'})
        r_errorbar(1:14, MeanAllClustCorrwithReward, AllClustCorrwithRewardSEM, black)
        r_errorbar(1:14, MeanNonClusteredCorrwithReward, NonClusteredCorrwithRewardSEM, dred)
        r_errorbar(1:14, MeanRewardRelatedClustersCorrwithReward, RewardRelatedClustersCorrwithRewardSEM, blue)
        r_errorbar(1:14, MeanRewardRelatedNonClusteredCorrwithReward, RewardRelatedNonClusteredCorrwithRewardSEM, gray)
        xlabel('Session', 'Fontsize', 14)
        set(gca, 'XTick', 0:15); set(gca, 'XTickLabel', 0:15)
        ylabel('Correlation', 'Fontsize', 14)
    subplot(2,5,6)
        plot(MeanAllCausalClustersCorrwithCue, 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanCausalNonClusteredCorrwithCue, 'Color', dred, 'Linewidth', 2);
        plot(MeanCausalCueRelatedClustersCorrwithCue, 'Color', blue, 'Linewidth', 2);
        plot(MeanCausalCueRelatedNonClusteredCorrwithCue, 'Color', gray, 'Linewidth', 2)
        xlim([0 15])
        title('Causal Events with Cue', 'Fontsize', 14)
        r_errorbar(1:14, MeanAllCausalClustersCorrwithCue, AllCausalClustersCorrwithCueSEM, black)
        r_errorbar(1:14, MeanCausalNonClusteredCorrwithCue, CausalNonClusteredCorrwithCue, dred)
        r_errorbar(1:14, MeanCausalCueRelatedClustersCorrwithCue, CausalCueRelatedClustersCorrwithCueSEM, blue)
        r_errorbar(1:14, MeanCausalCueRelatedNonClusteredCorrwithCue, CausalCueRelatedNonClusteredCorrwithCueSEM, gray)
        xlabel('Session', 'Fontsize', 14)
        set(gca, 'XTick', 0:15); set(gca, 'XTickLabel', 0:15)
        ylabel('Correlation', 'Fontsize', 14)
    subplot(2,5,7)
        plot(MeanAllCausalClustersCorrwithMDC, 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanCausalNonClusteredCorrwithMDC, 'Color', dred, 'Linewidth', 2);
        plot(MeanCausalMDCRelatedClustersCorrwithMDC, 'Color', blue, 'Linewidth', 2);
        plot(MeanCausalMDCRelatedNonClusteredCorrwithMDC, 'Color', gray, 'Linewidth', 2);
        title('Synaptic Events with Cue', 'Fontsize', 14)
        xlim([0 15])
        r_errorbar(1:14, MeanAllCausalClustersCorrwithMDC, AllCausalClustersCorrwithMDCSEM, black)
        r_errorbar(1:14, MeanCausalNonClusteredCorrwithMDC, CausalNonClusteredCorrwithMDCSEM, dred)
        r_errorbar(1:14, MeanCausalMDCRelatedClustersCorrwithMDC, CausalMDCRelatedClustersCorrwithMDCSEM, blue)
        r_errorbar(1:14, MeanCausalMDCRelatedNonClusteredCorrwithMDC, CausalMDCRelatedNonClusteredCorrwithMDCSEM, gray)
        xlabel('Session', 'Fontsize', 14)
        set(gca, 'XTick', 0:15); set(gca, 'XTickLabel', 0:15)
        ylabel('Correlation', 'Fontsize', 14)
    subplot(2,5,8)
        plot(MeanAllCausalClustersCorrwithMovement, 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanCausalNonClusteredCorrwithMovement, 'Color', dred, 'Linewidth', 2);
        plot(MeanCausalMovementRelatedClustersCorrwithMovement, 'Color', blue, 'Linewidth', 2);
        plot(MeanCausalMovementRelatedNonClusteredCorrwithMovement, 'Color', gray, 'Linewidth', 2);
        xlim([0 15])
        title('Causal Events with Movement', 'Fontsize', 14)
        r_errorbar(1:14, MeanAllCausalClustersCorrwithMovement, AllCausalClustersCorrwithMovementSEM, black)
        r_errorbar(1:14, MeanCausalNonClusteredCorrwithMovement, CausalNonClusteredCorrwithMovementSEM, dred)
        r_errorbar(1:14, MeanCausalMovementRelatedClustersCorrwithMovement, CausalMovementRelatedClustersCorrwithMovementSEM, blue)
        r_errorbar(1:14, MeanCausalMovementRelatedNonClusteredCorrwithMovement, CausalMovementRelatedNonClusteredCorrwithMovementSEM, gray)
        xlabel('Session', 'Fontsize', 14)
        set(gca, 'XTick', 0:15); set(gca, 'XTickLabel', 0:15)
        ylabel('Correlation', 'Fontsize', 14)
    subplot(2,5,9)
        plot(MeanAllCausalClustersCorrwithSuccess, 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanCausalNonClusteredCorrwithSuccess, 'Color', dred, 'Linewidth', 2);
        plot(MeanCausalSuccessRelatedClustersCorrwithSuccess, 'Color', blue, 'Linewidth', 2);
        plot(MeanCausalSuccessRelatedNonClusteredCorrwithSuccess, 'Color', gray, 'Linewidth', 2);
        xlim([0 15])
        title([{'Causal Events with'}, {'Successful Movements'}], 'Fontsize', 14)
        r_errorbar(1:14, MeanAllCausalClustersCorrwithSuccess, AllCausalClustersCorrwithSuccessSEM, black)
        r_errorbar(1:14, MeanCausalNonClusteredCorrwithSuccess, CausalNonClusteredCorrwithSuccessSEM, dred)
        r_errorbar(1:14, MeanCausalSuccessRelatedClustersCorrwithSuccess, CausalSuccessRelatedClustersCorrwithSuccessSEM, blue)
        r_errorbar(1:14, MeanCausalSuccessRelatedNonClusteredCorrwithSuccess, CausalSuccessRelatedNonClusteredCorrwithSuccessSEM, gray);
        xlabel('Session', 'Fontsize', 14)
        set(gca, 'XTick', 0:15); set(gca, 'XTickLabel', 0:15)
        ylabel('Correlation', 'Fontsize', 14)
    subplot(2,5,10)
        plot(MeanAllCausalClustCorrwithReward, 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanCausalNonClusteredCorrwithReward, 'Color', dred, 'Linewidth', 2);
        plot(MeanCausalRewardRelatedClustersCorrwithReward, 'Color', blue, 'Linewidth', 2);
        plot(MeanCausalRewardRelatedNonClusteredCorrwithReward, 'Color', gray, 'Linewidth', 2)
        xlim([0 15])
        title('Causal Events with Reward', 'Fontsize', 14)
        r_errorbar(1:14, MeanAllCausalClustCorrwithReward, AllCausalClustCorrwithRewardSEM, black)
        r_errorbar(1:14, MeanCausalNonClusteredCorrwithReward, CausalNonClusteredCorrwithRewardSEM, dred)
        r_errorbar(1:14, MeanCausalRewardRelatedClustersCorrwithReward, CausalRewardRelatedClustersCorrwithRewardSEM, blue)
        r_errorbar(1:14, MeanCausalRewardRelatedNonClusteredCorrwithReward, CausalRewardRelatedNonClusteredCorrwithRewardSEM, gray)
        xlabel('Session', 'Fontsize', 14)
        set(gca, 'XTick', 0:15); set(gca, 'XTickLabel', 0:15)
        ylabel('Correlation', 'Fontsize', 14)
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Figure 2: Clustered vs. nonclustered frequency, amp, etc. and
    %%%           dendrite information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure('Position', scrsz);
    subplot(2,3,1);
        plot(MeanClusterFreq, 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanCausalClusterFreq, 'Color', bgreen, 'Linewidth', 2); hold on;
        plot(MeanNonClusteredFreq, 'Color', dred, 'Linewidth', 2);
        plot(MeanNonClusteredCausalFreq, 'Color', gray, 'Linewidth', 2);
        ylabel('Event Frequency', 'FontSize', 14)
        xlabel('Session', 'FontSize', 14)
        xlim([0 15])
        r_errorbar(1:14, MeanClusterFreq, ClusterFreqSEM, black)
        r_errorbar(1:14, MeanCausalClusterFreq, CausalClusterFreqSEM, bgreen)
        r_errorbar(1:14, MeanNonClusteredFreq, NonClusteredFreqSEM, dred)
        r_errorbar(1:14, MeanNonClusteredCausalFreq, NonClusteredCausalFreqSEM, bgreen)
    legend({'Clustered', 'Causal', 'Nonclustered', 'NonClust Caus'});
    subplot(2,3,2)
        plot(MeanClusteredSpineAmp, 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanCausalClusteredSpineAmp, 'Color', bgreen, 'Linewidth', 2);
        plot(MeanNonClusteredSpineAmp, 'Color', dred, 'Linewidth', 2)
        plot(MeanCausalNonClusteredSpineAmp, 'Color', gray, 'Linewidth', 2);
        legend({'Clustered', 'Causal clustered', 'Nonclustered', 'Causal nonclustered'})
        ylabel('Event Amp', 'FontSize', 14);
        xlabel('Session', 'FontSize', 14);
        xlim([0 15])
        r_errorbar(1:14, MeanClusteredSpineAmp, ClusteredSpineAmpSEM, black)
        r_errorbar(1:14, MeanCausalClusteredSpineAmp, CausalClusteredSpineAmpSEM, bgreen)
        r_errorbar(1:14, MeanNonClusteredSpineAmp, NonClusteredSpineAmpSEM, dred)
        r_errorbar(1:14, MeanCausalNonClusteredSpineAmp, CausalNonClusteredSpineAmpSEM, gray)
    subplot(2,3,3)
        plot(MeanClustDendFreq, 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanNonClustDendFreq, 'Color', dred, 'Linewidth', 2); 
        ylabel('Event Frequency', 'Fontsize', 14);
        xlabel('Session', 'Fontsize', 14)
        xlim([0 15])
        legend({'Dendrites with Clusters', 'Dendrites w/o Clusters'})
        r_errorbar(1:14, MeanClustDendFreq, ClustDendFreqSEM, black)
        r_errorbar(1:14, MeanNonClustDendFreq, NonClustDendFreqSEM, dred)
    subplot(2,3,4); 
        plot(MeanCueClusterFrequency, 'Color', lgreen, 'Linewidth', 2); hold on;
        plot(MeanMovementClusterFrequency, 'Color', black, 'Linewidth', 2);
        plot(MeanMovementDuringCueClusterFrequency, 'Color', green, 'Linewidth', 2);
        plot(MeanPreSuccessClusterFrequency, 'Color', bgreen, 'Linewidth', 2);
        plot(MeanSuccessClusterFrequency, 'Color', lblue, 'Linewidth', 2);
        plot(MeanRewardClusterFrequency, 'Color', purple, 'Linewidth', 2);
        title([{'Frequency of Functionally relevant'}, {'clustered spines'}])
        ylabel('Event Frequency', 'Fontsize', 14);
        xlabel('Session', 'Fontsize', 14);
        xlim([0 15])
        r_errorbar(1:14, MeanCueClusterFrequency, CueClusterFrequencySEM, lgreen)
        r_errorbar(1:14, MeanMovementClusterFrequency, MovementClusterFrequencySEM, black)
        r_errorbar(1:14, MeanMovementDuringCueClusterFrequency, MovementDuringCueClusterFrequencySEM, green)
        r_errorbar(1:14, MeanPreSuccessClusterFrequency, PreSuccessClusterFrequencySEM, bgreen)
        r_errorbar(1:14, MeanSuccessClusterFrequency, SuccessClusterFrequencySEM, lblue)
        r_errorbar(1:14, MeanRewardClusterFrequency, RewardClusterFrequencySEM, purple)
    subplot(2,3,5);
        plot(MeanClusteredCueSpineAmp, 'Color', lgreen, 'Linewidth', 2); hold on;
        plot(MeanClusteredMoveSpineAmp, 'Color', black, 'Linewidth', 2);
        plot(MeanClusteredMovDuringCueSpineAmp, 'Color', green, 'Linewidth', 2);
        plot(MeanClusteredPreSuccessSpineAmp, 'Color', bgreen, 'Linewidth', 2);
        plot(MeanClusteredSuccessSpineAmp, 'Color', lblue, 'Linewidth', 2);
        plot(MeanClusteredRewardSpineAmp, 'Color', purple, 'Linewidth', 2);
        title([{'Amp. of Functionally relevant'}, {'clustered spines'}])
        ylabel('Event Amp', 'FontSize', 14);
        xlabel('Session', 'FontSize', 14);
        xlim([0 15])
        r_errorbar(1:14, MeanClusteredCueSpineAmp, ClusteredCueSpineAmpSEM, lgreen)
        r_errorbar(1:14, MeanClusteredMoveSpineAmp, ClusteredMoveSpineAmpSEM, black)
        r_errorbar(1:14, MeanClusteredSuccessSpineAmp, ClusteredSuccessSpineAmpSEM, lblue)
        r_errorbar(1:14, MeanClusteredRewardSpineAmp, ClusteredRewardSpineAmpSEM, purple)
    subplot(2,3,6)
        plot(MeanCueClustDendFreq, 'Color', lgreen, 'Linewidth', 2); hold on;
        plot(MeanMovClustDendFreq, 'Color', black, 'Linewidth', 2);
        plot(MeanMovDuringCueClustDendFreq, 'Color', green, 'Linewidth', 2);
        plot(MeanPreSucClustDendFreq, 'Color', bgreen, 'Linewidth', 2);
        plot(MeanSucClustDendFreq, 'Color', lblue, 'Linewidth', 2);
        plot(MeanRewClustDendFreq, 'Color', purple, 'Linewidth', 2);
        plot(MeanNonMovClustDendFreq, 'Color', dred, 'Linewidth', 2);
        title([{'Frequency of dendrites with functionally'}, {'relevant clustered spines'}])
        ylabel('Event Frequency', 'Fontsize', 14);
        xlabel('Session', 'Fontsize', 14)
        xlim([0 15])
        legend({'Dends with CueClusts','Dends with MovClusts', 'Dends with SucClusts', 'Dends with RewClusts', 'Dends w/o MovClusts'})
        r_errorbar(1:14, MeanCueClustDendFreq, CueClustDendFreqSEM,lgreen)
        r_errorbar(1:14, MeanMovClustDendFreq, MovClustDendFreqSEM, black)
        r_errorbar(1:14, MeanMovDuringCueClustDendFreq, MovDuringCueClustDendFreqSEM, green)
        r_errorbar(1:14, MeanPreSucClustDendFreq, PreSucClustDendFreqSEM, bgreen)
        r_errorbar(1:14, MeanSucClustDendFreq, SucClustDendFreqSEM, lblue)
        r_errorbar(1:14, MeanRewClustDendFreq, RewClustDendFreqSEM, purple)
        r_errorbar(1:14, MeanNonMovClustDendFreq, NonMovClustDendFreqSEM, dred)

  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Figure 3: Num of Mvmnt-related spines over time, in different
    %%%           categories
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure('Position', scrsz); 
    sub1 = 3;
    sub2 = 3;
    subplot(sub1,sub2,1)
        plot(MeanNumCueRelSpines,'Color', lgreen, 'Linewidth', 2); hold on;
        plot(MeanNumMovRelSpines,'Color', black, 'Linewidth', 2);
        plot(MeanNumCueORMovRelSpines, 'Color', red, 'Linewidth', 2);
        plot(MeanNumPreSucRelSpines, 'Color', bgreen, 'Linewidth', 2);
        plot(MeanNumSucRelSpines,'Color', lblue, 'Linewidth', 2);
        plot(MeanNumMovDuringCueRelSpines, 'Color', green, 'Linewidth', 2);
        plot(MeanNumRewRelSpines,'Color', purple, 'Linewidth', 2);
        legend({'All Cue Spines', 'All Mvmt Spines', 'Cue OR Mov Spines', 'Pre Success Spines', 'All Success Spines', 'Mov. during Cue Spines', 'All Reward Spines'})
        r_errorbar(1:14, MeanNumCueRelSpines, NumCueRelSpinesSEM, lgreen)
        r_errorbar(1:14, MeanNumMovRelSpines, NumMovRelSpinesSEM, black)
        r_errorbar(1:14, MeanNumCueORMovRelSpines, NumCueORMovRelSpinesSEM, red)
        r_errorbar(1:14, MeanNumPreSucRelSpines, NumPreSucRelSpinesSEM, bgreen)
        r_errorbar(1:14, MeanNumSucRelSpines, NumSucRelSpinesSEM, lblue)
        r_errorbar(1:14, MeanNumMovDuringCueRelSpines, NumMovDuringCueRelSpinesSEM, green)
        r_errorbar(1:14, MeanNumRewRelSpines, NumRewRelSpinesSEM, purple)
        xlim([0 15])
        xlabel('Session', 'FontSize', 14)
        ylabel('Fraction of total spines', 'FontSize', 14)
        title('Classes of Spines', 'Fontsize', 14)
    
    subplot(sub1,sub2,2)
        plot(MeanNumClustSpines, 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanNumCausClustSpines, 'Color', dred, 'Linewidth', 2)
        plot(MeanNumFarClustSpines, 'Color', gray, 'Linewidth', 2)
        legend({'Clustered', 'Causal Clustered', 'Far'})
        r_errorbar(1:14, MeanNumClustSpines, NumClustSpinesSEM, black)
        r_errorbar(1:14, MeanNumCausClustSpines, NumCausClustSpinesSEM, dred)
        r_errorbar(1:14, MeanNumFarClustSpines, NumFarClustSpinesSEM, gray)
        xlabel('Session','FontSize', 14)
        xlim([0 15])
        ylabel('Fraction of total spines','FontSize', 14)
    
    subplot(sub1,sub2,3)
        plot(MeanNumberofClusters, '-', 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanNumberofCausalClusters, '-', 'Color', bgreen, 'Linewidth', 2)
        plot(MeanNumberofSpinesinEachCluster, '-', 'Color', gray, 'Linewidth', 2)
        plot(MeanNumberofSpinesinEachCausalCluster, '-', 'Color', dorange, 'Linewidth', 2)
        legend({'Number of Clusters', 'Number of Causal Clusters', 'Spines in each cluster', 'Spines in each causal cluster'})
        r_errorbar(1:14, MeanNumberofClusters, NumberofClustersSEM, black)
        r_errorbar(1:14, MeanNumberofCausalClusters, NumberofCausalClustersSEM, bgreen)
        r_errorbar(1:14, MeanNumberofSpinesinEachCluster, NumberofSpinesinEachClusterSEM, gray)
        r_errorbar(1:14, MeanNumberofSpinesinEachCausalCluster, NumberofSpinesinEachCausalClusterSEM, dorange)
        xlabel('Session', 'FontSize', 14)
        xlim([0 15])
        ylabel('Raw Number', 'FontSize', 14)
        
    subplot(sub1,sub2,4)
        plot(MeanPercentofCueRelatedDendrites, 'Color', lgreen, 'Linewidth', 2); hold on;
        plot(MeanPercentofMovementRelatedDendrites, 'Color', black, 'Linewidth', 2); 
        plot(MeanPercentofPreSuccessRelatedDendrites, 'Color', bgreen, 'Linewidth', 2);
        plot(MeanPercentofSuccessRelatedDendrites, 'Color', lblue, 'Linewidth', 2); 
        plot(MeanPercentofMovementDuringCueRelatedDendrites, 'Color', green, 'Linewidth', 2);
        plot(MeanPercentofRewardRelatedDendrites, 'Color', purple, 'Linewidth', 2);
        xlim([0 15])
        legend({'Cue Dends', 'Mov Dends', 'PreSuc Dends','Suc Dends', 'Mov During Cue Dends', 'Rew Dends'})
        r_errorbar(1:14, MeanPercentofCueRelatedDendrites, PercentofCueRelatedDendritesSEM, lgreen)
        r_errorbar(1:14, MeanPercentofMovementRelatedDendrites, PercentofMovementRelatedDendritesSEM, black)
        r_errorbar(1:14, MeanPercentofPreSuccessRelatedDendrites, PercentofPreSuccessRelatedDendritesSEM, bgreen)
        r_errorbar(1:14, MeanPercentofSuccessRelatedDendrites, PercentofSuccessRelatedDendritesSEM, lblue)
        r_errorbar(1:14, MeanPercentofMovementDuringCueRelatedDendrites, PercentofMovementDuringCueRelatedDendritesSEM, green)
        r_errorbar(1:14, MeanPercentofRewardRelatedDendrites, PercentofRewardRelatedDendritesSEM, purple)
        xlabel('Session', 'FontSize', 14)
        ylabel('Fraction of dendrites', 'FontSize', 14)
        
    subplot(sub1,sub2,5)
        plot(MeanNumClustCueSpines, 'Color', lgreen, 'Linewidth', 2); hold on;
            plot(MeanNumCausClustCueSpines, '--', 'Color', lgreen, 'Linewidth', 2);
        plot(MeanNumClustMovSpines, 'Color', black, 'Linewidth', 2)
            plot(MeanNumCausClustMovSpines, '--', 'Color', black, 'Linewidth', 2);
        plot(MeanNumClustMixSpines, 'Color', red, 'Linewidth', 2)
        plot(MeanNumClustPreSucSpines, 'Color', bgreen, 'Linewidth', 2)
            plot(MeanNumCausClustPreSucSpines, '--', 'Color', bgreen, 'Linewidth', 2);
        plot(MeanNumClustSucSpines, 'Color', lblue, 'Linewidth', 2)
            plot(MeanNumCausClustSucSpines, '--', 'Color', lblue, 'Linewidth', 2);
        plot(MeanNumClustMovDuringCueSpines, 'Color', green, 'Linewidth', 2)
            plot(MeanNumCausClustMovDuringCueSpines, '--', 'Color', green, 'Linewidth', 2);
        plot(MeanNumClustRewSpines, 'Color', purple, 'Linewidth', 2)
            plot(MeanNumCausClustRewSpines, '--', 'Color', purple, 'Linewidth', 2);
        xlim([0 15])
        legend({'Clust. Cue Spines','Clust Caus Cue', 'Clust. Mov. Spines','Clust Caus Mov', 'Clust. Mixed Spines', 'Clust Pre suc.','Clust Caus Presuc', 'Clust. Suc. Spines','Clust Caus Suc', 'Clust Mov during Cue','Clust Caus MDC', 'Clust Rew. Spines', 'Clust Caus Rew'})
        r_errorbar(1:14, MeanNumClustCueSpines, NumClustCueSpinesSEM, lgreen)
            r_errorbar(1:14, MeanNumCausClustCueSpines, NumCausClustCueSpinesSEM, lgreen)
        r_errorbar(1:14, MeanNumClustMovSpines, NumClustMovSpinesSEM, black)
            r_errorbar(1:14, MeanNumCausClustMovSpines, NumCausClustMovSpinesSEM, black)
        r_errorbar(1:14, MeanNumClustMixSpines, NumClustMixSpinesSEM, red)
        r_errorbar(1:14, MeanNumClustPreSucSpines, NumClustPreSucSpinesSEM, bgreen)
            r_errorbar(1:14, MeanNumCausClustPreSucSpines, NumCausClustPreSucSpinesSEM, bgreen)
        r_errorbar(1:14, MeanNumClustSucSpines, NumClustSucSpinesSEM, lblue)
            r_errorbar(1:14, MeanNumCausClustSucSpines, NumCausClustSucSpinesSEM, lblue)
        r_errorbar(1:14, MeanNumClustMovDuringCueSpines, NumClustMovDuringCueSpinesSEM, green) 
            r_errorbar(1:14, MeanNumCausClustMovDuringCueSpines, NumCausClustMovDuringCueSpinesSEM, green) 
        r_errorbar(1:14, MeanNumClustRewSpines, NumClustRewSpinesSEM, purple)
            r_errorbar(1:14, MeanNumCausClustRewSpines, NumCausClustRewSpinesSEM, purple)
        xlabel('Session', 'FontSize', 14)
        ylabel('Fraction of total spines', 'FontSize', 14)
        title('Clustered function-related spines')
        
    subplot(sub1,sub2,6)
        plot(MeanNumberofMovClusters, '-', 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanNumberofSpinesinEachMovCluster, '-', 'Color', gray, 'Linewidth', 2)
        legend({'Number of Mov Clusters', 'Spines in each mov cluster'})
        r_errorbar(1:14, MeanNumberofMovClusters, NumberofMovClustersSEM, black)
        r_errorbar(1:14, MeanNumberofSpinesinEachMovCluster, NumberofSpinesinEachMovClusterSEM, gray)
        xlabel('Session', 'FontSize', 14)
        xlim([0 15])
        ylabel('Raw Number', 'FontSize', 14)
        
    subplot(sub1,sub2,7)
        plot(MeanNumFarClustCueSpines, 'Color', lgreen, 'Linewidth', 2); hold on;
        plot(MeanNumFarClustMovSpines, 'Color', black, 'Linewidth', 2)
        plot(MeanNumFarClustPreSucSpines, 'Color', bgreen, 'Linewidth', 2)
        plot(MeanNumFarClustSucSpines, 'Color', lblue, 'Linewidth', 2)
        plot(MeanNumFarClustMovDuringCueSpines, 'Color', green, 'Linewidth', 2)
        plot(MeanNumFarClustRewSpines, 'Color', purple, 'Linewidth', 2)
        xlim([0 15])
        legend({'Clust. Cue Spines' 'Clust. Mov. Spines','Clust Pre suc.','Clust. Suc. Spines','Clust Mov during Cue', 'Clust Rew. Spines'})
        r_errorbar(1:14, MeanNumFarClustCueSpines, NumFarClustCueSpinesSEM, lgreen)
        r_errorbar(1:14, MeanNumFarClustMovSpines, NumFarClustMovSpinesSEM, black)
        r_errorbar(1:14, MeanNumFarClustPreSucSpines, NumFarClustPreSucSpinesSEM, bgreen)
        r_errorbar(1:14, MeanNumFarClustSucSpines, NumFarClustSucSpinesSEM, lblue)
        r_errorbar(1:14, MeanNumFarClustMovDuringCueSpines, NumFarClustMovDuringCueSpinesSEM, green) 
        r_errorbar(1:14, MeanNumFarClustRewSpines, NumFarClustRewSpinesSEM, purple)
        xlabel('Session', 'FontSize', 14)
        ylabel('Fraction of total spines', 'FontSize', 14)
        title('Correlated spines on sep. dendrites')

    subplot(sub1,sub2,8)
        plot(MeanFractionofCueSpinesThatAreClustered, 'Color', lgreen, 'Linewidth', 2); hold on;
        plot(MeanFractionofMovementSpinesThatAreClustered, 'Color', black, 'Linewidth', 2);
        plot(MeanFractionofPreSuccessSpinesThatAreClustered, 'Color', bgreen, 'Linewidth', 2);
        plot(MeanFractionofSuccessSpinesThatAreClustered, 'Color', lblue, 'Linewidth', 2)
        plot(MeanFractionofMovementDuringCueSpinesThatAreClustered, 'Color', green, 'Linewidth', 2)
        plot(MeanFractionofRewardSpinesThatAreClustered, 'Color', purple, 'Linewidth', 2)
        xlim([0 15])
        xlabel('Session', 'Fontsize', 14)
        ylabel('Fraction of Function-related Spines', 'Fontsize', 14)
        title([{'Fraction of (function) spines'},{'that are clustered'}], 'Fontsize', 14)
        legend({'Cue related','Movement related', 'Pre Success', 'Success related', 'MovDuringCue', 'Reward related'})
        r_errorbar(1:14, MeanFractionofCueSpinesThatAreClustered, FractionofCueSpinesThatAreClusteredSEM, lgreen)
        r_errorbar(1:14, MeanFractionofMovementSpinesThatAreClustered, FractionofMovementSpinesThatAreClusteredSEM, black)
        r_errorbar(1:14, MeanFractionofPreSuccessSpinesThatAreClustered, FractionofPreSuccessSpinesThatAreClusteredSEM, bgreen)
        r_errorbar(1:14, MeanFractionofSuccessSpinesThatAreClustered, FractionofSuccessSpinesThatAreClusteredSEM, lblue)
        r_errorbar(1:14, MeanFractionofMovementDuringCueSpinesThatAreClustered, FractionofMovementDuringCueSpinesThatAreClusteredSEM, green)
        r_errorbar(1:14, MeanFractionofRewardSpinesThatAreClustered, FractionofRewardSpinesThatAreClusteredSEM, purple)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Figure 4: Spatial extent of clusters
    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure('Position', scrsz);
    subplot(2,3,1)
        plot(MeanAllClusterLength, 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanAllCausalClusterLength, 'Color', bgreen, 'Linewidth', 2)
        legend({'Clustered spines', 'Caus. clust.'})
        r_errorbar(1:14, MeanAllClusterLength, AllClusterLengthSEM, black)
        r_errorbar(1:14, MeanAllCausalClusterLength, AllCausalClusterLengthSEM, bgreen);
        xlabel('Session', 'FontSize', 14)
        xlim([0 15])
        ylabel('Ave. Dist. between spines', 'FontSize', 14)
        title('Mean Spatial Extent of Clusters', 'Fontsize', 14)
    subplot(2,3,2)
        plot(MeanAllClusterMax, 'Color', black, 'Linewidth', 2); hold on;
        plot(MeanAllCausalClusterMax, 'Color', bgreen, 'Linewidth', 2)
        legend({'Clustered spines', 'Caus. clust.'})
        r_errorbar(1:14, MeanAllClusterMax, AllClusterMaxSEM, black)
        r_errorbar(1:14, MeanAllCausalClusterMax, AllCausalClusterMaxSEM, bgreen);
        xlabel('Session', 'FontSize', 14)
        xlim([0 15])
        ylabel('Max Dist. between spines', 'FontSize', 14)
        title('MAX Spatial Extent of Clusters', 'Fontsize', 14)
    subplot(2,3,3)
        plot(MeanDistanceBetweenCueSpines, 'Color', lgreen, 'Linewidth', 2); hold on;
        plot(MeanDistanceBetweenMovementSpines, 'Color', black, 'Linewidth', 2)
        plot(MeanDistanceBetweenSuccessSpines, 'Color', lblue, 'Linewidth', 2)
        plot(MeanDistanceBetweenRewardSpines, 'Color', purple, 'Linewidth', 2);
        xlabel('Session', 'FontSize', 14)
        ylabel('Distance (um)', 'FontSize', 14)
        ylim([0 30])
        xlim([0 15])
        legend({'Cue Spines', 'Movement Spines', 'Success Spines', 'Reward Spines'})
        r_errorbar(1:14, MeanDistanceBetweenCueSpines, DistanceBetweenCueSpinesSEM, lgreen)
        r_errorbar(1:14, MeanDistanceBetweenMovementSpines, DistanceBetweenMovementSpinesSEM, black)
        r_errorbar(1:14, MeanDistanceBetweenSuccessSpines, DistanceBetweenSuccessSpinesSEM, lblue)
        r_errorbar(1:14, MeanDistanceBetweenRewardSpines, DistanceBetweenRewardSpinesSEM, purple)

    subplot(2,3,4)
        plot(MeanCueClusterLength, 'Color', lgreen, 'LineWidth', 2); hold on;
        plot(MeanMovClusterLength, 'Color', black, 'LineWidth', 2)
        plot(MeanMixClusterLength, 'Color', red, 'Linewidth', 2)
        plot(MeanPreSucClusterLength, 'Color', bgreen, 'Linewidth', 2)
        plot(MeanSucClusterLength, 'Color', lblue, 'Linewidth', 2)
        plot(MeanMovDuringCueClusterLength, 'Color', green, 'Linewidth', 2)
        plot(MeanRewClusterLength, 'Color', purple, 'Linewidth', 2)
        legend({'Cue Clusters', 'Mov clusters', 'Mix Clusters', 'PreSuc', 'Suc Clusters', 'MovDuringCue', 'Rew Clusters'});
        r_errorbar(1:14, MeanCueClusterLength, CueClusterLengthSEM, lgreen)
        r_errorbar(1:14, MeanMovClusterLength, MovClusterLengthSEM, black)
        r_errorbar(1:14, MeanMixClusterLength, MixClusterLengthSEM, red)
        r_errorbar(1:14, MeanPreSucClusterLength, PreSucClusterLengthSEM, bgreen)
        r_errorbar(1:14, MeanSucClusterLength, SucClusterLengthSEM, lblue)
        r_errorbar(1:14, MeanMovDuringCueClusterLength, MovDuringCueClusterLengthSEM, green)
        r_errorbar(1:14, MeanRewClusterLength, RewClusterLengthSEM, purple)
        xlabel('Session', 'FontSize', 14)
        ylabel('Mean spatial extent of clusters', 'FontSize', 14)
        ylim([0 30])
        xlim([0 15])
    subplot(2,3,5)
        plot(MeanCueClusterMax, 'Color', lgreen, 'LineWidth', 2); hold on;
        plot(MeanMovClusterMax, 'Color', black, 'LineWidth', 2)
        plot(MeanMixClusterMax, 'Color', red, 'Linewidth', 2)
        plot(MeanSucClusterMax, 'Color', lblue, 'Linewidth', 2)
        plot(MeanRewClusterMax, 'Color', purple, 'Linewidth', 2)
        legend({'Cue Clusters', 'Mov clusters', 'Mix Clusters', 'Suc Clusters', 'Rew Clusters'});
        r_errorbar(1:14, MeanCueClusterMax, CueClusterMaxSEM, lgreen)
        r_errorbar(1:14, MeanMovClusterMax, MovClusterMaxSEM, black)
        r_errorbar(1:14, MeanMixClusterMax, MixClusterMaxSEM, red)
        r_errorbar(1:14, MeanSucClusterMax, SucClusterMaxSEM, lblue)
        r_errorbar(1:14, MeanRewClusterMax, RewClusterMaxSEM, purple)
        xlabel('Session', 'FontSize', 14)
        ylabel('Max spatial extent of clusters', 'FontSize', 14)
        ylim([0 30])
        xlim([0 15])
        
    subplot(2,3,6)
        plot(MeanFarCueClusterLength, 'Color', lgreen, 'Linewidth', 2); hold on; 
        plot(MeanFarMovClusterLength, 'Color', black, 'Linewidth', 2);
        plot(MeanFarMixClusterLength, 'Color', red, 'Linewidth', 2);
        plot(MeanFarSucClusterLength, 'Color', lblue, 'Linewidth', 2);
        plot(MeanFarRewClusterLength, 'Color', purple, 'Linewidth', 2);
        plot(MeanAllFarClusterLength, 'Color', gray, 'Linewidth', 2);
        legend({'Far Cue', 'Far Mov', 'Far Mix', 'Far Suc', 'Far Rew', 'Far All'})
        r_errorbar(1:14, MeanFarCueClusterLength, FarCueClusterLengthSEM, lgreen)
        r_errorbar(1:14, MeanFarMovClusterLength, FarMovClusterLengthSEM, black)
        r_errorbar(1:14, MeanFarMixClusterLength, FarMixClusterLengthSEM, red)
        r_errorbar(1:14, MeanFarSucClusterLength, FarSucClusterLengthSEM, lblue)
        r_errorbar(1:14, MeanFarRewClusterLength, FarRewClusterLengthSEM, purple)
        r_errorbar(1:14, MeanAllFarClusterLength, AllFarClusterLengthSEM, gray)
        
    %%%
    %%% Figure 5: Correlation with Dendrite
    %%%
    
    figure('Position', scrsz);
    
    subplot(3,2,1)
        plot(MeanClusteredSpines_CorrwithDend, 'Color',black, 'LineWidth', 2); hold on;
        plot(MeanFilteredClusteredSpines_CorrwithDend, '--','Color', orange, 'Linewidth', 2);
        plot(MeanNonClusteredSpines_CorrwithDend, '--','Color',dred, 'LineWidth', 2);

        legend({'Clust','Filt Clust','Non Clust'})
    
        r_errorbar(1:14, MeanClusteredSpines_CorrwithDend, ClusteredSpines_CorrwithDendSEM, black)
        r_errorbar(1:14, MeanFilteredClusteredSpines_CorrwithDend, FilteredClusteredSpines_CorrwithDendSEM, orange)
        r_errorbar(1:14, MeanNonClusteredSpines_CorrwithDend, NonClusteredSpines_CorrwithDendSEM, dred)
    
        ylabel('Correlation with dendrite', 'Fontsize', 14)
        xlabel('Session', 'Fontsize', 14)
        title('Clustered Spines', 'Fontsize', 14)
        xlim([0 15])

    subplot(3,2,2)
        plot(MeanCausalClusteredSpines_CorrwithDend, 'Color',black, 'Linewidth', 2); hold on;
        plot(MeanFilteredCausalClusteredSpines_CorrwithDend, 'Color', orange, 'Linewidth', 2);
        plot(MeanNonCausalClusteredSpines_CorrwithDend, 'Color',dred, 'Linewidth', 2)
    
        legend({'Caus Clust','Filt Caus Clust', 'Caus Non Clust'})

        r_errorbar(1:14, MeanCausalClusteredSpines_CorrwithDend, CausalClusteredSpines_CorrwithDendSEM, 'k')
        r_errorbar(1:14, MeanFilteredCausalClusteredSpines_CorrwithDend, FilteredCausalClusteredSpines_CorrwithDendSEM, orange)
        r_errorbar(1:14, MeanNonCausalClusteredSpines_CorrwithDend, NonCausalClusteredSpines_CorrwithDendSEM, 'r')

        ylabel('Correlation with dendrite', 'Fontsize', 14)
        xlabel('Session', 'Fontsize', 14)
        title('Causal Clustered Spines', 'Fontsize', 14)
        xlim([0 15])
        
    subplot(3,2,3)
        plot(MeanCorrelationofClusters, 'Color', black, 'Linewidth', 2); hold on;
        r_errorbar(1:14, MeanCorrelationofClusters, CorrelationofClustersSEM, black)
        
        ylabel('Correlation', 'Fontsize', 14)
        xlabel('Session', 'Fontsize', 14)
        title('Mean Correlation of Clustered Spines', 'Fontsize', 14)
        xlim([0 15])
        
    subplot(3,2,4)
        plot(nanmean(MeanCorrelationBetweenMovementSpines,1), 'Color', black, 'Linewidth', 2); hold on;
        plot(AllMeanCorrelationBetweenFarMovementSpines, 'Color', gray)
        legend({'Mov spines same dend', 'Mov spines sep dends'})
        r_errorbar(1:14, nanmean(MeanCorrelationBetweenMovementSpines,1), MeanCorrelationBetweenMovementSpinesSEM, black)
        r_errorbar(1:14, AllMeanCorrelationBetweenFarMovementSpines, MeanCorrelationBetweenFarMovementSpinesSEM, gray)
        xlim([0 15])
        ylabel('Correlation', 'Fontsize', 14)
        xlabel('Session', 'Fontsize', 14)
        title('Correlation between mov spines', 'Fontsize', 14);

    subplot(3,2,5)
        plot(MeanCueRelClusteredSpines_CorrwithDend, 'Color', lgreen, 'Linewidth', 2); hold on;
        plot(MeanMovRelClusteredSpines_CorrwithDend, 'Color', black, 'LineWidth', 2);
        plot(MeanSucRelClusteredSpines_CorrwithDend, 'Color', lblue, 'Linewidth', 2)
        plot(MeanRewRelClusteredSpines_CorrwithDend, 'Color', purple, 'Linewidth', 2)
        
        legend({'Cue rel clusters', 'Mov-rel clusters', 'Suc-rel clusters', 'Rew-rel clusters'})
        
        r_errorbar(1:14, MeanCueRelClusteredSpines_CorrwithDend, CueRelClusteredSpines_CorrwithDendSEM, lgreen)
        r_errorbar(1:14, MeanMovRelClusteredSpines_CorrwithDend, MovRelClusteredSpines_CorrwithDendSEM, black)
        r_errorbar(1:14, MeanSucRelClusteredSpines_CorrwithDend, SucRelClusteredSpines_CorrwithDendSEM, lblue)
        r_errorbar(1:14, MeanRewRelClusteredSpines_CorrwithDend, RewRelClusteredSpines_CorrwithDendSEM, purple)
        
        xlim([0 15])
        ylabel('Correlation with Dendrite')
        xlabel('Session')
        title('Functional Clusters')
       
    subplot(3,2,6)
        plot(MeanCueRelCausalClusteredSpines_CorrwithDend, 'Color', lgreen, 'Linewidth', 2); hold on;
        plot(MeanMovRelCausalClusteredSpines_CorrwithDend, 'Color', black, 'LineWidth', 2);
        plot(MeanSucRelCausalClusteredSpines_CorrwithDend, 'Color', lblue, 'Linewidth', 2)
        plot(MeanRewRelCausalClusteredSpines_CorrwithDend, 'Color', purple, 'Linewidth', 2)
        
        legend({'Cue rel clusters', 'Mov-rel clusters', 'Suc-rel clusters', 'Rew-rel clusters'})
        
        r_errorbar(1:14, MeanCueRelCausalClusteredSpines_CorrwithDend, CueRelCausalClusteredSpines_CorrwithDendSEM, lgreen)
        r_errorbar(1:14, MeanMovRelCausalClusteredSpines_CorrwithDend, MovRelCausalClusteredSpines_CorrwithDendSEM, black)
        r_errorbar(1:14, MeanSucRelCausalClusteredSpines_CorrwithDend, SucRelCausalClusteredSpines_CorrwithDendSEM, lblue)
        r_errorbar(1:14, MeanRewRelCausalClusteredSpines_CorrwithDend, RewRelCausalClusteredSpines_CorrwithDendSEM, purple)
        
        xlim([0 15])
        title('Causal Functional Clusters')
        xlabel('Session')
        ylabel('Correlation with Dendrite')
        
    %%%
    %%% Figure 6: Fraction of clusters that are movement related
    %%%
    
%     figure('Position', scrsz); hold on;
%     plot(MeanFractionofClusterThatsMovementRelated,'Color',black, 'LineWidth', 2)
%     plot(MeanFractionofCausalClusterThatsMovementRelated, 'Color', bgreen, 'Linewidth',2)
%     
%     legend({'Synapse only clusters', 'Causal Clusters'}, 'Location', 'SouthEast')
%     r_errorbar(1:14, MeanFractionofClusterThatsMovementRelated, FractionofClusterThatsMovementRelatedSEM, 'k')
%     r_errorbar(1:14, MeanFractionofCausalClusterThatsMovementRelated, FractionofCausalClusterThatsMovementRelatedSEM, bgreen)
%     
%     xlabel('Session')
%     ylabel('Fraction of Cluster That is Movement Related', 'Fontsize', 14)
%     ylim([0 1.2])
    
    
    %%%
    %%% Figure 7: Spectral graph analysis of clusters
    %%%
                
    figure('Position', scrsz); hold on;
    subplot(2,4,1)
    plot(1:14, MeanSpatialDegree, 'Color', black, 'LineWidth', 2); hold on;
    plot(1:14, MeanTemporalDegree, 'Color',red, 'Linewidth', 2);
    plot(1:14, MeanSpatioTemporalDegree, 'Color',green, 'Linewidth', 2); 

    legend({'Spatial Degree', 'Temporal Degree', 'Spatiotemporal Degree'});

    r_errorbar(1:14, MeanSpatioTemporalDegree, SpatioTemporalDegreeSEM, green)
    r_errorbar(1:14, MeanSpatialDegree, SpatialDegreeSEM, black)
    r_errorbar(1:14, MeanTemporalDegree, TemporalDegreeSEM, red)
    ylabel('Mean Degree', 'Fontsize', 14)
    xlabel('Session', 'Fontsize', 14)
    xlim([0 15])
    set(gca, 'XTick', 0:15); set(gca, 'XTickLabel', 0:15)
    
    subplot(2,4,2)
    plot(1:14, MeanSpatialMovementCorrelation, 'Color',black, 'Linewidth', 2); hold on;
    plot(1:14, MeanTemporalMovementCorrelation, 'Color',red, 'Linewidth', 2);
    plot(1:14, MeanSpatioTemporalMovementCorrelation, 'Color',green, 'Linewidth', 2)
        r_errorbar(1:14, MeanSpatialMovementCorrelation, SpatioTemporalDegreeSEM, black)
        r_errorbar(1:14, MeanTemporalMovementCorrelation, TemporalMovementCorrelationSEM, red)
        r_errorbar(1:14, MeanSpatioTemporalMovementCorrelation, SpatioTemporalMovementCorrelationSEM, green)
    ylabel([{'Mean Correlation of Spatiotemporal'}, {'Degree with Movement'}],'Fontsize', 14)
    xlabel('Session', 'Fontsize', 14);
    xlim([0 15])
    set(gca, 'XTick', 0:15); set(gca, 'XTickLabel', 0:15)
    
    subplot(2,4,3)
    xlim([0 15])
    plot(1:14,MeanDendSpatialFiedlerValue, 'Color',black, 'Linewidth', 2); hold on;
    plot(1:14,MeanDendTemporalFiedlerValue, 'Color',red, 'Linewidth', 2);
    plot(1:14,MeanDendSpatioTemporalFiedlerValue, 'Color',green, 'Linewidth', 2);
    legend({'Spatial Fiedler', 'Temporal Fiedler', 'Spatiotemporal Fiedler'});
    r_errorbar(1:14, MeanDendSpatialFiedlerValue, DendSpatialFiedlerValueSEM, black);
    r_errorbar(1:14, MeanDendTemporalFiedlerValue, DendTemporalFiedlerValueSEM, red);
    r_errorbar(1:14, MeanDendSpatioTemporalFiedlerValue, DendSpatioTemporalFiedlerValueSEM, green);
    
    ylabel('Mean Algebraic Connectivity of Dendrites (Clustering)')
    xlabel('Session')
    xlim([0 15])
    set(gca, 'XTick', 0:15); set(gca, 'XTickLabel', 0:15)
    
    subplot(2,4,4)
    plot(MeanSpatioTemporalOverlap, 'Color',black, 'Linewidth', 2)
    r_errorbar(1:14, MeanSpatioTemporalOverlap, SpatioTemporalOverlapSEM, 'k')
    ylabel('Mean Correlation of Spatial and Temporal Eigenvectors')
    xlabel('Session')
    xlim([0 15])
    set(gca, 'XTick', 0:15); set(gca, 'XTickLabel', 0:15)
    
    subplot(2,4,5)
    plot(MeanSpatialDegreeofCueSpines, 'Color', lgreen, 'Linewidth', 2); hold on;
    plot(MeanSpatialDegreeofMovementSpines, 'Color', black, 'Linewidth', 2)
    plot(MeanSpatialDegreeofMovementDuringCueSpines, 'Color', green , 'Linewidth', 2)
    plot(MeanSpatialDegreeofPreSuccessSpines, 'Color', bgreen, 'Linewidth', 2)
    plot(MeanSpatialDegreeofSuccessSpines, 'Color', lblue, 'Linewidth', 2)
    plot(MeanSpatialDegreeofRewardSpines, 'Color', purple, 'Linewidth', 2)
    
    xlim([0 15])
    xlabel('Session', 'Fontsize', 14)
    ylabel('Mean Degree')
    legend({'Cue spines', 'Movement Spines', 'Success Spines', 'Reward Spines'})
    title([{'Mean Spatial Degree of'},{'feature-related spines'}], 'Fontsize', 14)
    
    r_errorbar(1:14, MeanSpatialDegreeofCueSpines, SpatialDegreeofCueSpinesSEM, lgreen)
    r_errorbar(1:14, MeanSpatialDegreeofMovementSpines, SpatialDegreeofMovementSpinesSEM, black)
    r_errorbar(1:14, MeanSpatialDegreeofMovementDuringCueSpines, SpatialDegreeofMovementDuringCueSpinesSEM, green)
    r_errorbar(1:14, MeanSpatialDegreeofPreSuccessSpines, SpatialDegreeofPreSuccessSpinesSEM, bgreen)
    r_errorbar(1:14, MeanSpatialDegreeofSuccessSpines, SpatialDegreeofSuccessSpinesSEM, lblue)
    r_errorbar(1:14, MeanSpatialDegreeofRewardSpines, SpatialDegreeofRewardSpinesSEM, purple)
    
    subplot(2,4,6)
    plot(MeanTemporalDegreeofCueSpines, 'Color', lgreen, 'Linewidth', 2); hold on;
    plot(MeanTemporalDegreeofMovementSpines, 'Color', black, 'Linewidth', 2)
    plot(MeanTemporalDegreeofMovementDuringCueSpines, 'Color', green, 'Linewidth', 2)
    plot(MeanTemporalDegreeofPreSuccessSpines, 'Color', bgreen, 'Linewidth', 2)
    plot(MeanTemporalDegreeofSuccessSpines, 'Color', lblue, 'Linewidth', 2)
    plot(MeanTemporalDegreeofRewardSpines, 'Color', purple, 'Linewidth', 2)    
    xlim([0 15])
    xlabel('Session', 'Fontsize', 14)
    ylabel('Mean Degree')
    legend({'Cue spines', 'Movement Spines', 'Success Spines', 'Reward Spines'})
    title([{'Mean Temporal Degree of'},{'feature-related spines'}], 'Fontsize', 14)
    
    r_errorbar(1:14, MeanTemporalDegreeofCueSpines, TemporalDegreeofCueSpinesSEM, lgreen)
    r_errorbar(1:14, MeanTemporalDegreeofMovementSpines, TemporalDegreeofMovementSpinesSEM, black)
    r_errorbar(1:14, MeanTemporalDegreeofMovementDuringCueSpines, TemporalDegreeofMovementDuringCueSpinesSEM, green)
    r_errorbar(1:14, MeanTemporalDegreeofPreSuccessSpines, TemporalDegreeofPreSuccessSpinesSEM, bgreen)
    r_errorbar(1:14, MeanTemporalDegreeofSuccessSpines, TemporalDegreeofSuccessSpinesSEM, lblue)
    r_errorbar(1:14, MeanTemporalDegreeofRewardSpines, TemporalDegreeofRewardSpinesSEM, purple)
    
    subplot(2,4,7)
    plot(MeanSpatioTemporalDegreeofCueSpines, 'Color', lgreen, 'Linewidth', 2); hold on;
    plot(MeanSpatioTemporalDegreeofMovementSpines, 'Color', black, 'Linewidth', 2)
    plot(MeanSpatioTemporalDegreeofMovementDuringCueSpines, 'Color', green, 'Linewidth', 2)
    plot(MeanSpatioTemporalDegreeofPreSuccessSpines, 'Color', bgreen, 'Linewidth', 2)
    plot(MeanSpatioTemporalDegreeofSuccessSpines, 'Color', lblue, 'Linewidth', 2)
    plot(MeanSpatioTemporalDegreeofRewardSpines, 'Color', purple, 'Linewidth', 2)
    
    xlim([0 15])
    xlabel('Session', 'Fontsize', 14)
    ylabel('Mean Degree')
    legend({'Cue spines', 'Movement Spines', 'Success Spines', 'Reward Spines'})
    title([{'Mean SpatioTemporal Degree of'},{'feature-related spines'}], 'Fontsize', 14)
    
    r_errorbar(1:14, MeanSpatioTemporalDegreeofCueSpines, SpatioTemporalDegreeofCueSpinesSEM, lgreen)
    r_errorbar(1:14, MeanSpatioTemporalDegreeofMovementSpines, SpatioTemporalDegreeofMovementSpinesSEM, black)
    r_errorbar(1:14, MeanSpatioTemporalDegreeofMovementDuringCueSpines, SpatioTemporalDegreeofMovementDuringCueSpinesSEM, green)
    r_errorbar(1:14, MeanSpatioTemporalDegreeofPreSuccessSpines, SpatioTemporalDegreeofPreSuccessSpinesSEM, bgreen)
    r_errorbar(1:14, MeanSpatioTemporalDegreeofSuccessSpines, SpatioTemporalDegreeofSuccessSpinesSEM, lblue)
    r_errorbar(1:14, MeanSpatioTemporalDegreeofRewardSpines, SpatioTemporalDegreeofRewardSpinesSEM, purple)
    
    %%%
    %%% Figure 8: Correlation vs. Distance Distributions
    %%%
    
    figure('Position', scrsz)
    currentplot = 1;
    subplot(2,4,1)
        xdata = [AllDistancesBetweenMovementSpines{1},AllDistancesBetweenMovementSpines{2}]';
        ydata = [CorrelationBetweenMovementSpinesMovePeriods{1}, CorrelationBetweenMovementSpinesMovePeriods{2}]';
        ydata(isnan(ydata)) = 0;
        %%% K-means clustering 
        X = [(xdata-nanmean(xdata))/nanstd(xdata), (ydata-nanmean(ydata))/nanstd(ydata)];   %%% Standardized data!!!!
        [idx, C] = kmeans(X,2);
%         x1 = min(X(:,1)):0.01:max(X(:,1));
%         x2 = min(X(:,2)):0.01:max(X(:,2));
%         [x1G,x2G] = meshgrid(x1,x2);
%         XGrid = [x1G(:),x2G(:)]; % Defines a fine grid on the plot
%         col1 = orange; col2 = lblue;
%         idx2Region = kmeans(XGrid,2,'MaxIter',1,'Start',C);
%         gscatter(XGrid(:,1),XGrid(:,2),idx2Region,...
%         [col1; col2],'..');hold on;
%         % Assigns each node in the grid to the closest centroid
%         %%%
        plot(xdata(idx==1), ydata(idx==1), '.', 'Color', lgreen); hold on;
        plot(xdata(idx==2), ydata(idx==2), '.k')
%         decay = fit(xdata, ydata, 'exp1'); 
%             fline = plot(decay); 
%             set(fline, 'Color', 'k')
%             legend off
%     %         plot(-1/decay.b,0.368*decay.a, 'ok', 'MarkerFaceColor', 'k') %%% 0.368 corresponds to the decay constant, tau
        xlim([0 100])
        ylim([-0.05 1])
        xlabel('Distance (\mum)', 'FontSize', 14)
        ylabel('Correlation', 'FontSize', 14)
        title('Movement Spines, Sessions 1-2',  'FontSize', 14)
        corratbin = cell(1,8);
        highcorratbin = cell(1,8);
        binstep = 5;
        bincount = 1;
        ydata2 = ydata(ydata>0.5);
        xdata2 = xdata(ydata>0.5);
        for i = 1:binstep:50
            corratbin{currentplot}(1,bincount) = nanmean(ydata(find(xdata>=(i-1) & xdata<(i+binstep))));
            highcorratbin{currentplot}(1,bincount) = nanmean(ydata2(find(xdata2>=(i-1) & xdata2<(i+binstep))));
            bincount = bincount+1;
        end
        pos = get(gca,'Position');
        axes('Position', [pos(1)+0.7*pos(3), pos(2)+0.7*pos(4), 0.35*pos(3), 0.25*pos(4)]);
        bar(highcorratbin{currentplot}, 'FaceColor', lgreen, 'EdgeColor', 'k'); hold on;
        bar(corratbin{currentplot}, 'FaceColor', 'k', 'EdgeColor', gray)
        xlim([-1 11])
        ylim([0 1])
    currentplot = 2;
    subplot(2,4,currentplot)
        xdata = [AllDistancesBetweenFarMovementSpines{1}, AllDistancesBetweenFarMovementSpines{2}]';
        ydata = [CorrelationBetweenFarMovementSpines{1}, CorrelationBetweenFarMovementSpines{2}]';
        ydata(isnan(ydata)) = 0;
        try
            plot(xdata(ydata>=0.5), ydata(ydata>=0.5), '.', 'Color', lgreen); hold on;
        catch
        end
        plot(xdata(ydata<0.5), ydata(ydata<0.5), '.', 'Color', gray)
%         decay = fit(xdata, ydata, 'exp1'); 
%             fline = plot(decay); 
%             set(fline, 'Color', 'k')
%             legend off
%     %         plot(-1/decay.b,0.368*decay.a, 'ok', 'MarkerFaceColor', 'k') %%% 0.368 corresponds to the decay constant, tau
        xlim([0 100])
        ylim([-0.05 1])
        xlabel('Distance (\mum)', 'FontSize', 14)
        ylabel('Correlation', 'FontSize', 14)
        title('Movement Spines, Sessions 1-2',  'FontSize', 14)
        bincount = 1;
        ydata2 = ydata(ydata>0.5);
        xdata2 = xdata(ydata>0.5);
        for i = 1:binstep:50
            corratbin{currentplot}(1,bincount) = nanmean(ydata(find(xdata>=(i-1) & xdata<(i+binstep))));
            highcorratbin{currentplot}(1,bincount) = nanmean(ydata2(find(xdata2>=(i-1) & xdata2<(i+binstep))));
            bincount = bincount+1;
        end
        pos = get(gca,'Position');
        axes('Position', [pos(1)+0.7*pos(3), pos(2)+0.7*pos(4), 0.35*pos(3), 0.25*pos(4)]);
        bar(highcorratbin{currentplot}, 'FaceColor', lgreen, 'EdgeColor', 'k'); hold on;
        bar(corratbin{currentplot}, 'FaceColor', 'k', 'EdgeColor', gray)
        xlim([-1 11])
        ylim([0 1])
    currentplot = 3;
    subplot(2,4,currentplot)
        xdata = [AllDistancesBetweenMovementSpines{10},AllDistancesBetweenMovementSpines{11}]';
        ydata = [CorrelationBetweenMovementSpinesMovePeriods{10},CorrelationBetweenMovementSpinesMovePeriods{11}]';
        ydata(isnan(ydata)) = 0;
        X = [(xdata-nanmean(xdata))/nanstd(xdata), (ydata-nanmean(ydata))/nanstd(ydata)];   %%% Standardized data!!!!
        [idx, C] = kmeans(X,2);
%         x1 = min(X(:,1)):0.01:max(X(:,1));
%         x2 = min(X(:,2)):0.01:max(X(:,2));
%         [x1G,x2G] = meshgrid(x1,x2);
%         XGrid = [x1G(:),x2G(:)]; % Defines a fine grid on the plot
%         col1 = orange; col2 = lblue;
%         idx2Region = kmeans(XGrid,2,'MaxIter',1,'Start',C);
%         gscatter(XGrid(:,1),XGrid(:,2),idx2Region,...
%         [col1; col2],'..');hold on;
%         % Assigns each node in the grid to the closest centroid
%         %%%
        plot(xdata(idx==1), ydata(idx==1), '.', 'Color', lgreen); hold on;
        plot(xdata(idx==2), ydata(idx==2), '.k')
%             decay = fit(xdata, ydata, 'exp1'); 
%             fline = plot(decay); 
%             set(fline, 'Color', 'k')
%             legend off
    %         plot(-1/decay.b,0.368*decay.a, 'ok', 'MarkerFaceColor', 'k') %%% 0.368 corresponds to the decay constant, tau
        xlim([0 100])
        ylim([-0.05 1])
        xlabel('Distance (\mum)', 'FontSize', 14)
        ylabel('Correlation', 'FontSize', 14)
        title('Movement Spines, Sessions 10-11', 'FontSize', 14)
        bincount = 1;
        ydata2 = ydata(ydata>0.5);
        xdata2 = xdata(ydata>0.5);
        for i = 1:binstep:50
            corratbin{currentplot}(1,bincount) = nanmean(ydata(find(xdata>=(i-1) & xdata<(i+binstep))));
            highcorratbin{currentplot}(1,bincount) = nanmean(ydata2(find(xdata2>=(i-1) & xdata2<(i+binstep))));
            bincount = bincount+1;
        end
        pos = get(gca,'Position');
        axes('Position', [pos(1)+0.7*pos(3), pos(2)+0.7*pos(4), 0.35*pos(3), 0.25*pos(4)]);
        bar(highcorratbin{currentplot}, 'FaceColor', lgreen, 'EdgeColor', 'k'); hold on;
        bar(corratbin{currentplot}, 'FaceColor', 'k', 'EdgeColor', gray)
        xlim([-1 11])
        ylim([0 1])
    currentplot = 4;
    subplot(2,4,currentplot)
        xdata = [AllDistancesBetweenFarMovementSpines{10},AllDistancesBetweenFarMovementSpines{11}]';
        ydata = [CorrelationBetweenFarMovementSpines{10},CorrelationBetweenFarMovementSpines{11}]';
        ydata(isnan(ydata)) = 0;
        try
            plot(xdata(ydata>=0.5), ydata(ydata>=0.5), '.', 'Color', lgreen); hold on;
        catch
        end
        plot(xdata(ydata<0.5), ydata(ydata<0.5), '.', 'Color', gray)
%             decay = fit(xdata, ydata, 'exp1'); 
%             fline = plot(decay); 
%             set(fline, 'Color', 'k')
%             legend off
%     %         plot(-1/decay.b,0.368*decay.a, 'ok', 'MarkerFaceColor', 'k') %%% 0.368 corresponds to the decay constant, tau
        xlim([0 100])
        ylim([-0.05 1])
        xlabel('Distance (\mum)', 'FontSize', 14)
        ylabel('Correlation', 'FontSize', 14)
        title('Movement Spines, Sessions 10-11', 'FontSize', 14)
        bincount = 1;
        ydata2 = ydata(ydata>0.5);
        xdata2 = xdata(ydata>0.5);
        for i = 1:binstep:50
            corratbin{currentplot}(1,bincount) = nanmean(ydata(find(xdata>=(i-1) & xdata<(i+binstep))));
            highcorratbin{currentplot}(1,bincount) = nanmean(ydata2(find(xdata2>=(i-1) & xdata2<(i+binstep))));
            bincount = bincount+1;
        end
        pos = get(gca,'Position');
        axes('Position', [pos(1)+0.7*pos(3), pos(2)+0.7*pos(4), 0.35*pos(3), 0.25*pos(4)]);
        bar(highcorratbin{currentplot}, 'FaceColor', lgreen, 'EdgeColor', 'k'); hold on;
        bar(corratbin{currentplot}, 'FaceColor', 'k', 'EdgeColor', gray)
        xlim([-1 11])
        ylim([0 1])
        
    currentplot = 5;
    subplot(2,4,currentplot)
        xdata = [AllDistancesBetweenAllSpines{1}, AllDistancesBetweenAllSpines{2}]';
        ydata = [CorrelationBetweenAllSpinesMovePeriods{1}, CorrelationBetweenAllSpinesMovePeriods{2}]';
        ydata(isnan(ydata)) = 0;
        X = [(xdata-nanmean(xdata))/nanstd(xdata), (ydata-nanmean(ydata))/nanstd(ydata)];   %%% Standardized data!!!!
        [idx, C] = kmeans(X,2);
%         x1 = min(X(:,1)):0.01:max(X(:,1));
%         x2 = min(X(:,2)):0.01:max(X(:,2));
%         [x1G,x2G] = meshgrid(x1,x2);
%         XGrid = [x1G(:),x2G(:)]; % Defines a fine grid on the plot
%         col1 = orange; col2 = lblue;
%         idx2Region = kmeans(XGrid,2,'MaxIter',1,'Start',C);
%         gscatter(XGrid(:,1),XGrid(:,2),idx2Region,...
%         [col1; col2],'..');hold on;
%         % Assigns each node in the grid to the closest centroid
%         %%%
        plot(xdata(idx==1), ydata(idx==1), '.', 'Color', dred); hold on;
        plot(xdata(idx==2), ydata(idx==2), '.k')
            decay = fit(xdata, ydata, 'exp1'); 
            fline = plot(decay); 
            set(fline, 'Color', 'k')
            legend off
    %         plot(-1/decay.b,0.368*decay.a, 'ok', 'MarkerFaceColor', 'k') %%% 0.368 corresponds to the decay constant, tau
        xlim([0 100])
        ylim([-0.05 1])
        xlabel('Distance (\mum)', 'FontSize', 14)
        ylabel('Correlation', 'FontSize', 14)
        title('All Spines, Sessions 1-2', 'FontSize', 14)
        
        bincount = 1;
        ydata2 = ydata(ydata>0.5);
        xdata2 = xdata(ydata>0.5);
        for i = 1:binstep:50
            corratbin{currentplot}(1,bincount) = nanmean(ydata(find(xdata>=(i-1) & xdata<(i+binstep))));
            highcorratbin{currentplot}(1,bincount) = nanmean(ydata2(find(xdata2>=(i-1) & xdata2<(i+binstep))));
            bincount = bincount+1;
        end
        pos = get(gca,'Position');
        axes('Position', [pos(1)+0.7*pos(3), pos(2)+0.7*pos(4), 0.35*pos(3), 0.25*pos(4)]);
        bar(highcorratbin{currentplot}, 'FaceColor', dred, 'EdgeColor', 'k'); hold on;
        bar(corratbin{currentplot}, 'FaceColor', 'k', 'EdgeColor', gray)
        xlim([-1 11])
        ylim([0 1])
        
    currentplot = 6;
    subplot(2,4,currentplot)
        xdata = [AllDistancesBetweenFarSpines{1}; AllDistancesBetweenFarSpines{2}]';
        ydata = [CorrelationBetweenFarSpines{1}; CorrelationBetweenFarSpines{2}]';
        ydata(isnan(ydata)) = 0;
        try
            plot(xdata(ydata>=0.5), ydata(ydata>=0.5), '.', 'Color', dred); hold on;
        catch
        end
        plot(xdata(ydata<0.5), ydata(ydata<0.5), '.', 'Color', gray)
%             decay = fit(xdata, ydata, 'exp1'); 
%             fline = plot(decay); 
%             set(fline, 'Color', 'k')
%             legend off
    %         plot(-1/decay.b,0.368*decay.a, 'ok', 'MarkerFaceColor', 'k') %%% 0.368 corresponds to the decay constant, tau
        xlim([0 100])
        ylim([-0.05 1])
        xlabel('Distance (\mum)', 'FontSize', 14)
        ylabel('Correlation', 'FontSize', 14)
        title('All Spines, Sessions 1-2', 'FontSize', 14)
        bincount = 1;
        ydata2 = ydata(ydata>0.5);
        xdata2 = xdata(ydata>0.5);
        for i = 1:binstep:50
            corratbin{currentplot}(1,bincount) = nanmean(ydata(find(xdata>=(i-1) & xdata<(i+binstep))));
            highcorratbin{currentplot}(1,bincount) = nanmean(ydata2(find(xdata2>=(i-1) & xdata2<(i+binstep))));
            bincount = bincount+1;
        end
        pos = get(gca,'Position');
        axes('Position', [pos(1)+0.7*pos(3), pos(2)+0.7*pos(4), 0.35*pos(3), 0.25*pos(4)]);
        bar(highcorratbin{currentplot}, 'FaceColor', dred, 'EdgeColor', 'k'); hold on;
        bar(corratbin{currentplot}, 'FaceColor', 'k', 'EdgeColor', gray)
        xlim([-1 11])
        ylim([0 1])
    currentplot = 7;
    subplot(2,4,currentplot)
        xdata = [AllDistancesBetweenAllSpines{10}, AllDistancesBetweenAllSpines{11}]';
        ydata = [CorrelationBetweenAllSpinesMovePeriods{10}, CorrelationBetweenAllSpinesMovePeriods{11}]';
        ydata(isnan(ydata)) = 0;
        X = [(xdata-nanmean(xdata))/nanstd(xdata), (ydata-nanmean(ydata))/nanstd(ydata)];   %%% Standardized data!!!!
        [idx, C] = kmeans(X,2);
%         x1 = min(X(:,1)):0.01:max(X(:,1));
%         x2 = min(X(:,2)):0.01:max(X(:,2));
%         [x1G,x2G] = meshgrid(x1,x2);
%         XGrid = [x1G(:),x2G(:)]; % Defines a fine grid on the plot
%         col1 = orange; col2 = lblue;
%         idx2Region = kmeans(XGrid,2,'MaxIter',1,'Start',C);
%         gscatter(XGrid(:,1),XGrid(:,2),idx2Region,...
%         [col1; col2],'..');hold on;
%         % Assigns each node in the grid to the closest centroid
%         %%%
        plot(xdata(idx==1), ydata(idx==1), '.', 'Color', dred); hold on;
        plot(xdata(idx==2), ydata(idx==2), '.k')
            decay = fit(xdata, ydata, 'exp1'); 
            fline = plot(decay); 
            set(fline, 'Color', 'k')
            legend off
            plot(-1/decay.b,0.368*decay.a, 'ok', 'MarkerFaceColor', 'k') %%% 0.368 corresponds to the decay constant, tau
        xlim([0 100])
        ylim([-0.05 1])
        xlabel('Distance (\mum)', 'FontSize', 14)
        ylabel('Correlation', 'FontSize', 14)
        title('All Spines, Sessions 10-11', 'Fontsize', 14)
        bincount = 1;
        ydata2 = ydata(ydata>0.5);
        xdata2 = xdata(ydata>0.5);
        for i = 1:binstep:50
            corratbin{currentplot}(1,bincount) = nanmean(ydata(find(xdata>=(i-1) & xdata<(i+binstep))));
            highcorratbin{currentplot}(1,bincount) = nanmean(ydata2(find(xdata2>=(i-1) & xdata2<(i+binstep))));
            bincount = bincount+1;
        end
        pos = get(gca,'Position');
        axes('Position', [pos(1)+0.7*pos(3), pos(2)+0.7*pos(4), 0.35*pos(3), 0.25*pos(4)]);
        bar(highcorratbin{currentplot}, 'FaceColor', dred, 'EdgeColor', 'k'); hold on;
        bar(corratbin{currentplot}, 'FaceColor', 'k', 'EdgeColor', gray)
        xlim([-1 11])
        ylim([0 1])
    currentplot = 8;
    subplot(2,4,currentplot)
        xdata = [AllDistancesBetweenFarSpines{10}; AllDistancesBetweenFarSpines{11}]';
        ydata = [CorrelationBetweenFarSpines{10}; CorrelationBetweenFarSpines{11}]';
        ydata(isnan(ydata)) = 0;
        try
            plot(xdata(ydata>=0.5), ydata(ydata>=0.5), '.', 'Color', dred); hold on;
        catch
        end
        plot(xdata(ydata<0.5), ydata(ydata<0.5), '.', 'Color', gray)
%             decay = fit(xdata, ydata, 'exp1'); 
%             fline = plot(decay); 
%             set(fline, 'Color', 'k')
%             legend off
    %         plot(-1/decay.b,0.368*decay.a, 'ok', 'MarkerFaceColor', 'k') %%% 0.368 corresponds to the decay constant, tau
        xlim([0 100])
        ylim([-0.05 1])
        xlabel('Distance (\mum)', 'FontSize', 14)
        ylabel('Correlation', 'FontSize', 14)
        title('All Spines, Sessions 10-11', 'Fontsize', 14)
        bincount = 1;
        ydata2 = ydata(ydata>0.5);
        xdata2 = xdata(ydata>0.5);
        for i = 1:binstep:50
            corratbin{currentplot}(1,bincount) = nanmean(ydata(find(xdata>=(i-1) & xdata<(i+binstep))));
            highcorratbin{currentplot}(1,bincount) = nanmean(ydata2(find(xdata2>=(i-1) & xdata2<(i+binstep))));
            bincount = bincount+1;
        end
        pos = get(gca,'Position');
        axes('Position', [pos(1)+0.7*pos(3), pos(2)+0.7*pos(4), 0.35*pos(3), 0.25*pos(4)]);
        bar(highcorratbin{currentplot}, 'FaceColor', dred, 'EdgeColor', 'k'); hold on;
        bar(corratbin{currentplot}, 'FaceColor', 'k', 'EdgeColor', gray)
        xlim([-1 11])
        ylim([0 1])
        

        
    %%% Figure 9

    figure('Position', scrsz, 'Name', 'Clustering Validation'); hold on;
    %%% Shuffled data plots: uncomment the following and set stop point

    corrtouse{1} = [CorrelationBetweenAllSpinesMovePeriods{1}, CorrelationBetweenAllSpinesMovePeriods{2}];
    disttouse{1} = [AllDistancesBetweenAllSpines{1}, AllDistancesBetweenAllSpines{2}];
    %         farcorrtouse = CorrelationBetweenFarSpines{sessiontouse};

    corrtouse{2} = [CorrelationBetweenAllSpinesMovePeriods{10}, CorrelationBetweenAllSpinesMovePeriods{11}];
    disttouse{2} = [AllDistancesBetweenAllSpines{10}, AllDistancesBetweenAllSpines{11}];

    corrtouse{3} = [CorrelationBetweenMovementSpinesMovePeriods{1}, CorrelationBetweenMovementSpinesMovePeriods{2}];
    disttouse{3} = [AllDistancesBetweenMovementSpines{1}, AllDistancesBetweenMovementSpines{2}];
    %         farcorrtouse = CorrelationBetweenFarSpines{sessiontouse};

    corrtouse{4} = [CorrelationBetweenMovementSpinesMovePeriods{10}, CorrelationBetweenMovementSpinesMovePeriods{11}];
    disttouse{4} = [AllDistancesBetweenMovementSpines{10}, AllDistancesBetweenMovementSpines{11}];
    
    usenorm = 1;

    for i = 1:4
        if i<3
            subplot(2,4,i); hold on;
            color = dred;
        else
            subplot(2,4,i+2); hold on;
            color = lgreen;
        end
        if i == 1 || i == 3
            sessiontag = '1-2';
        elseif i == 2 || i ==4 
            sessiontag = '10-11';
        end
        [sortedDistances, sortedDistIndices] = sort(disttouse{i});
    %         [sortedFarDistances, sortedFarDistInd] = sort(AllDistancesBetweenFarSpines{sessiontouse});
        shuffled = [];
        sortedshuffled = [];
        for j = 1:1000
            shuffled(1:length(corrtouse{i}),j) = corrtouse{i}(randperm(length(corrtouse{i})));
            shuffled(isnan(shuffled(:,j)),:) = 0;
            sortedshuffled(:,j) = shuffled(sortedDistIndices,j)./nansum(shuffled(sortedDistIndices,j));
            sortedshuffled(isnan(sortedshuffled(:,j)),:) = 0;
            plot(sortedDistances, cumsum(sortedshuffled(:,j)),'color', [0.5 0.5 0.5]);
        %             farshuffled(1:length(farcorrtouse),j) = farcorrtouse(randperm(length(farcorrtouse)));
        %             sortedFarshuffled(:,j) = farshuffled(sortedFarDistInd,j)./nansum(farshuffled(sortedFarDistInd,j));
        %             plot(sortedFarDistances,cumsum(sortedFarshuffled(:,j)),'color', [0.5 0.5 0.5]);
        end
        plot(sortedDistances, cumsum(nanmean(sortedshuffled,2)), 'k', 'Linewidth', 2)
        Correlations_fraction = corrtouse{i}(sortedDistIndices)./nansum(corrtouse{i}(sortedDistIndices));
        Correlations_fraction(isnan(Correlations_fraction))= 0;
        plot(sortedDistances,cumsum(Correlations_fraction),'color',color, 'LineWidth', 3)
        plot(sortedDistances,cumsum(Correlations_fraction)-cumsum(nanmean(sortedshuffled,2))', 'color', green, 'Linewidth', 2)
        ylabel('Cumulative Correlation', 'Fontsize', 14)
        xlabel('Distance (\mum)')
        title(['Session(s) ', sessiontag], 'Fontsize', 14)
        xlim([0 100])
        ylim([0 1])
    end

    subplot(2,4,3); hold on;
        neardist = [NearestMovementRelatedSpine{1},NearestMovementRelatedSpine{2}];
        n = hist(neardist,round(max(neardist)/5));
        if usenorm
            n = n/sum(n);
        else
        end
        bar(n, 'FaceColor', lpurple); hold on;
        nextdist = [NextClosestMovementRelatedSpine{1}, NextClosestMovementRelatedSpine{2}];
        nx = hist(nextdist, round(max(nextdist)/5));
        if usenorm
            nx = nx/sum(nx);
        else
        end
        bar(nx, 'FaceColor', orange);
        thirddist = [ThirdClosestMovementRelatedSpine{1}, ThirdClosestMovementRelatedSpine{2}];
        nt = hist(thirddist, round(max(thirddist)/5));
        if usenorm
            nt = nt/sum(nt);
        else
        end
        bar(nt, 'FaceColor', lgreen);
        fourthdist = [FourthClosestMovementRelatedSpine{1}, FourthClosestMovementRelatedSpine{2}];
        nf = hist(fourthdist, round(max(fourthdist)/5));
        if usenorm
            nf = nf/sum(nf);
            ylim([0 1])
        else
        end
        bar(nf, 'FaceColor', blue)
        legend({'Nearest MR spine', 'Next Closest', 'Third Closest', 'Fourth'}, 'Location', 'Northwest')
        title('Session 1-2', 'Fontsize',14)
        set(gca, 'XTick', 0:30, 'XTickLabel', mat2cell(num2str([0:5:150]')), 'Fontsize', 6)
        xlabel('Distance bins', 'Fontsize', 14)
        ylabel('Fraction', 'Fontsize', 14)
        
        text(1.5, n(1), ['mean = ', num2str(nanmean(neardist))], 'Color', lpurple);
        text(2.5, nx(1), ['mean = ', num2str(nanmean(nextdist))], 'Color', orange);
        text(3.5, nt(1), ['mean = ', num2str(nanmean(thirddist))], 'Color', lgreen);
        text(4.5, nf(1), ['mean = ', num2str(nanmean(fourthdist))], 'Color', blue);
        
        pos = get(gca,'Position');
        axes('Position', [pos(1)+0.7*pos(3), pos(2)+0.7*pos(4), 0.25*pos(3), 0.25*pos(4)], 'Fontsize', 6);
        nearcorr = [CorrelationwithNearestMovementRelatedSpine{1}, CorrelationwithNearestMovementRelatedSpine{2}]; hold on;
        plot(nanmean(nearcorr)*ones(1,max(hist(nearcorr))), 1:max(hist(nearcorr)), ':k')
        hist(nearcorr); set(findobj(gca, 'Facecolor', 'flat'), 'FaceColor', dred)
        ylabel('Count')
        xlabel('Correlation with nearest MRS')
        
        axes('Position', [pos(1)+0.7*pos(3), pos(2)+0.3*pos(4), 0.25*pos(3), 0.25*pos(4)]);
        neardist = [NearestHighlyCorrelatedMovementRelatedSpine{1}, NearestHighlyCorrelatedMovementRelatedSpine{2}];
        nextdist = [NextClosestHighlyCorrelatedMovementRelatedSpine{1}, NextClosestHighlyCorrelatedMovementRelatedSpine{2}];
        thirddist = [ThirdClosestHighlyCorrelatedMovementRelatedSpine{1}, ThirdClosestHighlyCorrelatedMovementRelatedSpine{2}];
        nd = hist(neardist, round(max(neardist)/5));
        nx = hist(nextdist, round(max(nextdist)/5));
        nt = hist(thirddist, round(max(thirddist)/5));
        if usenorm
            nd = nd/sum(nd);
            nx = nx/sum(nx);
            nt = nt/sum(nt);
        else
        end
        bar(nd, 'FaceColor', purple); hold on;
        bar(nx, 'FaceColor', orange)
        bar(nt, 'FaceColor', lgreen)
        if usenorm
            ylim([0 1])
        else
        end
        set(gca, 'XTick', [0:30], 'XTickLabel', mat2cell(num2str([0:5:150]')), 'Fontsize', 6)

        
    subplot(2,4,4); hold on;
        h = hist(disttouse{1}, round(max(disttouse{1}))/5);
        if usenorm
            h = h/sum(h);
        else
        end
        bar(h, 'FaceColor', dred)
        
        m = hist(disttouse{3}, round(max(disttouse{3}))/5);
        if usenorm
            m = m/sum(m);
        else
        end
        bar(m, 'FaceColor', lgreen)
        if usenorm
            ylim([-0.1 1])
        else
        end
        bar(m-h(1:length(m)), 'FaceColor', blue);
        legend({'All spine pairs', 'MR spines', 'Diff'})
        title('Session 1-2', 'Fontsize', 14)
        set(gca, 'XTick', [0:30], 'XTickLabel', mat2cell(num2str([0:5:150]')), 'Fontsize', 6)
        xlabel('Distance bins')
        ylabel('Fraction of Distances Measured', 'Fontsize', 14)
        
%         maxlength = max([length(n),length(h),length(m)]);
%         h(length(h)+1:maxlength) = 0;
%         m(length(m)+1:maxlength) = 0;
%         n(length(n)+1:maxlength) = 0;

        
    subplot(2,4,7); hold on;
        neardist = [NearestMovementRelatedSpine{10},NearestMovementRelatedSpine{11}];
        n = hist(neardist,round(max(neardist)/5));
        if usenorm
            n = n/sum(n);
        else
        end
        bar(n, 'FaceColor', lpurple); hold on;
        nextdist = [NextClosestMovementRelatedSpine{10}, NextClosestMovementRelatedSpine{11}];
        nx = hist(nextdist, round(max(nextdist)/5));
        if usenorm
            nx = nx/sum(nx);
        else
        end
        bar(nx, 'FaceColor', orange); 
        thirddist = [ThirdClosestMovementRelatedSpine{10}, ThirdClosestMovementRelatedSpine{11}];
        nt = hist(thirddist, round(max(thirddist)/5));
        if usenorm
            nt = nt/sum(nt);
        else
        end
        bar(nt, 'FaceColor', lgreen);
        fourthdist = [FourthClosestMovementRelatedSpine{10}, FourthClosestMovementRelatedSpine{11}];
        nf = hist(fourthdist, round(max(fourthdist)/5));
        if usenorm
            nf = nf/sum(nf);
            ylim([0 1])
        else
        end
        bar(nf, 'FaceColor', blue)
        legend({'Nearest MR spine', 'Next Closest', 'Third Closest', 'Fourth'}, 'Location', 'Northwest')
        title('Session 10-11', 'Fontsize', 14)
        set(gca, 'XTick', [0:30], 'XTickLabel', mat2cell(num2str([0:5:150]')), 'Fontsize', 6);
        xlabel('Distance bins', 'Fontsize', 14)
        ylabel('Fraction', 'Fontsize', 14)
        
        text(1.5, max(n), ['mean = ', num2str(nanmean(neardist))], 'Color', lpurple);
        text(2.5, max(nx), ['mean = ', num2str(nanmean(nextdist))], 'Color', orange);
        text(3.5, max(nt), ['mean = ', num2str(nanmean(thirddist))], 'Color', lgreen);
        text(4.5, max(nf), ['mean = ', num2str(nanmean(fourthdist))], 'Color', blue);
        
        pos = get(gca,'Position');
        axes('Position', [pos(1)+0.7*pos(3), pos(2)+0.7*pos(4), 0.25*pos(3), 0.25*pos(4)], 'Fontsize', 6);
        nearcorr = [CorrelationwithNearestMovementRelatedSpine{10}, CorrelationwithNearestMovementRelatedSpine{11}];
        hist(nearcorr); set(findobj(gca, 'Facecolor', 'flat'), 'FaceColor', dred); hold on;
        plot(nanmean(nearcorr)*ones(1,max(hist(nearcorr))), 1:max(hist(nearcorr)), ':k')
        ylabel('Count')
        xlabel('Correlation with nearest MRS')
        
        axes('Position', [pos(1)+0.7*pos(3), pos(2)+0.3*pos(4), 0.25*pos(3), 0.25*pos(4)]);
        neardist = [NearestHighlyCorrelatedMovementRelatedSpine{10}, NearestHighlyCorrelatedMovementRelatedSpine{11}];
        nextdist = [NextClosestHighlyCorrelatedMovementRelatedSpine{10}, NextClosestHighlyCorrelatedMovementRelatedSpine{11}];
        thirddist = [ThirdClosestHighlyCorrelatedMovementRelatedSpine{10}, ThirdClosestHighlyCorrelatedMovementRelatedSpine{11}];
        nd = hist(neardist, round(max(neardist)/5));
        nx = hist(nextdist, round(max(nextdist)/5));
        nt = hist(thirddist, round(max(thirddist)/5));
        if usenorm
            nd = nd/sum(nd);
            nx = nx/sum(nx);
            nt = nt/sum(nt);
        else
        end
        bar(nd, 'FaceColor', purple); hold on;
        bar(nd, 'FaceColor', orange)
        bar(nt, 'FaceColor', lgreen)
        if usenorm
            ylim([0 1])
        else
        end
        
        set(gca, 'XTick', [0:30], 'XTickLabel', mat2cell(num2str([0:5:150]')), 'Fontsize', 6)
        
        
    subplot(2,4,8); hold on;
        h = hist(disttouse{2}, round(max(disttouse{2}))/5);
        h = h/sum(h);
        bar(h, 'FaceColor', dred)
        
        m = hist(disttouse{4}, round(max(disttouse{4}))/5);
        if usenorm
            m = m/sum(m);
        else
        end
        bar(m, 'FaceColor', lgreen)
        if usenorm
            ylim([-0.1 1])
        else
        end
        bar(m-h(1:length(m)), 'FaceColor', blue)
        legend({'All spine pairs', 'MR spines', 'Diff'})
        title('Session 10-11', 'Fontsize', 14)
        set(gca, 'XTick', [0:30], 'XTickLabel', mat2cell(num2str([0:5:150]')), 'Fontsize', 6)
        xlabel('Distance bins')
        ylabel('Fraction of Distances Measured', 'Fontsize', 14)
        
%         maxlength = max([length(n),length(h),length(m)]);
%         h(length(h)+1:maxlength) = 0;
%         m(length(m)+1:maxlength) = 0;
%         n(length(n)+1:maxlength) = 0;
end



