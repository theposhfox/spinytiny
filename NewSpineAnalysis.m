function NewSpineAnalysis(varargin)

experimentnames = varargin;

if length(experimentnames) == 1
    experimentnames = varargin{1}; 
    %%%%%%%%%%%% Load Spine Dynamics Registry for a given animal

    if strcmpi(getenv('computername'), 'Nathan-Lab-PC')
        cd(['C:\Users\Komiyama\Desktop\Output Data', filesep, experimentnames, ' New Spine Analysis'])
    end

    fieldsource = fastdir(cd, 'Field');

    filecount = 1;
    for f = 1:length(fieldsource)
        load(fieldsource{f})
        fieldnumber = regexp(fieldsource{f}, '\d+.Spine');
        eval(['FieldData{', num2str(filecount), '} = SpineRegistry;']);
        clear SpineRegistry
        filecount = filecount+1;
    end

    NumFields = length(FieldData);

    for f = 1:NumFields
        FieldChanges{f} = diff(FieldData{f}.Data,1,2);
    end

    %%%%%%%%%%%% Load calcium imaging data for the animal

    if strcmpi(getenv('computername'), 'Nathan-Lab-PC')
        cd('C:\Users\Komiyama\Desktop\ActivitySummary_UsingRawData')
    end

    activitydata = fastdir(cd, [experimentnames, '.+_Summary']);

    for f = 1:length(activitydata)
        load(activitydata{f})
    end

    %%%%%%%%%%%% Match the loaded data with the session numbers from the spine
    %%%%%%%%%%%% registry data

    wrkspc = whos;
    for f = 1:NumFields
        for j = 1:length(FieldData{f}.DatesAcquired)
            locate =(regexp(who, FieldData{f}.DatesAcquired{j}));
            FieldData{f}.CalciumData{j} = eval(wrkspc(~cell2mat(cellfun(@isempty, locate, 'uni',false))).name);
        end
    end

    for f = 1:length(activitydata)
        clear(activitydata{f})
    end

    %%%%%%%%%%%% Separate the spine dynamics arrays into dendrites

    for f = 1:NumFields
        DendriteDynamics{f} = cellfun(@(x) FieldChanges{f}(x),FieldData{f}.CalciumData{1}.SpineDendriteGrouping,'uni', false);  %%% Calculate the CHANGE in spines (-1 is a lost spine, +1 is a new spine) between sessions for each dendrite
    end

    %%%%%%%%%%%% Load Statistical classification data

    if strcmpi(getenv('computername'), 'Nathan-Lab-PC')
        cd('C:\Users\Komiyama\Desktop\Output Data')
    end

    statdata = fastdir(cd, [experimentnames, '_StatClassified']);
    if ~isempty(statdata)
        load(statdata{1});
    else
        disp(['Cannot load stat data for animal', experimentnames]);
    end

    eval(['statclass = ', experimentnames, '_StatClassified;'])

    for f = 1:NumFields
        for s = 1:length(FieldData{f}.DatesAcquired)
            FieldData{f}.StatClass{s} = statclass{FieldData{f}.CalciumData{s}.Session};
        end
    end

    %%%%%%%%%%%%

    for f = 1:NumFields
        IsDendriteUsed{f} = sum([FieldData{f}.StatClass{1}.MovementDends, FieldData{f}.StatClass{2}.MovementDends],2);
        DendriteFunctionChange{f} = diff([FieldData{f}.StatClass{1}.MovementDends, FieldData{f}.StatClass{2}.MovementDends],1,2);
    end

    NumberofImagedDendrites = sum(cell2mat(cellfun(@length, DendriteDynamics, 'uni', false)));
    NumberofDynamicDendrites = 0;
    NumberofAdditionDendrites = 0;
    NumberofEliminationDendrites = 0;
    NumberofAdditionandEliminationDendrites = 0;
    NumberofStaticDendrites = 0;
    NumberofDynamicDendritesUsedForMovement = 0;
    NumberofAdditionDendritesUsedForMovement = 0;
    NumberofEliminationDendritesUsedForMovement = 0;
    NumberofStaticDendritesUsedForMovement = 0;
    for f = 1:NumFields
        for d = 1:length(DendriteDynamics{f})
            if sum(abs(DendriteDynamics{f}{d}))
                NumberofDynamicDendrites = NumberofDynamicDendrites+1;
                if ~isempty(find(DendriteDynamics{f}{d}>0,1))
                    NumberofAdditionDendrites = NumberofAdditionDendrites+1;
                end
                if ~isempty(find(DendriteDynamics{f}{d}<0,1))
                    NumberofEliminationDendrites = NumberofEliminationDendrites + 1;
                end
                if ~isempty(find(DendriteDynamics{f}{d}>0,1)) && ~isempty(find(DendriteDynamics{f}{d}<0,1))
                    NumberofAdditionandEliminationDendrites = NumberofAdditionandEliminationDendrites + 1;
                end
                if IsDendriteUsed{f}(d)
                    NumberofAdditionDendritesUsedForMovement = NumberofAdditionDendritesUsedForMovement+1;
                    NumberofEliminationDendritesUsedForMovement = NumberofEliminationDendritesUsedForMovement+1;
                    NumberofDynamicDendritesUsedForMovement = NumberofDynamicDendritesUsedForMovement+1;
                end
            elseif ~sum(abs(DendriteDynamics{f}{d}))
                NumberofStaticDendrites = NumberofStaticDendrites+1;
                if IsDendriteUsed{f}(d)
                    NumberofStaticDendritesUsedForMovement = NumberofStaticDendritesUsedForMovement+1;
                end
            end
        end
    end
    NumberofDendritesThatAreEverMovementRelated = sum(cell2mat(cellfun(@sum, IsDendriteUsed, 'uni', false)));
    FractionofDendritesThatAreEverMovementRelated = NumberofDendritesThatAreEverMovementRelated/NumberofImagedDendrites; 
    FractionofDendritesThatAreDynamic = NumberofDynamicDendrites/NumberofImagedDendrites;
    FractionofDendriteswithAddition = NumberofAdditionDendrites/NumberofImagedDendrites;
    FractionofDendriteswithElimination = NumberofEliminationDendrites/NumberofImagedDendrites;
    FractionofDynamicDendritesUsedForMovement = NumberofDynamicDendritesUsedForMovement/NumberofDynamicDendrites;
    FractionofAdditionDendritesUsedForMovement = NumberofAdditionDendritesUsedForMovement/NumberofAdditionDendrites;
    FractionofEliminationDendritesUsedForMovement = NumberofEliminationDendritesUsedForMovement/NumberofEliminationDendrites;
    FractionofStaticDendritesUsedForMovement = NumberofStaticDendritesUsedForMovement/NumberofStaticDendrites;

    a.SpineDynamics = FieldChanges;
    a.DendriteDynamics = DendriteDynamics;
    a.FractionofDendritesThatAreDynamic = FractionofDendritesThatAreDynamic;
    a.FractionofDendriteswithAddition = FractionofDendriteswithAddition;
    a.FractionofDendriteswithElimination = FractionofDendriteswithElimination; 
    a.NumberofDendritesThatAreEverMovementRelated = NumberofDendritesThatAreEverMovementRelated;
    a.FractionofDendritesThatAreEverMovementRelated = FractionofDendritesThatAreEverMovementRelated;
    a.NumberofImagedDendrites = NumberofImagedDendrites;
    a.NumberofDynamicDendrites = NumberofDynamicDendrites;
    a.NumberofAdditionDendrites = NumberofAdditionDendrites;
    a.NumberofEliminationDendrites = NumberofEliminationDendrites;
    a.NumberofAdditionandEliminationDendrites = NumberofAdditionandEliminationDendrites;
    a.NumberofStaticDendrites = NumberofStaticDendrites;
    a.IsDendriteEverMovementRelated = IsDendriteUsed;
    a.NumberofDynamicDendritesUsedForMovement = NumberofDynamicDendritesUsedForMovement;
    a.NumberofAdditionDendritesUsedForMovement = NumberofAdditionDendritesUsedForMovement;
    a.NumberofEliminationDendritesUsedForMovement = NumberofEliminationDendritesUsedForMovement;
    a.NumberofStaticDendritesUsedForMovement = NumberofStaticDendritesUsedForMovement;
    a.FractionofDynamicDendritesUsedForMovement = FractionofDynamicDendritesUsedForMovement;
    a.FractionofAdditionDendritesUsedForMovement = FractionofAdditionDendritesUsedForMovement;
    a.FractionofEliminationDendritesUsedForMovement = FractionofEliminationDendritesUsedForMovement;
    a.FractionofStaticDendritesUsedForMovement = FractionofStaticDendritesUsedForMovement;

    eval([experimentnames, '_SpineDynamicsSummary = a'])
    fname = [experimentnames, '_SpineDynamicsSummary'];
    save(fname, fname)
else
    if strcmpi(getenv('computername'), 'Nathan-Lab-PC')
        cd('C:\Users\Komiyama\Desktop\Output Data')
    end
    
    for i = 1:length(experimentnames)
        targetfile = [experimentnames{i}, '_SpineDynamicsSummary'];
        load(targetfile)
        eval(['currentdata = ',targetfile])
        SpineDynamics{i} = currentdata.SpineDynamics;
        DendriteDynamics{i} =  currentdata.DendriteDynamics;
        FractionofDendritesThatAreDynamic(1,i) = currentdata.FractionofDendritesThatAreDynamic;
        FractionofDendriteswithAddition(1,i) = currentdata.FractionofDendriteswithAddition;
        FractionofDendriteswithElimination(1,i) = currentdata.FractionofDendriteswithElimination;
        FractionofDendritesThatAreEverMovementRelated(1,i) = currentdata.FractionofDendritesThatAreEverMovementRelated;
        NumberofImagedDendrites(1,i) = currentdata.NumberofImagedDendrites;
        NumberofDynamicDendrites(1,i) = currentdata.NumberofDynamicDendrites;
        NumberofAdditionDendrites(1,i) = currentdata.NumberofAdditionDendrites;
        NumberofEliminationDendrites(1,i) = currentdata.NumberofEliminationDendrites;
        NumberofAdditionandEliminationDendrites(1,i) = currentdata.NumberofAdditionandEliminationDendrites;
        NumberofStaticDendrites(1,i) = currentdata.NumberofStaticDendrites;
        NumberofDendritesThatAreEverMovementRelated(1,i) = currentdata.NumberofDendritesThatAreEverMovementRelated;
        NumberofDynamicDendritesUsedForMovement(1,i) = currentdata.NumberofDynamicDendritesUsedForMovement;
        NumberofAdditionDendritesUsedForMovement(1,i) = currentdata.NumberofAdditionDendritesUsedForMovement;
        NumberofEliminationDendritesUsedForMovement(1,i) = currentdata.NumberofEliminationDendritesUsedForMovement;
        NumberofStaticDendritesUsedForMovement(1,i) = currentdata.NumberofStaticDendritesUsedForMovement;
        FractionofDynamicDendritesUsedForMovement(1,i) = currentdata.FractionofDynamicDendritesUsedForMovement;
        NumberofAdditionDendritesUsedForMovement(1,i) = currentdata.FractionofAdditionDendritesUsedForMovement;
        NumberofEliminationDendritesUsedForMovement(1,i) = currentdata.FractionofEliminationDendritesUsedForMovement;
        NumberofStaticDendritesUsedForMovement(1,i) = currentdata.FractionofStaticDendritesUsedForMovement;
        clear currentdata
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% Color Information %%
    %%%%%%%%%%%%%%%%%%%%%%%%

    lgray = [0.50 0.51 0.52];       brown = [0.28 0.22 0.14];
    gray = [0.50 0.51 0.52];        lbrown = [0.59 0.45 0.28];
    yellow = [1.00 0.76 0.05];      orange = [0.95 0.40 0.13];
    lgreen = [0.45 0.8 0.35];       green = [0.00 0.43 0.23];
    lblue = [0.30 0.65 0.94];       blue = [0.00 0.33 0.65];
    magenta = [0.93 0.22 0.55];     purple = [0.57 0.15 0.56];
    pink = [0.9 0.6 0.6];           lpurple  = [0.7 0.15 1];
    red = [0.85 0.11 0.14];         black = [0.1 0.1 0.15];
    dred = [0.6 0 0];               dorange = [0.8 0.3 0.03];
    bgreen = [0 0.6 0.7];
    colorj = {red,lblue,green,lgreen,gray,brown,yellow,blue,purple,lpurple,magenta,pink,orange,brown,lbrown};
    rnbo = {dred, red, dorange, orange, yellow, lgreen, green, bgreen, blue, lblue, purple, magenta, lpurple, pink}; 

    %%%%%%%%%%%%%%%%%%%%%%%

    scrsz = get(0, 'ScreenSize');

    figure; 
%     allmat = [nanmean(FractionofDendritesThatAreDynamic); nanmean(FractionofDendriteswithAddition); nanmean(FractionofDendriteswithElimination)];
%     allerror = [nanstd(FractionofDendritesThatAreDynamic)/sqrt(length(FractionofDendritesThatAreDynamic)); nanstd(FractionofDendriteswithAddition)/sqrt(length(FractionofDendriteswithAddition)); nanstd(FractionofDendriteswithElimination)/sqrt(length(FractionofDendriteswithElimination))];
    allmat = [NumberofDynamicDendrites/NumberofImagedDendrites; NumberofAdditionDendrites/NumberofImagedDendrites; NumberofEliminationDendrites/NumberofImagedDendrites; NumberofAdditionandEliminationDendrites/NumberofImagedDendrites];
    bar(allmat, 'FaceColor', lgreen)
%     r_errorbar(1:3, allmat, allerror, 'k')
    ylabel({'Fraction of Dendrites'; 'with Dynamic Spines'}, 'Fontsize', 12)
    set(gca, 'XTick', 1:4, 'XTickLabel',{'AorE', 'A', 'E', 'A&E'})
    title('Prevalence of Spine Dynamics on Imaged Dendrites')
    ylim([0 1])


    figure;
%     allmat = [nanmean(FractionofDendritesThatAreEverMovementRelated), nanmean(FractionofDynamicDendritesUsedForMovement),nanmean(FractionofAdditionDendritesUsedForMovement),nanmean(FractionofEliminationDendritesUsedForMovement),nanmean(FractionofStaticDendritesUsedForMovement)];
%     allerror = [nanstd(FractionofDendritesThatAreEverMovementRelated)/sqrt(length(FractionofDendritesThatAreEverMovementRelated)); nanstd(FractionofDynamicDendritesUsedForMovement)/sqrt(length(FractionofDynamicDendritesUsedForMovement)); nanstd(FractionofAdditionDendritesUsedForMovement)/sqrt(length(FractionofAdditionDendritesUsedForMovement));nanstd(FractionofEliminationDendritesUsedForMovement)/sqrt(length(FractionofEliminationDendritesUsedForMovement));nanstd(FractionofStaticDendritesUsedForMovement)/sqrt(length(FractionofStaticDendritesUsedForMovement))];
    allmat = [nansum(NumberofDendritesThatAreEverMovementRelated)/nansum(NumberofImagedDendrites), nansum(NumberofDynamicDendritesUsedForMovement)/nansum(NumberofImagedDendrites), nansum(NumberofAdditionDendritesUsedForMovement)/nansum(NumberofImagedDendrites), nansum(NumberofEliminationDendritesUsedForMovement)/nansum(NumberofImagedDendrites), nansum(NumberofStaticDendritesUsedForMovement)/nansum(NumberofImagedDendrites)];
    bar(allmat, 'FaceColor', blue)
%     r_errorbar(1:5, allmat, allerror, 'k')
    ylabel({'Fraction of Dendrites'; 'That Are Movement Related'}, 'Fontsize', 12)
    set(gca, 'XTick', 1:5, 'XTickLabel',{'All Dendrites','A&E','A','E','St'})
    title('Likelihood of Movement Relatedness')
    text(1,0.9, ['Total Dendrites = ', num2str(nansum(NumberofImagedDendrites))])
    ylim([0 1])

end
end