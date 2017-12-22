function a = AdjustFrequency_SynapseOnly(File, currentsession, showFig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Color Information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Find the file being called %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isstruct(File)
    user = File.Filename(1:2);
    folder = regexp(File.Filename, [user, '000?\d+'], 'match');
    folder = folder{1};
    Date = regexp(File.Filename, '_\d+_');
    Date = Date{1};
    cd(['Z:\People\Nathan\Data\', folder, '\16', Date(2:end-1), '\summed'])
else
    user = regexp(File, '[ABCDEFGHIJKLMNOPQRSTUVWXYZ]{2}', 'match');
    user = user{1};
%     folder = regexp(File, [user, '00\d+_'], 'match');
    folder = regexp(File, [user, '\d+[^_]'], 'match');
    folder = folder{1};
    Date = regexp(File, '\d{6}', 'match');
    Date = Date{1};
%     Date = '1713';
    try
        if strcmpi(user, 'NH')
            cd(['Z:\People\Nathan\Data\', folder, '\', Date, '\summed'])
            load([folder, '_', Date(3:end), '_001_001_summed_50_Analyzed'])
        elseif strcmpi(user, 'XR')
            cd(['Z:\People\Nathan\Xiangyu\', folder, '\', Date, '\summed'])
            load([folder, '_', Date(3:end), '_001_001_summed_50_Analyzed'])
        elseif strcmpi(user, 'SC')
            cd(['Z:\People\Simon\Data\Simon\', folder, '\', Date, '\summed'])
            load([folder, '_', Date, '_001_001_summed_50_Analyzed'])
        end
    catch      %%% File naming gets wonked up sometimes; change whatever you need to make the program read the file
        if strcmpi(user, 'NH') 
%             folder = 'NH005';
%             Date = '160713'
            cd(['Z:\People\Nathan\Data','\', folder, '\', Date, '\summed'])
%             folder = 'NH0005';
%             Date = '160330';
            load([folder, '_', Date(3:end), '_001_001_summed_50_Analyzed'])
        elseif strcmpi(user, 'XR')                                              
            cd(['Z:\People\Nathan\Xiangyu\', folder, '\', Date, '\summed'])
            load(['XR0002_', Date(3:end), '_001_001_summed_50_Analyzed'])
        elseif strcmpi(user, 'SC')
            cd(['Z:\People\Simon\Data\Simon\', folder, '\', Date, '\summed'])
            Date = Date(3:end);
            load([folder, '_', Date, '_001_001_summed_50_Analyzed'])
        end
    end
    try
        eval(['File =' folder, '_', Date, '_001_001_summed_50_Analyzed;'])
    catch
        temp = who(['*', user, '*']);
        eval(['File =', temp{1}, ';']);
    end
end

filename = regexp(File.Filename, '.tif', 'split');
filename = filename{1};
File.Filename = [folder, '_', Date(3:end), '_001_001_summed_50_Analyzed'];

a = File;
Scrsz = get(0, 'Screensize');

%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%% Controlled variables
%%%%%%%%%%%%%%%%%%%%%%%%

a.UsePreviousPreferences = 0;

% vardir = '/usr/local/lab/People/Nathan/Data/ActivitySummary';
vardir = 'C:\Users\Komiyama\Desktop\ActivitySummary';

if a.UsePreviousPreferences
    cd(vardir); %%% Variable directory
    load([folder, '_', Date, '_Summary']);
    eval(['SummaryFile = ', folder, '_', Date, '_Summary;'])
    spinethreshmultiplier = SummaryFile.spinethresholdmultiplier;
    spinevalueslimit = SummaryFile.spinevalueslimit;
    spinebaselinesmoothwindow = SummaryFile.spinebaselinesmoothwindow;
    spinesmoothwindow = SummaryFile.spinesmoothwindow;
    Dendthreshmultiplier = SummaryFile.Dendthreshmultiplier;
    Dendvalueslimit = SummaryFile.Dendvalueslimit;
    dendbaselinesmoothwindow = SummaryFile.dendbaselinesmoothwindow;
    dendsmoothwindow = SummaryFile.dendsmoothwindow;  
%     ClusterThresh = SummaryFile.ClusterThresh;
    ClusterThresh = 0.6;
    SpectralLengthConstant = SummaryFile.SpectralLengthConstant;
%     SpectralLengthConstant = 10;
    currentsession = SummaryFile.Session;
else
    cd(vardir); %%% Variable directory
    load([folder, '_', Date, '_Summary']);
    eval(['SummaryFile = ', folder, '_', Date, '_Summary;'])

    spinethreshmultiplier = 0.5;
    spinevalueslimit = 2;
    spinebaselinesmoothwindow = 200;
    spinesmoothwindow = 30;
    Dendthreshmultiplier = .3;
    Dendvalueslimit = 1;
    dendbaselinesmoothwindow = 200;
    dendsmoothwindow = 30;

    ClusterThresh = 0.5;
    SpectralLengthConstant = 10;
    
    currentsession = SummaryFile.Session;
    
end

a.spinethresholdmultiplier = spinethreshmultiplier;
a.spinevalueslimit = spinevalueslimit; 
a.spinebaselinesmoothwindow = spinebaselinesmoothwindow;
a.spinesmoothwindow = spinesmoothwindow;
a.Dendthreshmultiplier = Dendthreshmultiplier;
a.Dendvalueslimit = Dendvalueslimit;
a.dendbaselinesmoothwindow = dendbaselinesmoothwindow;
a.dendsmoothwindow = dendsmoothwindow;
a.ClusterThresh = ClusterThresh;
a.SpectralLengthConstant = SpectralLengthConstant;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Select which spine to plot %%%

% SpineNo = randi(File.NumberofSpines,1);
SpineNo = 1;


DendNum = File.NumberofDendrites;
% DendriteChoice = randi(DendNum,1);
DendriteChoice = 4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Select which spine data to use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Options:    1) File.Fluorescence_Measurement %%% Raw data
%%%             2) File.deltaF                   %%% Baseline-subtracted
%%%             3) File.dF_over_F                %%% Baseline-sub and div.
%%%             4) File.SynapticEvents           %%% All above + dend-subtract

spinetraceoption = 4;

if spinetraceoption == 1
    spinedatatouse = File.Fluorescence_Measurement;
    correspondingnewdata = a.Fluorescence_Measurement; %%% The "new" data is sometimes changed in parallel, and this should always be accounted for
elseif spinetraceoption == 2
    spinedatatouse = File.deltaF;
    correspondingnewdata = a.deltaF;
elseif spinetraceoption == 3
    spinedatatouse = File.dF_over_F;
    correspondingnewdata = a.dF_over_F;
elseif spinetraceoption ==4
    spinedatatouse = File.SynapticEvents;
    correspondingnewdata = a.SynapticEvents;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Describe the basic shape of each calcium trace


%%%%%%%%%%%%%%%%
% spinethreshmultiplier = 0.5;
%%%%%%%%%%%%%%%%

rawspread = zeros(length(File.deltaF),1);
ave = zeros(length(File.deltaF),length(File.deltaF{1}));
blsub = zeros(length(File.deltaF),length(File.deltaF{1}));
smoo = zeros(length(File.deltaF),length(File.deltaF{1}));
all = zeros(length(File.deltaF),length(File.deltaF{1}));
med = zeros(length(File.deltaF),1);
spread = zeros(length(File.deltaF),1);
binarized = zeros(length(File.deltaF),length(File.deltaF{1}));
amp = zeros(1,length(File.deltaF));
freq = zeros(1,length(File.deltaF),1);

for i = 1:length(File.deltaF)
    if any(isnan(spinedatatouse{i}))
        try
            spinedatatouse{i}(isnan(spinedatatouse{i})) = nanmean([spinedatatouse(find(isnan(spinedatatouse{i}))-1),spinedatatouse{i}(find(isnan(spinedatatouse))+1)]);
        catch
            spinedatatouse{i}(isnan(spinedatatouse{i})) = 0;
        end
    end
    
    spinedatatouse{i}(1:10) = nanmedian(spinedatatouse{i}(1:1000));
    correspondingnewdata{i}(1:10) = nanmedian(spinedatatouse{i}(1:1000));
    spinedatatouse{i}(end-10:end) = nanmedian(spinedatatouse{i}(end-1000:end));
    correspondingnewdata{i}(end-10:end) = nanmedian(spinedatatouse{i}(end-1000:end));
    raw = spinedatatouse{i};

    rawmed = nanmedian(raw);
    rawspread(i,1) = nanstd(raw);
    raw(raw>rawmed+spinevalueslimit*rawspread(i,1)) = rawmed; %%% Cap off large and small values to pinch the data towards the true baseline
    raw(raw<rawmed-spinevalueslimit*rawspread(i,1)) = rawmed; %%%
    ave(i,:) = smooth(raw,spinebaselinesmoothwindow);         %%% Baseline value
    blsub(i,:) = spinedatatouse{i}-ave(i,:);             %%% Baseline-subtracted value
    smoothed = smooth(blsub(i,:),spinesmoothwindow);
    smoothed(smoothed<(median(smoothed)-nanstd(smoothed))) = median(smoothed)-nanstd(smoothed);   %%% This is the dendrite-subtracted data, so it's possible that there are very large negative events, which are artificial, and should be reduced to ~the level of noise
    smoo(i,:) = smoothed;
    all(i,:) = smoo(i,:);
    med(i,1) = nanmedian(smoo(i,:));
    spread(i,1) = rawspread(i,1);
end

% a.EventNumber = freq;
% a.ActivityMap = binarized;
% a.MeanEventAmp = amp;

numberofSpines = size(binarized,1);
DendNum = File.NumberofDendrites;

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

%%% Figure 1 %%%

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

showFig = 1;

Scrsz = get(0, 'Screensize');
if showFig == 1
    k = zeros(1,length(File.Time));
    k(1:end) = med(SpineNo,1)+spinethreshmultiplier*spread(SpineNo,1);
    m = ones(1,length(File.Time))*med(SpineNo,1);
    fig1 = figure('Position', [10, Scrsz(4)/2.5,Scrsz(3)/2,Scrsz(4)/2]); 
    rawplot = subplot(2,2,1:2);
    plot(File.Time, spinedatatouse{SpineNo}, 'Color', [0.2 0.2 0.2]); hold on;
    axpos = get(rawplot, 'Position');
    plot(File.Time, ave(SpineNo, :), 'Color', red, 'Linewidth', 3)
    ylabel(['Dend-subtracted dF/F for spine no. ', num2str(SpineNo)])
    title('Example of event detection')
    legend({'Raw Data', 'Baseline'}, 'Location', 'SouthEastOutside');
else
end

%%% Change all values below peak detection to zero
floored = [];
riderthresh = 1.5;
riderthresh2 = 2;

for i = 1:numberofSpines
    temp = all(i,:);
    temp(temp<(med(i,1)+(spinethreshmultiplier*spread(i,1)))) = 0;
    floored(i,:) = temp;
    temp(temp<(med(i,1)+(riderthresh*spinethreshmultiplier*spread(i,1)))) = nan;
%     riders(i,:) = temp;
%     temp(temp<(med(i,1)+(riderthresh2*spinethreshmultiplier*spread(i,1)))) = nan;
    tamp = temp;
    tamp(isnan(tamp)) = 0;
    tamp = smooth(tamp,20);
    dtamp = diff(tamp);     %%% first derivative of the binarized data
    dtamp = [0;dtamp];
    dtamp(dtamp>0) = 1; dtamp(dtamp<0) = -1;
%     dtamp = smooth(dtamp,20);
    d2tamp = diff(dtamp);
    d2tamp = [0;d2tamp];    %%% Second derivative of the binarized data (concavity)
    d2tamp(d2tamp>0) = 1; d2tamp(d2tamp<0) = -1;
%     d2tamp = smooth(d2tamp, 20);
    temp(d2tamp>0) = nan; %% For plateau spikes, when the 2nd derivative is positive (concave up, corresponding to dips), punch a 'hole' in the data, so that multiple peaks will be counted
    riderstop(i,:) = temp;
end

if showFig == 1
    procplot = subplot(2,2,3:4);
    hold on; plot(File.Time, smoo(SpineNo,:), 'Color',[0.2 0.2 0.2], 'LineWidth', 1);
    plot(File.Time, k, '--', 'Color', lgreen, 'LineWidth', 2)
    plot(File.Time, m, '--', 'Color', purple)
    abovethresh = floored(SpineNo,:);
    abovethresh(abovethresh == 0) = nan;
%     plot(1:length(File.Time),abovethresh, 'Color',orange, 'LineWidth',2);
%     plateauspikes = riders(SpineNo, :);
%     plot(1:length(File.Time), plateauspikes, 'Color',blue, 'LineWidth', 2);
    topspikes = riderstop(SpineNo,:);
    plot(1:length(File.Time), topspikes, 'Color', dred, 'LineWidth', 2);
    xlabel('Frames')
    ylabel(['Smoothed dF/F for spine no. ', num2str(SpineNo)])
    % plot(File.Time(peak_loc{SpineNo}), spine_peaks{SpineNo}+0.1, 'kv', 'markerfacecolor', 'g');
else
end

%%% Change all events to one, making square pulses corresponding to
%%% activity
square = [];
% riders(~isnan(riders)) = riderthresh-1;
% riders(isnan(riders)) = 0;
riderstop(~isnan(riderstop)) = riderthresh2-1.5;
riderstop(isnan(riderstop)) = 0;

for i = 1:numberofSpines
    temp = floored(i,:);   %%% This value will eventually be used to define "synapse only" events, which only requires knowledge of when spines are above a threshold (e.g. spikes riding on top of activity need not be considered)
    temp(temp~=0)= 1;
    square(i,:) = temp;
    temp = [];
%     temp = square(i,:)+riders(i,:)+riderstop(i,:); %% Can remove 'riderstop to get rid of plateau spike summing
    temp = square(i,:)+riderstop(i,:); %% Can remove 'riderstop to get rid of plateau spike summing
    both(i,:) = temp;
    temp2 = (diff(temp)>0.1)>0;
    temp3 = [0, temp2];          %%% Any use of 'diff' shortens the vector by 1
    smeared = smooth(temp3, 5); %%% Smoothing factor is taken from the reported decay constant of GCaMP6f (~150ms), converted to frames 
    smeared(smeared>0) = 1;
    trueeventcount(i,:) = smeared;
end

if showFig == 1
    plot(File.Time, both(SpineNo,:)*(med(SpineNo,1)+spinethreshmultiplier*spread(SpineNo,1)), 'Color',yellow, 'LineWidth', 2)
    plot(File.Time, trueeventcount(SpineNo,:)*(med(SpineNo,1)+spinethreshmultiplier*spread(SpineNo,1)),'Color', lblue, 'LineWidth', 2);
    legend({'Smoothed Data', 'Threshold', 'Baseline', 'Events above thresh', 'Ternarized', 'Counted Events'}, 'Location', 'SouthEastOutside')
    axpos2 = get(procplot, 'Position');
    set(rawplot, 'Position', [axpos2(1), axpos(2), axpos2(3), axpos(4)])
    set(procplot, 'Box', 'on')
    drawnow;
else
end


for i = 1:numberofSpines
    frequency(i,1) = (nnz(diff(trueeventcount(i,:)>0.5)>0)/((length(File.Time)/30.49)/60))';
end


for i = 1:size(floored,1)
    [peaks, loc] = findpeaks(smooth(floored(i,:),10), 'MinPeakDistance', 5);   %%% The "floored" variable is used to count events, and so should be used to find the corresponding amplitude
    spine_peaks{i} = peaks;
    peak_loc{i} = loc;
    bin = zeros(1,length(spinedatatouse{i}));
    bin(loc) = 1;
    binarized(i,:) = bin;
    amp(1,i) = mean(peaks);
    freq(1,i) = length(peaks);
end

if ~isfield(File, 'SpineDendriteGrouping')
    if File.NumberofDendrites == 1
        File.SpineDendriteGrouping{1} = 1:numberofSpines;
    else
        disp('SpineDendriteGrouping field is absent!!!')
%        File.SpineDendriteGrouping{1} = 1:4;
 %       File.SpineDendriteGrouping{2} = 5:13;
      %  File.SpineDendriteGrouping{3} = 14:16;
       % File.SpineDendriteGrouping{4} = 17:18;
%         File.SpineDendriteGrouping{5} = 17;
%         File.SpineDendriteGrouping{6} = 18:19;
%         File.SpineDendriteGrouping{7}= 20:22;
    end
end

a.ActivityMap = trueeventcount;
a.MeanEventAmp = amp;
a.AllSpineAmpData = spine_peaks;

a.EventNumber = frequency;
a.Frequency = frequency;
newamp = amp';

%%%%%%%%%%%%%%%%
%%%%%%%%%%%
% Dendthreshmultiplier = 0.5;
%%%%%%%%%%%
%%%%%%%%%%%%%%%%

Drawspread = zeros(DendNum,1);
Dave = zeros(DendNum,length(File.Dendrite_dFoF(1,:)));
Dblsub = zeros(DendNum,length(File.Dendrite_dFoF(1,:)));
Dsmoo = zeros(DendNum,length(File.Dendrite_dFoF(1,:)));
Dall = zeros(DendNum,length(File.Dendrite_dFoF(1,:)));
Dmed = zeros(DendNum,1);
Dspread = zeros(DendNum,1);
Dbinarized = zeros(DendNum,length(File.Dendrite_dFoF(1,:)));
Damp = zeros(DendNum,1);
Dfreq = zeros(DendNum,1);

for i = 1:DendNum
    if any(isnan(File.Dendrite_dFoF(i,:)))
        temp = File.Dendrite_dFoF(i,:);
        try
            temp(isnan(temp)) = nanmean([temp(find(isnan(temp))-1),temp(find(isnan(temp))+1)]);
        catch
            temp(isnan(temp)) = 0;
        end
        File.Dendrite_dFoF(i,:) = temp;
    end
    File.Dendrite_dFoF(i,1:10) = 0;
    a.Dendrite_dFoF(i,1:10) = 0;
    File.Dendrite_dFoF(i,end-10:end) = 0;
    a.Dendrite_dFoF(i,end-10:end) = 0;
    Draw = File.Dendrite_dFoF(i,:);
    Drawmed = nanmedian(Draw);
    Drawspread = nanstd(Draw);
    Draw(Draw>Drawmed+Dendvalueslimit*Drawspread) = Drawmed;
    Draw(Draw<Drawmed-Dendvalueslimit*Drawspread) = Drawmed;
    Dave(i,:) = smooth(Draw,dendbaselinesmoothwindow); 
    blDend(i,:) = File.Dendrite_dFoF(i,:)-Dave(i,:);
    Dsmoo(i,:) = smooth(blDend(i,:),dendsmoothwindow);
    Dmed(i,1) = nanmedian(Dsmoo(i,:));
    Dspread(i,1) = nanstd(Dsmoo(i,:));
    Dtemp = Dsmoo(i,:);
%     Dtemp(Dtemp>Dmed(i,1)+Dspread(i,1)) = Dmed(i,1)+Dspread(i,1);
    Dmed(i,1) = nanmedian(Dtemp);
    Dspread(i,1) = Drawspread;
%     Dtemp = Dsmoo(i,:);
%     Dtemp(Dtemp>Dmed(i,1)+Dspread(i,1)) = Dmed(i,1)+Dspread(i,1);
%     Dtemp(Dtemp<Dmed(i,1)-Dspread(i,1)) = Dmed(i,1)-Dspread(i,1);
%     Dave(i,:) = smooth(Dtemp,1000);
%     Dsmoo(i,:) = Dsmoo(i,:)-Dave(i,:);
    [Dpeaks, Dloc] = findpeaks(Dsmoo(i,:), 'MinPeakHeight', Dmed(i,1)+(Dendthreshmultiplier*Dspread(i,1)), 'MinPeakDistance', 20, 'Threshold', 0);
    Dend_Peaks{i} = Dpeaks;
    Dend_Locations{i} = Dloc;
    Damp(i,1) = mean(Dpeaks);
    Dfreq(i,1) = length(Dpeaks)/((length(a.Time)/30.49)/60);
end


for i = 1:File.NumberofDendrites
    temp = Dsmoo(i,:);
    temp(temp<(Dmed(i,1)+(Dendthreshmultiplier*Dspread(i,1)))) = 0;
    floored_Dend(i,:) = temp;
    temp(temp<(Dmed(i,1)+(Dendthreshmultiplier*Dspread(i,1)))) = nan;
%     Driders(i,:) = temp;
%     temp(temp<(med(i,1)+(riderthresh2*Dendthreshmultiplier*spread(i,1)))) = nan;
    tamp = temp;
    tamp(isnan(tamp)) = 0;
%     tamp = smooth(tamp,5);
    dtamp = diff(tamp);
%     dtamp = [0;dtamp];
%     dtamp = smooth(dtamp,10);
    dtamp(dtamp>0) = 1; dtamp(dtamp<0) = -1;
    d2tamp = diff(dtamp);
%     d2tamp = [0;d2tamp];
%     d2tamp = smooth(d2tamp, 10);
    d2tamp(d2tamp>0) = 1; d2tamp(d2tamp<0) = -1;
    temp(d2tamp>0) = nan; %% For plateau spikes, when the 1st derivative is negative (val. decreasing) and 2nd derivative is positive (concave up, corresponding to dips), punch a 'hole' in the data, so that multiple peaks will be counted
    Driders(i,:) = temp;
end

Driders(~isnan(Driders)) = riderthresh-1;
Driders(isnan(Driders)) = 0;
% Driderstop(~isnan(Driderstop)) = riderthresh2-1.5;
% Driderstop(isnan(Driderstop)) = 0;

for i = 1:File.NumberofDendrites
    temp = floored_Dend(i,:);
    temp(temp~=0) = 1;
    square_Dend(i,:) = temp;
    temp = [];
    temp = square_Dend(i,:)+Driders(i,:);
    Dboth(i,:) = temp;
    temp2 = (diff(temp)>0.1)>0;
    temp3 = [0, temp2];          %%% Any use of 'diff' shortens the vector by 1; correct for this
    Dsmeared = smooth(temp3, 5); %%% Smoothing factor is taken from the reported decay constant of GCaMP6f (~150ms), converted to frames 
    Dsmeared(Dsmeared>0) = 1;
    Dendtrueeventcount(i,:) = Dsmeared;
end


%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

%%% Figure 2 %%%

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%


cofires = sum(square);
cofires(cofires == 1) = 0; %%% values of '1' correspond to a single event, meaning it's not actually a 'cofire' event

if showFig == 1
    figure('Position', [Scrsz(3)/2, Scrsz(4)/2.5,Scrsz(3)/2,Scrsz(4)/2]); hold on;

    subplot(3,3, 1:6); hold on; 
    for i = 1:numberofSpines
        for j = 1:size(square, 2)
            if square(i,j) == 1
                line([j j],[i-0.5:i+0.5],'color', [0.2 0.2 0.2]);
            end
            if trueeventcount(i,j) == 2
                line([j j], [i,i+1], 'color', red);
            end
        end
    end

    title('Activity Map');
    ylabel('Spine No.');
    xlabel('Frames');
    xlim([0 size(binarized,2)+1]);
    ylim([0 numberofSpines+1]);
    set(gca, 'YTick', [1:numberofSpines]);
    subplot(3,3,7:9); plot(cofires, 'b')
    xlim([0 size(binarized,2)+1]);
else
end

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

%%% Figure 3 %%%

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%


if showFig == 1
    figure('Position', [ 10, 50 ,Scrsz(3)/2,Scrsz(4)/2]); hold on;
else
end

causal = zeros(numberofSpines,size(all,2));
back = 1;   %%% Number of frames prior to a dendritic event that a spine must be active to be considered "causal"

for i = 1:numberofSpines
    onDend = 1;
    for j = 1:DendNum
        if ~isempty(find(File.SpineDendriteGrouping{j} == i))
            onDend = j;
            synapseOnly(i,:) = square(i,:)-square_Dend(onDend,:);
            synapseOnly(i,synapseOnly(i,:)<1) = 0;
            synapseOnlyFreq(i,1) = (nnz(diff(synapseOnly(i,:)>0.5)>0)/((length(a.Time)/30.49)/60))';
            withAPs(i,:) = square(i,:)+square_Dend(onDend,:);   %%% Add binarized spine data to dendrite data to illustrate when dendrite and spines are co-firing
            withAPsFreq(i,1) = (nnz(diff(withAPs(i,:)>1.5)>0)/((length(a.Time)/30.49)/60))';
            withAPs(i,(withAPs(i,:)==1)) = 0;
            withAPs(i,(withAPs(i,:)==2)) = 1;
%             dendOnly(i,:) = square_Dend(onDend,:)-square(i,:);
            APstart(i,:) = [0,diff(withAPs(i,:))]; %%% Find where the derivative of the AP-paired events == 2 (i.e. where it is first increasing, making it the start of an AP)
            for k = 1:size(binarized,2)
                if showFig == 1
                    if synapseOnly(i,k) == 1 && withAPs(i,k) ~=1
                        line([k k],[i-0.5:i+0.5],'color', [0.2 0.2 0.2]);
                    end
                    if withAPs(i,k) == 1
                        line([k k],[i-0.5:i+0.5],'color', red);
                    end
    %                 if dendOnly(i,k) == 1
    %                     line([k k],[i:i+1],'color', 'r'); hold on;
    %                 end
                    if k > back %%% can't index at zero (see below)
                        if APstart(i,k) == 1 && synapseOnly(i,k-back) == 1 %%% If a spine event precedes a putative AP
                            causal(i,(find(diff(synapseOnly(i,1:k))==1, 1,'last'))+1:k) = 1;
                            line([k-back k-back], i:0.5:(i+0.5), 'color', lblue, 'LineWidth', 2); hold on;
                            line([k-1 k-1], i:0.5:(i+0.5), 'color', lblue, 'LineWidth', 2); hold on;
                        end
                    end
                else
                    if k > back %%% can't index at zero (see below)
                        if APstart(i,k) == 1 && synapseOnly(i,k-2) == 1 %%% If a spine event precedes a putative AP
                            causal(i,(find(diff(synapseOnly(i,1:k))==1, 1,'last'))+1:k) = 1;
                        end
                    end
                end
            end 
        else
        end
    end
end

drawnow;

% withAPs(withAPs == 1) = 0; 
% withAPs(withAPs ==2 ) = 1;

for i = 1:DendNum
    Dfrequency(i,1) = (nnz(diff(Dendtrueeventcount(i,:)>0.5)>0)/((length(a.Time)/30.49)/60))';
end

a.Dendritic_Frequency = Dfrequency;
a.Dendritic_Amp = Damp;
a.SynapseOnlyFreq = synapseOnlyFreq;
a.SpikeTimedEvents = withAPsFreq;
a.CausalBinarized = causal;
%%%%%%%%%%%%

if showFig == 1
    xlabel('Frames');
    ylabel('Spine number');
    xlim([0 size(binarized,2)+1]);
    ylim([0 numberofSpines+1]);
    title('Synapse only vs. AP-paired events');
else
end
%%%%%%%%%%%%%%%%

Spine_No = numberofSpines;
Time = size(binarized,2)/30.49;
Time = Time/60;
cofire_index = sum(cofires)/(Spine_No*Time);

cofires(cofires ==0) = NaN;
coactive_percentage = 100*(max(cofires))/numberofSpines
spec = 'Approximately %2.0f percent of spines showed co-active firing\n';
fprintf(spec, coactive_percentage);


%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

%%% Figure 4 %%%

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%


DendriteChoice = randi(DendNum,1);

if showFig == 1
    k = zeros(1,length(File.Time));
    k(1:end) = Dmed(DendriteChoice,1)+Dendthreshmultiplier*Dspread(DendriteChoice,1);
    m = ones(1,length(File.Time)) * Dmed(DendriteChoice,1);

    figure('Position', [Scrsz(3)/2, 50 ,Scrsz(3)/2,Scrsz(4)/2]); 
    subplot(2,2,1:2)
    plot(File.Time, File.Dendrite_dFoF(DendriteChoice,:), 'Color', [0.2 0.2 0.2]);
    hold on; 
    plot(File.Time, Dave(DendriteChoice,:), 'r', 'LineWidth', 2);
    
    xlabel('Frames')
    ylabel(['Events for Dendrite ', num2str(DendNum)])
    
    subplot(2,2,3:4)
    plot(File.Time, Dsmoo(DendriteChoice,:), 'Color', blue, 'LineWidth', 2.5); hold on;
%     plot(File.Time(Dend_Locations{DendriteChoice}), Dend_Peaks{DendriteChoice}+0.05, 'kv', 'markerfacecolor', lgreen);
    plot(File.Time, Dboth(DendriteChoice,:)*(Dmed(DendriteChoice,1)+Dendthreshmultiplier*Dspread(DendriteChoice,1)), 'Color', red, 'Linewidth', 2)
    plot(File.Time, Dendtrueeventcount(DendriteChoice,:)*(Dmed(DendriteChoice,1)+Dendthreshmultiplier*Dspread(DendriteChoice,1)),'Color', orange, 'LineWidth', 2);
    plot(File.Time, k, '--c')
    plot(File.Time, m, '--m')

    xlabel('Frames')
    ylabel(['Events for Dendrite ', num2str(DendNum)])
    title('Example of event detection')
    drawnow
else
end

Dboth(Dboth>1) = 1;

a.Dendrite_Binarized = Dboth;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 5 : Spatial Analysis %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pixpermicron = 6.6667;
SpineToSpineDistance = nan(numberofSpines,numberofSpines);
OverallCorrelation = nan(numberofSpines, numberofSpines);
SpineToSpineCorrelation = nan(numberofSpines,numberofSpines);
SpineToSpinePValue = nan(numberofSpines, numberofSpines);
SpwAPCorrelation = nan(numberofSpines, numberofSpines);
SpwAP_PValue = nan(numberofSpines, numberofSpines);
CausalCorrelation = nan(numberofSpines, numberofSpines);
CausalPValue = nan(numberofSpines,numberofSpines);


%%% Correlate spines on the same dendrite as a function of dendritic
%%% distance between them

if DendNum ~= length(File.DendritePolyPointNumber)
    error('The polyline ROIs are not binned correctly; this file needs to be re-analyzed!')
end
counter = 1;
for i = 1:File.NumberofDendrites
    Branch{i}(1,1) = 0;
%     zeropoint = [(File.PolyLinePos{counter}(1)+File.PolyLinePos{counter}(3)/2), (File.PolyLinePos{counter}(2)+File.PolyLinePos{counter}(4)/2)]; %%% zeropoint(x,y) = center in x,y coordinates of the beginning of a dendrite
    PolyX_center{i}(1,1) = File.PolyLinePos{counter}(1)+File.PolyLinePos{counter}(3)/2;
    PolyY_center{i}(1,1) = File.PolyLinePos{counter}(2)+File.PolyLinePos{counter}(4)/2;
    Pix_Dist{i}(1,1) = 0;
    Mic_Dist{i}(1,1) = 0;
    for j = 2:File.DendritePolyPointNumber(i)
        counter = counter+1;
        PolyX_center{i}(1,j) = File.PolyLinePos{counter}(1)+File.PolyLinePos{counter}(3)/2;
        PolyY_center{i}(1,j) = File.PolyLinePos{counter}(2)+File.PolyLinePos{counter}(4)/2;
        Pix_Dist{i}(1,j) = sqrt((PolyX_center{i}(1,j)-PolyX_center{i}(j-1)).^2 + (PolyY_center{i}(j)-PolyY_center{i}(j-1)).^2);
        Mic_Dist{i}(1,j) = Pix_Dist{i}(1,j)/pixpermicron;
    end
    counter = counter+1;
    for j = File.SpineDendriteGrouping{i}(1):File.SpineDendriteGrouping{i}(end)
        spine_pos{j} = [File.ROIPosition{j+1}(1)+File.ROIPosition{j+1}(3)/2, File.ROIPosition{j+1}(2)+File.ROIPosition{j+1}(4)/2]; %%% Don't forget that position 1 in this cell is actually ROI0/background ROI!!!! 
        [distance, index] = min(sqrt(((PolyX_center{i}-spine_pos{j}(1)).^2)+(PolyY_center{i}-spine_pos{j}(2)).^2)); %%% Find the closest ROI along the dendrite (usually spaced evenly and regularly enough that it should be right at the base of the spine, more or less)
%         spine_address{j} = [PolyX_center{i}(1,index), PolyY_center{i}(1,index)]; %%% Set a spine's "address" as that point along the dendrite, found above
        spine_address{j}.Dendrite = i;
        spine_address{j}.Index = index;
    end
    if length(File.SpineDendriteGrouping{i})>1
        for j = File.SpineDendriteGrouping{i}(1):File.SpineDendriteGrouping{i}(end-1)
            for k = (j+1):File.SpineDendriteGrouping{i}(end)
                SpineToSpineDistance(j,k) = abs(sum(Mic_Dist{spine_address{j}.Dendrite}(spine_address{j}.Index:spine_address{k}.Index))-Mic_Dist{spine_address{j}.Dendrite}(spine_address{j}.Index));  %%% Find the sum of linear distances from the current point to the nearby spine
                if sum(square(j,:))>0 && sum(square(k,:))>0 %%% In the case of flat lines, correlations are NaN, so make it zero instead
                    [r_all, p_all] = corrcoef(square(j,:)', square(k,:)');
                    OverallCorrelation(j,k) = r_all(1,2);
                    OverallPValue(j,k) = p_all(1,2);
                else
                    OverallCorrelation(j,k) = NaN;
                    OverallPValue(j,k) = 1;
                end
                if sum(synapseOnly(j,:))>0 && sum(synapseOnly(k,:))>0 %%% In the case of flat lines, correlations are NaN, so make it zero instead
                    [r, p] = corrcoef(synapseOnly(j,:)', synapseOnly(k,:)');
                    SpineToSpineCorrelation(j,k) = r(1,2);
                    SpineToSpinePValue(j,k) = p(1,2);
                else
                    SpineToSpineCorrelation(j,k) = NaN;
                    SpineToSpinePValue(j,k) = 1;
                end
                if sum(withAPs(j,:))>0 && sum(withAPs(k,:))>0 %%% In the case of flat lines, correlations are NaN, so make it zero instead
                    [r_AP, p_AP] = corrcoef(withAPs(j,:)', withAPs(k,:)');
                    SpwAPCorrelation(j,k) = r_AP(1,2);
                    SpwAP_PValue(j,k) = p_AP(1,2);
                else
                    SpwAPCorrelation(j,k) = NaN;
                    SpwAP_PValue(j,k) = 1;
                end
                if sum(causal(j,:))>0 && sum(causal(k,:))>0
                    [r_causal, p_causal] = corrcoef(causal(j,:)', causal(k,:)');
                    CausalCorrelation(j,k) = r_causal(1,2);
                    CausalPValue(j,k) = p_causal(1,2);
                else
                    CausalCorrelation(j,k) = NaN;
                    CausalPValue(j,k) = 1;
                end
            end
        end 
    else
    end
end


nonnan = find(~isnan(SpineToSpineDistance)); %% Find the indices for  non-NaN values
Correlations = SpineToSpineCorrelation(nonnan);
pValues = SpineToSpinePValue(nonnan);
wAPCorrelations = SpwAPCorrelation(nonnan);
SpwAP_PValue = SpwAP_PValue(nonnan);
Distances = SpineToSpineDistance(nonnan);
CausalCorrelations = CausalCorrelation(nonnan);
CausalPValues = CausalPValue(nonnan);

a.SpineToSpineDistance = Distances;
a.OverallCorrelationsHeatMap = OverallCorrelation;
a.OverallCorrelation = OverallCorrelation(nonnan);
a.SpineToSpineCorrelation = Correlations;
a.CorrelationHeatMap = SpineToSpineCorrelation;
a.PValueHeatMap = SpineToSpinePValue;
a.CausalHeatMap = CausalCorrelation;
a.DistanceHeatMap = SpineToSpineDistance;
a.CausalPValueHeatMap = CausalPValue;
a.SpinewithAP_Correlation = wAPCorrelations;
a.SpineToSpine_PValues = pValues;
a.SpinewithAP_PValues = SpwAP_PValue;
a.CausalCorrelations = CausalCorrelations;
a.CausalPValues = CausalPValues;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Analysis of individual clusters %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ClusterThresh = 0.5;

Clustnum = 1;
%%% Find all values that are greater than the cluster threshold
[row, col] = find(SpineToSpineCorrelation>=ClusterThresh);
[Crow, Ccol] = find(CausalCorrelation>=ClusterThresh);
Dendind = [];
CDendind = [];
addresses = cell(1,DendNum);
Caddresses = cell(1,DendNum);

%%% Make full correlation matrix %%%
tempA = SpineToSpineCorrelation;
tempB = SpineToSpineCorrelation';
tempA(isnan(tempA) & ~isnan(tempB)) = tempB(isnan(tempA) & ~isnan(tempB));
fullmat = tempA;
tempC = SpineToSpineDistance;
tempD = SpineToSpineDistance';
tempC(isnan(tempC) & ~isnan(tempD)) = tempD(isnan(tempC) & ~isnan(tempD));
fullDist = tempC;

%%% 'Synapse only' clusters
for i = 1:length(row)
    for j = 1:DendNum
        if ~isempty(find(File.SpineDendriteGrouping{j} == row(i)))
            Dendind = [Dendind; j];
            addresses{j} = [addresses{j}; row(i), col(i)];
        end
    end
end
usedDend = unique(Dendind);

for i = 1:length(addresses)
    if ~isempty(addresses{i})
        if size(addresses{i},1)>1
            spines = unique(addresses{i});
            clust = spines(1);
            nonclust = [];
            for j = 2:length(spines)
                if sum(fullmat(clust, spines(j))>ClusterThresh)==length(clust)  %%% If all the indices in the 'clust' array yield correlation values > ClustThresh, then it should return all logical == 1, so the sum should be the same as the length of the 'clust' array
                    clust = [clust; spines(j)];
                else
                    nonclust = [nonclust; spines(j)];
                end
            end
            Clustered{Clustnum} = clust;
            Clustnum = Clustnum+1;
            while length(clust)<length(spines)
                for j = 1:length(clust)
                    spines = spines(spines~=clust(j));  %%% Replace 'spines' array with only the ones that haven't been used yet
                end
%                 if length(spines) == 1
%                     continue
%                 end
                clust = spines(1);
                options = unique(addresses{i});
                for j = 1:length(options)
                    if sum(fullmat(clust, options(j))>ClusterThresh)==length(clust)
                        clust = [clust; options(j)];
                    end
                end
                if length(clust)>1
                    Clustered{Clustnum} = clust;
                    Clustnum = Clustnum+1;
                end
            end
        else
            spines = unique(addresses{i});
            Clustered{Clustnum} = spines';
            Clustnum = Clustnum+1;
        end
    else
        Clustered{Clustnum} = [];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find the largest distance covered by the spines in a each cluster

for i = 1:length(Clustered)
    if length(Clustered{i})>1
        combinations = [];
        dist = [];
        combinations = nchoosek(Clustered{i},2); %%% Find all combinations of spines (two at a time) in a given cluster
        for j = 1:size(combinations, 1);
            dist(j) = fullDist(combinations(j,1), combinations(j,2));
        end
        ClustLength{i} = nanmean(dist);
    else
        ClustLength{i} = [];
    end
end

a.Clustered_Spines = Clustered;
a.Cluster_Length = ClustLength;

%%% Repeat the above for causal clusters
%%% Make full correlation matrix %%%
tempA = CausalCorrelation;
tempB = CausalCorrelation';
tempA(isnan(tempA) & ~isnan(tempB)) = tempB(isnan(tempA) & ~isnan(tempB));
fullmat = tempA;

%%% Causal clusters

Clustnum = 1;

for i = 1:length(Crow)
    for j = 1:DendNum
        if ~isempty(find(File.SpineDendriteGrouping{j} == Crow(i)))
            CDendind = [CDendind; j];
            Caddresses{j} = [Caddresses{j}; Crow(i), Ccol(i)];
        end
    end
end

CusedDend = unique(CDendind);


for i = 1:length(Caddresses)
    if ~isempty(Caddresses{i})
        if size(Caddresses{i},1)>1
            spines = unique(Caddresses{i});
            clust = spines(1);
            nonclust = [];
            for j = 2:length(spines)
                if sum(fullmat(clust, spines(j))>ClusterThresh)==length(clust)  %%% If all the indices in the 'clust' array yield correlation values > 0.5, then it should return all logical == 1, so the sum should be the same as the length of the 'clust' array
                    clust = [clust; spines(j)];
                else
                    nonclust = [nonclust; spines(j)];
                end
            end
            CausalClustered{Clustnum} = clust;
            Clustnum = Clustnum+1;
            while length(clust)<length(spines)
                for j = 1:length(clust)
                    spines = spines(spines~=clust(j));  %%% Replace 'spines' array with only the ones that haven't been used yet
                end
%                 if length(spines) == 1
%                     continue
%                 end
                clust = spines(1);
                options = unique(Caddresses{i});
                for j = 1:length(options)
                    if sum(fullmat(clust, options(j))>ClusterThresh)==length(clust)
                        clust = [clust; options(j)];
                    end
                end
                if length(clust)>1
                    CausalClustered{Clustnum} = clust;
                    Clustnum = Clustnum+1;
                end
            end
        else
            spines = unique(Caddresses{i});
            CausalClustered{Clustnum} = spines';
            Clustnum = Clustnum+1;
        end
    else
        CausalClustered{Clustnum} = [];
    end
end


for i = 1:length(CausalClustered)
    if length(CausalClustered{i})>1
        combinations = [];
        dist = [];
        combinations = nchoosek(CausalClustered{i},2); %%% Find all combinations of spines (two at a time) in a given cluster
        for j = 1:size(combinations, 1);
            dist(j) = fullDist(combinations(j,1), combinations(j,2));
        end
        CausalClustLength{i} = nanmean(dist);
    else
        CausalClustLength{i} = [];
    end
end

a.CausalClustered_Spines = CausalClustered;
a.CausalCluster_Length = CausalClustLength;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Spectral graph analysis of degree of clustering for each spine and
%%% dendrite

DendClustering = zeros(1,DendNum);
weightedCluster = cell(1,DendNum);
DistanceEigenVectors = cell(1,DendNum);
DistanceEigenValues = cell(1,DendNum);

for i = 1:DendNum
    firstspine = File.SpineDendriteGrouping{i}(1);
    lastspine = File.SpineDendriteGrouping{i}(end);
    if firstspine ~= lastspine
        A{i} = fullDist(firstspine:lastspine, firstspine:lastspine);
%         if max(max(A{i}))>0
%             A{i} = A{i}/max(max(A{i}));     %%% Normalize to dendrite length;
%         else
%         end
%         A{i}(A{i}<10) = 1;
%         A{i}(A{i}>10) = 0;
        A{i}(A{i}<1) = 1;
        A{i} = 1./exp(A{i}./SpectralLengthConstant);          %%% Adjacency matrix --> 1/e^x, where x = distance
        A{i}(isnan(A{i})) = 0;                                %%% Since the diagonal of the laplacian == the degree, set NaNs in A to be zero to maintain this identity;
        degs = sum(A{i},2);
        degs(degs==0) = eps;
        D{i} = sparse(1:size(A{i},1),1:size(A{i},2),degs);    %%% Degree matrix
        L{i} = D{i}-A{i};                                     %%% Laplacian matrix
        Dinv{i} = inv(D{i});                                  %%% Determine the inverse Degree matrix
        nL{i}= Dinv{i} * L{i};
        [eVecs eVals] = eig(nL{i});                           %%% Find the eigenvectors and eigenvalues for the Laplacian
        DistanceEigenVectors{i} = eVecs;
        DistanceEigenValues{i} = eVals;
        e = diag(eVals);
        DendClustering(1,i) = min(e(~ismember(e,min(e))));   %%% Finds the SECOND smallest eigenvalue (the Fiedler value or algebraic connectivity) which corresponds to the extent of clustering for the whole dendrit
        weightedCluster{i} = [];
    else
        A{i} = [];
        D{i} = [];
        L{i} = [];
        nL{i} = [];
        DendClustering(1,i) = NaN;
        weightedCluster{i} = [];
    end
end

a.AdjacencyMatrix = A;
a.DegreeMatrix = D;
a.LaplacianMatrix = L;
a.NormalizedLaplacian = nL;
a.DistanceEigenVectors = DistanceEigenVectors;
a.DistanceEigenValues = DistanceEigenValues; 
a.DendriteClusterDegree = DendClustering;
a.WeightedClustering = weightedCluster;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%%% Correlate distant spines on separate dendrites
if File.NumberofDendrites > 1 
    nofSp = File.NumberofSpines;
    nofD = File.NumberofDendrites;
    farCorr = nan(nofSp, nofSp);
    farDist = nan(nofSp, nofSp);     
    for DendCount = 1:nofD-1                                                                        %%% For each dendrite
        sptarg = 1; %%% spine that is being compared to all others (in 1v2, 1v3, etc. setup)
        for j = File.SpineDendriteGrouping{DendCount}(1):File.SpineDendriteGrouping{DendCount}(end) %%% For each spine on the current dendrite
            spcomp = 1; %%% comparator spine
            for k = File.SpineDendriteGrouping{DendCount+1}(1):File.SpineDendriteGrouping{end}(end) %%% Compare to each spine on a different dendrite
                FarSpineDist(j,k) = sqrt((spine_pos{j}(1)-spine_pos{k}(1)).^2 + (spine_pos{j}(2)-spine_pos{k}(2)).^2); 
                FarSpineDist(j,k) = FarSpineDist(j,k)/pixpermicron;
                if sum(synapseOnly(j,:))>0 && sum(synapseOnly(k,:))>0
                    [r, p] = corrcoef(synapseOnly(j,:)', synapseOnly(k,:)');
                    FarSpineCorrelation(j,k) = r(1,2);
                    FarPValue(j,k) = p(1,2);
                else
                    FarSpineCorrelation(j,k) = 0;
                    FarPValue(j,k) = 1;
                end
                spcom = spcomp+1;
            end
            sptarg = sptarg + 1;
        end
    end
else
    FarSpineDist = 0;
    FarSpineCorrelation = 0;
end

nonzero = find(FarSpineCorrelation);
FarCorrelations = FarSpineCorrelation(nonzero);
FarDistances = FarSpineDist(nonzero);


a.FarSpineToSpineDistance = FarDistances;
a.FarSpineToSpineCorrelation = FarCorrelations;
a.SynapseOnlyBinarized = synapseOnly;
a.OverallSpineActivity = square;


options = statset('MaxIter', 1000); 

X = [];
X(:,1) = Distances;
X(:,2) = Correlations;
try
    obj = []; obj = fitgmdist(X,3, 'Options', options);
    idx = []; idx = cluster(obj,X);
    cluster1 = X(idx ==1,:); cluster2 = X(idx ==2,:); cluster3 = X(idx ==3,:);
    group = [mean(cluster1(:,2)) mean(cluster2(:,2)) mean(cluster3(:,2))];
    poss = [1 2 3];
    [val, ind] = max(group);
    eval(['TopGroup = cluster', num2str(ind), ';']); poss = poss(poss~=ind);
    [val, ind] = min(group);
    eval(['BottomGroup = cluster', num2str(ind), ';']); poss = poss(poss~=ind);
    eval(['MiddleGroup = cluster', num2str(poss), ';']); 
catch
    TopGroup(:,1) = Distances;
    TopGroup(:,2) = Correlations;
    MiddleGroup(:,1) = Distances;
    MiddleGroup(:,2) = Correlations;
    BottomGroup(:,1) = Distances;
    BottomGroup(:,2) = Correlations;
end


X = [];
X(:,1) = Distances;
X(:,2) = CausalCorrelations;
try
    obj = []; obj = fitgmdist(X,3, 'Options', options);
end
try
    idx = []; idx = cluster(obj,X);
    causalcluster1 = X(idx ==1,:); causalcluster2 = X(idx ==2,:); causalcluster3 = X(idx ==3,:);
    causalgroup = [mean(causalcluster1(:,2)) mean(causalcluster2(:,2)) mean(causalcluster3(:,2))];
    poss = [1 2 3];
    [val ind] = max(causalgroup);
    eval(['TopCausalGroup = causalcluster', num2str(ind), ';']); poss = poss(poss~=ind);
    [val ind] = min(causalgroup);
    eval(['BottomCausalGroup = causalcluster', num2str(ind), ';']); poss = poss(poss~=ind);
    eval(['MiddleCausalGroup = causalcluster', num2str(poss), ';']); 
catch
    TopCausalGroup(:,1) = Distances;
    TopCausalGroup(:,2) = CausalCorrelations;
    MiddleCausalGroup(:,1) = Distances;
    MiddleCausalGroup(:,2) = CausalCorrelations;
    BottomCausalGroup(:,1) = Distances;
    BottomCausalGroup(:,2) = CausalCorrelations;
end

a.TopCluster = TopGroup;
a.MiddleCluster = MiddleGroup;
a.BottomCluster = BottomGroup;
a.TopCausalCluster = TopCausalGroup;
a.MiddleCausalCluster = MiddleCausalGroup;
a.BottomCausalCluster = BottomCausalGroup;

%%%%
a.Session = currentsession;
cd(vardir); %%% Variable directory
%%%%

savefile = [folder, '_' Date, '_Summary'];
eval([savefile, '= a;']);
save(savefile, savefile)

if showFig == 1
    figure('Position', [(Scrsz(3)/2)-(Scrsz(3)/2)/2, Scrsz(2)+(Scrsz(4)/2)/2 ,Scrsz(3)/2,Scrsz(4)/2]); 
    subplot(1,3,1)
%     plot(Distances, Correlations, 'ok')
    plot(BottomGroup(:,1), BottomGroup(:,2), 'o', 'Color', [0.2 0.2 0.2]); hold on;
    plot(MiddleGroup(:,1), MiddleGroup(:,2), 'o', 'Color', red)
    plot(TopGroup(:,1), TopGroup(:,2), 'o', 'Color', lblue)
    xlabel('Proximity of Spines (um)')
    ylabel('Correlation')
    ylim([-0.05 1])
    title('Spine on the same dendrite')
    subplot(1,3,2)
    plot(FarDistances, FarCorrelations, 'or')
    ylim([-0.05 1])
    xlabel('Proximity of Spines (um)')
    ylabel('Correlation')
    title('Spines on separate dendrites')
    subplot(1,3,3)
%     plot(Distances, CausalCorrelations, 'oc')
    plot(BottomCausalGroup(:,1), BottomCausalGroup(:,2), 'o', 'Color', [0.2 0.2 0.2]); hold on;
    plot(MiddleCausalGroup(:,1), MiddleCausalGroup(:,2), 'o', 'Color', red)
    plot(TopCausalGroup(:,1), TopCausalGroup(:,2), 'o', 'Color', lblue)
    xlabel('Proximity of Spines (um)')
    ylabel('Correlation')
    title('Events Causing APs')
    ylim([-0.05 1])
else
end

% spatialfit = robustfit(SpineToSpineDistance, SpineToSpineCorrelation);
% fitcurve = spatialfit(2)*([0:0.1:max(SpineToSpineDistance)])+spatialfit(1);







