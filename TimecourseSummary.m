function TimecourseSummary(varargin)

ns = length(varargin);

Frequency = nan(ns,14);
Amp = nan(ns,14);
SynapseOnlyFreq = nan(ns,14);
SpikeTimedEvents = nan(ns,14);
Dendritic_Freq = nan(ns,14);
Dendritic_Amp = nan(ns,14);

for j = 1:ns
    Frequency(j,1:14) = varargin{j}.Frequency;
    Amp(j,1:14) = varargin{j}.Amplitude;
    SynapseOnlyFreq(j,1:14) = varargin{j}.SynapseOnlyFreq;
    SpikeTimedEvents(j,1:14) = varargin{j}.SpikeTimedEvents;
    Dendritic_Freq(j,1:14) = varargin{j}.Dendritic_Frequency;
    Dendritic_Amp(j,1:14) = varargin{j}.Dendritic_Amp;
end


meanF = nanmean(Frequency,1);
semF = nanstd(Frequency,0,1)./sqrt(sum(~isnan(Frequency),1));
meanA = nanmean(Amp,1);
semA = nanstd(Amp,0,1)./sqrt(sum(~isnan(Amp),1));
meanSOF = nanmean(SynapseOnlyFreq,1);
semSOF = nanstd(SynapseOnlyFreq,0,1)./sqrt(sum(~isnan(SynapseOnlyFreq),1));
meanSTE = nanmean(SpikeTimedEvents,1);
semSTE = nanstd(SpikeTimedEvents,1)./sqrt(sum(~isnan(SpikeTimedEvents),1));
meanDF = nanmean(Dendritic_Freq,1);
semDF = nanstd(Dendritic_Freq,0,1)./sqrt(sum(~isnan(Dendritic_Freq),1));
meanDA = nanmean(Dendritic_Amp,1);
semDA = nanstd(Dendritic_Amp,0,1)./sqrt(sum(~isnan(Dendritic_Amp),1));


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
    
scrsz = get(0, 'ScreenSize');
figure('Position', scrsz)

subplot(2,3,1); plot(1:14,meanF,'k')
r_errorbar(1:14,meanF,semF,'k')
title('Event Frequency')
xlabel('Session')
ylabel('Events/min')

xlim([0 15])
subplot(2,3,2); plot(1:14,meanA,'k')
r_errorbar(1:14,meanA,semA,'k')
title('Event Amp')
xlim([0 15])
xlabel('Session')
ylabel('dF/F')

subplot(2,3,3); plot(1:14, meanSOF,'k')
r_errorbar(1:14,meanSOF,semSOF,'k')
title('Synapse Only Freq')
xlim([0 15])
xlabel('Session')
ylabel('Events/min')

subplot(2,3,4); plot(1:14, meanSTE,'k')
r_errorbar(1:14,meanSTE,semSTE,'k')
title('Spike Timed Freq')
xlim([0 15])
xlabel('Session')
ylabel('Events/min')

subplot(2,3,5); plot(1:14, meanDF,'k')
r_errorbar(1:14,meanDF,semDF,'k')
title('Dendritic Freq')
xlim([0 15])
xlabel('Session')
ylabel('Events/min')

subplot(2,3,6); plot(1:14, meanDA,'k')
r_errorbar(1:14,meanDA,semDA,'k')
title('Dendritic Amp')
xlim([0 15])
xlabel('Session')
ylabel('dF/F')
