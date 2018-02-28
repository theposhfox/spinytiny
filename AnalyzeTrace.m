function [Threshold, DriftBaseline, ProcessedData] = AnalyzeTrace(Data, Options)
    
    driftbaselinesmoothwindow = Options.DriftBaselineSmoothWindow;
    baselinesmoothwindow = Options.BaselineSmoothWindow;
    smoothwindow = Options.SmoothWindow;
    valueslimitfornoise = Options.ValuesLimitforNoise;
    valueslimitforbaseline = Options.ValuesLimitforBaseline;
    traceoption = Options.TraceOption;
    BeingAnalyzed = Options.BeingAnalyzed;
    
   %%% Data with NaN cannot be smoothed well, so find and fix any NaNs
    if any(isnan(Data))
        try
            Data(isnan(Data)) = nanmean([Data(find(isnan(Data))-1),Data(find(isnan(Data))+1)]);
        catch
            Data(isnan(Data)) = 0;
        end
    end
    
    %%% Values at the ends can mess up smoothing; set the first few to the
    %%% median of the first 1000 frames
    
    Data(1:10) = nanmedian(Data(1:100));
    Data(end-10:end) = nanmedian(Data(end-100:end));
    raw = Data;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Smoothing
    smoothed = smooth(raw,smoothwindow);     
        smoothed = reshape(smoothed, 1, length(smoothed));
        %%% Remove large negative values 
        if traceoption ==4
            smoothed(smoothed<(median(smoothed)-nanstd(smoothed))) = median(smoothed)-nanstd(smoothed);   %%% If using dendrite-subtracted data, it's possible that there are very large negative events, which are artificial, and should be reduced to ~the level of noise
        end
        smoothed_forbaseline = smoothed;
    rawmed = nanmedian(smoothed);
    rawspread = nanstd(smoothed);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Baseline %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Estimate baseline by performing x successful rounds of capping off large and small values,
    %%% then re-calculating the median 
    
%     roundstodo = 50;
%     roundnum = 1;    
%     
%     while roundnum<=roundstodo 
%         smoothed_forbaseline(smoothed_forbaseline>rawmed+valueslimitforbaseline*rawspread) = rawmed+valueslimitforbaseline*rawspread;      %%% Cap off large and small values to pinch the data towards the true baseline
%         smoothed_forbaseline(smoothed_forbaseline<rawmed-valueslimitforbaseline*rawspread) = rawmed-valueslimitforbaseline*rawspread;      %%%
%             rawspread = nanstd(smoothed_forbaseline);
%             rawmed = nanmedian(smoothed_forbaseline);
%         roundnum = roundnum+1;
%     end
% 
%     DriftBaseline = reshape(smooth(smoothed_forbaseline,driftbaselinesmoothwindow), 1, length(smoothed));
%     
%     driftsub = (smoothed-DriftBaseline)+nanmedian(DriftBaseline);
%     
%     roundnum = 1;    
%     
%     smoothed_forbaseline = driftsub;
%     rawmed = nanmedian(driftsub);
%     rawspread = nanstd(driftsub);
% 
%     while roundnum<=roundstodo 
%         smoothed_forbaseline(smoothed_forbaseline>rawmed+valueslimitforbaseline*rawspread) = rawmed+valueslimitforbaseline*rawspread;      %%% Cap off large and small values to pinch the data towards the true baseline
%         smoothed_forbaseline(smoothed_forbaseline<rawmed-valueslimitforbaseline*rawspread) = rawmed-valueslimitforbaseline*rawspread;      %%%
%             rawspread = nanstd(smoothed_forbaseline);
%             rawmed = nanmedian(smoothed_forbaseline);
%         roundnum = roundnum+1;
%     end

%     truebaseline = smooth(smoothed_forbaseline, baselinesmoothwindow)';


    starterdata = smoothed;

    %%% Kernel Density Estimation (Aki's method) %%%
    truebaseline = baseline_kde(starterdata',100,50,20);
    DriftBaseline = truebaseline;
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Baseline Subtraction
    
    if nanmedian(truebaseline)<0
        starterdata = starterdata+abs(min(truebaseline));
        truebaseline = truebaseline+abs(min(truebaseline));
    end
    
%     blsub = driftsub-truebaseline;                                                             %%% Baseline-subtracted value
    blsub = starterdata-truebaseline';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Baseline division (if using raw traces)
    
    if traceoption == 1 
        if nanmedian(truebaseline)~= 0
            dFoF = blsub/nanmedian(truebaseline);
%             blsub(blsub<0) = 0;
            rawmed = nanmedian(dFoF);
            rawspread = nanstd(dFoF);
        else
            blsub = blsub+1;
            truebaseline = truebaseline+1;
            dFoF = blsub/nanmedian(truebaseline);
%             blsub(blsub<0) = 0;
            rawmed = nanmedian(dFoF);
            rawspread = nanstd(dFoF);
        end
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Final variable
    
             processed_dFoF = smooth(dFoF, smoothwindow);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        fornoise = processed_dFoF;
        rawmed = nanmedian(processed_dFoF);
        rawspread = nanstd(processed_dFoF);
        
            
        roundnum = 1;
        roundstodo = 50;
        while roundnum<=roundstodo 
            fornoise(fornoise>rawmed+(valueslimitfornoise*rawspread)) = rawmed+(valueslimitfornoise*rawspread);      %%% Cap off large and small values to pinch the data towards the true baseline
            fornoise(fornoise<rawmed-(valueslimitfornoise*rawspread)) = rawmed-(valueslimitfornoise*rawspread);      %%%
                rawspread = nanstd(fornoise);
                rawmed = nanmedian(fornoise);
            roundnum = roundnum+1;
        end
    end
             
%     spread = rawmed+spinevalueslimitfornoise*nanstd(fornoise);
    spread = nanmax(fornoise);
    
    med = nanmedian(fornoise);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Noise Estimation %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Adjust spine threshold multiplier based on processed data for each
    %%% spine
    
    %%% Estimate signal by finding peaks
    pks = findpeaks(processed_dFoF, 'MinPeakHeight', spread, 'MinPeakDistance', 200,'sortstr', 'descend');
    
%     [f, xi] = ksdensity(blsub);
%     lowerlimit = prctile(xi,75);
    
    if isempty(pks)
        thresh = 1;
    else
        switch BeingAnalyzed
            case 'Spine';
                if spread < 0.5
                    thresh = 0.5;
                else
                    thresh = spread;
                end
            case 'Poly';
                if spread < 0.5
                    thresh = 0.5;
                else
                    thresh = spread;
                end
            case 'Dendrite'
                if spread < 0.5
                    thresh = 0.5;
                else
                    thresh = spread;
                end
        end

        signal = nanmean(pks);
        noise = spread;
        StN = signal/noise;
%         maxaddition = 1;
    end
    
    Threshold = thresh;
    ProcessedData = processed_dFoF;