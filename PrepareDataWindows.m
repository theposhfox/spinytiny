function X = PrepareDataWindows(File, spine)

spinethreshmultiplier = 1.5;
raw = File.SynapticEvents{spine};
rawmed = nanmedian(raw);
rawspread = nanstd(raw);
raw(raw>rawmed+2*spinethreshmultiplier*rawspread) = rawmed+2*spinethreshmultiplier*rawspread; %%% Cap off large and small values to pinch the data towards the true baseline
raw(raw<rawmed-2*spinethreshmultiplier*rawspread) = rawmed-2*spinethreshmultiplier*rawspread; %%%
ave(spine,:) = smooth(raw,1000); %%% Baseline value
smoothed = smooth(File.SynapticEvents{spine},20);
smoothed = smoothed-ave(spine,:)';
smoothed(smoothed<median(smoothed)-spinethreshmultiplier*std(smoothed)) = median(smoothed)-spinethreshmultiplier*std(smoothed);   %%% This is the dendrite-subtracted data, so it's possible that there are very large negative events, which are artificial, and should be reduced to ~the level of noise
smoo = smoothed;
trace = smoo;
med = nanmedian(smoo);
spread = nanstd(smoo);

start = 1;
fin = start+100;
sample = 1;


figure; hold on; 

while fin<length(trace)
    X = trace(start:fin);
    [Y,~,Af] = ClassifyCaTransients(X');
    if Y <0.5
        start = start+1;
        fin = start+100;
        sample = sample+1;
    else
        start = fin+1;
        fin = start+100;
        sample = sample+1
        plot(1:size(X), X, 'k')
    end
end

