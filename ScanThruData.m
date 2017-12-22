function [X y] = ScanThruData(File, spine)

spinethreshmultiplier = 1.5;
raw = File.SynapticEvents{spine};
rawmed = nanmedian(raw);
rawspread = nanstd(raw);
raw(raw>rawmed+2*spinethreshmultiplier*rawspread) = rawmed+2*spinethreshmultiplier*rawspread; %%% Cap off large and small values to pinch the data towards the true baseline
raw(raw<rawmed-2*spinethreshmultiplier*rawspread) = rawmed-2*spinethreshmultiplier*rawspread; %%%
ave = smooth(raw,1000); %%% Baseline value
smoothed = smooth(File.SynapticEvents{spine},20);
smoothed = smoothed-ave;
smoothed(smoothed<median(smoothed)-spinethreshmultiplier*std(smoothed)) = median(smoothed)-spinethreshmultiplier*std(smoothed);   %%% This is the dendrite-subtracted data, so it's possible that there are very large negative events, which are artificial, and should be reduced to ~the level of noise
smoo = smoothed;
trace = smoo;
med = nanmedian(smoo);
spread = nanstd(smoo);


scrsz = get(0, 'ScreenSize');
h1 = figure('Position', [((scrsz(3)/3)-(scrsz(3)/6)), scrsz(4)/2.2, scrsz(3)/3, scrsz(4)/2.2]);
a1 = axes; hold on;
wholecurve = plot(trace, 'k');


h2 = figure('Position', [((scrsz(3)/3)+(scrsz(3)/6)), scrsz(4)/2.2, scrsz(3)/3, scrsz(4)/2.2], 'KeyPressFcn', @(obj, evt) 0);
a2 = axes; hold on;
set(gcf, 'CurrentCharacter', '2');
counter = 0;
start = 1;
fin = 1;
divider = 1000;
example = 0;

while fin < length(trace)
    fin = start+100;
    if fin > length(trace)
        fin = length(trace)
    end
    axes(a1); cla; plot(trace(1:end), 'k'); plot(start:fin, trace(start:fin), 'r', 'LineWidth', 2);
    axes(a2); cla; plot(trace(start:fin), 'k')
    try
        waitfor(gcf, 'CurrentCharacter')
    catch
        continue
    end
    inpt = get(gcf, 'CurrentCharacter');
    while strcmpi(inpt, 'a') || strcmpi(inpt, 'd') || strcmpi(inpt, 'b') %%% the character b corresponds to an intermediate choice, since 'waitfor' only recognizes changes in the character
        axes(a2)
%         waitfor(gcf, 'CurrentCharacter')
        newinpt = get(gcf, 'CurrentCharacter');
        if strcmpi(newinpt, 'd')
            start = start+1;
            fin = fin+1;
            axes(a2); cla; plot(trace(start:fin), 'k')
            axes(a1); cla; plot(trace, 'k'); hold on; plot(start:fin, trace(start:fin), 'r', 'LineWidth', 2)
            set(h2, 'CurrentCharacter', 'b')
        elseif strcmpi(newinpt, 'a')
            if start~=1
                start = start-1;
            else
                start = start;
                fprintf('Cannot shift to the left')
            end
            fin = fin-1;
            axes(a2); cla; plot(trace(start:fin), 'k')
            axes(a1); cla; plot(trace, 'k'); hold on; plot(start:fin, trace(start:fin), 'r', 'LineWidth', 2)
            set(h2, 'CurrentCharacter', 'b')
        elseif strcmpi(newinpt, '0') || strcmpi(newinpt, '1')
            inpt = newinpt;
            break
        end
    end
    if strcmpi(inpt, '0') || strcmpi(inpt, '1')
        example = example+1;
        decision(example,1) = inpt;
        if strcmpi(inpt, '1')  
            fprintf('Example %d marked as positive \n', example);
        end
        if fin == length(trace)
            sweep(example,:) = zeros(1,length(sweep(1,:)));
        else
            sweep(example, :) = trace(start:fin);
        end
    end
    set(gcf, 'CurrentCharacter', '2');
    counter = counter+1;
    start = fin+1;
end

y = str2num(decision);
X = sweep;

close(h1); close(h2);



