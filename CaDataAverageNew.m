function CaDataAverage(Subject_File, Condition, Excluder, Red);

load(Subject_File);

global gui_CaAnalysis

WS_Data = who(['*',Condition, '*']);
WS_Data2 = [];

Compare_cond = strfind(Condition, 'Compare');

selected_files = 1;

axes(gui_CaAnalysis.figure.handles.Graph); hold on;

fig_axes = get(gca);

offset = fig_axes.Position;

plot_empty_check = fig_axes.Children;

if isempty(plot_empty_check)
    gui_CaAnalysis.color_counter = 1;
    gui_CaAnalysis.UTA =[];
    gui_CaAnalysis.VolumeData = [];
else
    gui_CaAnalysis.color_counter = gui_CaAnalysis.color_counter+1;
end

trace_num = gui_CaAnalysis.color_counter;



for i = 1:length(WS_Data);
    if ~isempty(WS_Data{i});
        a = strfind(WS_Data{i}, Condition);
        if isempty(Compare_cond)     ;           %%% If the desired condition includes the term 'Compare', include files that meet these criteria
            b = strfind(WS_Data{i}, ['Compare', Condition]);  %%% Otherwise, exclude files that contain the 'Compare' modifier (as in the paired control for a drug condition)
        else
            b = [];                           
        end
        if ~isempty(Excluder);
            for j = 1:length(Excluder);
                c{j} = strfind(WS_Data{i}, Excluder{j});
            end
        else
            c = 1;
        end
        if ~isempty(a) && isempty(b)
            d = 0;
            for j = 1:length(Excluder);
                if isempty(c{j})
                    d = d + 0;
                else
                    d = d + 1;
                end
            end
            if d == 0
                WS_Data2{selected_files} = WS_Data{i};
                selected_files = selected_files+1;
            end
        end
    end
end

for i = 1:length(WS_Data2)
    Selected_Samples(i,1) = WS_Data2(1,i);
end

Selected_Samples
n = length(Selected_Samples)

file_num = 0;

for i = 1:length(WS_Data2)
    file_num = file_num+1;
    evalc(['va{', num2str(file_num), '}=', WS_Data2{i}]);
end

% for i = 1:length(va)
%     spine_data(:,i) = va{i}.deltaF{1}(:);
%     for j = 2:length(va{i}.deltaF)
%         nearby_spine_data{j}(:,i) = va{i}.deltaF{j}(:);
%     end
%     for k = 1:length(va{i}.Poly_deltaF)
%         dendrite_data{k}(:,i) = va{i}.Poly_deltaF{j}(:);
%     end
% end


%%%%% Uncomment to calculate deltaF/F and/or red values %%%%%%%
gui_CaAnalysis.RawTrace{trace_num} = [];

for i = 1:length(va)
    spine_delta(:,i) = va{i}.Fluorescence_Measurement{1}- mean(va{i}.Fluorescence_Measurement{1}(1:32));
    red_delta(:,i) = va{i}.Red_Measurement{1}./mean(va{i}.Red_Measurement{1}(1:32));
    spine_dFoF(:,i) = (spine_delta(:,i))./(mean(va{i}.Fluorescence_Measurement{1}(1:32)));
        delta_ratio(:,i) = spine_delta(:,i)./(imfilter(red_delta(:,i), ones(16,1)/16, 'replicate'));
        baseline_ratio(:,i) = va{i}.Fluorescence_Measurement{1}(1:32)./va{i}.Red_Measurement{1}(1:32);
        Red_norm(:,i) = delta_ratio(:,i)./mean(baseline_ratio(:,i));
        
    red_data(:,i) = imfilter(red_delta(:,i), ones(16,1)/16, 'replicate');  %%% uncomment to see the red change instead of the green (for spine volume)
%     spine_data(:,i) = (spine_dFoF(:,i)+1)./(red_data(:,i)+1);  %%% uncomment to see the red change instead of the green (for spine volume)
%     spine_data(:,i) = spine_dFoF(:,i)./va{i}.Red_Measurement{1}';
    spine_data(:,i) = spine_dFoF(:,i)+1;
    
    for j = 2:length(va{i}.deltaF)
        nearby_spine_data{j}(:,i) = va{i}.deltaF{j}(:);
    end
    for k = 1:length(va{i}.Poly_deltaF)
        dendrite_data{k}(:,i) = va{i}.Poly_deltaF{j}(:);
    end
    for j = 0:29
        UTA(j+1, 1:17) = spine_dFoF(32 + j*16 : 32 + (j+1) * 16,i) - spine_dFoF(32 + j*16,i);
    end
    gui_CaAnalysis.UTA{trace_num}(i,:) = mean(UTA,1);
end

gui_CaAnalysis.VolumeData{trace_num} = red_data;

if Red == 0;
    gui_CaAnalysis.RawTrace{trace_num} = spine_dFoF';
%     gui_CaAnalysis.RawTrace{trace_num} = spine_delta;
else
    gui_CaAnalysis.RawTrace{trace_num} = Red_norm';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stim_spine_average = mean(spine_data, 2);
stim_spine_SEM = (std(spine_data,0,2))/sqrt(size(spine_data,2));
red_spine_average = mean(red_data,2);
red_spine_SEM = (std(red_data,0,2))/sqrt(size(red_data,2));
vol_corr_average = stim_spine_average./red_spine_average;
Red_norm_average = mean(Red_norm,2);
Red_norm_SEM = (std(Red_norm,0,2))/sqrt(size(Red_norm,2));

for i = 1:length(nearby_spine_data)
    nearby_spine_average{i} = nanmean(nearby_spine_data{i},2);
    nearby_spine_SEM{i} = (nanstd(nearby_spine_data{i},0,2))/sqrt(size(nearby_spine_data{i},2));
end

for i = 1:length(dendrite_data)
    dendrite_average{i} = nanmean(dendrite_data{i},2);
    dendrite_SEM{i} = (nanstd(dendrite_data{i},0,2))/sqrt(size(dendrite_data{i},2));
end

spine_SEM_upper = stim_spine_average+stim_spine_SEM;
spine_SEM_lower = stim_spine_average-stim_spine_SEM;
red_SEM_upper = red_spine_average+red_spine_SEM;
red_SEM_lower = red_spine_average-red_spine_SEM;
dendrite_SEM_upper = dendrite_average{1}+dendrite_SEM{1};
dendrite_SEM_lower = dendrite_average{1}-dendrite_SEM{1};
norm_SEM_upper = Red_norm_average+Red_norm_SEM;
norm_SEM_lower = Red_norm_average-Red_norm_SEM;

Time = 0+0.125:0.125:70;
Time = Time-4;
Time = Time(1:560);


%%%% Graph %%%%



    %%%%%%%%%%%%%%%%%%%%%%
    %%Color Information%%%
    %%%%%%%%%%%%%%%%%%%%%%

    lgray = [0.50 0.51 0.52];       brown = [0.28 0.22 0.14];
    gray = [0.30 0.31 0.32];        lbrown = [0.59 0.45 0.28];
    yellow = [1.00 0.76 0.05];      orange = [0.95 0.40 0.13];
                                    lorange  = [0.75 0.20 0.06];
    green = [0.00 0.43 0.23];       blue = [0.00 0.33 0.65];
    lgreen = [0.55 0.78 0.25];      lblue = [0.00 0.68 0.94];     
    magenta = [0.93 0.22 0.55];     
    lmagenta = [0.93 0.44 0.75];
    purple = [0.57 0.15 0.56];      
    lpurple = [0.85 0.15 0.86];   
    red = [0.93 0.11 0.14];         black = [0 0 0];
    pink = [1.00 0.6 0.6];
    
    colorj = {red,   blue,  green,  gray,  brown,  orange,  purple, magenta,  black, yellow};
    colork = {pink, lblue, lgreen, lgray, lbrown, lorange, lpurple, lmagenta, gray, lgreen};
    

% pdend = patch([Time fliplr(Time)], [dendrite_SEM_upper', fliplr(dendrite_SEM_lower')], colork{gui_CaAnalysis.color_counter+1});
% set(pdend, 'EdgeColor', 'none');


if Red == 0
    p = patch([Time fliplr(Time)], [spine_SEM_upper', fliplr(spine_SEM_lower')], colork{gui_CaAnalysis.color_counter});
    set(p, 'EdgeColor', 'none');

    % Patch_alpha = findobj('Type', 'patch');

    % set(Patch_alpha, 'FaceAlpha', 0.5);

    hold on; gui_CaAnalysis.Curve(trace_num) = plot(Time,stim_spine_average, 'color', colorj{gui_CaAnalysis.color_counter}, 'Linewidth', 2);

    kcheck_ones = mean(spine_data(1:32),2)-1;
    
elseif Red == 1
    %%% Red data corrections, etc. %%%
    %     figure; hold on; plot(Time,vol_corr_average, 'color', colorj{gui_CaAnalysis.color_counter}, 'Linewidth', 2);
    p = patch([Time fliplr(Time)], [norm_SEM_upper', fliplr(norm_SEM_lower')], colork{gui_CaAnalysis.color_counter});
    set(p, 'EdgeColor', 'none');
    hold on; gui_CaAnalysis.Curve(trace_num) = plot(Time, Red_norm_average, 'color', colorj{gui_CaAnalysis.color_counter}, 'Linewidth',2)
        % set(p_red, 'EdgeColor', 'none');
        % plot(Time, dendrite_average{1}, 'color', colorj{gui_CaAnalysis.color_counter+1}, 'LineWidth', 2);
    kcheck_ones = mean(Red_norm_average(1:32),1);
end



if kcheck_ones < -0.5 
    k = zeros(1,length(Time));
elseif kcheck_ones > -0.5
    k = ones(1,length(Time));
end

plot(Time, k, '--k')
set(gca, 'XTickMode', 'auto');
set(gca, 'YTickMode', 'auto');


%%%% Volume inset %%%%%%

if isfield(gui_CaAnalysis, 'vol_axes')
    axes(gui_CaAnalysis.vol_axes)
else
    gui_CaAnalysis.vol_axes = axes('Position', [0.4 0.75 0.1 0.1]);
end
hold on; 
p_red = patch([Time fliplr(Time)], [red_SEM_upper', fliplr(red_SEM_lower')], colork{gui_CaAnalysis.color_counter});
set(p_red, 'EdgeColor', 'none');
gui_CaAnalysis.VolCurve(trace_num) = plot(Time, red_spine_average, 'color', colorj{gui_CaAnalysis.color_counter});
xlim([-5 65]);
ylim([0.75 2.5])

%%%%%%%%%%%%%%%%%%%%%%%

if isfield(gui_CaAnalysis, 'legend')
    gui_CaAnalysis.legend{trace_num} = [Condition, '; n = ', num2str(length(va))];
    l1 = legend(gui_CaAnalysis.Curve, gui_CaAnalysis.legend, 'Location', 'Northwest', 'boxoff');
    set(l1, 'Interpreter', 'none');
else
    gui_CaAnalysis.legend{trace_num} = [Condition, '; n = ', num2str(length(va))];
    l1 = legend(gui_CaAnalysis.Curve(trace_num), [Condition, ', n = ', num2str(length(va))], 'Location', 'Northwest', 'boxoff');
    set(l1, 'Interpreter', 'none');
end

if isfield(gui_CaAnalysis, 'File_n')
    gui_CaAnalysis.File_n{trace_num} = length(va);
else
    gui_CaAnalysis.File_n{trace_num} = length(va);
end

if get(gui_CaAnalysis.figure.handles.AutoFigure_CheckBox, 'Value')
    if isfield(gui_CaAnalysis, 'AutoFig')
        if ~isempty(gui_CaAnalysis.AutoFig)
            AutoAx = get(gui_CaAnalysis.AutoFig, 'Children');
            axes(AutoAx);
            p = patch([Time fliplr(Time)], [spine_SEM_upper', fliplr(spine_SEM_lower')], colork{gui_CaAnalysis.color_counter});
            set(p, 'EdgeColor', 'none');
%             set(p, 'FaceAlpha', 0.5);
            plot(Time,stim_spine_average, 'color', colorj{gui_CaAnalysis.color_counter}, 'Linewidth', 2);
            k = ones(1,length(Time));
            plot(Time, k, '--k');
        else
            gui_CaAnalysis.AutoFig = figure; hold on;
            p = patch([Time fliplr(Time)], [spine_SEM_upper', fliplr(spine_SEM_lower')], colork{gui_CaAnalysis.color_counter});
            set(p, 'EdgeColor', 'none');
%             set(p, 'FaceAlpha', 0.5);
            plot(Time,stim_spine_average, 'color', colorj{gui_CaAnalysis.color_counter}, 'Linewidth', 2);
            k = ones(1,length(Time));
            plot(Time, k, '--k');
        end
    else
        gui_CaAnalysis.AutoFig = figure; hold on;
        p = patch([Time fliplr(Time)], [spine_SEM_upper', fliplr(spine_SEM_lower')], colork{gui_CaAnalysis.color_counter});
        set(p, 'EdgeColor', 'none');
%         set(p, 'FaceAlpha', 0.5);
        plot(Time,stim_spine_average, 'color', colorj{gui_CaAnalysis.color_counter}, 'Linewidth', 2);
        k = ones(1,length(Time));
        plot(Time, k, '--k');
    end
end


% plot(Time,spine_SEM_upper)
% plot(Time,spine_SEM_lower)
xlabel('Time(sec)');
ylabel('\DeltaF');