function CaDataAverage(Subject_File, Condition, Excluder);

load(Subject_File);

global gui_CaAnalysis

WS_Data = who(['*',Condition, '*']);
WS_Data2 = [];

Compare_cond = strfind(Condition, 'Compare');

selected_files = 1;

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
    
file_num = 0;

for i = 1:length(WS_Data2)
    file_num = file_num+1;
    evalc(['va{', num2str(file_num), '}=', WS_Data2{i}]);
end

for i = 1:length(va)
    spine_data(:,i) = (va{i}.raw(:,4)/mean(va{i}.raw(1:32,4)))./va{i}.red(:,4);
end

stim_spine_average = mean(spine_data, 2);

stim_spine_SEM = (std(spine_data,0,2))/sqrt(size(spine_data,2));

spine_SEM_upper = stim_spine_average+stim_spine_SEM;
spine_SEM_lower = stim_spine_average-stim_spine_SEM;

Time = 0:0.125:70;
Time = Time-4;
Time = Time(1:560);


%%%% Graph %%%%

axes(gui_CaAnalysis.figure.handles.Graph); hold on;

fig_axes = get(gca);

    %%%%%%%%%%%%%%%%%%%%%%
    %%Color Information%%%
    %%%%%%%%%%%%%%%%%%%%%%

    lgray = [0.50 0.51 0.52]; brown = [0.28 0.22 0.14];
    gray = [0.50 0.51 0.52]; brown = [0.59 0.45 0.28];
    yellow = [1.00 0.76 0.05]; orange = [0.95 0.40 0.13];
    lgreen = [0.55 0.78 0.25]; green = [0.00 0.43 0.23];
    lblue = [0.00 0.68 0.94]; blue = [0.00 0.33 0.65];
    magenta = [0.93 0.22 0.55]; purple = [0.57 0.15 0.56];
    red = [0.93 0.11 0.14]; black = [0 0 0];
    colorj = {red,lblue,green,lgreen,gray,brown,yellow,blue,purple,magenta,black};
    
plot_empty_check = fig_axes.Children;

if isempty(plot_empty_check)
    gui_CaAnalysis.color_counter = 1;
else
    gui_CaAnalysis.color_counter = gui_CaAnalysis.color_counter+1;
end

patch([Time fliplr(Time)], [spine_SEM_upper', fliplr(spine_SEM_lower')], colorj{gui_CaAnalysis.color_counter});

Patch_alpha = findobj('Type', 'patch');

set(Patch_alpha, 'FaceAlpha', 0.2);

hold on; plot(Time,stim_spine_average, 'color', colorj{gui_CaAnalysis.color_counter}, 'Linewidth', 2);


% plot(Time,spine_SEM_upper)
% plot(Time,spine_SEM_lower)
xlabel('Time(sec)');
ylabel('\DeltaF/F');