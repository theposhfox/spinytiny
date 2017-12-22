function CaAnalysisBarGraph

global gui_CaAnalysis;

Graph_Data = get(gui_CaAnalysis.Curve);

num_samples = length(gui_CaAnalysis.Curve);
counter = 1;

for i = 1:length(Graph_Data)      
    Yvalue{counter} = (Graph_Data(i).YData);
    counter = counter+1;
end

for i = 1:num_samples
    Y1(1,i) = mean(Yvalue{i}(475:500));
    YSEM(1,i) = (std(Yvalue{i}(475:500)))/sqrt(gui_CaAnalysis.File_n{i});
end

figure; barweb(Y1, YSEM);