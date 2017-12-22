function PlaceROIs(ROI_stamp, coordinates);

%%% Place ROIs %%%

% This script will place an ROI on the newly formed image when parameters
% of the image are changed. The positions of the ROIs are stored, the image
% redrawn anew with the new parameters, and the old ROIs are drawn on the
% new image using their stored positions.

program = get(gcf);

running = program.FileName;

if ~isempty(regexp(running, 'CaImageViewer'))
    global gui_CaImageViewer
    glovar = gui_CaImageViewer;
    axes1 = glovar.figure.handles.GreenGraph;
    axes2 = glovar.figure.handles.RedGraph;
elseif ~isempty(regexp(running, 'FluorescenceSuite'));
    global gui_FluorescenceSuite
    glovar = gui_FluorescenceSuite;
    axes1 = glovar.figure.handles.Green_Axes;
    axes2 = glovar.figure.handles.Red_Axes;
end

twochannels = get(glovar.figure.handles.TwoChannels_CheckBox, 'Value');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Place all existing ROIs on the new frame %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Delete all regular ROIs that aren't for PolyLines drawing
a = findobj('Type', 'Rect', '-and', '-not', {'Tag', 'PolyROI', '-or', 'Tag', 'RedPolyROI'});
for i = 1:length(a)-length(ROI_stamp)
    delete(a(i));
end
b = findobj('Type', 'Text');
for i = 1:length(b)-length(ROI_stamp)
    delete(b(i));
end

for a = 1:length(ROI_stamp)
    if ~isempty(ROI_stamp{a})
        ROInum = a-1;
        axes(axes1);
        glovar.ROI(a) = rectangle('Position', ROI_stamp{a}, 'EdgeColor', 'green', 'Curvature', [1 1],'Tag', ['ROI', num2str(ROInum)], 'ButtonDownFcn', {@DragROI, ROInum});
        glovar.Roitext(a) = text(ROI_stamp{a}(1)-2, ROI_stamp{a}(2)-2, num2str(a-1), 'color', 'white', 'Tag', ['ROI', num2str(a-1), ' Text'],'ButtonDownFcn', 'DeleteROI');
        if twochannels == 1
            axes(axes2);
            glovar.ROIred(a) = rectangle('Position', ROI_stamp{a}, 'EdgeColor', 'red', 'Curvature', [1 1],'Tag', ['ROIred', num2str(ROInum)],'ButtonDownFcn', {@DragROI, ROInum});
            glovar.ROIredtext(a) = text(ROI_stamp{a}(1)-2, ROI_stamp{a}(2)-2, num2str(a-1), 'color', 'white', 'Tag', ['ROI', num2str(a-1), ' Text'],'ButtonDownFcn', 'DeleteROI');
        else
        end
    end
end

radius = 5;
x = [];
x_red = [];
y = [];
y_red = [];


DendNum = glovar.Dendrite_Number;

if isfield(glovar, 'PolyROI') && ~isempty(coordinates)
    PPsperDend = glovar.DendritePolyPointNumber;
    if ~isempty(coordinates{1})
    axes(axes1);
    currDend = 1;
        for i = 1:length(coordinates)
            glovar.PolyLinePos{i} = [coordinates{i}(1), coordinates{i}(2), radius*2, radius*2];
            glovar.PolyROI{i} = rectangle('Position', glovar.PolyLinePos{i}, 'EdgeColor', 'red', 'Tag', ['Dendrite ', num2str(currDend), ' PolyROI ', num2str(i)], 'Curvature', [1 1], 'ButtonDownFcn', 'Drag_Poly');
            x = [x,coordinates{i}(1)+radius];
            y = [y,coordinates{i}(2)+radius];
            if i < sum(PPsperDend(1:currDend))
                currDend = 1;
            else
                currDend = currDend+1;
            end
        end
        if DendNum == 1
            for i = 1:2:length(glovar.PolyLinePos)
                hold on;
                plot(x,y, 'color', 'cyan', 'Tag', 'PolyLine');
            end
        else
            counter = 1;
            for i = 1:DendNum
                hold on;
                plot(x(counter:(counter+PPsperDend(i)-1)),y(counter:(counter+PPsperDend(i)-1)), 'color', 'cyan', 'Tag', 'PolyLine');
                counter = sum(PPsperDend(1:i))+1;
            end
        end
    if twochannels == 1
        axes(axes2)
        for i = 1:length(glovar.RedPolyLinePos)
            glovar.RedPolyLinePos{i} = [coordinates{i}(1)-radius, coordinates{i}(2)-radius, radius*2, radius*2];
            glovar.RedPolyROI{i} = rectangle('Position', glovar.RedPolyLinePos{i}, 'EdgeColor', 'cyan', 'Tag', ['Dendrite ', num2str(DendriteNum), ' RedPolyROI', num2str(i)], 'Curvature', [1 1], 'ButtonDownFcn', 'Drag_Poly');
            x_red = [x_red,coordinates{i}(1)];
            y_red = [y_red,coordinates{i}(2)];
        end
        for i = 1:2:length(glovar.RedPolyLinePos)
            hold on;
            plot(x_red,y_red, 'color', 'cyan', 'Tag', 'PolyLine');
        end
    else
    end
    end
end

 %%% Overwrite the previous existing global workspace with the newly imprinted one
 
if ~isempty(regexp(running, 'CaImageViewer'))
    gui_CaImageViewer = glovar;
elseif ~isempty(regexp(running, 'FluorescenceSuite'));
    gui_FluorescenceSuite = glovar;
end

