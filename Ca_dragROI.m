function DragROI(hObject, eventdata, ROInum)

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

point1 = get(gca, 'CurrentPoint');  %%% Button down detected and position parameters stored (x, y, width, height)
point1 = point1(1,1:2)   ;          %%% Extract x and y

RoiRect = get(gco, 'Position') ;    %%% Get placement of current object on screen

set(gcf, 'Units', 'pixels');
rectFig = get(gcf, 'Position');     %%% Get position of current figure on the screen

set(gca, 'Units', 'pixels');
rectAx = get(gca, 'Position');      %%% Get position of current axes on the screenddd

scrnsz = get(0,'ScreenSize');

xmag = (rectFig(3)*rectAx(3))/scrnsz(3);  %pixel/screenpixel
xoffset =rectAx(1)*rectFig(3);
ymag = (rectFig(4)*rectAx(4))/scrnsz(2);
yoffset = rectAx(2)*rectFig(4);
rect1 = [xmag*RoiRect(1)+xoffset+0.5, ymag*(scrnsz(2)-RoiRect(2)-RoiRect(4))+yoffset+.5, xmag*RoiRect(3), ymag*RoiRect(4)];

if point1(1)>(RoiRect(1)+2*RoiRect(3)/3) && point1(2)> (RoiRect(2)+2*RoiRect(4)/3) %%%resize
%     fixedpoint = [RoiRect(1)+rectAx(1), rectAx(2)+rectAx(4)-RoiRect(2)]
%     rbbox([RoiRect(1)+rectAx(1), rectAx(2)+rectAx(4)-RoiRect(2)+RoiRect(4), RoiRect(3),RoiRect(4)], fixedpoint);
    fixedpoint = [rect1(1), rect1(2)+rect1(4)];
    rbbox(rect1, fixedpoint);
    point2 = get(gca, 'CurrentPoint');
    point2 = point2(1,1:2);
    RoiRect(3) = point2(1)-RoiRect(1);
    RoiRect(4) = point2(2)-RoiRect(2);
    RoiRect_final = [RoiRect(1), RoiRect(2), point2(1)-RoiRect(1), point2(2)-RoiRect(2)];
else %%% drag
    rbbox;      %%% initiate function to detect transition from button press to button release
    point2 = get(gca, 'CurrentPoint');
    point2 = point2(1,1:2);
    RoiRect_new = [point2(1)-(RoiRect(3)/2), point2(2)-RoiRect(4), RoiRect(3), RoiRect(4)];
    RoiRect_final = dragrect(RoiRect_new);
end

Roi_Rect_tag = get(gco, 'Tag');
RoiN = Roi_Rect_tag(end);

ROItext = findobj('Type', 'Text', 'Tag', ['ROI', num2str(ROInum), ' Text']);
ROIredtext = findobj('Type', 'Text', 'Tag', ['ROIred', num2str(ROInum), ' Text']);
set(ROItext, 'Position', [RoiRect_final(1), RoiRect_final(2), 0]);
% set(ROItext, 'Position', [RoiRect_final(1), RoiRect_final(2), 0]);

set(gco, 'Position', RoiRect_final);

if ~isempty(strfind(Roi_Rect_tag, 'red'))
    other_graph = regexp(Roi_Rect_tag, 'red', 'split');
    other_graph = [other_graph{1}, other_graph{2}];
    Other_Rect_tag = other_graph;
    Other_RoiN_tag = [other_graph, ' Text'];
elseif isempty(strfind(Roi_Rect_tag, 'red'))
    other_graph = [Roi_Rect_tag(1:3), 'red', Roi_Rect_tag(4)];
    Other_Rect_tag = other_graph;
    Other_RoiN_tag = [other_graph, ' Text'];
end

Translate_other = RoiRect_final - RoiRect;

OtherROI = findobj('Type', 'Rectangle', 'Tag', Other_Rect_tag);
if length(OtherROI)>1 
    OtherROI = OtherROI(1);
end
OtherRect = get(OtherROI, 'Position');

Other_Texts = findobj('Type', 'Text', 'Tag', Other_RoiN_tag);
if length(Other_Texts)>1;
    Other_Texts = Other_Texts(1);
end

OtherFinal = [OtherRect(1)+Translate_other(1), OtherRect(2)+Translate_other(2), RoiRect(3), RoiRect(4)];
set(OtherROI, 'Position', OtherFinal);
set(Other_Texts, 'Position', [OtherFinal(1)-2, OtherFinal(2)-2]);
