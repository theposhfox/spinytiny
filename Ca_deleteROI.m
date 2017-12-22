function Ca_deleteROI

global gui_CaImageViewer

twochannels = get(gui_CaImageViewer.figure.handles.TwoChannels_CheckBox, 'Value');

RoiN_tag = get(gco, 'Tag');
RoiN = RoiN_tag(4);

Roi_Rect_tag = regexp(RoiN_tag, ' Text', 'split');
if isempty(Roi_Rect_tag{2})
    Roi_Rect_tag = Roi_Rect_tag(1);
end

if twochannels == 1
    if ~isempty(strfind(Roi_Rect_tag{1}, 'red'))
        other_graph = regexp(Roi_Rect_tag{1}, 'red', 'split');
        other_graph = [other_graph{1}, other_graph{2}];
        Roi_Rect_tag{2} = other_graph;
        Other_RoiN_tag = [other_graph, ' Text'];
    elseif isempty(strfind(Roi_Rect_tag{1}, 'red'))
        other_graph = [Roi_Rect_tag{1}(1:3), 'red', Roi_Rect_tag{1}(4)];
        Roi_Rect_tag{2} = other_graph;
        Other_RoiN_tag = [other_graph, ' Text'];
    end
else
    other_graph = 'No red graph';
    Other_RoiN_tag = 'No red channel';
end

for i = 1:length(Roi_Rect_tag)
    ROIs(i) = findobj('Tag', Roi_Rect_tag{i});
end

Texts  = findobj('Tag', RoiN_tag);
Other_Texts = findobj('Tag', Other_RoiN_tag);

if ishandle(gui_CaImageViewer.ROI(str2num(RoiN)+1));
    ROI_data = get(gui_CaImageViewer.ROI(str2num(RoiN)+1));
    if twochannels == 1
        ROIred_data = get(gui_CaImageViewer.ROIred(str2num(RoiN)+1));
    else
    end
else
    ROI_data = 'empty';
    if twochannels == 1
        ROIred_data = 'empty';
    else
    end
end

for i = 1:length(Texts)
    delete(Texts(i));
    if twochannels == 1
        delete(Other_Texts(i));
    else
    end
end

for i = 1:length(ROIs)
    delete(ROIs(i));
end

if ~strcmpi(ROI_data, 'empty')
    gui_CaImageViewer.ROI = []
%     if twochannels == 1
%          delete(ROIred_data);
%     else
%     end
end

