function AdjustROI
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


global gui_Microglia

set(gca, 'Units', 'Pixels');
set(gcf, 'Units', 'Pixels');
set(0, 'Units', 'Pixels');

axes_position = get(gca, 'Position')

point1 = get(gca, 'CurrentPoint')

ROI = get(gco, 'Position')

ROI_upperleft = [ROI(1), ROI(2)+ROI(4)];
fixed_point = ROI_upperleft;

newbox = rbbox([ROI(1)+axes_position(1), ROI(2)+axes_position(2), ROI(3), ROI(4)], fixed_point)

point2 = get(gca, 'CurrentPoint')


point1 = point1(1,1:2);
point2 = point2(1,1:2);

point = min(fixed_point,point2)
offset = abs(fixed_point-point2)

newROI = round([point,offset]);

finalbox = [point1(1)-axes_position(1), point2(2)-axes_position(2), point2(1)-point1(1), point2(2)-point1(2)]

program = get(gcf);

oldROIs = findobj(program.Children, 'Tag', 'ROI');
delete(oldROIs);

axes(gui_Microglia.handles.MicrogliaImage_Axes);

gui_Microglia.FrameofInt = rectangle('Position', newROI, 'Tag', 'ROI', 'EdgeColor', 'Cyan', 'ButtonDownFcn', 'AdjustROI')

end

