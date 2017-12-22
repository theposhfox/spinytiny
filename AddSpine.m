function AddSpine(eventdata)

global gui_SpineSegmentation

clicktype = get(gcbf, 'SelectionType');

if strcmp(clicktype, 'normal')
    cent_click = [];
    cent_click = get(gca, 'CurrentPoint');
    cent_click = cent_click(1,1:2);
    cent_click = round(cent_click);
    rectangle('Position', [cent_click, 2, 2], 'Curvature', [1 1], 'EdgeColor', 'b', 'Tag', ['Branchpoint ', num2str(size(gui_SpineSegmentation.Spines,1)+1)], 'ButtonDownFcn', 'DeleteROI')
    gui_SpineSegmentation.Spines = [gui_SpineSegmentation.Spines; cent_click]
else
end


