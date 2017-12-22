function DeleteROI

global gui_SpineSegmentation

clicktype = get(gcbf, 'SelectionType')

if strmp(clicktype, 'normal')
    return
else
    BP_tag = get(gco, 'Tag')
    BP_num = BP_tag(13:end);
    ROI = findobj('Tag', BP_tag);
    delete(ROI)
end