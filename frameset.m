function frameset(hObject, eventdata)

if strcmpi(eventdata.Key, 'return')
    
    global gui_CaImageViewer
    
    aveproj = get(gui_CaImageViewer.figure.handles.AveProjection_CheckBox, 'Value');
    maxproj = get(gui_CaImageViewer.figure.handles.MaxProjection_CheckBox, 'Value');
% 
    ImageNum = str2num(get(gui_CaImageViewer.figure.handles.Frame_EditableText, 'String'));
% 
    set(gui_CaImageViewer.figure.handles.ImageSlider_Slider, 'Value', ImageNum);

    if aveproj || maxproj
        ch1image = gui_CaImageViewer.ch1image;
    else
        ch1image = gui_CaImageViewer.GCaMP_Image{ImageNum};
    end
    
    PlaceImages(ch1image, [], 'Slider');
else
end