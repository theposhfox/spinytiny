function AdjustLUT(hObject, eventdata)

if strcmpi(eventdata.Key, 'return')
    global gui_CaImageViewer
    
    pause(0.1);
    
    NewUpper = str2num(get(gui_CaImageViewer.figure.handles.UpperLUT_EditableText, 'String'));
    NewLower = str2num(get(gui_CaImageViewer.figure.handles.LowerLUT_EditableText, 'String'));
    
    
    ImageNum = str2num(get(gui_CaImageViewer.figure.handles.Frame_EditableText, 'String'));

    CaImageSlider(ImageNum);