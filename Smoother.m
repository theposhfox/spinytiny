function Smoother(hObject, eventdata, channel1image, channel2image, CommandSource);

global gui_CaImageViewer

set(gui_CaImageViewer.figure.handles.MaxProjection_CheckBox, 'Value', 0)

%%%% Create ROI stamps if ROIs exist %%%%

[ROI_stamp, coordinates] = CaCreateROIstamps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi('CommandSource', 'Loader')

    ImageNum = get(gui_CaImageViewer.figure.handles.ImageSlider_Slider, 'Value');
    if ~isinteger(ImageNum)
        ImageNum = ceil(ImageNum);
    end
    filterwindow = str2num(get(gui_CaImageViewer.figure.handles.SmoothingFactor_EditableText, 'String'));

    smoothing_green = filter2(ones(filterwindow, filterwindow)/filterwindow^2, channel1image(ImageNum));
    channel1 = smoothing_green;
    smoothing_red = filter2(ones(filterwindow, filterwindow)/filterwindow^2, channel2image(ImageNum));
    channel2 = smoothing_red;

    if ~isinteger(ImageNum)
        ImageNum = ceil(ImageNum);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Place Images with modifications (if any) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    CommandSource = 'Smoother';
    PlaceImages(channel1,channel2, CommandSource);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Place all existing ROIs on the new frame %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    PlaceROIs(ROI_stamp, coordinates);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

try eventdata.Key;    
    if strcmpi(eventdata.Key, 'return')
        ImageNum = get(gui_CaImageViewer.figure.handles.ImageSlider_Slider, 'Value');
        if ~isinteger(ImageNum)
            ImageNum = ceil(ImageNum);
        end
        filterwindow = str2num(get(gui_CaImageViewer.figure.handles.SmoothingFactor_EditableText, 'String'));

        if ~isinteger(ImageNum)
            ImageNum = ceil(ImageNum);
        end
        
        smoothing_green = filter2(ones(filterwindow, filterwindow)/filterwindow^2, gui_CaImageViewer.GCaMP_Image{ImageNum});
        channel1 = smoothing_green;
        smoothing_red = filter2(ones(filterwindow, filterwindow)/filterwindow^2, gui_CaImageViewer.Red_Image{ImageNum});
        channel2 = smoothing_red;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Place Images with modifications (if any) %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CommandSource = 'Smoother';
        PlaceImages(channel1,channel2, CommandSource);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%

        %%% Place all existing ROIs on the new frame %%%

        PlaceROIs(ROI_stamp, coordinates);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
catch
    return
end
