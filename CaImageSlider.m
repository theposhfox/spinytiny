function CaImageSlider(ImageNum)

global gui_CaImageViewer

ImageNum = get(gui_CaImageViewer.figure.handles.ImageSlider_Slider, 'Value');

twochannels = get(gui_CaImageViewer.figure.handles.TwoChannels_CheckBox, 'Value');

set(gui_CaImageViewer.figure.handles.MaxProjection_CheckBox, 'Value', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if ImageNum > gui_CaImageViewer.imageserieslength
    ImageNum = gui_CaImageViewer.imageserieslength;
    set(gui_CaImageViewer.figure.handles.ImageSlider_Slider, 'Value', gui_CaImageViewer.imageserieslength);
elseif ImageNum < 1
    ImageNum = 1;
    set(gui_CaImageViewer.figure.handles.ImageSlider_Slider, 'Value', 1);
end


%%%% Create ROI stamps if ROIs exist %%%%

% [ROI_stamp, coordinates] = CaCreateROIstamps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Modify and filter the new frame like the previous one(s) %%%

filterwindow = str2num(get(gui_CaImageViewer.figure.handles.SmoothingFactor_EditableText, 'String'));
if ~isnumeric(filterwindow)
    filterwindow = 1;
end


if ~isinteger(ImageNum)
    ImageNum = ceil(ImageNum);
end



if filterwindow == 1;
    
    channel1 = gui_CaImageViewer.GCaMP_Image{ImageNum};
    if twochannels == 1
        channel2 = gui_CaImageViewer.Red_Image{ImageNum};
    else
        channel2 = [];
    end
    
    CommandSource = 'Slider';
    
    %%%%%%%%%
    PlaceImages(channel1,channel2, CommandSource);
    %%%%%%%%%
    
else
    smoothing_green = filter2(ones(filterwindow, filterwindow)/filterwindow^2, gui_CaImageViewer.GCaMP_Image{ImageNum});
    channel1 = smoothing_green;
    if twochannels == 1
        smoothing_red = filter2(ones(filterwindow, filterwindow)/filterwindow^2, gui_CaImageViewer.Red_Image{ImageNum});
        channel2 = smoothing_red;
    else
        channel2 = [];
    end

    CommandSource = 'Slider';

    %%%%%%%%%
    PlaceImages(channel1,channel2, CommandSource);
    %%%%%%%%%
end


%%% Place all existing ROIs on the new frame %%%

% PlaceROIs(ROI_stamp, coordinates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gui_CaImageViewer.figure.handles.Frame_EditableText, 'String', num2str(ImageNum));
set(gui_CaImageViewer.figure.handles.Frame_EditableText, 'String', ImageNum);
