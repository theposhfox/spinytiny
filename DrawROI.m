function DrawROI(hObject, eventdata, ROInum, Router)

%%% This function is designed to work through the "CaImageViewer" main
%%% function. After initiating the drawing procedure (clicking in the image
%%% window after selecting "Add Spine"), this function will first draw a
%%% box in the general area of choice, then bring up a new window
%%% ("Fine-Select") that allows careful drawing of the elliptical ROI at an
%%% automatically zoomed-in scale (note: this does not occur for the
%%% background ROI). The ROI can be dragged by clicking and dragging the
%%% bottom of the ellipse, and can be resized by selecting the bottom-right
%%% corner of the encasing (white-dashed) rectangle. 
%%% The resulting window also gives you the option of drawing a "surround
%%% background" ROI, or a concentric ring around the actual ROI that
%%% subtracts overlaying contamination from other non-targeted signals. 
%%% When ROI drawing is complete in the "Fine-Select" window, right-click
%%% on the ellipitical ROI to signal termination, at which point the ROI
%%% coordinates will be passed to the function "DragROI", where the final
%%% drawing into the image axes is handled. 

program = get(gcf);

running = program.FileName;
% 
% if isempty(running)
%     running = 'CaImageViewer';
% end

if ~isempty(regexp(running, 'CaImageViewer'))
    global gui_CaImageViewer
    glovar = gui_CaImageViewer;
    axes1 = glovar.figure.handles.GreenGraph;
    axes2 = glovar.figure.handles.RedGraph;
    twochannels = glovar.figure.handles.TwoChannels_CheckBox;
    zoom = str2num(get(glovar.figure.handles.Zoom_EditableText, 'String'));
    if zoom > 10 
        radius = 30;
    else
        radius = 20;
    end
elseif ~isempty(regexp(running, 'FluorescenceSuite'));
    global gui_FluorescenceSuite
    glovar = gui_FluorescenceSuite;
    axes1 = glovar.figure.handles.Green_Axes;
    axes2 = glovar.figure.handles.Red_Axes;
    twochannels = [];
    radius = 15;
else
    global sideline
    global gui_CaImageViewer;
    glovar = gui_CaImageViewer;
    axes1 = gca;
    axes2 = [];
    twochannels = [];
    radius = 10;
end

if ~isnumeric(ROInum)
    ROInum = str2double(ROInum);
end

if strcmpi(Router, 'FineSelect')
    edgewidth = 2;
    textsize = 12;
    set(gcf, 'WindowButtonDownFcn', [])
else
    edgewidth = 1;
    textsize = 10;
end

cmap = glovar.CurrentCMap; 

if strcmpi(cmap, 'RGB')
    linecolor = 'b';
    cmap = [zeros(length([0:0.01:1]),1),[0:0.01:1]', zeros(length([0:0.01:1]),1)];
elseif strcmpi(cmap, 'Jet')
    linecolor = 'w';
elseif strcmpi(cmap, 'Hot')
    linecolor = 'r';
elseif strcmpi(cmap, 'Fire')
    linecolor = 'g'; %when this linecolor changed, the linecolor of zoomin window changed, ZL  comment
end

% Backgroundtoggle = get(glovar.figure.handles.BackgroundROI_ToggleButton, 'Value');
% StimSpinetoggle = get(glovar.figure.handles.SpineROI_ToggleButton, 'Value');
% Nearbytoggle = get(glovar.figure.handles.NearbySpine_ToggleButton, 'Value');

% waitforbuttonpress;
pointer_location1 = get(gca, 'CurrentPoint');

% Use the following to allow single-click selection of ROIs
point = pointer_location1(1,1:2);
Fl_ROI(1) = point(1)-radius/2; Fl_ROI(2) = point(2)-radius/2;
Fl_ROI(3) = radius; Fl_ROI(4) = radius;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Draw ROI (method depends on invoking function) %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmpi(Router, 'FineSelect') && ROInum == 0 %%% For drawing the general background ROI
    axes(axes1)
    delete(findobj('Tag', 'ROI confine'));
    delete(findobj('Tag', 'ROI0'));
    delete(findobj('Tag', 'BackgroundROI0'))
    delete(findobj('Tag', 'ROI0 Text'));
    glovar.ROI(ROInum+1) = rectangle('Position', Fl_ROI, 'EdgeColor', linecolor, 'Curvature', [1 1], 'Tag', ['ROI', num2str(ROInum)], 'ButtonDownFcn', {@DragROI, ROInum, 'ZoomWindow'}, 'Linewidth', edgewidth);
    glovar.ROItext(ROInum+1) = text(Fl_ROI(1)-2, Fl_ROI(2)-2, num2str(ROInum), 'color', 'white', 'Tag', ['ROI', num2str(ROInum), ' Text'], 'ButtonDownFcn', @DeleteROI, 'FontSize', textsize);
    glovar.BackgroundROI(ROInum+1) = NaN;
    if twochannels == 1;
        axes(axes2)
        glovar.ROIred(ROInum+1) = rectangle('Position', Fl_ROI, 'EdgeColor', 'red', 'Curvature', [1 1], 'Tag', ['ROIred', num2str(ROInum)],'ButtonDownFcn', {@DragROI, ROInum});
        glovar.ROIredtext(ROInum+1) = text(Fl_ROI(1)-2, Fl_ROI(2)-2, num2str(ROInum), 'color', 'white', 'Tag', ['ROIred', num2str(ROInum), ' Text'],'ButtonDownFcn', 'Ca_deleteROI');
    else
    end
elseif strcmpi(Router, 'Spine') %%% First round drawing, just temporary while the fine-select mode is initiated
    axes(axes1)
    InsertOpt = get(glovar.figure.handles.InsertSpine_ToggleButton,'Value');
    if InsertOpt
        ROInum = glovar.InsertPoint;
    end
    glovar.ROI(ROInum+1) = rectangle('Position', Fl_ROI, 'EdgeColor', linecolor, 'Curvature', [0 0], 'Tag', ['ROI', num2str(ROInum), ' Starter'], 'Linewidth', edgewidth);
    glovar.ROItext(ROInum+1) = text(Fl_ROI(1)-2, Fl_ROI(2)-2, num2str(ROInum), 'color', 'white', 'Tag', ['ROI', num2str(ROInum), ' Text Starter'], 'ButtonDownFcn', 'DeleteROI', 'FontSize', textsize);
    
    if twochannels == 1;
        axes(axes2)
        glovar.ROIred(ROInum+1) = rectangle('Position', Fl_ROI, 'EdgeColor', 'red', 'Curvature', [1 1], 'Tag', ['ROIred', num2str(ROInum)]);
        glovar.ROIredtext(ROInum+1) = text(Fl_ROI(1)-2, Fl_ROI(2)-2, num2str(ROInum), 'color', 'white', 'Tag', ['ROIred', num2str(ROInum), ' Text'],'ButtonDownFcn', 'DeleteROI');
    else
    end
elseif strcmpi(Router, 'FineSelect') && ROInum ~=0
        axes(axes1)
    if glovar.UsingSurroundBackground
        surroundoffset = glovar.SurroundBackgroundBuffer;
        glovar.BackgroundROI(ROInum+1) = rectangle('Position', [Fl_ROI(1)-surroundoffset/2, Fl_ROI(2)-surroundoffset/2, Fl_ROI(3)+surroundoffset, Fl_ROI(4)+surroundoffset], 'EdgeColor', 'w', 'Curvature', [1 1], 'Tag', ['BackgroundROI', num2str(ROInum)], 'Linewidth', 0.75);
    else
        glovar.BackgroundROI(ROInum+1) = NaN;
    end
    glovar.ROI(ROInum+1) = rectangle('Position', Fl_ROI, 'EdgeColor', linecolor, 'Curvature', [1 1], 'Tag', ['ROI', num2str(ROInum)], 'ButtonDownFcn', {@DragROI, ROInum, 'ZoomWindow'}, 'Linewidth', edgewidth);
    rectangle('Position', Fl_ROI, 'EdgeColor', 'w', 'Curvature', [0 0], 'Tag', 'ROI confine','Linewidth', 1, 'LineStyle', ':');
    uistack(glovar.ROI(ROInum+1), 'top');
    glovar.ROItext(ROInum+1) = text(Fl_ROI(1)-2, Fl_ROI(2)-2, num2str(ROInum), 'color', 'white', 'Tag', ['ROI', num2str(ROInum), ' Text'], 'ButtonDownFcn', @DeleteROI, 'FontSize', textsize);

    if twochannels == 1;
        axes(axes2)
        glovar.ROIred(ROInum+1) = rectangle('Position', Fl_ROI, 'EdgeColor', 'red', 'Curvature', [1 1], 'Tag', ['ROIred', num2str(ROInum)],'ButtonDownFcn', {@DragROI, ROInum, 'ZoomWindow'});
        glovar.ROIredtext(ROInum+1) = text(Fl_ROI(1), Fl_ROI(2), num2str(ROInum), 'color', 'white', 'Tag', ['ROIred', num2str(ROInum), ' Text'],'ButtonDownFcn', 'Ca_deleteROI');
    else
    end
elseif strcmpi(Router, 'Other')      
    set(glovar.figure.handles.OtherFeature_ToggleButton, 'Value', 0)
    set(glovar.figure.handles.output, 'WindowButtonDownFcn', [])
    oROI = impoly(axes1);
    selectedPos = getPosition(oROI);
    glovar.ROI(ROInum+1) = NaN;
    glovar.OtherFeature_ROIpos{ROInum+1} = selectedPos;
    selectedPos = [selectedPos; selectedPos(1,:)];
    delete(oROI)
    glovar.OtherFeature_ROI = line(selectedPos(:,1), selectedPos(:,2), 'Tag', ['ROI', num2str(ROInum)]);
    glovar.ROItext(ROInum+1) = text(selectedPos(1,1)-2, selectedPos(1,2)-2, num2str(ROInum), 'color', 'white', 'Tag', ['ROI', num2str(ROInum), ' Text'], 'ButtonDownFcn', @DeleteROI, 'FontSize', textsize);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%
% Create temporary window for automatic zoom into ROI area
%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(Router, 'Background')
    set(glovar.figure.handles.BackgroundROI_ToggleButton, 'Value', 0)
elseif strcmpi(Router, 'Spine')
    set(glovar.figure.handles.SpineROI_ToggleButton, 'Value', 0)
    
    r = round([Fl_ROI(2), Fl_ROI(2), Fl_ROI(2)+Fl_ROI(4), Fl_ROI(2)+Fl_ROI(4)]); %%% ROI Y values, to be used as row index for sub-image (indices are: 1) top-left corner, 2) top-right corner, 3) bottom-right, 4)bottom-left
    c = round([Fl_ROI(1), Fl_ROI(1)+Fl_ROI(3), Fl_ROI(1)+Fl_ROI(3), Fl_ROI(1)]); %%% ROI X values
    
    if ~glovar.NewSpineAnalysis
        im = gui_CaImageViewer.GCaMP_Image;
        im = cat(3, im{:});
        immax = max(im, [], 3);
    else
        immax = gui_CaImageViewer.ch1image;
    end
    poschange = 0;
    if r(1)<=0
        r(1) = 1;
        r(2) = 1;
        poschange = 1;
    elseif r(1)>length(immax)
        r(1) = length(immax);
        r(2) = length(immax);
        poschange = 1;
    else
    end
    if r(3)<=0
        r(3) = 1;
        r(4) = 1;
        poschange = 1;
    elseif r(3)>length(immax)
        r(3) = length(immax);
        r(4) = length(immax);
        poschange = 1;
    else
    end
    if c(1)<=0
        c(1) = 1;
        c(4) = 1;
        poschange = 1;
    elseif c(1)>length(immax)
        c(1) = length(immax);
        c(4) = length(immax);
        poschange = 1;
    else
    end
    if c(2) <=0
        c(2) = 1;
        c(3) = 1;
        poschange = 1;
    elseif c(2)>length(immax)
        c(2) = length(immax);
        c(3) = length(immax);
        poschange = 1;
    else
    end
    if poschange 
        set(glovar.ROI(ROInum+1), 'Position', [c(1), r(1), abs(c(3)-c(4)), abs(r(4)-r(1))]); %%% 
        set(glovar.ROItext(ROInum+1),'color', 'red');
    end
    
    zoomim = immax(r(1):r(3), c(1):c(2));
    filterwindow = str2num(get(gui_CaImageViewer.figure.handles.SmoothingFactor_EditableText, 'String'));
    padvalue = 5;
    zoomim = filter2(ones(filterwindow, filterwindow)/filterwindow^2, padarray(zoomim, [padvalue padvalue], nanmedian(nanmedian(zoomim))));
    zoomim = zoomim(padvalue+1:end-padvalue, padvalue+1:end-padvalue);
    
    zi = figure('Name', 'Auto zoom window', 'WindowButtonDownFcn', {@DrawROI, ROInum, 'FineSelect'}, 'NumberTitle', 'off');
    bgcol = get(zi, 'Color');
    usingsurround = gui_CaImageViewer.UsingSurroundBackground;
    uicontrol('Parent', zi, 'Style', 'Checkbox', 'Position', [15 15 200 30], 'String', 'Draw Surround Background?', 'BackgroundColor',bgcol,'Value', usingsurround, 'Callback', {@SurroundToggle, ROInum});
    uicontrol('Parent', zi, 'Style', 'text', 'Position', [350 10 80 30], 'String', 'Added Diameter:', 'BackgroundColor', bgcol)
    uicontrol('Parent', zi, 'Style', 'edit', 'Position', [475 15 50 30], 'String', num2str(gui_CaImageViewer.SurroundBackgroundBuffer), 'KeyPressFcn', {@SurroundStretch, ROInum})
    imagesc(zoomim); colormap(cmap)
    title(['Fine outline of ROI ', num2str(ROInum)])
    set(gca, 'XTick', [], 'YTick', [])
    
    if InsertOpt
        glovar.Spine_Number = sum(cell2mat(cellfun(@length, gui_CaImageViewer.SpineDendriteGrouping, 'uni', false)));
    else
        glovar.Spine_Number = ROInum;
    end
    
    global sideline
    
elseif strcmpi(Router, 'Other')
    
elseif strcmpi(Router, 'FineSelect')
    
end


if ~isempty(regexp(running, 'CaImageViewer'))
    gui_CaImageViewer = glovar;
    set(gui_CaImageViewer.figure.handles.output, 'WindowButtonDownFcn', []);
elseif ~isempty(regexp(running, 'FluorescenceSuite'));
    gui_FluorescenceSuite = glovar;
    set(gui_FluorescenceSuite.figure.handles.output, 'WindowButtonDownFcn', []);
end



function SurroundToggle(hObject, eventdata, ROInum)

global gui_CaImageViewer

toggleval = get(hObject, 'Value');
isdrawn = ~isempty(findobj(gcf, 'Type', 'rectangle', '-and', {'-regexp', 'Tag', 'ROI'}));

if toggleval
    gui_CaImageViewer.UsingSurroundBackground = 1;
    if isdrawn
        surroundoffset = gui_CaImageViewer.SurroundBackgroundBuffer;
        Fl_ROI = get(findobj(gcf,'Type', 'Rectangle', '-and', {'-regexp', 'Tag', 'ROI'}, '-and', 'Curvature', [1 1]), 'Position');
        gui_CaImageViewer.BackgroundROI(ROInum+1) = rectangle('Position', [Fl_ROI(1)-surroundoffset/2, Fl_ROI(2)-surroundoffset/2, Fl_ROI(3)+surroundoffset, Fl_ROI(4)+surroundoffset], 'EdgeColor', 'w', 'Curvature', [1 1], 'Tag', ['BackgroundROI', num2str(ROInum)], 'Linewidth', 0.75);
    else
    end
else
%     delete(zi, findobj('Type', 'rectangle', '-and', {'-regexp', 'Tag', 'Background'})
    gui_CaImageViewer.UsingSurroundBackground = 0;
    if isdrawn
        gui_CaImageViewer.BackgroundROI(ROInum+1) = NaN;
        delete(findobj(gcf, 'Type', 'Rectangle', '-and', {'-regexp', 'Tag', 'Background'}))
    else
    end
end


function SurroundStretch(hObject, eventdata, ROInum)

if strcmpi(eventdata.Key, 'return')
    global gui_CaImageViewer
    isdrawn = ~isempty(findobj(gcf, 'Type', 'rectangle', '-and', {'-regexp', 'Tag', 'ROI'}));
    gui_CaImageViewer.SurroundBackgroundBuffer = str2num(get(hObject, 'String'));
    if isdrawn
        Fl_ROI = get(findobj(gcf,'Type', 'Rectangle', '-and', 'Tag', ['ROI', num2str(ROInum)], '-and', 'Curvature', [1 1]), 'Position');    %%% Based your modification of the surrounding background ROI on the original ROI so that errors don't accumulate
        oldBG = findobj(gcf, 'Type', 'rectangle', '-and', 'Tag', ['BackgroundROI', num2str(ROInum)]);
        surroundoffset = gui_CaImageViewer.SurroundBackgroundBuffer;
        set(oldBG, 'Position', [Fl_ROI(1)-surroundoffset/2, Fl_ROI(2)-surroundoffset/2, Fl_ROI(3)+surroundoffset, Fl_ROI(4)+surroundoffset]);
    end
end