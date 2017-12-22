function DragROI(hObject, eventdata, ROInum, CurrentWindow)

%%% This function is generally designed to manipulate the eliptical ROIs
%%% created using the function DrawROI. 

program = get(gcf);

running = program.FileName;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize: set parameters based on the parent function

    global gui_CaImageViewer
    glovar = gui_CaImageViewer;
    axes1 = glovar.figure.handles.GreenGraph;
    axes2 = glovar.figure.handles.RedGraph;
    twochannels = get(glovar.figure.handles.TwoChannels_CheckBox, 'Value');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Determine desired behavior based on click type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clicktype = get(gcbf, 'SelectionType');

if strcmpi(clicktype, 'alt') %%% This is the terminal call for the "Fine Select" portion of the draw ROI program
    if strcmp(CurrentWindow, 'ZoomWindow')
        newROI = findobj(gcf, 'Type', 'Rectangle', 'Tag', ['ROI', num2str(ROInum)]);
        newROIpos = round(get(newROI, 'Position'));
        oldROIpos = round(get(gui_CaImageViewer.ROI(ROInum+1), 'Position'));

        %%% Correction Factor
        correctionfactor = str2num(get(gui_CaImageViewer.figure.handles.ROIoffset_EditableText, 'String'));

        adjustedpos = [oldROIpos(1)+newROIpos(1)+correctionfactor, oldROIpos(2)+newROIpos(2)+correctionfactor, newROIpos(3), newROIpos(4)];
        try
            close('Auto zoom window');
        catch
        end
        axes(gui_CaImageViewer.figure.handles.GreenGraph)
        delete(findobj('Type', 'rectangle', '-and', 'Tag', ['ROI', num2str(ROInum)]))
        delete(findobj('Type', 'text', '-and', 'Tag', ['ROI', num2str(ROInum), ' Text']))

        cmap = glovar.CurrentCMap; 

        if strcmpi(cmap, 'RGB')
            linecolor = 'b';
        elseif strcmpi(cmap, 'Jet')
            linecolor = 'w';
        elseif strcmpi(cmap, 'Hot')
            linecolor = 'c';
        elseif strcmpi(cmap, 'Fire')
            linecolor = 'g'; %ZL comment, change the color of ROI during drawing
        end

        if gui_CaImageViewer.NewSpineAnalysis
            c = uicontextmenu;
            uimenu(c, 'Label', 'Set as eliminated', 'Callback', @CategorizeSpines);
            uimenu(c, 'Label', 'Set as active', 'Callback', @CategorizeSpines);
            gui_CaImageViewer.ROI(ROInum+1) = rectangle('Position', adjustedpos, 'EdgeColor', linecolor, 'ButtonDownFcn', {@DragROI, ROInum, 'HomeWindow'}, 'Curvature', [1 1],'Tag', ['ROI', num2str(ROInum)], 'UIContextMenu', c);
            gui_CaImageViewer.NewSpineAnalysisInfo.SpineList = [gui_CaImageViewer.NewSpineAnalysisInfo.SpineList; 1];
        else
            gui_CaImageViewer.ROI(ROInum+1) = rectangle('Position', adjustedpos, 'EdgeColor', linecolor, 'Curvature', [1 1],'Tag', ['ROI', num2str(ROInum)]);
        end
        gui_CaImageViewer.ROItext(ROInum+1) = text(adjustedpos(1)-4, adjustedpos(2)-3, num2str(ROInum), 'color', 'white', 'Tag', ['ROI', num2str(ROInum), ' Text'],'ButtonDownFcn', 'DeleteROI', 'Fontsize', 6);
    else
        
    end
else
    point1 = get(gca, 'CurrentPoint');  %%% Button down detected and position parameters stored (x, y, width, height)
    point1 = point1(1,1:2)   ;          %%% Extract x and y

    RoiRect = get(gco, 'Position') ;    %%% Get placement of current object in REFERENCE TO THE AXES

    set(gcf, 'Units', 'normalized');
    rectFig = get(gcf, 'Position');     %%% Get position of current figure on the screen

    set(gca, 'Units', 'normalized');
    rectAx = get(gca, 'Position');      %%% Get position of current axes in REFERENCE TO THE FIGURE

    set(0, 'Units', 'Pixels');
    scrnsz = get(0,'ScreenSize');

    xmag = rectAx(3)/128;                %%% screen pixel/image pixel
    xoffset =rectAx(1)*rectFig(3);
    ymag = rectAx(4)/128;
    yoffset = rectAx(2)*rectFig(4);
    rect1 = [xmag*RoiRect(1)+xoffset+0.5, ymag*(scrnsz(2)-RoiRect(2)-RoiRect(4))+yoffset+.5, xmag*RoiRect(3), ymag*RoiRect(4)];

    if point1(1)>(RoiRect(1)+2*RoiRect(3)/3) && point1(2)> (RoiRect(2)+2*RoiRect(4)/3) %%%resize
        fixedpoint = [(rectAx(1)+(xmag*RoiRect(1))), rectAx(2)+rectAx(4)-(ymag*RoiRect(2)+RoiRect(4))+(ymag*RoiRect(4))];
        rbbox([rectAx(1)+(xmag*RoiRect(1)), rectAx(2)+rectAx(4)-(ymag*RoiRect(2)+RoiRect(4))+(ymag*RoiRect(4)), xmag*RoiRect(3),ymag*RoiRect(4)], fixedpoint);
        point2 = get(gca, 'CurrentPoint');
        point2 = point2(1,1:2);
        RoiRect(3) = point2(1)-RoiRect(1);
        RoiRect(4) = point2(2)-RoiRect(2);
        RoiRect_final = [RoiRect(1), RoiRect(2), point2(1)-RoiRect(1), point2(2)-RoiRect(2)];
        oldbox = findobj('Tag', 'ROI confine');
        delete(oldbox);
    else %%% drag
        rbbox;      %%% initiate function to detect transition from button press to button release
        point2 = get(gca, 'CurrentPoint');
        point2 = point2(1,1:2);
        RoiRect_new = [point2(1)-(RoiRect(3)/2), point2(2)-RoiRect(4), RoiRect(3), RoiRect(4)];
        RoiRect_final = dragrect(RoiRect_new);
        oldbox = findobj('Tag', 'ROI confine');
        delete(oldbox);
    end

    Roi_Rect_tag = get(gco, 'Tag');
    RoiN = Roi_Rect_tag(end);

    ROItext = findobj(gca,'Type', 'Text', 'Tag', ['ROI', num2str(ROInum), ' Text']);
    if twochannels == 1
        ROIredtext = findobj(gca, 'Type', 'Text', 'Tag', ['ROIred', num2str(ROInum), ' Text']);
    end
    set(ROItext, 'Position', [RoiRect_final(1)-2, RoiRect_final(2)-2, 0]);
    % set(ROItext, 'Position', [RoiRect_final(1), RoiRect_final(2), 0]);

    actualROI = gco;
    set(actualROI, 'Position', RoiRect_final);
    rectangle('Position', RoiRect_final, 'EdgeColor', 'w', 'Curvature', [0 0], 'Tag', 'ROI confine', 'Linewidth', 1, 'LineStyle', ':');
    uistack(actualROI, 'top');
        
    if twochannels == 1
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
    end
end
