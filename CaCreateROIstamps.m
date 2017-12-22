function [ROI_stamp, coordinates] = CreateROIstamps

ROI_stamp = [];
coordinates = [];

program = get(gcf);

running = program.FileName;

if ~isempty(regexp(running, 'CaImageViewer'))
    global gui_CaImageViewer
    glovar = gui_CaImageViewer;
elseif ~isempty(regexp(running, 'FluorescenceSuite'));
    global gui_FluorescenceSuite
    glovar = gui_FluorescenceSuite;
end

if isfield (glovar, 'ROI')
    if ~isempty(glovar.ROI)
        for i = 1:length(glovar.ROI)
            if ishandle(glovar.ROI(i))
                existing_ROI{i} = get(glovar.ROI(i));
            else
                existing_ROI{i} = [];
            end
        end
    else existing_ROI = [];
    end
else existing_ROI = [];
end

if isfield(glovar, 'PolyROI')
    for i = 1:length(glovar.PolyROI)
        if ishandle(glovar.PolyROI{i})
            coordinates{i} = glovar.PolyLinePos{i};
        else
            coordinates{i} = [];
        end
    end
end

for i = 1:length(existing_ROI)
    if ~isempty(existing_ROI{i})
        ROI_stamp{i} = get(glovar.ROI(i), 'Position');
    end
end
