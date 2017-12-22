function GroupImagingWindows(~,~)

global gui_CaImageViewer

selectedaxes = findobj(gcf, 'XColor', [0 1 0]);     %%% Finds the selected axes based on the color set to 'XColor' in function HighLightAxis
imagingfieldlabel = get(findobj(gcf,'Type', 'uicontrol', 'Style', 'edit'), 'String');

for i = 1:length(selectedaxes)
    dates(i,:) = get(get(selectedaxes(i), 'Title'), 'String');
    axes(selectedaxes(i))
    xlabel(['Imaging field ', imagingfieldlabel], 'Color', 'k');
end

%%%%% Move to parent folder

fullpath = gui_CaImageViewer.save_directory;
allseps = strfind(fullpath, '\');
stepsup = 2;
newpath = fullpath(1:allseps(end-stepsup)-1); %%% move two steps up in the path directory to get bath to the main animal folder (e.g. Z:/People/Nathan/Data/NH004 instead of Z:/People/Nathan/Data/NH004/160316/summed)
cd(newpath)

%%%% Construct features for building a table to show spine lifetimes over
%%%% each imaging field

SpineRegistry = struct;
SpineRegistry.ColumnNames = cellfun(@(x) ['Day ', num2str(x)], mat2cell(1:length(selectedaxes), 1, ones(1,length(selectedaxes))), 'Uni', false);
SpineRegistry.Data = [];
SpineRegistry.ColumnFormat = repmat({'logical'}, 1, length(selectedaxes));
SpineRegistry.ColumnEditable = true(1,length(selectedaxes));
SpineRegistry.RowName = [];
SpineRegistry.DatesAcquired = flipud(mat2cell(dates,ones(1,length(selectedaxes)), 6));

save(['Imaging Field ', imagingfieldlabel, ' Spine Registry'], 'SpineRegistry')

% cd(fullpath)

nextimagingwindow = num2str(str2num(imagingfieldlabel)+1);
set(findobj(gcf, 'Type', 'uicontrol', 'Style', 'edit'), 'String', nextimagingwindow);


