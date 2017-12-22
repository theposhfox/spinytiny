function DeleteROI

global gui_SpineSegmentation

clicktype = get(gcbf, 'SelectionType');

if strcmp(clicktype, 'alt')
    BP_tag = get(gco, 'Tag');
    BP_num = str2num(BP_tag(13:end));
    ROI = findobj('Tag', BP_tag);
    delete(ROI);
    oldspines = gui_SpineSegmentation.Spines;
    AllROIs = findobj('Type', 'rectangle');         %%% ***** NOTE: this captures the ROIs in reverse order (i.e. top--> bottom is arranged as last-->first 
    gui_SpineSegmentation.Spines = [];              %%% Need to restructure the spine data array to account for the deletion
    if BP_num == 1
        gui_SpineSegmentation.Spines = oldspines(2:end,:);
        for i = 1:size(oldspines,1)-1
            set(AllROIs(size(oldspines,1)-i), 'Tag', ['Branchpoint ', num2str(i-1)]) 
        end
    elseif BP_num > 1 && BP_num < size(oldspines,1)
        gui_SpineSegmentation.Spines(1:BP_num-1, :) = oldspines(1:BP_num-1,:);                  %%% captures all the spines from the first to the one before the one being deleted
        gui_SpineSegmentation.Spines(BP_num:size(oldspines,1)-1,:) = oldspines((BP_num+1):size(oldspines,1),:) %%% captures all the spines after the one being deleted and shifts it to fill the place of the one being deleted
        for i = 1:((size(oldspines,1))-BP_num)
            set(AllROIs(i), 'Tag', ['Branchpoint ', num2str(size(oldspines,1)-i)])
        end
    elseif BP_num == size(oldspines,1)
        gui_SpineSegmentation.Spines(1:BP_num-1,:) = oldspines(1:BP_num-1,:);
        set(AllROIs(1), 'Tag', ['Branchpoint ', num2str(BP_num-1)])
    end
end