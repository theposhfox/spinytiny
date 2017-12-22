function SetEndFrame(~,~,timecourse_length)

global gui_CaImageViewer

newlength = inputdlg('Stop analysis at frame:', 'Choose end frame', 1, {num2str(timecourse_length)});

gui_CaImageViewer.SelectedStopFrame = str2num(newlength{1});