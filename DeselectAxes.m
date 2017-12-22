function DeselectAxes(~,~)

global gui_CaImageViewer

selectedaxes = findobj(gcf, 'XColor', [0 1 0]);

for i = 1:length(selectedaxes)
    set(selectedaxes(i), 'XColor', 'k')
    set(selectedaxes(i), 'YColor', 'k')
    set(selectedaxes(i), 'Linewidth', 0.5) 
end