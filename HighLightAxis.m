function HighLightAxis(~,~)

currentwidth = get(gca, 'Linewidth');

ActImageMat = get(gcf, 'UserData');
session = get(gco, 'UserData');

if currentwidth == 0.5
    set(gca, 'Linewidth', 3)
    set(gca, 'Box', 'on')
    set(gca, 'XColor', 'g')
    set(gca, 'YColor', 'g')
    ActImageMat(1,session) = 1;
else
    set(gca, 'Linewidth', 0.5)
    set(gca, 'Box', 'off')
    set(gca, 'XColor', 'k')
    set(gca, 'YColor', 'k')
    ActImageMat(1,session) = 0;
end

set(gcf, 'UserData', ActImageMat);