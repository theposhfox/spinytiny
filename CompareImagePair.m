function CompareImagePair(~,~)

global gui_CaImageViewer

selectedaxes = findobj(gcf, 'XColor', [0 1 0]);     %%% Finds the selected axes based on the color set to 'XColor' in function HighLightAxis)


%%% If there are more than two axes selected, ignore 
if length(selectedaxes)>2
    selectedaxes = selectedaxes(1:2);
end

%%% Pull images from selected axes
for i = 1:length(selectedaxes)
    im{i} = get(get(selectedaxes(i), 'Children'), 'CData');
end

alignchoice = get(findobj('Tag', 'Alignment_CheckBox'), 'Value');

if alignchoice 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[RESULTS, WARP, WARPEDiMAGE] = ECC(IMAGE, TEMPLATE, LEVELS, NOI, TRANSFORM, DELTA_P_INIT)

    delta_p_init = zeros(2,3); delta_p_init(1,1) = 1; delta_p_init(2,2) = 1;
    [results, warp, comparatorimage] = ecc(im{2}, im{1}, 3,50, 'affine', delta_p_init);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    comparatorimage = uint8(im{2});
end


figure('Name', 'Side-by-side Comparison','Position', [336,285,1217,522], 'NumberTitle', 'off');  
im1 = axes('Units', 'Normalized','Position', [0.05 0.05 0.4 0.9], 'XTick', [], 'YTick', []); 
imshow(uint8(im{1}),[]); colormap(fire);
title('Image 1')
im2 = axes('Units', 'Normalized','Position', [0.55 0.05 0.4 0.9], 'XTick', [], 'YTick', []);
imshow(comparatorimage, []); colormap(fire);
title('Image 2')

figure('Name', 'Side-by-side Comparison','Position', [336,255,1217,522], 'NumberTitle', 'off');  
im3 = axes('Units', 'Normalized','Position', [0.05 0.05 0.4 0.9], 'XTick', [], 'YTick', []); 
imshow(uint8(im{2}-im{1}),[]); colormap(jet);
title('Image difference pre-correction', 'Fontsize', 12)
im4 = axes('Units', 'Normalized','Position', [0.55 0.05 0.4 0.9], 'XTick', [], 'YTick', []); 
imshow(comparatorimage-uint8(im{1}),[]); colormap(jet);
title('Image difference post-correction', 'Fontsize', 12)
linkaxes([im1,im2,im3,im4], 'xy')




