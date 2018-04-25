function warpmatrix = CompareImagePair(~,~)

global gui_CaImageViewer

selectedaxes = findobj(gcf, 'XColor', [0 1 0]);     %%% Finds the selected axes based on the color set to 'XColor' in function HighLightAxis)


%%% If there are more than two axes selected, ignore 
if length(selectedaxes)>2
    selectedaxes = selectedaxes(1:2);
end

selectedaxes = flipud(selectedaxes);

%%% Pull images from selected axes
for i = 1:length(selectedaxes)
    im{i} = uint16(get(get(selectedaxes(i), 'Children'), 'CData'));
end

%%%%%%%%%%%%%%%%%%%%%%%
centeredimage = im{2};
mobileimage = im{1};
%%%%%%%%%%%%%%%%%%%%%%%

alignchoice = get(findobj('Tag', 'Alignment_CheckBox'), 'Value');

if alignchoice 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[RESULTS, WARP, WARPEDiMAGE] = ECC(IMAGE, TEMPLATE, LEVELS, NOI, TRANSFORM, DELTA_P_INIT)

    delta_p_init = zeros(2,3); delta_p_init(1,1) = 1; delta_p_init(2,2) = 1;
    [results, warpmatrix, shiftedimage] = ecc(mobileimage, centeredimage, 5,25, 'affine', delta_p_init);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
else
    shiftedimage = uint16(mobileimage);
end


figure('Name', 'Side-by-side Comparison','Position', [336,285,1217,522], 'NumberTitle', 'off');  
im1 = axes('Units', 'Normalized','Position', [0.05 0.05 0.4 0.9], 'XTick', [], 'YTick', []); 
imshow(centeredimage,[]); colormap(fire);
title('Image 1')
im2 = axes('Units', 'Normalized','Position', [0.55 0.05 0.4 0.9], 'XTick', [], 'YTick', []);
imshow(shiftedimage, []); colormap(fire);
title('Image 2')

figure('Name', 'Side-by-side Comparison','Position', [336,255,1217,522], 'NumberTitle', 'off');  
im3 = axes('Units', 'Normalized','Position', [0.05 0.05 0.4 0.9], 'XTick', [], 'YTick', []); 
imshow(im2double(mobileimage)/max(im2double(mobileimage(:)))-im2double(centeredimage)/max(im2double(centeredimage(:))),[]); 
title('Image difference pre-correction', 'Fontsize', 12)
im4 = axes('Units', 'Normalized','Position', [0.55 0.05 0.4 0.9], 'XTick', [], 'YTick', []); 
imshow(im2double(shiftedimage)/max(im2double(shiftedimage(:)))-im2double(centeredimage)/max(im2double(centeredimage(:))),[]);
title('Image difference post-correction', 'Fontsize', 12)

figure('Name', 'RGB Overlap','Position', [336,255,1217,522], 'NumberTitle', 'off');  
im1pre = im2double(centeredimage)/max(im2double(centeredimage(:)));
im1pre3 = zeros(size(centeredimage,1),size(centeredimage,1), 3);
im1pre3(1:end,1:end,2) = im1pre;
im2post = im2double(shiftedimage)/max(im2double(shiftedimage(:)));
im2post3 = zeros(size(centeredimage,1),size(centeredimage,1), 3);
im2post3(1:end,1:end,1) = im2post;
im5 = axes('Units', 'Normalized','XTick', [], 'YTick', []); 
imshow(im1pre3+im2post3)
title('Color-coded alignment', 'Fontsize', 12)
linkaxes([im1,im2,im3,im4,im5, gui_CaImageViewer.figure.handles.GreenGraph], 'xy')

warpmatrix




