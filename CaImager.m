function CaImager

fname = '/Users/theposhfox/Desktop/NH002_4_24_14/NH003g001.tif'

CaImage_File_info = imfinfo(fname);

timecourse_image_number = numel(CaImage_File_info);

Green_Frame = 1;
Red_Frame = 1;

for i = 1:2:timecourse_image_number
    [GCaMP_Image{1,Green_Frame}, mapg] = imread(fname, i);
    Green_Frame = Green_Frame+1;
end
for i = 2:2:timecourse_image_number
    [Red_Image{1,Red_Frame}, mapr] = imread(fname, i);
    Red_Frame = Red_Frame+1;
end
% 
% Green_Figure = subplot(1,2,1);
% imshow(GCaMP_Image{1}, mapg);
% 
% Red_Figure = subplot(1,2,2);
% imshow(Red_Image{1}, mapr);
% 
% axes(Green_Figure);
% 
% ROI_1 = imellipse;
% 
% vertices = wait(ROI_1);
% 
% copyobj(ROI_1, Red_Figure);


% Green_ROI = roipoly(Green_Figure);

%figure; imshow(Green_ROI);
% 



