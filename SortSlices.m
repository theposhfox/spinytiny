clc, clear, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% User define Variable (Change this part if necessary)  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Brain section number             %%%%%%%%%%%%%%%%%%%
animal_num = 2;  %% For the second brain, change it to 2

%%% Input folder name and output stack name    %%%%%%%%%
%%% For Windows System, you may need to change "/" to "\" manually
cd('E:')
folder_pre = ['E:', filesep, 'CavCre',num2str(animal_num),'.'];
folder_post = '.vsi.Collection';
layer_name = [filesep,'Layer'];
file_name = [filesep, 'Layer.btf'];

output_name = 'stack';

%%% Loop info %%% 
max_collection = 8;  %%% Loop from 2.1 to 2.8
max_layer = 0;       %%% Don't need to care!
min_layer = 0;       %%% Start the loop from 0 to max_layer

height = 1800;    
width = 1200;

rotate_clock_wise_90 = 1;   %%% Clock wise rotate 90 degree; 0 to turn off
flip_ud = 1;                %%% Top - Down flip; 0 to turn off
flip_lr = 0;                %%% Left - right flip ; 1 to turn on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compression and Save to stack    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('E:\')
slicemap = fastdir(cd, ['Animal ', num2str(animal_num), ' Slice Mapping']);
load(slicemap{1})
current_slice_map = SliceMapping{1};
reordered_layer = find(current_slice_map == min_layer + 1)-1;


f_dir = [folder_pre,num2str(1),folder_post ...
        ,layer_name,num2str(reordered_layer),file_name];
img = imread(f_dir);
img_r = imresize(img,[height, width]);

if rotate_clock_wise_90 == 1
    img_r = permute(img_r,[2,1,3]);
end

if flip_ud == 1
    img_r = flipdim(img_r, 1);
end

if flip_lr == 1
    img_r = flipdim(img_r, 2);
end

imwrite(img_r,[output_name,num2str(animal_num),'.tif'],...
        'Compression','packbits');
    
sortfig = figure('Position', [456 94 1035 401]);
keep_button = uicontrol('style', 'pushbutton', 'Units', 'Normalize','Position', [0.5703 0.025 0.16 0.08], 'callback', @flipdecision, 'String', 'Keep as-is');
flip_button = uicontrol('style', 'pushbutton', 'Units', 'Normalize','Position', [0.745 0.025 0.16 0.08], 'callback', @flipdecision, 'String', 'Flip L/R');
currentimage = ['Slide ', num2str(1), ' Slice ', num2str(min_layer+1), ' (Layer ', num2str(reordered_layer), ')'];

for i = 1:max_collection
    figure('Position', [456 582 1035 401])
    max_layer = size(fastdir([folder_pre, num2str(i), folder_post], 'Layer'),1)-1;
    current_slice_map = SliceMapping{i};
    if isempty(fastdir([folder_pre, num2str(i), folder_post, filesep],'Overview_Downsampled.btf'))
        cd([folder_pre, num2str(i), folder_post, filesep])
        OV = imread([folder_pre, num2str(i), folder_post, filesep,'Overview.btf']);
        OV_r = imresize(OV, [height, width]);
        imwrite(OV_r,'Overview_Downsampled.tif', 'Compression', 'packbits')
        OV_r = permute(OV_r,[2,1,3]); OV_r = flipdim(OV_r,1);
    else
        imread([folder_pre, num2str(i), folder_post, filesep,'Overview_Downsampled.btf'])
    end
    imagesc(OV_r); set(gca, 'XTick', []); set(gca, 'YTick', []);
    clear('OV'); clear('OV_r');
    if i==1
        min_layer = 1;
    else
        min_layer = 0;
    end
    for j = min_layer : max_layer
        if i == 0 && j == 0
            continue
        else
            reordered_layer = find(current_slice_map == j+1)-1;
            f_dir = [folder_pre,num2str(i)  ...
                    ,folder_post,layer_name,num2str(reordered_layer),file_name];
            figure(sortfig)
            subplot(1,2,1); cla;
            imagesc(img_r); set(gca, 'XTick', []); set(gca, 'YTick', []);
            title(currentimage)
            img = imread(f_dir);
            img_r = imresize(img, [height, width]);
        
            if rotate_clock_wise_90 == 1
                img_r = permute(img_r,[2,1,3]);
            end
            
            if flip_ud == 1
                img_r = flipdim(img_r, 1);
            end
            
            currentimage = ['Slide ', num2str(i), ' Slice ', num2str(j+1), ' (Layer ', num2str(reordered_layer), ')'];
            subplot(1,2,2);
            imagesc(img_r); set(gca, 'XTick', []); set(gca, 'YTick', []);
            title(currentimage)
            uiwait
            
            flipchoice = get(gcf, 'UserData');
            if strcmpi(flipchoice, 'Flip L/R');
                flip_lr = 1;
            else
                flip_lr = 0;
            end
            
            if flip_lr == 1
                img_r = flipdim(img_r, 2);
            end
            cla; imagesc(img_r);set(gca, 'XTick', []); set(gca, 'YTick', [])
            cd('E:\')
            imwrite(img_r, [output_name,num2str(animal_num), ...
                '.tif'], 'WriteMode', 'append', ...
                'Compression','packbits');
        end
    end
end