function varargout = SpineSegmentationGUI(varargin)
% SPINESEGMENTATIONGUI MATLAB code for SpineSegmentationGUI.fig
%      SPINESEGMENTATIONGUI, by itself, creates a new SPINESEGMENTATIONGUI or raises the existing
%      singleton*.
%
%      H = SPINESEGMENTATIONGUI returns the handle to a new SPINESEGMENTATIONGUI or the handle to
%      the existing singleton*.
%
%      SPINESEGMENTATIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPINESEGMENTATIONGUI.M with the given input arguments.
%
%      SPINESEGMENTATIONGUI('Property','Value',...) creates a new SPINESEGMENTATIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SpineSegmentationGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SpineSegmentationGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SpineSegmentationGUI

% Last Modified by GUIDE v2.5 09-Feb-2016 12:27:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SpineSegmentationGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SpineSegmentationGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before SpineSegmentationGUI is made visible.
function SpineSegmentationGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SpineSegmentationGUI (see VARARGIN)

% Choose default command line output for SpineSegmentationGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global gui_SpineSegmentation

gui_SpineSegmentation.handles = handles;

set(handles.Image_Axes, 'YTick', []);
set(handles.Image_Axes, 'XTick', []);
set(handles.Image_Axes, 'ButtonDownFcn', []);
set(handles.Image_Axes, 'Visible', 'off');
set(handles.TimeCourse_Slider, 'Visible', 'off');
gui_SpineSegmentation.LoadedImage = 0;
gui_SpineSegmentation.TC_Length = [];
gui_SpineSegmentation.Arbitrary_Shaping = 0;
gui_SpineSegmentation.Isolated_Dendrite = [];
gui_SpineSegmentation.Dendrite_Bridges = [];
gui_SpineSegmentation.Spines = [];
gui_SpineSegmentation.ProcessedView = 0;


% UIWAIT makes SpineSegmentationGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SpineSegmentationGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function LoadFile_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to LoadFile_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LoadFile_EditableText as text
%        str2double(get(hObject,'String')) returns contents of LoadFile_EditableText as a double


% --- Executes during object creation, after setting all properties.
function LoadFile_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LoadFile_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadFile_PushButton.
function LoadFile_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadFile_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_SpineSegmentation

gui_SpineSegmentation.Segmented = 0;

mode = gui_SpineSegmentation.Mode;


if strcmpi(mode, 'Cluster')

    [fname pname] = uigetfile('.tif');

    Image_file = [pname, fname];

    cd(pname);

    set(handles.LoadFile_EditableText, 'String', Image_file);

    gui_SpineSegmentation.Image_File = Image_file;
    gui_SpineSegmentation.Resized_Image = [];
    img_inf = imfinfo(fname);
    frames = numel(img_inf);
    gui_SpineSegmentation.TC_Length = frames;
    
    for i = 1:frames
        Im(:,:,i) = imread(Image_file,i);
    end

    gui_SpineSegmentation.Image = Im;
    gui_SpineSegmentation.Unprocessed_Image = Im;
    
    set(handles.TimeCourse_Slider, 'Value', 1);
    if frames > 1
        set(handles.TimeCourse_Slider, 'Max', frames);
        set(handles.TimeCourse_Slider, 'SliderStep', [1/(frames-1) 5/(frames-1)]);
    else
        set(handles.TimeCourse_Slider, 'Max', 2);
        set(handles.TimeCourse_Slider, 'SliderStep', [0 0]);
    end
    set(handles.TimeCourse_Slider, 'Min', 1);


    [counts, x] = imhist(Im);

    [val pos] = max(counts);

    axes(handles.Image_Axes);

    min1 = min(min(Im));
    max1 = max(max(Im));

    imshow(Im, [pos max1/4]);
    gui_SpineSegmentation.LoadedImage = 1;
    
elseif strcmpi(mode, 'Tracking')
    [fname pname] = uigetfile('.tif');

    Image_file = [pname, fname];

    cd(pname);

    set(handles.LoadFile_EditableText, 'String', Image_file);

    gui_SpineSegmentation.Image_File = Image_file;
    gui_SpineSegmentation.Resized_Image = [];
    img_inf = imfinfo(fname);
    frames = numel(img_inf);
    gui_SpineSegmentation.TC_Length = frames;
    
    for i = 1:frames
        Im(:,:,i) = imread(Image_file,i);
    end
    
    set(handles.TimeCourse_Slider, 'Value', 1);
    if frames > 1
        set(handles.TimeCourse_Slider, 'Max', frames);
        set(handles.TimeCourse_Slider, 'SliderStep', [1/(frames-1) 5/(frames-1)]);
    else
        set(handles.TimeCourse_Slider, 'Max', 2);
        set(handles.TimeCourse_Slider, 'SliderStep', [0 0]);
    end
    set(handles.TimeCourse_Slider, 'Min', 1);

    [counts, x] = imhist(Im(:,:,1));

    [val pos] = max(counts);

    axes(handles.Image_Axes);

    min1 = min(min(Im(:,:,1)));
    max1 = max(max(Im(:,:,1)));

    imshow(Im(:,:,1), [pos max1/4]);
    
    gui_SpineSegmentation.LoadedImage = 1;  
end

    gui_SpineSegmentation.Image = Im;
    gui_SpineSegmentation.Unprocessed_Image = Im;



% --- Executes on button press in DrawBox_ToggleButton.
function DrawBox_ToggleButton_Callback(hObject, eventdata, handles)
% hObject    handle to DrawBox_ToggleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_SpineSegmentation
Box = get(handles.DrawBox_ToggleButton, 'Value');

if Box == 1
    set(handles.ArbShape_ToggleButton, 'Value', 0);
    gui_SpineSegmentation.Arbitary_Shaping = 0;

    axes(handles.Image_Axes);

    oldROIs = findobj('Tag', 'ROI');
    delete(oldROIs);

    waitforbuttonpress;
    pointer_location1 = get(gca, 'CurrentPoint');
    rbbox;
    pointer_location2 = get(gca, 'CurrentPoint');
    point1 = pointer_location1(1,1:2);
    point2 = pointer_location2(1,1:2);

    point = min(point1, point2);

    offset = abs(point1-point2);

    ROI = round([point, offset]);

    axes(gui_SpineSegmentation.handles.Image_Axes);

    restricted_Im = [];
    
    %%% Check if the ROI rbbox goes outside of the image boundaries
    if ROI(2)+ROI(4) > size(gui_SpineSegmentation.Image,1)
        ROI(4) = size(gui_SpineSegmentation.Image,1)-ROI(2);
    else
    end
    
    if ROI(1)+ROI(3) > size(gui_SpineSegmentation.Image,2)
        ROI(3) = size(gui_SpineSegmentation.Image,2)-ROI(1);
    else
    end
    
    for i = 1 : gui_SpineSegmentation.TC_Length
        restricted_Im(:,:,i) = gui_SpineSegmentation.Image(ROI(2):ROI(2)+ROI(4), ROI(1):ROI(1)+ROI(3),i);
    end

    restricted_Im = uint16(restricted_Im);

    gui_SpineSegmentation.Image = restricted_Im;
    gui_SpineSegmentation.Resized_Image = restricted_Im;

    frame_num = get(handles.TimeCourse_Slider, 'Value');
    frame_num = floor(frame_num);

    [counts, x] = imhist(restricted_Im(:,:,frame_num));

    [val pos] = max(counts);

    min1 = min(min(restricted_Im(:,:,frame_num)));
    max1 = max(max(restricted_Im(:,:,frame_num)));

    imshow(restricted_Im(:,:,frame_num), [pos, max1/2]);
    set(handles.DrawBox_ToggleButton, 'Value', 0);
end


% --- Executes on button press in RevertImage_PushButton.
function RevertImage_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to RevertImage_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)\

global gui_SpineSegmentation

gui_SpineSegmentation.Segmented = 0;
gui_SpineSegmentation.Resized_Image = [];

if gui_SpineSegmentation.LoadedImage == 0;
    return
end

for i = 1:gui_SpineSegmentation.TC_Length
    Im(:,:,i) = imread(gui_SpineSegmentation.Image_File,i);
end

gui_SpineSegmentation.Image = Im;
gui_SpineSegmentation.Unprocessed_Image = Im;

frame = get(gui_SpineSegmentation.handles.TimeCourse_Slider, 'Value');

[counts, x] = imhist(Im(:,:,frame));

[val pos] = max(counts);

axes(handles.Image_Axes);

min1 = min(min(Im(:,:,frame)));
max1 = max(max(Im(:,:,frame)));

imshow(Im(:,:,frame), [pos max1/4]);


% --- Executes on button press in ProcessImage_PushButton.
function ProcessImage_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to ProcessImage_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_SpineSegmentation

if ~isempty(gui_SpineSegmentation.Resized_Image)
    Im = gui_SpineSegmentation.Resized_Image;
else
    Im = gui_SpineSegmentation.Unprocessed_Image;
end

gui_SpineSegmentation.Unprocessed_Image = Im;

frames = size(Im,3);

if frames > 1
    currentframe = get(gui_SpineSegmentation.handles.TimeCourse_Slider, 'Value');
else
    currentframe = 1;
end

PSF = get(handles.EstimatedPSF_EditableText, 'String');
PSF = str2num(PSF);
PSF = fspecial('gaussian', PSF, 10);
PSF = ones(size(PSF));
routed = 1;
GaussFilt_size = str2num(get(handles.GaussFilt_EditableText, 'String'));
MedFilt_size = str2num(get(handles.MedFilt_EditableText, 'String'));
MinBlob = str2num(get(handles.MinBlob_EditableText, 'String'));
MaxBlob = str2num(get(handles.MaxBlob_EditableText, 'String'));

[Processed_Image, Label_Matrix] = NHimagedeconvolve(Im, PSF, routed, GaussFilt_size, MedFilt_size, MinBlob, MaxBlob);

axes(handles.Image_Axes); cla; pause(0.1)

Im = Processed_Image;

[counts, x] = imhist(Im(:,:,currentframe));

[val pos] = max(counts);

axes(handles.Image_Axes); cla

max1 = max(max(Im(:,:,currentframe)));

imshow(Im(:,:,currentframe), [pos max1/4]);

gui_SpineSegmentation.Processed_Image = Processed_Image;
gui_SpineSegmentation.Label_Matrix = Label_Matrix;
gui_SpineSegmentation.Image = Processed_Image;
gui_SpineSegmentation.Segmented = 0;
gui_SpineSegmentation.ProcessedView = 1;



function EstimatedPSF_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to EstimatedPSF_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EstimatedPSF_EditableText as text
%        str2double(get(hObject,'String')) returns contents of EstimatedPSF_EditableText as a double


% --- Executes during object creation, after setting all properties.
function EstimatedPSF_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EstimatedPSF_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Segment_PushButton.
function Segment_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Segment_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_SpineSegmentation

Processed_Image = gui_SpineSegmentation.Label_Matrix;

frames = size(Processed_Image,3);
currentframe = get(handles.TimeCourse_Slider, 'Value');

%% Detached spine segmentation

%%% see 'Automatic dendritic spine analysis in two-photon laser scanning
%%% microscopy images


%%% First, perform a 2d interpolation to connect potentially disjointed
%%% parts of the dendrite

Processed_Image = double(Processed_Image);

%%%

MinBlob = str2num(get(handles.MinSpine_EditableText, 'String'));
MaxBlob = str2num(get(handles.MaxSpine_EditableText, 'String'));

h_blob = vision.BlobAnalysis('AreaOutputPort', true, 'BoundingBoxOutputPort', true, 'OutputDataType', 'double', 'MinimumBlobArea', MinBlob, 'MaximumBlobArea', MaxBlob, 'LabelMatrixOutputPort', true);


Processed_Image = logical(Processed_Image);


for i = 1:frames
    [~, ~, ~, LabelMatrix(:,:,i)] = step(h_blob, Processed_Image(:,:,i));
end

%%% Create a red image out of the LabelMatrix (a copy of the identified
%%% blobs)

if frames == 1
    red_Im = cast(cat(3,LabelMatrix,zeros(size(LabelMatrix)), zeros(size(LabelMatrix))), 'double');
    white_im = cast(cat(3,Processed_Image, Processed_Image, Processed_Image), class(red_Im));
    white_im(:,:,2) = white_im(:,:,2)-red_Im(:,:,1);
    white_im(:,:,3) = white_im(:,:,3)-red_Im(:,:,1);
else
    for i = 1:frames
        framedim = size(LabelMatrix);
        red_Im{i} = cast(cat(3,LabelMatrix(:,:,i),zeros(framedim(1:2)), zeros(framedim(1:2))), 'double');
        white_im{i} = cast(cat(3,Processed_Image(:,:,i), Processed_Image(:,:,i), Processed_Image(:,:,i)), 'double');
        overlay_im{i} = white_im{i};
        overlay_im{i}(:,:,2) = overlay_im{i}(:,:,2)-red_Im{i}(:,:,1);
        overlay_im{i}(:,:,3) = overlay_im{i}(:,:,3)-red_Im{i}(:,:,1);
    end
end

axes(handles.Image_Axes);

if frames == 1
    imshow(overlay_im);
else
    imshow(overlay_im{currentframe})
end

%% Attached segmentation with parallel thinning 

%%% see Parallel thinning with two subiteration algorithms

if frames == 1
    dendrite_im = overlay_im;
    dendrite_im(:,:,1) = overlay_im(:,:,1)-red_Im(:,:,1); %%% reduces to only dendrites and attached spines

    dendrite_im2 = dendrite_im(:,:,1); 
    dendrite_im2(dendrite_im2<0) = 0;
    dendrite_im2 = logical(dendrite_im2);

    skeleton_im = bwmorph(dendrite_im2, 'thin', Inf);
    gui_SpineSegmentation.Skeleton = skeleton_im;
    
    %%%%%%%%
    SpurGUI
    %%%%%%%%
    
    SpurNum = gui_SpineSegmentation.SpurNum;

    cleaned = bwmorph(skeleton_im, 'spur', SpurNum); %%% de-spur the image 'n' times; this should get rid of branches for analysis of the main dendrite only
    % branchie_ps = bwmorph(rem_branches, 'branchpoints');
    % brokenLines = rem_branches-branchie_ps;

    cleaned = double(cleaned);
    cleaned = imdilate(cleaned, strel('disk',2));
    cleaned = bwmorph(cleaned, 'thin', Inf);
    cleaned = logical(cleaned);

    skeleton_frame = 0.5*skeleton_im; %%% turns values gray instead of white

    skeleton_frame = cat(3,skeleton_frame, skeleton_frame, skeleton_frame); %%% Convert to RGB by dimensional triplicate

    attached_segmentation_im = overlay_im;

    attached_segmentation_im(:,:,1) = overlay_im(:,:,1)-skeleton_frame(:,:,1);
    attached_segmentation_im(:,:,3) = attached_segmentation_im(:,:,3)-skeleton_frame(:,:,1);

    gui_SpineSegmentation.Final_Image_Data = attached_segmentation_im;
    gui_SpineSegmentation.Final_Image = imshow(attached_segmentation_im);
    set(gui_SpineSegmentation.Final_Image, 'ButtonDownFcn', 'AddSpine');
    set(gui_SpineSegmentation.Final_Image, 'HitTest', 'on');

    [i, j] = find(bwmorph(skeleton_im, 'branchpoints'));

    branches = [j i];
    branches = branches - 1; %%% corrects for apparent 'shift' when drawing markers
    gui_SpineSegmentation.Spines = branches;
    
    for i = 1:size(branches,1)
        rectangle('Position', [branches(i,1), branches(i,2), 2,2], 'Curvature', [1 1], 'EdgeColor', 'b', 'Tag', ['Branchpoint ', num2str(i)], 'ButtonDownFcn', 'DeleteROI')
    end

else
    dendrite_im = overlay_im;
    for i = 1:frames
        dendrite_im{i}(:,:,1) = overlay_im{i}(:,:,1)-red_Im{i}(:,:,1); %%% reduces to only dendrites and attached spines
        dendrite_im2{i} = dendrite_im{i}(:,:,1); 
        dendrite_im2{i}(dendrite_im2{i}<0) = 0;
        dendrite_im2{i} = logical(dendrite_im2{i});
        skeleton_im{i} = bwmorph(dendrite_im2{i}, 'thin', Inf);
    end
    
    gui_SpineSegmentation.Skeleton = skeleton_im;
    gui_SpineSegmentation.SpurMessage = 'First, remove all branchpoints to get just the dendrite';
    
    %%%%%%%%
    SpurGUI
    %%%%%%%%
    
    SpurNum = gui_SpineSegmentation.SpurNum;
    
    for i = 1:frames
        cleaned{i} = bwmorph(skeleton_im{i}, 'spur', SpurNum); %%% de-spur the image 'n' times; this should get rid of branches for analysis of the main dendrite only
        cleaned{i} = double(cleaned{i});
        cleaned{i} = imdilate(cleaned{i}, strel('disk',2));
        cleaned{i} = bwmorph(cleaned{i}, 'thin', Inf);
        cleaned{i} = logical(cleaned{i});
        skeleton_im{i} = cleaned{i};
    end
    
%     gui_SpineSegmentation.Skeleton = skeleton_im;
    % branchie_ps = bwmorph(rem_branches, 'branchpoints');
    % brokenLines = rem_branches-branchie_ps;
    
    h = waitbar(0, 'Calculating lengths of each segment...');
    
    for k = 1:frames
        B{k} = bwmorph(skeleton_im{k}, 'branchpoints');
        E{k} = bwmorph(skeleton_im{k}, 'endpoints');
        [y x] = find(E{k});
        B_loc = find(B{k});
        Dmask = false(size(skeleton_im{1}));
        for i = 1:numel(x)
            D = bwdistgeodesic(skeleton_im{1},x(i),y(i));
            distanceToBranchPt(1,i) = min(D(B_loc));
            if distanceToBranchPt(1,i) > 20
                continue
                else
            end
            Dmask(D < distanceToBranchPt(1,i)) = true;
        end
        skelD{k} = skeleton_im{k}-Dmask;
        waitbar(k/frames,h,['Measuring frame ', num2str(k)])
    end
    
    close(h);
    
    skeleton_im = skelD;
    
    for i = 1:frames
        skeleton_frame{i} = 0.5*skeleton_im{i};                  %%% turns values gray instead of white
        skeleton_frame{i} = cat(3,skeleton_im{i}, skeleton_im{i}, skeleton_im{i}); %%% Convert to RGB by dimensional triplicate
    end
        
    DilateNum = gui_SpineSegmentation.DilateNum;
    IsoDend = cell(frames,1);

    for i = 1:frames
        IsoDend{i} = imdilate(skeleton_im{i}, strel('disk', DilateNum));
    end
    
    AttachedSpines = cell(frames,1);
    AllSpines = cell(frames,1);
    for i = 1:frames
        AttachedSpines{i} = dendrite_im2{i}-IsoDend{i};
        AllSpines{i} = overlay_im{i};
        AllSpines{i}(:,:,1) = AllSpines{i}(:,:,1)-AttachedSpines{i};
        AllSpines{i}(:,:,3) = AllSpines{i}(:,:,3)-AttachedSpines{i};
    end
    
    imshow(AllSpines{currentframe});
    

            
    attached_segmentation_im = overlay_im;
    
    for i = 1:frames   %%%% Subtract the skeleton from the R and G channels to make a blue-ish trace of the skeleton on the white dendrite
        attached_segmentation_im{i}(:,:,1) = overlay_im{i}(:,:,1)-skeleton_frame{i}(:,:,1);
        attached_segmentation_im{i}(:,:,3) = attached_segmentation_im{i}(:,:,3)-skeleton_frame{i}(:,:,1);
    end
    
    gui_SpineSegmentation.Final_Image_Data = AllSpines;
    gui_SpineSegmentation.Final_Image = imshow(AllSpines{currentframe});
    set(gui_SpineSegmentation.Final_Image, 'ButtonDownFcn', 'AddSpine');
    set(gui_SpineSegmentation.Final_Image, 'HitTest', 'on');
    gui_SpineSegmentation.Image = gui_SpineSegmentation.Final_Image_Data;

    x = [];
    y = [];
    for i = 1:frames
        [y{i} x{i}] = find(bwmorph(skeleton_im{i}, 'branchpoints'));
        branches{i} = [x{i} y{i}];
        branches{i} = branches{i} -1;
    end
    
    gui_SpineSegmentation.Spines = branches;
    
%     for i = 1:size(branches{currentframe},1)
%         rectangle('Position', [branches{currentframe}(i,1), branches{currentframe}(i,2), 2,2], 'Curvature', [1 1], 'EdgeColor', 'b', 'Tag', ['Branchpoint ', num2str(i)], 'ButtonDownFcn', 'DeleteROI')
%     end
    
    gui_SpineSegmentation.Segmented = 1;
    gui_SpineSegmentation.ProcessedView = 1;
    gui_SpineSegmentation.Image = gui_SpineSegmentation.Final_Image_Data;
end


%%% Playgroung for attempting cross correlation

% for i = 1:frames
%     fft_frame = fft2(Image(:,:,i));
%     prod = fft_ref.*conj(fft_frame);
%     cc = ifft2(prod)
%     [maxy maxx] = find(fftshift(cc)==max(max(cc)));
%     shifty(i) = maxy(1)-centery;
%     shiftx(i) = maxx(1)-centerx;
%     if i > 1
%         if abs(shifty(i)-shifty(i-1))>(size(Image(:,:,1),1))/2;
%             shifty(i) = shifty(i)-sign(shifty(i)-shift(i-1))*(size(Image(:,:,1),1));
%         end
%         if abs(shiftx(i)-shiftx(i-1))>(size(Image(:,:,1),2))/2;
%             shiftx(i)= shiftx(i)-sign(shiftx(i)-shiftx(i-1))*(size(Image(:,:,1),2));
%         end
%     end
% end

% map = insertMarker(attached_segmentation_im, branches,'o', 'color', 'g');
% imshow(map)



% --- Executes on button press in Analyze_PushButton.
function Analyze_PushButton_Callback(~, ~, ~)
% hObject    handle to Analyze_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_SpineSegmentation

%% Trace main dendrite

skeleton_im = gui_SpineSegmentation.Skeleton;

%%%% De-branch the main dendrite by using an interactive spurring function
%%%% (scroll through the slider until you get rid of the branches, then hit
%%%% 'Accept')

    SpurGUI;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Now, perform the actual spurring based on the values you chose

SpurNum = gui_SpineSegmentation.SpurNum;

rem_branches = bwmorph(skeleton_im, 'spur', SpurNum); %%% de-spur the image 'n' times; this should get rid of branches for analysis of the main dendrite only
% branchie_ps = bwmorph(rem_branches, 'branchpoints');
% brokenLines = rem_branches-branchie_ps;

rem_branches = double(rem_branches);
rem_branches = imdilate(rem_branches, strel('disk',2));
rem_branches = bwmorph(rem_branches, 'thin', Inf);
rem_branches = logical(rem_branches);

%%% Perform another round of blob analysis to get rid of small/artifactual
%%% objects
h_blob = vision.BlobAnalysis('AreaOutputPort', true, 'BoundingBoxOutputPort', true, 'OutputDataType', 'double', 'MinimumBlobArea', 50, 'MaximumBlobArea', 10000, 'LabelMatrixOutputPort', true);
[~, ~, ~, LabelMatrix] = step(h_blob, rem_branches);
rem_branches = logical(LabelMatrix);

endpoints = bwmorph(rem_branches, 'endpoints');

%%% Progress notes: you have calculated the closest point on the dendrite
%%% to each spine; now find a way to bin these into measured distances
%%% along the dendrite. The value 'D' should have this from the
%%% bwdistgeodesic function...

[i j] = find(endpoints);
endpoints(i(1),j(1)) = 0;   %%% The geodesic distance calculation below works best for this purposes when there is only ONE ENDPOINT

D = bwdistgeodesic(rem_branches, endpoints, 'quasi-euclidean');
D(isnan(D)) = 0;
[i j] = find(D);
c = [j i];  %%% Coordinates of all points along the line found by thining and removing spurs (the main dendrites)


branches = gui_SpineSegmentation.Spines;


%%% Calculate the pixel distance from the branch points to the nearest
%%% point on the skeleton

for j = 1:size(branches,1)
    for i = 1:size(c,1)
        dist_vecs{j}(i,1:2) = branches(j,:)-c(i,:);
    end
end

for j = 1:size(branches,1)
    dist{j} = sqrt(abs((dist_vecs{j}(:,1).^2)+dist_vecs{j}(:,2).^2));
    [val pos] = min(dist{j});
    address{j} = c(pos,1:2);
    position_bin(j) = D(address{j}(2), address{j}(1));
end

figure; hist(position_bin, 10); title('Histogram of Spine Positions'); xlabel('Address'); ylabel('Number of Putative Spines')

%%% Label the dendrite with the addressing scheme used to build the above
%%% histogram

figure; imshow(gui_SpineSegmentation.Processed_Image); title('Numeric Addresses for Dendrite'); 
for i = 1:10:size(c,1)
    text(c(i,1),c(i,2),num2str(D(c(i,2), c(i,1))), 'color', 'g')
end

%%% Calculate clustering index

for i = 1:length(position_bin)
    for j = 2:length(position_bin)
        comparison_mat(i,j-1) = abs(position_bin(i)-position_bin(j));
    end
end

nearest_neighbor = sort(comparison_mat,2);
nearest_neighbor = nearest_neighbor(:,2:end); %%% Need to take into account self-comparison always equal to zero, and sorting puts this in the first posiiton

% figure;

AnalyzedData = [];
AnalyzedData.PositionsAlongDendrite = position_bin;
AnalyzedData.SpineDistribution = nearest_neighbor;
AnalyzedData.DendriteImage = gui_SpineSegmentation.Final_Image_Data;

save('AnalyzedData', 'AnalyzedData');


% --- Executes on button press in TrackingMode_PushButton.
function TrackingMode_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to TrackingMode_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_SpineSegmentation

gui_SpineSegmentation.Mode = 'Tracking';

set(handles.TrackingMode_PushButton, 'Visible', 'off');
set(handles.ClusterMode_PushButton, 'Visible', 'off');
set(handles.Image_Axes, 'Visible', 'on');
set(handles.Image_Axes, 'Box', 'on');
set(handles.TimeCourse_Slider, 'Visible', 'on');

% --- Executes on button press in ClusterMode_PushButton.
function ClusterMode_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClusterMode_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_SpineSegmentation

gui_SpineSegmentation.Mode = 'Cluster';

set(handles.TrackingMode_PushButton, 'Visible', 'off');
set(handles.ClusterMode_PushButton, 'Visible', 'off');
set(handles.Image_Axes, 'Visible', 'on');
set(handles.Image_Axes, 'Box', 'on');
set(handles.Image_Axes, 'YTick', []);
set(handles.Image_Axes, 'XTick', []);




% --- Executes on button press in ChangeModes_PushButton.
function ChangeModes_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to ChangeModes_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_SpineSegmentation

axes(handles.Image_Axes); cla;

set(handles.TrackingMode_PushButton, 'Visible', 'on');
set(handles.ClusterMode_PushButton, 'Visible', 'on');
set(handles.Image_Axes, 'Visible', 'off');
set(handles.TimeCourse_Slider, 'Visible', 'off');
gui_SpineSegmentation.LoadedImage = 0;
gui_SpineSegmentation.TC_Length = [];
gui_SpineSegmentation.Segmented = 0;


% --- Executes on slider movement.
function TimeCourse_Slider_Callback(hObject, eventdata, handles)
% hObject    handle to TimeCourse_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global gui_SpineSegmentation

frame = get(handles.TimeCourse_Slider, 'Value');
currentframe = floor(frame);

Im = gui_SpineSegmentation.Image;

axes(handles.Image_Axes);

branches = gui_SpineSegmentation.Spines;

if gui_SpineSegmentation.Segmented == 0 || gui_SpineSegmentation.ProcessedView == 0
    [counts, x] = imhist(Im(:,:,currentframe));

    [val pos] = max(counts);

    axes(handles.Image_Axes);

    min1 = min(min(Im(:,:,currentframe)));
    max1 = max(max(Im(:,:,currentframe)));

    imshow(Im(:,:,currentframe), [pos max1/4]);
elseif gui_SpineSegmentation.Segmented == 1
    imshow(Im{currentframe}(:,:,:));
%     for i = 1:size(branches{currentframe},1)
%         rectangle('Position', [branches{currentframe}(i,1), branches{currentframe}(i,2), 2,2], 'Curvature', [1 1], 'EdgeColor', 'b', 'Tag', ['Branchpoint ', num2str(i)], 'ButtonDownFcn', 'DeleteROI')
%     end
end



% --- Executes during object creation, after setting all properties.
function TimeCourse_Slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeCourse_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

global gui_SpineSegmentation



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

global gui_SpineSegmentation

gui_SpineSegmentation.Image_File = [];
gui_SpineSegmentation.Resized_Image = [];
gui_SpineSegmentation.Image = [];
gui_SpineSegmentation.Processed_Image = [];
gui_SpineSegmentation.Skeleton = [];
gui_SpineSegmentation.Final_Image_Data = [];
gui_SpineSegmentation.Final_Image = [];
gui_SpineSegmentation.Spines = [];
gui_SpineSegmentation.SpurNum = [];
gui_SpineSegmentation.Mode = [];
gui_SpineSegmentation.LoadedImage = 0;
gui_SpineSegmentation.TC_Length = [];
gui_SpineSegmentation.ProcessedView = 0;
gui_SpineSegmentation.Segmented = 1;

delete(hObject);


% --- Executes on button press in Unprocess_PushButton.
function Unprocess_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Unprocess_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_SpineSegmentation

gui_SpineSegmentation.Spines = [];

frame = get(handles.TimeCourse_Slider, 'Value');

gui_SpineSegmentation.ProcessedView = 0;


if isempty(gui_SpineSegmentation.Resized_Image)
    Im = gui_SpineSegmentation.Unprocessed_Image;

    gui_SpineSegmentation.Image = Im;

    [counts, x] = imhist(Im(:,:,frame));

    [val pos] = max(counts);

    axes(handles.Image_Axes); cla

    max1 = max(max(Im(:,:,frame)));

    imshow(Im(:,:,frame), [pos max1/4]);
else
    Im = gui_SpineSegmentation.Resized_Image;

    gui_SpineSegmentation.Image = Im;

    [counts, x] = imhist(Im(:,:,frame));

    [val pos] = max(counts);

    axes(handles.Image_Axes); cla

    max1 = max(max(Im(:,:,frame)));

    imshow(Im(:,:,frame), [pos max1/4]);
end



% --- Executes on button press in Reprocess_Pushbutton.
function Reprocess_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Reprocess_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_SpineSegmentation

currentframe = get(handles.TimeCourse_Slider, 'Value');

gui_SpineSegmentation.ProcessedView = 1;
gui_SpineSegmentation.Segmented = 0;


    Im = gui_SpineSegmentation.Processed_Image;
    gui_SpineSegmentation.Image = Im;

    [counts, x] = imhist(Im(:,:,currentframe));

    [val pos] = max(counts);

    axes(handles.Image_Axes);

    min1 = min(min(Im(:,:,currentframe)));
    max1 = max(max(Im(:,:,currentframe)));

    imshow(Im(:,:,currentframe), [pos max1/4]);



function GaussFilt_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to GaussFilt_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GaussFilt_EditableText as text
%        str2double(get(hObject,'String')) returns contents of GaussFilt_EditableText as a double


% --- Executes during object creation, after setting all properties.
function GaussFilt_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GaussFilt_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MedFilt_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to MedFilt_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MedFilt_EditableText as text
%        str2double(get(hObject,'String')) returns contents of MedFilt_EditableText as a double


% --- Executes during object creation, after setting all properties.
function MedFilt_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MedFilt_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MinBlob_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to MinBlob_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinBlob_EditableText as text
%        str2double(get(hObject,'String')) returns contents of MinBlob_EditableText as a double


% --- Executes during object creation, after setting all properties.
function MinBlob_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinBlob_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MaxBlob_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to MaxBlob_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxBlob_EditableText as text
%        str2double(get(hObject,'String')) returns contents of MaxBlob_EditableText as a double


% --- Executes during object creation, after setting all properties.
function MaxBlob_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxBlob_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MinSpine_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to MinSpine_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinSpine_EditableText as text
%        str2double(get(hObject,'String')) returns contents of MinSpine_EditableText as a double


% --- Executes during object creation, after setting all properties.
function MinSpine_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinSpine_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MaxSpine_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to MaxSpine_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxSpine_EditableText as text
%        str2double(get(hObject,'String')) returns contents of MaxSpine_EditableText as a double


% --- Executes during object creation, after setting all properties.
function MaxSpine_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxSpine_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ArbShape_ToggleButton.
function ArbShape_ToggleButton_Callback(hObject, eventdata, handles)
% hObject    handle to ArbShape_ToggleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_SpineSegmentation
axes(handles.Image_Axes);
Arb = get(handles.ArbShape_ToggleButton, 'Value');

if Arb == 1
    set(handles.DrawBox_ToggleButton, 'Value', 0);
    Im = gui_SpineSegmentation.Image;
    Im = double(Im);

    frames = size(Im,3);

    if frames == 1
        [mask xi yi] = roipoly(Im);
        mask = uint16(mask);
        restricted_Im = Im.*mask;
    else
    %     for i = 1:frames
    %         [mask(:,:,i) xi{i} yi{i}] = roipoly
    %         mask(:,:,i) = uint16(mask(:,:,i));
    %         restricted_Im{i} = Im(:,:,i).*mask(:,:,i)
    %     end
            [mask(:,:,1) xi{1} yi{1}] = roipoly;
            mask(:,:,1) = uint16(mask(:,:,1));
            restricted_Im(:,:,1) = Im(:,:,1).*mask(:,:,1);
    end

    ystart = round(min(yi{1}));
    yend = round(max(yi{1}));
    xstart = round(min(xi{1}));
    xend = round(max(xi{1}));

    restricted_Im = restricted_Im(ystart:yend,xstart:xend);
    imshow(restricted_Im, []);

    gui_SpineSegmentation.Resized_Image = restricted_Im;
    gui_SpineSegmentation.Arbitrary_Shaping = 1;
    set(handles.ArbShape_ToggleButton, 'Value', 0);
end


% --- Executes on button press in ExcludeRegion_PushButton.
function ExcludeRegion_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to ExcludeRegion_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_SpineSegmentation

axes(handles.Image_Axes)

BlackoutMask = roipoly;
BlackoutMask = imcomplement(BlackoutMask);

Im = gui_SpineSegmentation.Image;
LabMat = gui_SpineSegmentation.Label_Matrix;
frames = gui_SpineSegmentation.TC_Length;
currentframe = get(handles.TimeCourse_Slider, 'Value');

if gui_SpineSegmentation.Segmented == 0
    for i = 1:frames
        NewIm(:,:,i) = Im(:,:,i).*BlackoutMask;
        NewLabel(:,:,i) = LabMat(:,:,i).*BlackoutMask;
    end
    [counts, x] = imhist(NewIm(:,:,currentframe));

    [val pos] = max(counts);

    axes(handles.Image_Axes);

    max1 = max(max(NewIm(:,:,currentframe)));

    imshow(NewIm(:,:,currentframe),[pos max1/4]);
else
    for i = 1:3
        for j = 1:frames
            NewIm{j}(:,:,i) = Im{j}(:,:,i).*BlackoutMask;
        end
    end
    imshow(NewIm{currentframe});
end

if gui_SpineSegmentation.ProcessedView == 0
    gui_SpineSegmentation.Image = NewIm;
    gui_SpineSegmentation.Unprocessed_Image = NewIm;
else
    gui_SpineSegmentation.Image = NewIm;
    gui_SpineSegmentation.Processed_Image = NewIm;
    gui_SpineSegmentation.Label_Matrix = NewLabel;
end


% --- Executes on button press in Segmented_PushButton.
function Segmented_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Segmented_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_SpineSegmentation

currentframe = get(handles.TimeCourse_Slider, 'Value');

gui_SpineSegmentation.ProcessedView = 1;
gui_SpineSegmentation.Segmented = 1;

    Im = gui_SpineSegmentation.Final_Image_Data{currentframe};
    gui_SpineSegmentation.Image = gui_SpineSegmentation.Final_Image_Data;
    imshow(Im);


    


% --- Executes on button press in Align_PushButton.
function Align_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Align_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global gui_SpineSegmentation

frames = gui_SpineSegmentation.TC_Length;
Image = gui_SpineSegmentation.Image;
LabMat = gui_SpineSegmentation.Label_Matrix;
h = waitbar(0, 'Aligning...');
Im(:,:,1) = Image(:,:,1);
Lab(:,:,1) = LabMat(:,:,1);

for i = 2:frames
    A = Image(:,:,i-1);
    B = Image(:,:,i);
    ptsOriginal = detectSURFFeatures(A);
    ptsDistorted = detectSURFFeatures(B);
    [featuresOriginal, validPtsOriginal] = extractFeatures(A, ptsOriginal);
    [featuresDistorted, validPtsDistorted] = extractFeatures(B,ptsDistorted);
    indexPairs = matchFeatures(featuresOriginal, featuresDistorted);
    matchedOriginal = validPtsOriginal(indexPairs(:,1));
    matchedDistorted = validPtsDistorted(indexPairs(:,2));
    [tform, inlierDistorted, inlierOriginal] = estimateGeometricTransform(matchedDistorted, matchedOriginal, 'similarity');
    outputView = imref2d(size(A));
    recovered = imwarp(B, tform, 'OutputView', outputView);
    recoveredLabels = imwarp(LabMat(:,:,i), tform, 'OutputView', outputView);
    Im(:,:,i) = recovered;
    Lab(:,:,i) = recoveredLabels;
    waitbar(i/frames, h, ['Aligning frame ', num2str(i)])
end

close(h)
gui_SpineSegmentation.Image = Im;
gui_SpineSegmentation.Processed_Image = Image;
gui_SpineSegmentation.Label_Matrix = Lab;


