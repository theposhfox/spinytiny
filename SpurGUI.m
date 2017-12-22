function varargout = SpurGUI(varargin)
% SPURGUI MATLAB code for SpurGUI.fig
%      SPURGUI, by itself, creates a new SPURGUI or raises the existing
%      singleton*.
%
%      H = SPURGUI returns the handle to a new SPURGUI or the handle to
%      the existing singleton*.
%
%      SPURGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPURGUI.M with the given input arguments.
%
%      SPURGUI('Property','Value',...) creates a new SPURGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SpurGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SpurGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SpurGUI

% Last Modified by GUIDE v2.5 03-Feb-2016 15:17:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SpurGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SpurGUI_OutputFcn, ...
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


% --- Executes just before SpurGUI is made visible.
function SpurGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SpurGUI (see VARARGIN)

% Choose default command line output for SpurGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global gui_SpineSegmentation;

skeleton_im = gui_SpineSegmentation.Skeleton;

datatype = class(skeleton_im);

if strcmpi(datatype, 'cell')
    currentframe = get(gui_SpineSegmentation.handles.TimeCourse_Slider, 'Value');
    skeleton_im = skeleton_im{currentframe};
end

axes(handles.Skeleton_Axes);
imshow(skeleton_im);
title(gui_SpineSegmentation.SpurMessage)

% UIWAIT makes SpurGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
uiwait

% --- Outputs from this function are returned to the command line.
function varargout = SpurGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes on slider movement.
function SpurNum_Slider_Callback(hObject, eventdata, handles)
% hObject    handle to SpurNum_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global gui_SpineSegmentation

skeleton_im = gui_SpineSegmentation.Skeleton;

datatype = class(skeleton_im);

if strcmpi(datatype, 'cell')
    currentframe = get(gui_SpineSegmentation.handles.TimeCourse_Slider, 'Value');
    skeleton_im = skeleton_im{currentframe};
end

mode = get(handles.Dilate_ToggleButton, 'Value');

if mode == 0
    SpurNum = get(handles.SpurNum_Slider, 'Value');
    SpurNum = floor(SpurNum);
    rem_branches = bwmorph(skeleton_im, 'spur', SpurNum); %%% de-spur the image 'n' times; this should get rid of branches for analysis of the main dendrite only
    axes(handles.Skeleton_Axes); cla;
    imshow(rem_branches);
    title(gui_SpineSegmentation.SpurMessage)
else
    rem_branches = bwmorph(skeleton_im, 'spur', gui_SpineSegmentation.SpurNum);
    DilateNum = get(handles.SpurNum_Slider, 'Value');
    DilateNum = floor(DilateNum);
    fatDend = imdilate(rem_branches, strel('disk', DilateNum));
    imshow(fatDend);
    gui_SpineSegmentation.DilateNum = DilateNum;
end
    


% --- Executes during object creation, after setting all properties.
function SpurNum_Slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SpurNum_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in Accept_Pushbutton.
function Accept_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Accept_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_SpineSegmentation

SpurNum = get(handles.SpurNum_Slider, 'Value');

gui_SpineSegmentation.SpurNum = SpurNum;

if strcmpi(gui_SpineSegmentation.SpurMessage, 'First, remove all branchpoints to get just the dendrite')
    gui_SpineSegmentation.SpurMessage = 'Now, dilate the dendrite to get it back to the original size';
    title('Now, dilate the dendrite to get it back to the original size')
    set(handles.Dilate_ToggleButton, 'Value', 1);
    set(handles.SpurNum_Slider, 'Value', 0);
    gui_SpineSegmentation.Isolated_Dendrite = [];
    return
end 

uiresume;

figure1_CloseRequestFcn(handles.figure1, eventdata,handles);





% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

delete(hObject);


% --- Executes on button press in Dilate_ToggleButton.
function Dilate_ToggleButton_Callback(hObject, eventdata, handles)
% hObject    handle to Dilate_ToggleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Dilate_ToggleButton

set(handles.SpurNum_Slider, 'Value', 0)


% --- Executes on button press in ConnectDendrites_ToggleButton.
function ConnectDendrites_ToggleButton_Callback(hObject, eventdata, handles)
% hObject    handle to ConnectDendrites_ToggleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ConnectDendrites_ToggleButton

global gui_SpineSegmentation

skeleton_im = gui_SpineSegmentation.Skeleton;

SpurNum = gui_SpineSegmentation.SpurNum;

    SpurNum = floor(SpurNum);
    rem_branches = bwmorph(skeleton_im, 'spur', SpurNum); %%% de-spur the image 'n' times; this should get rid of branches for analysis of the main dendrite only
    axes(handles.Skeleton_Axes); cla;
    imshow(rem_branches);

shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'white')

[x y] = roipoly;
bridge = [x(1) y(1) x(2) y(2)];

if isempty(gui_SpineSegmentation.Dendrite_Bridges)
    gui_SpineSegmentation.Dendrite_Bridges = bridge;
else
    gui_SpineSegmentation.Dendrite_Bridges = [gui_SpineSegmentation.Dendrite_Bridges; bridge];
end

axes(handles.Skeleton_Axes);

step(skeleton_im, bridge);

uiwait


