function varargout = CaImageViewer(varargin)
% CAIMAGEVIEWER MATLAB code for CaImageViewer.fig
%      CAIMAGEVIEWER, by itself, creates a new CAIMAGEVIEWER or raises the existing
%      singleton*.
%
%      H = CAIMAGEVIEWER returns the handle to a new CAIMAGEVIEWER or the handle to
%      the existing singleton*.
%
%      CAIMAGEVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CAIMAGEVIEWER.M with the given input arguments.
%
%      CAIMAGEVIEWER('Property','Value',...) creates a new CAIMAGEVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CaImageViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CaImageViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CaImageViewer

% Last Modified by GUIDE v2.5 01-Dec-2017 10:17:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CaImageViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @CaImageViewer_OutputFcn, ...
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


% --- Executes just before CaImageViewer is made visible.
function CaImageViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CaImageViewer (see VARARGIN)

% Choose default command line output for CaImageViewer
global gui_CaImageViewer

handles.output = hObject;

gui_CaImageViewer.figure.handles = handles;

Users = {'Nathan', 'Zhongmin', 'Giulia', 'Mingyuan', 'Assaf'}; % Assaf was add as a user.

Scrsz = get(0, 'Screensize');
    dialogboxwidth = 315;
    dialogboxheight = 150;
    d = dialog('Position', [(Scrsz(3)/2)-dialogboxwidth/2 Scrsz(4)/2-dialogboxheight/2 dialogboxwidth dialogboxheight ], 'Name', 'User');
    txt = uicontrol('Parent', d, 'Style', 'text', 'Position', [15 100 270 30], 'String', 'Select User:');
    
    btn1 = uicontrol('Parent', d, 'Style', 'pushbutton', 'Position', [15 30 50 25], 'String', Users{1}, 'Callback', @UserName);
    btn2 = uicontrol('Parent', d, 'Style', 'pushbutton', 'Position', [65.5 30 70 25], 'String', Users{2}, 'Callback', @UserName);
    btn3 = uicontrol('Parent', d, 'Style', 'pushbutton', 'Position', [136 30 50 25], 'String', Users{3}, 'Callback', @UserName);
    btn4 = uicontrol('Parent', d, 'Style', 'pushbutton', 'Position', [186.5 30 70 25], 'String', Users{4}, 'Callback', @UserName);
    btn5 = uicontrol('Parent', d, 'Style', 'pushbutton', 'Position', [257 30 50 25], 'String', Users{5}, 'Callback', @UserName);
    
    uiwait(d)
    choice = get(d, 'UserData');
    set(gui_CaImageViewer.figure.handles.figure1, 'UserData', choice);
    delete(d);


% Update handles structure
guidata(hObject, handles);

%%% Set appearance of GUI
set(handles.GreenGraph, 'YTick', []);
set(handles.GreenGraph, 'XTick', []);
set(handles.GreenGraph, 'Box', 'on');
set(handles.RedGraph, 'YTick', []);
set(handles.RedGraph, 'XTick', []);
set(handles.RedGraph, 'Box', 'on');

%%% Initialize Key press functions for various editable text boxes
set(gui_CaImageViewer.figure.handles.Frame_EditableText, 'KeyPressFcn', @frameset);
set(gui_CaImageViewer.figure.handles.UpperLUT_EditableText, 'KeyPressFcn', @frameset);
set(gui_CaImageViewer.figure.handles.RedUpperLUT_EditableText, 'KeyPressFcn', @frameset);
set(gui_CaImageViewer.figure.handles.LowerLUT_EditableText, 'KeyPressFcn', @frameset);
set(gui_CaImageViewer.figure.handles.RedLowerLUT_EditableText, 'KeyPressFcn', @frameset);
set(gui_CaImageViewer.figure.handles.SmoothingFactor_EditableText, 'KeyPressFcn', @frameset);
set(gui_CaImageViewer.figure.handles.GreenGamma_EditableText, 'KeyPressFcn', @frameset);
set(gui_CaImageViewer.figure.handles.RedGamma_EditableText, 'KeyPressFcn', @frameset);


%%% Initialize Random Parameters
gui_CaImageViewer.CurrentCMap = 'RGB';
gui_CaImageViewer.NewSpineAnalysis = 0;
gui_CaImageViewer.NewSpineAnlysisInfo.CurrentDate = [];
gui_CaImageViewer.NewSpineAnalysisInfo.CurrentImagingField = [];
gui_CaImageViewer.NewSpineAnalysisInfo.SpineList = [];


%%%%% Clear All ROIs and Associated Labels %%%%%%

clearROIs = findobj('Type', 'rectangle');
for i = 1:length(clearROIs);
    delete(clearROIs(i));
end
clearTexts = findobj('Type', 'Text');

clearLines = findobj('Type', 'Line');
delete(clearLines)
    
gui_CaImageViewer.ROI = [];
gui_CaImageViewer.ROItext = [];
gui_CaImageViewer.PolyROI = [];
gui_CaImageViewer.PolyLinePos = [];
gui_CaImageViewer.PolyLineVertices = [];
gui_CaImageViewer.Spine_Number = 0;
gui_CaImageViewer.Dendrite_Number = 0;
gui_CaImageViewer.DendritePolyPointNumber = [];
gui_CaImageViewer.SpineDendriteGrouping = [];
gui_CaImageViewer.Dendrite_ROIs = 0;
gui_CaImageViewer.LoadedFile = 0;
gui_CaImageViewer.filename = [];
gui_CaImageViewer.GCaMP_Image = [];
gui_CaImageViewer.Red_Image = [];
gui_CaImageViewer.ch1image = [];
gui_CaImageViewer.imageserieslength = [];

gui_CaImageViewer.GreenGraph_loc = get(handles.GreenGraph, 'Position');
gui_CaImageViewer.RedGraph_loc = get(handles.RedGraph, 'Position');

function UserName(hObject, eventdata, ~)

button = get(hObject);

choice = button.String;

sourcewindow = button.Parent;

set(sourcewindow, 'UserData', choice);

uiresume



% UIWAIT makes CaImageViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CaImageViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on slider movement.
function ImageSlider_Slider_Callback(hObject, eventdata, handles)
% hObject    handle to ImageSlider_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global gui_CaImageViewer

ImageNum = str2num(get(gui_CaImageViewer.figure.handles.Frame_EditableText, 'String'));

CaImageSlider(ImageNum);

% --- Executes during object creation, after setting all properties.
function ImageSlider_Slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageSlider_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% figure_children = get(gui_CaImageViewer.figure.handles.output, 'Children');
% figure_axes = findobj(figure_children, 'Type', 'axes');
% set(findobj(figure_children, 'Type', 'axes'), 'ButtonDownFcn', @Ca_DrawROI);
% red_axes_children = get(figure_axes(1), 'Children');
% green_axes_children = get(figure_axes(2), 'Children');

% set(findobj(green_axes_children, 'Type', 'Image'), 'ButtonDownFcn', @Ca_DrawROI);
% set(findobj(red_axes_children, 'Type', 'Image'), 'ButtonDownFcn', @Ca_DrawROI);


function Frame_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to Frame_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Frame_EditableText as text
%        str2double(get(hObject,'String')) returns contents of Frame_EditableText as a double


% --- Executes during object creation, after setting all properties.
function Frame_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Frame_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global gui_CaImageViewer
    


% --- Executes on button press in CalcROIs_PushButton.
function CalcROIs_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to CalcROIs_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


CaCalculateROIs



function UpperLUT_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to UpperLUT_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UpperLUT_EditableText as text
%        str2double(get(hObject,'String')) returns contents of UpperLUT_EditableText as a double


% --- Executes during object creation, after setting all properties.
function UpperLUT_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UpperLUT_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LowerLUT_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to LowerLUT_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LowerLUT_EditableText as text
%        str2double(get(hObject,'String')) returns contents of LowerLUT_EditableText as a double


% --- Executes during object creation, after setting all properties.
function LowerLUT_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LowerLUT_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SmoothingFactor_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to SmoothingFactor_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SmoothingFactor_EditableText as text
%        str2double(get(hObject,'String')) returns contents of SmoothingFactor_EditableText as a double


% --- Executes during object creation, after setting all properties.
function SmoothingFactor_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SmoothingFactor_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ViewMovie_Checkbox.
function ViewMovie_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to ViewMovie_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ViewMovie_Checkbox
global gui_CaImageViewer


% --- Executes on button press in BackgroundROI_ToggleButton.
function BackgroundROI_ToggleButton_Callback(hObject, eventdata, handles)
% hObject    handle to BackgroundROI_ToggleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BackgroundROI_ToggleButton

global gui_CaImageViewer

% set(gui_CaImageViewer.figure.handles.ZoomIn_ToggleTool, 'state', 'off')
% set(gui_CaImageViewer.figure.handles.ZoomOut_ToggleTool, 'state', 'off')
% set(gui_CaImageViewer.figure.handles.Pan_ToggleTool, 'state', 'off')

BackgroundROI = get(gui_CaImageViewer.figure.handles.BackgroundROI_ToggleButton, 'Value');
SpineROI = get(gui_CaImageViewer.figure.handles.SpineROI_ToggleButton, 'Value');
NearbySpineROI = get(gui_CaImageViewer.figure.handles.NearbySpine_ToggleButton, 'Value');
Dendrite_PolyLines = get(gui_CaImageViewer.figure.handles.DendritePolyLines_ToggleButton, 'Value');
Router = 'Background';

if BackgroundROI == 1
    set(handles.SpineROI_ToggleButton, 'Value', 0);
    set(handles.NearbySpine_ToggleButton, 'Value', 0);
    set(handles.DendritePolyLines_ToggleButton, 'Value', 0);
    set(gui_CaImageViewer.figure.handles.output, 'WindowButtonDownFcn', {@DrawROI, 0, Router});
end

if BackgroundROI == 0 && SpineROI == 0 && NearbySpineROI == 0 && Dendrite_PolyLines == 0
    set(gui_CaImageViewer.figure.handles.output, 'WindowButtonDownFcn', []);
end


% --- Executes on button press in SpineROI_ToggleButton.
function SpineROI_ToggleButton_Callback(hObject, eventdata, handles)
% hObject    handle to SpineROI_ToggleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SpineROI_ToggleButton

global gui_CaImageViewer

BackgroundROI = get(gui_CaImageViewer.figure.handles.BackgroundROI_ToggleButton, 'Value');
SpineROI = get(gui_CaImageViewer.figure.handles.SpineROI_ToggleButton, 'Value');
NearbySpineROI = get(gui_CaImageViewer.figure.handles.NearbySpine_ToggleButton, 'Value');
Dendrite_PolyLines = get(gui_CaImageViewer.figure.handles.DendritePolyLines_ToggleButton, 'Value');
ROInum = gui_CaImageViewer.Spine_Number + 1;

% if ~isempty(gui_CaImageViewer.ROI)
%     for i = 1:ROInum-1
%         set(gui_CaImageViewer.ROI(i), 'ButtonDownFcn', [])
%     end
% else
%     ClearROIs
%     set(gui_CaImageViewer.figure.handles.SpineROI_ToggleButton, 'Value', 0)
%     error('Make sure to draw the background ROI first...')
% end

Router = 'Spine';

if SpineROI == 1
    set(handles.BackgroundROI_ToggleButton, 'Value', 0);
    set(handles.NearbySpine_ToggleButton, 'Value', 0);
    set(handles.DendritePolyLines_ToggleButton, 'Value', 0);
    set(gui_CaImageViewer.figure.handles.output, 'WindowButtonDownFcn', {@DrawROI, ROInum, Router})
end

if BackgroundROI == 0 && SpineROI == 0 && NearbySpineROI == 0 && Dendrite_PolyLines == 0
    set(gui_CaImageViewer.figure.handles.output, 'WindowButtonDownFcn', []);
end

% --- Executes on button press in ClearROIS_PushButton.
function ClearROIS_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClearROIS_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


ClearROIs


% --- Executes on button press in NearbySpine_ToggleButton.
function NearbySpine_ToggleButton_Callback(hObject, eventdata, handles)
% hObject    handle to NearbySpine_ToggleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NearbySpine_ToggleButton

global gui_CaImageViewer;

BackgroundROI = get(gui_CaImageViewer.figure.handles.BackgroundROI_ToggleButton, 'Value');
SpineROI = get(gui_CaImageViewer.figure.handles.SpineROI_ToggleButton, 'Value');
NearbySpineROI = get(gui_CaImageViewer.figure.handles.NearbySpine_ToggleButton, 'Value');
Dendrite_PolyLines = get(gui_CaImageViewer.figure.handles.DendritePolyLines_ToggleButton, 'Value');
ROInum = gui_CaImageViewer.Spine_Number + 1;
Router = 'Nearby';

if NearbySpineROI == 1
    set(handles.BackgroundROI_ToggleButton, 'Value', 0);
    set(handles.SpineROI_ToggleButton, 'Value', 0);
    set(handles.DendritePolyLines_ToggleButton, 'Value', 0);
    set(gui_CaImageViewer.figure.handles.output, 'WindowButtonDownFcn', {@DrawROI, ROInum, Router})
end

if BackgroundROI == 0 && SpineROI == 0 && NearbySpineROI == 0 && Dendrite_PolyLines == 0
    set(gui_CaImageViewer.figure.handles.output, 'WindowButtonDownFcn', []);
end

% --- Executes on button press in DendritePolyLines_ToggleButton.
function DendritePolyLines_ToggleButton_Callback(hObject, eventdata, handles)
% hObject    handle to DendritePolyLines_ToggleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DendritePolyLines_ToggleButton

global gui_CaImageViewer;


BackgroundROI = get(gui_CaImageViewer.figure.handles.BackgroundROI_ToggleButton, 'Value');
SpineROI = get(gui_CaImageViewer.figure.handles.SpineROI_ToggleButton, 'Value');
NearbySpineROI = get(gui_CaImageViewer.figure.handles.NearbySpine_ToggleButton, 'Value');
Dendrite_PolyLines = get(gui_CaImageViewer.figure.handles.DendritePolyLines_ToggleButton, 'Value');
CurrentDendNum = gui_CaImageViewer.Dendrite_Number+1;

if Dendrite_PolyLines == 1
    set(handles.BackgroundROI_ToggleButton, 'Value', 0);
    set(handles.SpineROI_ToggleButton, 'Value', 0);
    set(handles.NearbySpine_ToggleButton, 'Value', 0);
        DendriteNum = inputdlg({'Dendrite number:'}, 'Input', 1, {num2str(CurrentDendNum)});
    DendriteNum = str2num(DendriteNum{1});
    set(gui_CaImageViewer.figure.handles.output, 'WindowButtonDownFcn', {@DrawPolyLines, DendriteNum});
end

if BackgroundROI == 0 && SpineROI == 0 && NearbySpineROI == 0 && Dendrite_PolyLines == 0
    set(gui_CaImageViewer.figure.handles.output, 'WindowButtonDownFcn', []);
end



function RedUpperLUT_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to RedUpperLUT_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RedUpperLUT_EditableText as text
%        str2double(get(hObject,'String')) returns contents of RedUpperLUT_EditableText as a double


% --- Executes during object creation, after setting all properties.
function RedUpperLUT_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RedUpperLUT_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RedLowerLUT_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to RedLowerLUT_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RedLowerLUT_EditableText as text
%        str2double(get(hObject,'String')) returns contents of RedLowerLUT_EditableText as a double


% --- Executes during object creation, after setting all properties.
function RedLowerLUT_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RedLowerLUT_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Zoom_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to Zoom_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Zoom_EditableText as text
%        str2double(get(hObject,'String')) returns contents of Zoom_EditableText as a double


% --- Executes during object creation, after setting all properties.
function Zoom_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Zoom_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function CodeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to CodeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

edit CaImageViewer



function GreenGamma_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to GreenGamma_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GreenGamma_EditableText as text
%        str2double(get(hObject,'String')) returns contents of GreenGamma_EditableText as a double


% --- Executes during object creation, after setting all properties.
function GreenGamma_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GreenGamma_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RedGamma_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to RedGamma_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RedGamma_EditableText as text
%        str2double(get(hObject,'String')) returns contents of RedGamma_EditableText as a double


% --- Executes during object creation, after setting all properties.
function RedGamma_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RedGamma_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TwoChannels_CheckBox.
function TwoChannels_CheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to TwoChannels_CheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TwoChannels_CheckBox


% --- Executes on button press in RecoverROIs_PushButton.
function RecoverROIs_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to RecoverROIs_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

program = get(gcf);

running = program.FileName;

if ~isempty(regexp(running, 'CaImageViewer'))
    global gui_CaImageViewer
    glovar = gui_CaImageViewer;
    axes1 = glovar.figure.handles.GreenGraph;
    axes2 = glovar.figure.handles.RedGraph;
elseif ~isempty(regexp(running, 'FluorescenceSuite'));
    global gui_FluorescenceSuite
    glovar = gui_FluorescenceSuite;
    axes1 = glovar.figure.handles.Green_Axes;
    axes2 = glovar.figure.handles.Red_Axes;
end

file = gui_CaImageViewer.filename;
file = file(1:end-4);
experiment = regexp(gui_CaImageViewer.filename, '[A-Z]{2}\d+[_]\d+', 'match');
experiment = experiment{1};
animal = experiment(1:5);
date = experiment(7:end);

twochannels = get(glovar.figure.handles.TwoChannels_CheckBox, 'Value');
fname = [];

if ispc
    save_directory = gui_CaImageViewer.save_directory;
else
    nameparts = regexp(gui_CaImageViewer.save_directory, filesep, 'split');
    linuxstarter = '/usr/local/lab/';
    save_directory = [linuxstarter, nameparts{2}, filesep, nameparts{3}, filesep, nameparts{4}, filesep, nameparts{5}, filesep, nameparts{6}, filesep, nameparts{7},filesep];
end

try
    cd(save_directory)
    folder = dir(save_directory);
catch
    [fname pname] = uigetfile();
    cd(pname)
    folder = dir(pname);
    glovar.save_directory = pname;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find all files that hold ROI positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    count = 1;
    roifile = [];
    if gui_CaImageViewer.NewSpineAnalysis
        for i = 1:length(folder)
            match = regexp(folder(i).name, 'NewSpineAnalysisROIs');
            if ~isempty(match)
                roifile{count} = folder(i).name;
                count = count+1;
            end
        end
    else
        for i = 1:length(folder)
            match = regexp(folder(i).name, [experiment, '_SavedROIs']);
            if ~isempty(match)
                roifile{count} = folder(i).name;
                count = count+1;
            end
        end
    end

    if length(roifile)>1
        for i = 1:length(roifile)
            temp = regexp(roifile{i}, ['DrawnBy'], 'split');
            if length(temp)>1           %%% If the file doesn't contain the "DrawnBy" tag, then it will only return 1 answer, and therefore was made prior to qualifying the file according to username
                useroption{i} = temp{2}(1:end-4);
            else
                useroption{i} = 'undefined';
            end
        end
        choice = questdlg('Which color do you prefere?', ... % Add by Assaf for choosing colors (Red, Green or Blue) for the ROIs
	'ROIs color', ...
	'Red','Green', 'Blue', 'Blue');
% Handle response
switch choice
    case 'Red'
        color_by_user = [1 0 0];
    case 'Green'
        color_by_user = [0 1 0];
    case 'Green'
        color_by_user = [0 0 1];
end

        Scrsz = get(0, 'Screensize');
        d = dialog('Position', [(Scrsz(3)/2)-125 Scrsz(4)/2-75 250 150], 'Name', 'Found multiple ROI files...');
        txt = uicontrol('Parent', d, 'Style', 'text', 'Position', [10 100 230 30], 'String', 'Load ROIs drawn by whom?');
        for j = 1:length(roifile)
            uicontrol('Parent', d, 'Style', 'pushbutton', 'Position', [30+((j-1)*55) 30 50 25], 'String', useroption{j}, 'Callback', @UserName)
        end
        uiwait(d)
        choice = get(d, 'UserData');
        drawer = choice;
        delete(d);
        roifilename = [experiment, '_SavedROIs_DrawnBy', choice];
    else
        roifilename = roifile{1}(1:end-4);
    end
        load(roifilename)
        savedFile = roifilename;
        disp('Successfully pulled ROIs from saved ROI file')
catch
    try
        load(roifilename)
        savedFile = roifilename;
        disp('Successfully pulled ROIs from saved ROI file')
    catch
        try
            cd(['Z:\People\Nathan\Data\', animal])
        catch
            cd(['/usr/local/lab/People/Nathan/Data/', animal])
        end
        try
            load([file, '_Analyzed'])
            savedFile = [file, '_Analyzed'];
            disp('Couldn''t find ROI file, but could pull ROIs from analyzed file')
        catch
            try
                [fname pname] = uigetfile();
                file = fname(1:end-4);
                savedFile = file;
                cd(pname)
                load(file)
            catch
                disp('This file has not yet been analyzed, so no ROIs are saved...')
                return
            end
        end
    end
end

try
    eval(['savedFile =', savedFile]);
catch
    if ~isempty(savedFile)
        temp = who(['*',savedFile(1:5), '*']);
        savedFile = temp{1};
        disp('File name discrepancy; using closest available')
        eval(['savedFile =', savedFile]);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Draw loaded ROIs %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ROIs = savedFile.ROIPosition; % how to categorize newly added spines?-ZL

glovar.Spine_Number = length(ROIs)-1;

for a = 1:length(ROIs)
    if ~isempty(ROIs{a})
        ROInum = a-1;
        axes(axes1);
        if glovar.NewSpineAnalysis
            c = uicontextmenu;
            uimenu(c, 'Label', 'Set as eliminated', 'Callback', @CategorizeSpines);
            uimenu(c, 'Label', 'Set as active', 'Callback', @CategorizeSpines);
            glovar.ROI(a) = rectangle('Position', ROIs{a}, 'EdgeColor', color_by_user, 'Curvature', [1 1],'Tag', ['ROI', num2str(ROInum)], 'ButtonDownFcn', {@DragROI, ROInum, 'HomeWindow'}, 'Linewidth', 1, 'UIContextMenu', c); % Assaf changed [0.2 0.4 0.9] to the variable color_by_user, that is defined by the user choice of color (lines 681-602) 
        else
            glovar.ROI(a) = rectangle('Position', ROIs{a}, 'EdgeColor', color_by_user, 'Curvature', [1 1],'Tag', ['ROI', num2str(ROInum)], 'ButtonDownFcn', {@DragROI, ROInum, 'HomeWindow'}, 'Linewidth', 1);
        end
            glovar.ROItext(a) = text(ROIs{a}(1)-6, ROIs{a}(2)-4, num2str(a-1), 'color', 'white', 'Tag', ['ROI', num2str(a-1), ' Text'],'ButtonDownFcn', 'DeleteROI', 'Fontsize', 6);
        if twochannels == 1
            axes(axes2);
            glovar.ROIred(a) = rectangle('Position', ROIs{a}, 'EdgeColor', 'red', 'Curvature', [1 1],'Tag', ['ROIred', num2str(ROInum)],'ButtonDownFcn', {@DragROI, ROInum, 'HomeWindow'});
            glovar.ROIredtext(a) = text(ROIs{a}(1)-2, ROIs{a}(2)-2, num2str(a-1), 'color', 'white', 'Tag', ['ROI', num2str(a-1), ' Text'],'ButtonDownFcn', 'DeleteROI');
        else
        end
    end
end

%%% New spine analysis section

if glovar.NewSpineAnalysis
    terminus = regexp(save_directory, animal, 'end');
    targ_folder = save_directory(1:terminus);
    currentfield = glovar.NewSpineAnalysisInfo.CurrentImagingField;
    load([targ_folder, filesep,'Imaging Field ', num2str(currentfield), ' Spine Registry'])
    instanceofappearance = logical(strcmpi(SpineRegistry.DatesAcquired, gui_CaImageViewer.NewSpineAnalysisInfo.CurrentDate));
    glovar.NewSpineAnalysisInfo.SpineList = ones(length(ROIs)-1,1); %%% Don't forget the first ROI is always the background ROI!
    if size(SpineRegistry.Data,2)>=find(instanceofappearance) %% && find(instanceofappearance)~=1 %%% ZL commentm, it is possible need to set another category of spines specifying the "true new spines"
        r = find(SpineRegistry.Data(:,instanceofappearance)==0);
        for i = 1:length(r)
            set(findobj(glovar.figure.handles.GreenGraph, 'Type', 'rectangle', 'Tag', ['ROI', num2str(r(i))]), 'FaceColor', [1 0 0])
        end
        glovar.NewSpineAnalysisInfo.SpineList(r) = 0;
    else
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Handle dendritic ROIs %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DendNum = savedFile.NumberofDendrites; glovar.Dendrite_Number = DendNum;
coordinates = savedFile.PolyLinePosition;
PPnum = cumsum(savedFile.DendritePolyPointNumber);
glovar.SpineDendriteGrouping = savedFile.SpineDendriteGrouping;

radius = 3;
x = [];
x_red = [];
y = [];
y_red = [];

pause(0.05)

try
    if isfield(glovar, 'PolyROI') && ~isempty(coordinates)
        PPsperDend = savedFile.DendritePolyPointNumber;
        if PPsperDend == 0
            warning('Dendrite ROIs could not be returned and must be redrawn')
            glovar.Dendrite_Number = 0;
            if ~isempty(regexp(running, 'CaImageViewer'))
                gui_CaImageViewer = glovar;
            elseif ~isempty(regexp(running, 'FluorescenceSuite'));
                gui_FluorescenceSuite = glovar;
            end
            return
        else
            glovar.DendritePolyPointNumber = PPsperDend;
            Dendrite_ROIs = length(savedFile.PolyROI);
        end
        if ~isempty(coordinates{1})
        axes(axes1);
        currDend = 1;
            for i = 1:length(coordinates)
                glovar.PolyLinePos{i} = [coordinates{i}(1), coordinates{i}(2), radius*2, radius*2];
                glovar.PolyROI{i} = rectangle('Position', glovar.PolyLinePos{i}, 'EdgeColor', 'red', 'Tag', ['Dendrite ', num2str(currDend), ' PolyROI ', num2str(i)], 'Curvature', [1 1], 'ButtonDownFcn', 'Drag_Poly');
                x = [x,coordinates{i}(1)+radius];
                y = [y,coordinates{i}(2)+radius];
                if i < sum(PPsperDend(1:currDend))
                    currDend = currDend;
                else
                    currDend = currDend+1;
                end
            end
            if DendNum == 1
%                 for i = 1:2:length(glovar.PolyLinePos)
                    hold on;
                    plot(x,y, 'color', 'cyan', 'Tag', ['PolyLine 1']);
%                 end
            else
                counter = 1;
                for i = 1:DendNum
                    hold on;
                    plot(x(counter:(counter+PPsperDend(i)-1)),y(counter:(counter+PPsperDend(i)-1)), 'color', 'cyan', 'Tag', ['PolyLine ', num2str(i)]);
                    counter = sum(PPsperDend(1:i))+1;
                end
            end
            if twochannels == 1
                axes(axes2)
                for i = 1:length(glovar.RedPolyLinePos)
                    glovar.RedPolyLinePos{i} = [coordinates{i}(1)-radius, coordinates{i}(2)-radius, radius*2, radius*2];
                    glovar.RedPolyROI{i} = rectangle('Position', glovar.RedPolyLinePos{i}, 'EdgeColor', 'cyan', 'Tag', ['Dendrite ', num2str(DendriteNum), ' RedPolyROI', num2str(i)], 'Curvature', [1 1], 'ButtonDownFcn', 'Drag_Poly');
                    x_red = [x_red,coordinates{i}(1)];
                    y_red = [y_red,coordinates{i}(2)];
                end
                for i = 1:2:length(glovar.RedPolyLinePos)
                    hold on;
                    plot(x_red,y_red, 'color', 'cyan', 'Tag', 'PolyLine');
                end
            else
            end
        end
    end
catch
    disp('Dendrite ROIs could not be extracted, and need to be redrawn')
end
    
% PPseries{1} = 1:PPnum(1);
% for i = 2:length(PPnum)
%     PPseries{i} = PPnum(i-1):PPnum(i);
% end
% if length(glovar.PolyROI)<PPnum(end)
%     disp('Could not recover all dendrite information')
%     for i = 1:DendNum
%         try
%             cellfun(@(x) isempty(x), PolyROIcheck, 'UniformOutput', false)
%         catch
%         end
%     end
% end


 %%% Overwrite the previous existing global workspace with the newly imprinted one
 
if ~isempty(regexp(running, 'CaImageViewer'))
    gui_CaImageViewer = glovar;
elseif ~isempty(regexp(running, 'FluorescenceSuite'));
    gui_FluorescenceSuite = glovar;
end

function BaselineFrames_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to BaselineFrames_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BaselineFrames_EditableText as text
%        str2double(get(hObject,'String')) returns contents of BaselineFrames_EditableText as a double


% --- Executes during object creation, after setting all properties.
function BaselineFrames_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BaselineFrames_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function LoadFile_Callback(hObject, eventdata, handles)
% hObject    handle to LoadFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global gui_CaImageViewer

%%% Get file information %%%

try
    cd('Z:\People\Nathan\Data')
catch
    cd('/usr/local/lab/People/Nathan/Data')
end

%%% Initialize/reset parameters and settings when loading new file
set(gui_CaImageViewer.figure.handles.MaxProjection_CheckBox, 'Value', 0);
set(gui_CaImageViewer.figure.handles.AveProjection_CheckBox, 'Value', 0);
set(gui_CaImageViewer.figure.handles.ImageSlider_Slider, 'Enable', 'on');
gui_CaImageViewer.NewSpineAnalysis = 0;


[filename, pathname] = uigetfile('.tif');

if isnumeric(pathname) && isnumeric(filename)
    return
end

fname = [pathname, filename];

load_type = listdlg('PromptString', 'Analyze compressed image (fast), or full series (slow)? (note that only compressed image will be displayed)', 'SelectionMode', 'single','ListString', {'Compressed', 'Full'});

timecourse_image_number = [];

if load_type == 1
    gui_CaImageViewer.Load_Type = 'Compressed';
elseif load_type == 2
    gui_CaImageViewer.Load_Type = 'Full';
end

% fname = fname;
CaImage_File_info = imfinfo(fname);
timecourse_image_number = numel(CaImage_File_info);


gui_CaImageViewer.filename = filename;
gui_CaImageViewer.save_directory = pathname;
cd(pathname)
twochannels = get(gui_CaImageViewer.figure.handles.TwoChannels_CheckBox, 'Value');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Set Image Properties %%%

Green_Frame = 1;
Red_Frame = 1;

gui_CaImageViewer.GCaMP_Image = [];
gui_CaImageViewer.Red_Image = [];

h = waitbar(0, 'Loading Image ');
TifLink = Tiff(fname, 'r');

if twochannels
    if load_type == 1
        for i = 1:2:timecourse_image_number
            TifLink.setDirectory(i);
            gui_CaImageViewer.GCaMP_Image{1,Green_Frame} = TifLink.read();
            Green_Frame = Green_Frame+1;
            waitbar(Green_Frame/(timecourse_image_number/2),h,['Loading Image ', num2str(Green_Frame)]);
        end
        for i = 2:2:timecourse_image_number
            [gui_CaImageViewer.Red_Image{1,Red_Frame}, mapr{1,Red_Frame}] = imread(fname, i, 'Info', CaImage_File_info);
            Red_Frame = Red_Frame+1;
        end
    elseif load_type == 2
        for i = 1:2:timecourse_image_number
            TifLink.setDirectory(i);
            gui_CaImageViewer.GCaMP_Image{1,Green_Frame} = TifLink.read();
            Green_Frame = Green_Frame+1;
        end
        for i = 2:2:timecourse_image_number
            [gui_CaImageViewer.Red_Image{1,Red_Frame}, mapr{1,Red_Frame}] = imread(fname, i, 'Info', CaImage_File_info);
            Red_Frame = Red_Frame+1;
        end
    end
else
    if load_type == 1 %%% if loading compressed 
        for i = 1:timecourse_image_number
            TifLink.setDirectory(i);
            gui_CaImageViewer.GCaMP_Image{1,Green_Frame} = TifLink.read();
            Green_Frame = Green_Frame+1;
            waitbar(Green_Frame/timecourse_image_number,h,['Loading Image ', num2str(Green_Frame)]);
        end
    elseif load_type == 2 %%% if loading full
        for i = 1:timecourse_image_number
            TifLink.setDirectory(i);
            gui_CaImageViewer.GCaMP_Image{1,Green_Frame} = TifLink.read();
            Green_Frame = Green_Frame+1;
            waitbar(Green_Frame/timecourse_image_number,h,['Loading Image ', num2str(Green_Frame)]);
        end
    end
        if ~gui_CaImageViewer.LoadedFile
            Green_loc = gui_CaImageViewer.GreenGraph_loc;
            Red_loc = gui_CaImageViewer.RedGraph_loc;
            set(handles.RedGraph, 'Visible', 'off')
            set(handles.Channel2_StaticText, 'Visible', 'off')
            set(handles.RedUpperLUT_EditableText, 'Visible', 'off')
            set(handles.RedLowerLUT_EditableText, 'Visible', 'off')
            set(handles.RedGamma_EditableText, 'Visible', 'off')
            set(handles.RedGamma_StaticText, 'Visible', 'off')
            gui_CaImageViewer.GraphPlacement = [Green_loc(1), Green_loc(2), Green_loc(3)+(Red_loc(1)-(Green_loc(1)+Green_loc(3))+Red_loc(3)), Green_loc(4)];
            set(handles.GreenGraph, 'Units', 'normalized')
            figure(gui_CaImageViewer.figure.handles.figure1)
            axes(gui_CaImageViewer.figure.handles.GreenGraph);
%             set(handles.GreenGraph, 'Position', [Green_loc(1), Green_loc(2), Green_loc(3)+(Red_loc(1)-(Green_loc(1)+Green_loc(3))+Red_loc(3)), Green_loc(4)])
            intergraphdistance = Red_loc(1)-(Green_loc(1)+Green_loc(3));
            set(handles.GreenGraph, 'Position', [Green_loc(1), Green_loc(2), Green_loc(3)+Red_loc(3)+intergraphdistance, Green_loc(4)])
        else
        end
end

close(h)

channel1 = gui_CaImageViewer.GCaMP_Image;
channel2 = gui_CaImageViewer.Red_Image;

CommandSource = 'Loader';

[ch1image, ch2image] = PlaceImages(channel1, channel2, CommandSource);

imageserieslength = size(gui_CaImageViewer.GCaMP_Image, 2);
gui_CaImageViewer.imageserieslength = imageserieslength;

set(handles.ImageSlider_Slider, 'Value', 1);
set(handles.ImageSlider_Slider, 'Min', 1);
set(handles.ImageSlider_Slider, 'Max', imageserieslength);
set(handles.ImageSlider_Slider, 'SliderStep', [(1/(imageserieslength-1)) (32/(imageserieslength-1))]);  %%% The Slider Step values indicate the minor and major transitions, which should be represented by the desired transition as the numerator and the length of the series as the denominator
set(handles.Frame_EditableText, 'String', 1);
set(handles.SmoothingFactor_EditableText, 'String', '1');

% 
% Smoothing = str2num(get(handles.SmoothingFactor_EditableText, 'String'));
% 
% if Smoothing ~= 1
%     Smoother(hObject, eventdata, ch1image,ch2image, CommandSource)
% end

set(gui_CaImageViewer.figure.handles.output, 'WindowButtonDownFcn', [])

gui_CaImageViewer.LoadedFile = 1;


% --- Executes on button press in MaxProjection_CheckBox.
function MaxProjection_CheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to MaxProjection_CheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MaxProjection_CheckBox

global gui_CaImageViewer

val = get(handles.MaxProjection_CheckBox, 'Value');
frame = get(handles.Frame_EditableText, 'String');

ImageNum = str2num(get(gui_CaImageViewer.figure.handles.Frame_EditableText, 'String'));
twochannels = get(gui_CaImageViewer.figure.handles.TwoChannels_CheckBox, 'Value');
filterwindow = str2num(get(gui_CaImageViewer.figure.handles.SmoothingFactor_EditableText, 'String'));

if val
    set(handles.AveProjection_CheckBox, 'Value', 0);
    im = gui_CaImageViewer.GCaMP_Image;
    im = cat(3, im{:});
    immax = max(im, [], 3);

    filterwindow = str2num(get(gui_CaImageViewer.figure.handles.SmoothingFactor_EditableText, 'String'));
    
    if filterwindow == 1;
    
        channel1 = immax;
        if twochannels == 1
            channel2 = gui_CaImageViewer.Red_Image{ImageNum};
        else
            channel2 = [];
        end

        CommandSource = 'Slider';

        %%%%%%%%%
        PlaceImages(channel1,channel2, CommandSource);
        %%%%%%%%%
    
    else
        smoothing_green = filter2(ones(filterwindow, filterwindow)/filterwindow^2, immax);
        channel1 = smoothing_green;
        if twochannels == 1
            smoothing_red = filter2(ones(filterwindow, filterwindow)/filterwindow^2, gui_CaImageViewer.Red_Image{ImageNum});
            channel2 = smoothing_red;
        else
            channel2 = [];
        end

        CommandSource = 'Slider';

        %%%%%%%%%
        PlaceImages(channel1,channel2, CommandSource);
        %%%%%%%%%
    end
else
    channel1 = gui_CaImageViewer.GCaMP_Image{ImageNum};
    if twochannels == 1
        channel2 = gui_CaImageViewer.Red_Image{ImageNum};
    else
        channel2 = [];
    end
    
    PlaceImages(channel1, channel2, 'Slider');
    
    CaImageSlider(ImageNum);
end
    

% --- Executes on button press in SaveROIs_PushButton.
function SaveROIs_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveROIs_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and  data (see GUIDATA)

global gui_CaImageViewer

drawer = get(gui_CaImageViewer.figure.handles.figure1, 'UserData');

a.SpineROIs = gui_CaImageViewer.ROI;
a.SpineROItext = gui_CaImageViewer.ROItext;
a.PolyROI = gui_CaImageViewer.PolyROI;
a.PolyLines = gui_CaImageViewer.PolyLine;
a.PolyLinePosition = gui_CaImageViewer.PolyLinePos;
a.PolyROIPos = gui_CaImageViewer.PolyLinePos;
a.PolyLineVertices = gui_CaImageViewer.PolyLineVertices;
a.NumberofSpines = gui_CaImageViewer.Spine_Number;

if a.SpineROIs(1) == 0
    msgbox('Cannot save ROIs without drawing background!');
    return
end

for i = 1:length(a.SpineROIs)
    a.ROIPosition{i} = get(a.SpineROIs(i), 'Position');
end

a.NumberofDendrites = gui_CaImageViewer.Dendrite_Number;
a.DendritePolyPointNumber = gui_CaImageViewer.DendritePolyPointNumber;

DendNum = gui_CaImageViewer.Dendrite_Number;

if DendNum > 1
    prompt = cell(1,DendNum);
    for i = 1:length(prompt)
        prompt{i} = ['Spines on dendrite ', num2str(i)];
    end
    name = 'Spine Grouping';
    numlines = 1;
    defaultanswer = cell(1,DendNum);
    for i = 1:length(defaultanswer)
        if ~isempty(gui_CaImageViewer.SpineDendriteGrouping)
            try
                first = min(gui_CaImageViewer.SpineDendriteGrouping{i});
                last = max(gui_CaImageViewer.SpineDendriteGrouping{i});
                defaultanswer{i} = sprintf('%d:%d', first, last);
            catch
                defaultanswer{i} = '';
            end
        else
            defaultanswer{i} = '';
        end
    end
    s_d_grouping = inputdlg(prompt, name, numlines, defaultanswer);
    for i = 1:DendNum
        DendSpines{i} = str2num(s_d_grouping{i});
    end
else 
    DendSpines{1} = 1:length(a.SpineROIs)-1;
end

a.SpineDendriteGrouping = DendSpines;

gui_CaImageViewer.SpineDendriteGrouping = DendSpines;

if gui_CaImageViewer.NewSpineAnalysis
    animal = regexp(gui_CaImageViewer.filename, '[A-Z]{2,3}[0-9]*', 'match');
    animal = animal{1};
    date = regexp(gui_CaImageViewer.save_directory, '[0-9]{5,7}', 'match');
    experiment = [animal, '_', date{1}];
    fname = [experiment, '_NewSpineAnalysisROIs', '_DrawnBy', drawer];
    
        %%%%% Move to parent folder
        fullpath = gui_CaImageViewer.save_directory;
        allseps = strfind(fullpath, '\');
        stepsup = 2;
        newpath = fullpath(1:allseps(end-stepsup)-1); %%% move two steps up in the path directory to get bath to the main animal folder (e.g. Z:/People/Nathan/Data/NH004 instead of Z:/People/Nathan/Data/NH004/160316/summed)
        cd(newpath)
        
        %%%%% Identify imaging field number
        currentimagingfield = gui_CaImageViewer.NewSpineAnalysisInfo.CurrentImagingField;
        
        %%% Identify instance of appearance of this imaging field (i.e. is this the first time imaging here? the second? etc.) 
        prompt = 'What imaging instance (of this field) is this?';
        name = 'Designate imaging instance';
        numlines = 1;
        defaultanswer = {'1'};
        
        currentsession = inputdlg(prompt, name, numlines, defaultanswer);
        currentsession = str2num(currentsession{1});
        
        try
            load(['Imaging Field ', num2str(currentimagingfield), ' Spine Registry'])
        catch
        end
        
     %%   if currentsession == 1     
      %%%      SpineRegistry.Data(1:length(a.SpineROIs)-1,currentsession) = gui_CaImageViewer.NewSpineAnalysisInfo.SpineList;

      %%%      a.SpineStatusList = gui_CaImageViewer.NewSpineAnalysisInfo.SpineList;

      %%%      save(['Imaging Field ', num2str(currentimagingfield), ' Spine Registry'], 'SpineRegistry')
      %%%  else ZL comment: this part causing more issues in saving ROIs for session 1, use with caution
            SpineRegistry.Data(1:length(gui_CaImageViewer.NewSpineAnalysisInfo.SpineList),currentsession) = gui_CaImageViewer.NewSpineAnalysisInfo.SpineList;

            a.SpineStatusList = gui_CaImageViewer.NewSpineAnalysisInfo.SpineList;

            save(['Imaging Field ', num2str(currentimagingfield), ' Spine Registry'], 'SpineRegistry');
      %%%  end
else
    experiment = regexp(gui_CaImageViewer.filename, '[A-Z]{2}\d+[_]\d+', 'match');

    fname = [experiment{1}, '_SavedROIs', '_DrawnBy', drawer];
end

eval([fname,'= a'])

target_dir = gui_CaImageViewer.save_directory;
cd(target_dir);

save(fname, fname)


% --------------------------------------------------------------------1:14
function User_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to User_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ChangeUser_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to ChangeUser_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


id = inputdlg('Enter new user name:', 'New User', 1);

set(handles.figure1, 'UserData', id{1})


% --------------------------------------------------------------------
function ImageOptions_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to ImageOptions_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_CaImageViewer


% --------------------------------------------------------------------
function ColorMap_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to ColorMap_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_CaImageViewer


% --------------------------------------------------------------------
function RGB_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to RGB_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_CaImageViewer

gui_CaImageViewer.CurrentCMap = 'RGB';

if size(gui_CaImageViewer.ch1image,3)>1
    ch1image = gui_CaImageViewer.ch1image(:,:,2);
else
    ch1image = gui_CaImageViewer.ch1image;
end

PlaceImages(ch1image,[],'Slider');

% --------------------------------------------------------------------
function Jet_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to Jet_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_CaImageViewer

gui_CaImageViewer.CurrentCMap = 'Jet';

if size(gui_CaImageViewer.ch1image,3)>1
    ch1image = gui_CaImageViewer.ch1image(:,:,2);
else
    ch1image = gui_CaImageViewer.ch1image;
end

PlaceImages(ch1image,[],'Slider');

% --------------------------------------------------------------------
function Hot_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to Hot_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_CaImageViewer

gui_CaImageViewer.CurrentCMap = 'Hot';

if size(gui_CaImageViewer.ch1image,3)>1
    ch1image = gui_CaImageViewer.ch1image(:,:,2);
else
    ch1image = gui_CaImageViewer.ch1image;
end

PlaceImages(ch1image,[],'Slider');

% --------------------------------------------------------------------
function Fire_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to Fire_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_CaImageViewer

gui_CaImageViewer.CurrentCMap = 'Fire';

if size(gui_CaImageViewer.ch1image,3)>1
    ch1image = gui_CaImageViewer.ch1image(:,:,2);
else
    ch1image = gui_CaImageViewer.ch1image;
end

PlaceImages(ch1image,[],'Slider');



% --------------------------------------------------------------------
function GraphScaling_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to GraphScaling_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Stretched_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to Stretched_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_CaImageViewer

PlaceImages(gui_CaImageViewer.ch1image,[],'Stretcher');

% --------------------------------------------------------------------
function Square_Dropdown_Callback(hObject, eventdata, handles)
% hObject    handle to Square_Dropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_CaImageViewer

PlaceImages(gui_CaImageViewer.ch1image,[],'Square');



function ROIoffset_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to ROIoffset_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ROIoffset_EditableText as text
%        str2double(get(hObject,'String')) returns contents of ROIoffset_EditableText as a double


% --- Executes during object creation, after setting all properties.
function ROIoffset_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROIoffset_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AveProjection_CheckBox.
function AveProjection_CheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to AveProjection_CheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AveProjection_CheckBox

global gui_CaImageViewer

val = get(handles.AveProjection_CheckBox, 'Value');
frame = get(handles.Frame_EditableText, 'String');

ImageNum = str2num(get(gui_CaImageViewer.figure.handles.Frame_EditableText, 'String'));
twochannels = get(gui_CaImageViewer.figure.handles.TwoChannels_CheckBox, 'Value');
filterwindow = str2num(get(gui_CaImageViewer.figure.handles.SmoothingFactor_EditableText, 'String'));

if val
    set(handles.MaxProjection_CheckBox, 'Value', 0);
    im = gui_CaImageViewer.GCaMP_Image;
    im = cat(3, im{:});
    immean = mean(im,3);

    filterwindow = str2num(get(gui_CaImageViewer.figure.handles.SmoothingFactor_EditableText, 'String'));
    
    if filterwindow == 1;
    
        channel1 = immean;
        if twochannels == 1
            channel2 = gui_CaImageViewer.Red_Image{ImageNum};
        else
            channel2 = [];
        end

        CommandSource = 'Slider';

        %%%%%%%%%
        PlaceImages(channel1,channel2, CommandSource);
        %%%%%%%%%
    
    else
        smoothing_green = filter2(ones(filterwindow, filterwindow)/filterwindow^2, immean);
        channel1 = smoothing_green;
        if twochannels == 1
            smoothing_red = filter2(ones(filterwindow, filterwindow)/filterwindow^2, gui_CaImageViewer.Red_Image{ImageNum});
            channel2 = smoothing_red;
        else
            channel2 = [];
        end

        CommandSource = 'Slider';

        %%%%%%%%%
        PlaceImages(channel1,channel2, CommandSource);
        %%%%%%%%%
    end
else
    channel1 = gui_CaImageViewer.GCaMP_Image{ImageNum};
    if twochannels == 1
        channel2 = gui_CaImageViewer.Red_Image{ImageNum};
    else
        channel2 = [];
    end
    
    PlaceImages(channel1, channel2, 'Slider');
    
    CaImageSlider(ImageNum);
end


% --------------------------------------------------------------------
function MultipleSessions_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to MultipleSessions_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MultiSession_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to MultiSession_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_CaImageViewer

animal = regexp(gui_CaImageViewer.save_directory, '[A-Z]{2,3}0*[0-9]*', 'match');
fileparts = regexp(gui_CaImageViewer.save_directory, '[A-Z]{2,3}0*[0-9]*', 'split');
directory = [fileparts{1}, animal{1}];

exp_folder = dir(directory);

numsessions = length(exp_folder)-2;


scrsz = get(0, 'ScreenSize');
OverSessionsFigure = figure('Position', scrsz, 'Name', 'Multiple Sessions Analysis', 'NumberTitle', 'off');
set(OverSessionsFigure, 'UserData', zeros(1,14));
h1 = waitbar(0, 'Loading images for session 1');

for i = 3:length(exp_folder)
    cd(directory)
    if isdir(exp_folder(i).name) && isempty(regexp(exp_folder(i).name, '[A-Z]*'))
        path = [directory, '\', exp_folder(i).name, '\summed'];
        cd(path);
        a = dir(cd);
        for j = 3:length(a)
            if ~isempty(strfind(a(j).name, 'summed_50.tif'))
                imagefile = a(j).name;
                break
            end
        end
        fname = [path,'\',imagefile];
        CaImage_File_info = imfinfo(fname);
        timecourse_image_number = numel(CaImage_File_info);
        TifLink = Tiff(fname, 'r');
        h2 = waitbar(0, 'Loading image ');
        Green_Frame = 1;
        GCaMP_Image = {};
            for j = 1:timecourse_image_number
                TifLink.setDirectory(j);
                GCaMP_Image{1,Green_Frame} = TifLink.read();
                Green_Frame = Green_Frame+1;
                waitbar(Green_Frame/timecourse_image_number,h2,['Loading Image ', num2str(Green_Frame)]);
            end
            delete(h2)
        im = cat(3, GCaMP_Image{:});
        immean = mean(im,3);
        figure(OverSessionsFigure);
    %     subplot(2,round(numsessions/2), i-2)
        figpos = get(gcf, 'Position');
        xint = figpos(3)/7;
        yint = figpos(4)/2;
        if (i-2)/7 <= 1
            yrow = 1;
        else
            yrow = 0;
        end
        posmat = [1:7,1:7];
        axes('Position', [((posmat(i-2)-1)*xint+10)/figpos(3), (yint*(yrow)+150)/figpos(4), 250/figpos(3), 250/figpos(4)])
        A = imagesc(immean); colormap(fire); set(gca, 'XTick', [], 'YTick', [])
        set(A,'ButtonDownFcn', @HighLightAxis)
        set(A, 'UserData', (i-2));
        title(exp_folder(i).name)
        waitbar((i-2)/numsessions, h1, ['Loading images for session ', num2str(i-2)])
    else
    end
end

%%% Run these variables before manually adding a button!
%        figpos = get(gcf, 'Position');hui
%        xint = figpos(3)/7;
%        yint = figpos(4)/2;


uicontrol('Style', 'pushbutton', 'String', 'Project to Analysis Window', 'FontSize', 12, 'Units', 'Normalized','Position', [0.4 0.925 0.2 0.05], 'CallBack', @ProjectToAnalysisWindow)
uicontrol('Style', 'pushbutton', 'String', 'Group Selected Imaging Fields', 'Units', 'Normalized', 'Position', [((5)*xint+10)/figpos(3) 0.925 250/figpos(3) 0.05], 'CallBack', @GroupImagingWindows) %%% Position over the 6th image window 
uicontrol('Style', 'text', 'String', 'Image Field Count:', 'Units', 'Normalized', 'Position', [((6)*xint+10)/figpos(3), 0.925-0.02, 125/figpos(3), 0.05], 'BackgroundColor', get(gcf, 'Color'))
uicontrol('Style', 'edit', 'String', '1', 'Units', 'Normalized', 'Position', [((6)*xint+10)/figpos(3)+125/figpos(3), 0.925, 125/figpos(3), 0.05])
uicontrol('Style', 'pushbutton', 'String', 'Tabulate spine lifetimes of selected field', 'Units', 'Normalized', 'Position', [(10)/figpos(3), 0.925, 250/figpos(3), 0.05], 'Callback', @TabulateSpineLifetimes)
uicontrol('Style', 'pushbutton', 'String', 'Deselect All Axes', 'Fontsize', 12, 'Units', 'Normalized', 'Position', [((6)*xint+10)/figpos(3), 0.05, 250/figpos(3), 0.05], 'Callback', @DeselectAxes)
uicontrol('Style', 'pushbutton', 'String', 'Compare Image Pair', 'Fontsize', 12, 'Units', 'Normalized', 'Position', [(10)/figpos(3), 0.05, 250/figpos(3), 0.05], 'Callback', @CompareImagePair)
uicontrol('Style', 'checkbox', 'String', 'Try image alignment?', 'Tag', 'Alignment_CheckBox', 'Units', 'Normalized', 'Position', [((2-1)*xint+10)/figpos(3), 0.05, 250/figpos(3), 0.05], 'BackgroundColor', get(gcf, 'Color'))

delete(h1)


% --------------------------------------------------------------------
function ImageIntegrity_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to ImageIntegrity_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function CheckMotionCorrection_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to CheckMotionCorrection_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_CaImageViewer

for i = 2:length(gui_CaImageViewer.GCaMP_Image)
    pwc(1,i) = corr2(gui_CaImageViewer.GCaMP_Image{i},gui_CaImageViewer.GCaMP_Image{i-1});
end

figure('NumberTitle', 'off', 'Name', 'Pairwise Image Correlation Over Timecourse');
plot(pwc, 'k')
xlabel('Frame', 'FontSize', 14)
ylabel('Pairwise Correlation', 'FontSize', 14);