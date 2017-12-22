function varargout = CaAnalysis(varargin)
% CAANALYSIS MATLAB code for CaAnalysis.fig
%      CAANALYSIS, by itself, creates a new CAANALYSIS or raises the existing
%      singleton*.
%
%      H = CAANALYSIS returns the handle to a new CAANALYSIS or the handle to
%      the existing singleton*.
%
%      CAANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CAANALYSIS.M with the given input arguments.
%
%      CAANALYSIS('Property','Value',...) creates a new CAANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CaAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CaAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CaAnalysis

% Last Modified by GUIDE v2.5 06-Sep-2015 17:09:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CaAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @CaAnalysis_OutputFcn, ...
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


% --- Executes just before CaAnalysis is made visible.
function CaAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CaAnalysis (see VARARGIN)

% Choose default command line output for CaAnalysis

handles.output = hObject;

global gui_CaAnalysis

gui_CaAnalysis.figure.handles = handles;

% Update handles structure
guidata(hObject, handles);

set(handles.Graph, 'YTick', []);
set(handles.Graph, 'XTick', []);

% UIWAIT makes CaAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Executes during object creation, after setting all properties.
function Graph_CreateFcn(~, ~, ~)
% hObject    handle to Graph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate Graph


% --- Outputs from this function are returned to the command line.
function varargout = CaAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function EnterCondition_EditableText_Callback(~, ~, handles)
% hObject    handle to EnterCondition_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EnterCondition_EditableText as text
%        str2double(get(hObject,'String')) returns contents of EnterCondition_EditableText as a double


% --- Executes during object creation, after setting all properties.
function EnterCondition_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EnterCondition_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Excluder_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to Excluder_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Excluder_EditableText as text
%        str2double(get(hObject,'String')) returns contents of Excluder_EditableText as a double


% --- Executes during object creation, after setting all properties.
function Excluder_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Excluder_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Calcium_PushButton.
function Calcium_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Calcium_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gui_CaAnalysis

h1 = get(handles.File_EditableText, 'String');

Subject_File = h1;

h2 = get(handles.EnterCondition_EditableText, 'String');

Condition = h2;

h3 = get(handles.Excluder_EditableText, 'String');

Excluder = h3;

Excluder = regexp(Excluder, '\,* ', 'split');

Red = get(handles.DivideByRed_Checkbox, 'Value');

CaDataAverageNew(Subject_File, Condition, Excluder, Red);



function File_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to File_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of File_EditableText as text
%        str2double(get(hObject,'String')) returns contents of File_EditableText as a double


% --- Executes during object creation, after setting all properties.
function File_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FileSelect_PushButton.
function FileSelect_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to FileSelect_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filenameCa, pathnameCa, filterindexCa] = uigetfile();

File_for_Ca = [pathnameCa, filenameCa];

set(handles.File_EditableText, 'String', File_for_Ca)


% --- Executes on button press in GraphClear_PushButton.
function GraphClear_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to GraphClear_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gui_CaAnalysis

cla(gui_CaAnalysis.figure.handles.Graph)
leg = findobj(legend);
delete(leg);
if isfield(gui_CaAnalysis, 'legend')
    gui_CaAnalysis = rmfield(gui_CaAnalysis, 'legend');
end
if isfield(gui_CaAnalysis, 'Curve')
    gui_CaAnalysis = rmfield(gui_CaAnalysis, 'Curve');
end
if isfield(gui_CaAnalysis, 'AutoFig')
    gui_CaAnalysis = rmfield(gui_CaAnalysis, 'AutoFig');
end

if isfield(gui_CaAnalysis, 'File_n')
    gui_CaAnalysis = rmfield(gui_CaAnalysis, 'File_n');
end

if isfield(gui_CaAnalysis, 'vol_axes')
    cla(gui_CaAnalysis.vol_axes)
end
    



function xlim_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to xlim_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xlim_EditableText as text
%        str2double(get(hObject,'String')) returns contents of xlim_EditableText as a double


% --- Executes during object creation, after setting all properties.
function xlim_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xlim_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ylim_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to ylim_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ylim_EditableText as text
%        str2double(get(hObject,'String')) returns contents of ylim_EditableText as a double


% --- Executes during object creation, after setting all properties.
function ylim_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ylim_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AdjustAxes_PushButton.
function AdjustAxes_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to AdjustAxes_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gui_CaAnalysis

x_limits = get(gui_CaAnalysis.figure.handles.xlim_EditableText, 'String');
y_limits = get(gui_CaAnalysis.figure.handles.ylim_EditableText, 'String');

set(gui_CaAnalysis.figure.handles.Graph, 'XLim', str2num(x_limits));

if ~strcmpi(y_limits, '--');
    set(gui_CaAnalysis.figure.handles.Graph, 'YLim', str2num(y_limits));
end


% --- Executes on button press in MakeFigure_PushButton.
function MakeFigure_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to MakeFigure_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gui_CaAnalysis

h1 = get(handles.File_EditableText, 'String');

Subject_File = h1;

h2 = get(handles.EnterCondition_EditableText, 'String');

Condition = h2;

h3 = get(handles.Excluder_EditableText, 'String');

Excluder = h3;

Excluder = regexp(Excluder, '\,* ', 'split');

MainFigure = gui_CaAnalysis.figure.handles.Graph;
InsetFigure = gui_CaAnalysis.vol_axes

scrsz = get(0, 'ScreenSize');

h4 = figure('Position', [1, scrsz(4)/3, scrsz(3)/2, scrsz(4)/2]); 

h5 = copyobj(MainFigure, h4);

legend(h5, gui_CaAnalysis.Curve, gui_CaAnalysis.legend, 'Location', 'NorthWest')


copyobj(InsetFigure, h4);


% --- Executes on button press in UTA_PushButton.
function UTA_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to UTA_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Make_UncagingTriggeredAve


% --------------------------------------------------------------------
function Analyze_Callback(hObject, eventdata, handles)
% hObject    handle to Analyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Spatial_Profile_Callback(hObject, eventdata, handles)
% hObject    handle to Spatial_Profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global gui_CaAnalysis

h1 = get(handles.File_EditableText, 'String');

Subject_File = h1;

h2 = get(handles.EnterCondition_EditableText, 'String');

Condition = h2;

h3 = get(handles.Excluder_EditableText, 'String');

Excluder = h3;

Excluder = regexp(Excluder, '\,* ', 'split');

FastSpatialProfile(Subject_File, Condition, Excluder);


% --- Executes on button press in AutoFigure_CheckBox.
function AutoFigure_CheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to AutoFigure_CheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AutoFigure_CheckBox


% --------------------------------------------------------------------
function BarGraph_Callback(hObject, eventdata, handles)
% hObject    handle to BarGraph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CaAnalysisBarGraph


% --------------------------------------------------------------------
function CodeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to CodeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

edit CaAnalysis


% --- Executes on button press in DivideByRed_Checkbox.
function DivideByRed_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to DivideByRed_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DivideByRed_Checkbox
