function varargout = KomiyamaLabHub(varargin)
% KOMIYAMALABHUB MATLAB code for KomiyamaLabHub.fig
%      KOMIYAMALABHUB, by itself, creates a new KOMIYAMALABHUB or raises the existing
%      singleton*.
%
%      H = KOMIYAMALABHUB returns the handle to a new KOMIYAMALABHUB or the handle to
%      the existing singleton*.
%
%      KOMIYAMALABHUB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KOMIYAMALABHUB.M with the given input arguments.
%
%      KOMIYAMALABHUB('Property','Value',...) creates a new KOMIYAMALABHUB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before KomiyamaLabHub_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to KomiyamaLabHub_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help KomiyamaLabHub

% Last Modified by GUIDE v2.5 16-Dec-2017 07:48:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @KomiyamaLabHub_OpeningFcn, ...
                   'gui_OutputFcn',  @KomiyamaLabHub_OutputFcn, ...
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


% --- Executes just before KomiyamaLabHub is made visible.
function KomiyamaLabHub_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to KomiyamaLabHub (see VARARGIN)

% Choose default command line output for KomiyamaLabHub
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes KomiyamaLabHub wait for user response (see UIRESUME)
% uiwait(handles.figure1);

set(handles.AnalyzeActivity_ToggleButton, 'Value', 0);
set(handles.AnalyzeBehavior_ToggleButton, 'Value', 0);
set(handles.Timecourse_PushButton, 'Enable', 'off');
set(handles.Clustering_PushButton, 'Enable', 'off');
set(handles.BehaviorTimecourse_PushButton, 'Enable', 'off');
set(handles.AlignActivty_PushButton, 'Enable', 'off');

% setappdata(gcf, 'Folder', 'C:\Users\Komiyama\Desktop\ActivitySummary_UsingRawData');

Scrsz = get(0, 'Screensize');
    d = dialog('Position', [(Scrsz(3)/2)-125 Scrsz(4)/2-75 250 150], 'Name', 'User');
    txt = uicontrol('Parent', d, 'Style', 'text', 'Position', [10 100 230 30], 'String', 'Select User:');
    btn1 = uicontrol('Parent', d, 'Style', 'pushbutton', 'Position', [35 30 50 25], 'String', 'Nathan', 'Callback', @UserName);
    btn2 = uicontrol('Parent', d, 'Style', 'pushbutton', 'Position', [85.5 30 70 25], 'String', 'Zhongmin', 'Callback', @UserName);
    btn3 = uicontrol('Parent', d, 'Style', 'pushbutton', 'Position', [156 30 50 25], 'String', 'Sara', 'Callback', @UserName);
    uiwait(d)
    choice = get(d, 'UserData');
    set(handles.figure1, 'UserData', choice);
    delete(d);
    
global gui_KomiyamaLabHub
gui_KomiyamaLabHub.figure.handles = handles;

    
function UserName(hObject, eventdata, ~)

button = get(hObject);

choice = button.String;

sourcewindow = button.Parent;

set(sourcewindow, 'UserData', choice);

uiresume



% --- Outputs from this function are returned to the command line.
function varargout = KomiyamaLabHub_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in AnimalName_ListBox.
function AnimalName_ListBox_Callback(hObject, eventdata, handles)
% hObject    handle to AnimalName_ListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AnimalName_ListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AnimalName_ListBox


% --- Executes during object creation, after setting all properties.
function AnimalName_ListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AnimalName_ListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AnalyzeActivity_ToggleButton.
function AnalyzeActivity_ToggleButton_Callback(hObject, eventdata, handles)
% hObject    handle to AnalyzeActivity_ToggleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AnalyzeActivity_ToggleButton

Activity = get(handles.AnalyzeActivity_ToggleButton, 'Value');
Behavior = get(handles.AnalyzeBehavior_ToggleButton, 'Value');

if Activity == 1 && Behavior == 0
    set(handles.Timecourse_PushButton, 'Enable', 'on');
    set(handles.Clustering_PushButton, 'Enable', 'on');
    set(handles.BehaviorTimecourse_PushButton, 'Enable', 'off');
    set(handles.AlignActivty_PushButton, 'Enable', 'off');
elseif Activity == 0 && Behavior == 1
    set(handles.Timecourse_PushButton, 'Enable', 'off');
    set(handles.Clustering_PushButton, 'Enable', 'off');
    set(handles.BehaviorTimecourse_PushButton, 'Enable', 'on');
    set(handles.AlignActivty_PushButton, 'Enable', 'off');
elseif Activity == 1 && Behavior == 1
    set(handles.Timecourse_PushButton, 'Enable', 'off');
    set(handles.Clustering_PushButton, 'Enable', 'on');
    set(handles.BehaviorTimecourse_PushButton, 'Enable', 'off');
    set(handles.AlignActivty_PushButton, 'Enable', 'on');
elseif Activity == 0 && Behavior == 0
    set(handles.Timecourse_PushButton, 'Enable', 'off');
    set(handles.Clustering_PushButton, 'Enable', 'off');
    set(handles.BehaviorTimecourse_PushButton, 'Enable', 'off');
    set(handles.AlignActivty_PushButton, 'Enable', 'off');
end


% --- Executes on button press in AnalyzeBehavior_ToggleButton.
function AnalyzeBehavior_ToggleButton_Callback(hObject, eventdata, handles)
% hObject    handle to AnalyzeBehavior_ToggleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AnalyzeBehavior_ToggleButton

Activity = get(handles.AnalyzeActivity_ToggleButton, 'Value');
Behavior = get(handles.AnalyzeBehavior_ToggleButton, 'Value');

if Activity == 1 && Behavior == 0
    set(handles.Timecourse_PushButton, 'Enable', 'on');
    set(handles.Clustering_PushButton, 'Enable', 'on');
    set(handles.BehaviorTimecourse_PushButton, 'Enable', 'off');
    set(handles.AlignActivty_PushButton, 'Enable', 'off');
elseif Activity == 0 && Behavior == 1
    set(handles.Timecourse_PushButton, 'Enable', 'off');
    set(handles.Clustering_PushButton, 'Enable', 'off');
    set(handles.BehaviorTimecourse_PushButton, 'Enable', 'on');   
    set(handles.AlignActivty_PushButton, 'Enable', 'off');
elseif Activity == 1 && Behavior == 1
    set(handles.Timecourse_PushButton, 'Enable', 'off');
    set(handles.Clustering_PushButton, 'Enable', 'on');
    set(handles.BehaviorTimecourse_PushButton, 'Enable', 'off');
    set(handles.AlignActivty_PushButton, 'Enable', 'on');
elseif Activity == 0 && Behavior == 0
    set(handles.Timecourse_PushButton, 'Enable', 'off');
    set(handles.Clustering_PushButton, 'Enable', 'off');
    set(handles.BehaviorTimecourse_PushButton, 'Enable', 'off');  
    set(handles.AlignActivty_PushButton, 'Enable', 'off');
end


function Session_EditableText_Callback(hObject, eventdata, handles)
% hObject    handle to Session_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Session_EditableText as text
%        str2double(get(hObject,'String')) returns contents of Session_EditableText as a double


% --- Executes during object creation, after setting all properties.
function Session_EditableText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Session_EditableText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Timecourse_PushButton.
function Timecourse_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Timecourse_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

listpos = get(handles.AnimalName_ListBox, 'Value');
list = get(handles.AnimalName_ListBox, 'String');
h1 = waitbar(0, 'Initializing...');
% datafolder = getappdata(KomiyamaLabHub, 'Folder');
datafolder = 'C:\Users\Komiyama\Desktop\ActivitySummary_UsingRawData';

if length(listpos) == 1
    folder = dir(datafolder);
    cd(datafolder);
    filestoanalyze = [];
    animal = list{listpos};
    for i = 1:length(folder)
        ismatch = strfind(folder(i).name, animal);
        wrongfile = strfind(folder(i).name, 'Timecourse');
        if ~isempty(ismatch) && isempty(wrongfile)
            load(folder(i).name)
            filestoanalyze = [filestoanalyze, ',', folder(i).name(1:end-4)];
        end
        waitbar(i/length(folder), h1, 'Looking for files...')
    end

    close(h1)
    filestoanalyze = filestoanalyze(2:end);

    eval(['SummarizeActivity(', filestoanalyze, ');']);
else
    cd('C:\Users\Komiyama\Desktop\Output Data')
    folder = dir('C:\Users\Komiyama\Desktop\Output Data');
    filestoanalyze = [];
    for n = 1:length(listpos)
        animal{n} = list{listpos(n)};
        for i = 1:length(folder)
            ismatch = strfind(folder(i).name, animal{n});
            rightfile = strfind(folder(i).name, 'Timecourse_Summary');
            if ~isempty(ismatch) && ~ isempty(rightfile)
                load(folder(i).name)
                filestoanalyze = [filestoanalyze, ',', folder(i).name(1:end-4)];
            end
            waitbar(i/(length(folder)*length(listpos)), h1, 'Looking for files')
        end
    end
    
    close(h1)
    filestoanalyze = filestoanalyze(2:end);
    
    eval(['TimecourseSummary(', filestoanalyze, ');'])
end





% --- Executes on button press in Clustering_PushButton.
function Clustering_PushButton_Callback(~, ~, handles)
% hObject    handle to Clustering_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

listpos = get(handles.AnimalName_ListBox, 'Value');
list = get(handles.AnimalName_ListBox, 'String');
h1 = waitbar(0, 'Initializing...');
Activity = get(handles.AnalyzeActivity_ToggleButton, 'Value');
Behavior = get(handles.AnalyzeBehavior_ToggleButton, 'Value');
% datafolder = getappdata(KomiyamaLabHub, 'Folder');
datafolder = 'C:\Users\Komiyama\Desktop\ActivitySummary_UsingRawData';

if Activity == 1 && Behavior == 0
    folder = dir(datafolder);
    cd(datafolder);
    filestoanalyze = [];
    if length(listpos) ==1
        source = 'single';
        animal = list{listpos};
        for i = 1:length(folder)
            ismatch = strfind(folder(i).name, animal);
            wrongfile = strfind(folder(i).name, 'Poly');
            if ~isempty(ismatch) && isempty(wrongfile)
                load(folder(i).name)
                filestoanalyze = [filestoanalyze, ',', folder(i).name(1:end-4)];
            end
            waitbar(i/length(folder), h1, 'Looking for files...')
        end
        filestoanalyze = filestoanalyze(2:end);
        eval([animal, '_ClusteringProfile = NHanalyClusteringAnalysis(', filestoanalyze, ',source);']);
        cd('C:\Users\Komiyama\Desktop\Output Data');    
        close(h1)
        fname = [animal, '_ClusteringProfile'];
        save(fname, fname);

    else
        animal = list(listpos);
        folder = dir('C:\Users\Komiyama\Desktop\Output Data');
        cd('C:\Users\Komiyama\Desktop\Output Data');    
        source = 'multiple';
        for i = 1:length(folder)
            ismatch = strfind(folder(i).name, '_ClusteringProfile');
            if isempty(ismatch)
                continue
            else
                ismatch2 = strncmp(folder(i).name, animal, 5);
            end
            if sum(ismatch2) ~=0
                load(folder(i).name)
                if isempty(filestoanalyze)
                    filestoanalyze = folder(i).name(1:end-4);
                else
                    filestoanalyze = [filestoanalyze, ',', folder(i).name(1:end-4)];
                end
            end
            waitbar(i/length(folder), h1, 'Looking for files')
        end
        eval(['NHanalyClusteringAnalysis(', filestoanalyze,',source);']);
        close(h1)
    end
             
elseif Activity == 1 && Behavior == 1
    numfiles = 0;
    filestoanalyze = [];
    if length(listpos) == 1
        folder = dir(datafolder);
        cd(datafolder);
        animal = list{listpos};
        for i = 1:length(folder)
            ispart = strfind(folder(i).name, animal);
            wrongfile = strfind(folder(i).name, 'Timecourse');
            wrongfile2 = strfind(folder(i).name, 'Poly');
            if ~isempty(ispart) && isempty(wrongfile) && isempty(wrongfile2)
                load(folder(i).name)
                eval(['filesesh = ', folder(i).name(1:end-4), '.Session;']);
                currentanimal = regexp(folder(i).name, '[A-Z]{2,3}0+\d+', 'match');
                currentanimal = currentanimal{1};
                usesessions = blacklist(currentanimal);
                if ~ismember(filesesh,usesessions)
                    clear(folder(i).name)
                    continue
                end
                numfiles = numfiles+1;
                filestoanalyze = [filestoanalyze, ',', folder(i).name(1:end-4)];
            end
            waitbar(i/length(folder), h1, 'Looking for files')
        end
        folder = dir('C:\Users\Komiyama\Desktop\Output Data');
        cd('C:\Users\Komiyama\Desktop\Output Data');
        for i = 1:length(folder)
            ismatchCorrelations = strfind(folder(i).name, [animal, '_Correlations']);
            ismatchStats = strfind(folder(i).name, [animal, '_StatClassified']);
            if ~isempty(ismatchCorrelations) || ~isempty(ismatchStats)
                load(folder(i).name)
%                 eval([inputlist{n}, 'file = folder(', num2str(i), ').name(1:end-4);']);
                if ~isempty(ismatchCorrelations)
                    Correlationsfile = folder(i).name(1:end-4);
                elseif ~isempty(ismatchStats)
                    StatClassifiedfile = folder(i).name(1:end-4);
                end
            end
        end
        filestoanalyze = filestoanalyze(2:end);
        close(h1);
        eval(['ClusterBehaviorCorrelations(', ...
            Correlationsfile, ',', ...
            StatClassifiedfile ',', ...
            filestoanalyze, ');'])
    else
        folder = dir('C:\Users\Komiyama\Desktop\Output Data');
        cd('C:\Users\Komiyama\Desktop\Output Data');
        animals = list(listpos);
        for i = 1:length(folder)
            ismatch = strncmp(folder(i).name, animals, 5);
            isfile = strfind(folder(i).name, 'SpineCorrelationTimecourse');
            if sum(ismatch) ~=0 && ~isempty(isfile)
                load(folder(i).name)
%                 currentanimal = regexp(folder(i).name, '[A-Z]{2,3}0+\d+', 'match');
%                 currentanimal = currentanimal{1};
%                 [usesessions] = blacklist(currentanimal);
%                 eval([folder(i).name(1:end-4), ' = structfun(@(x) x([', num2str(usesessions), ']), ', folder(i).name(1:end-4), ', ''uni'', false)'])
                filestoanalyze = [filestoanalyze, ',', folder(i).name(1:end-4)];
            end
            waitbar(i/length(folder), h1, 'Looking for files')
        end
        filestoanalyze = filestoanalyze(2:end);
        close(h1)
        eval(['ClusterBehaviorCorrelations(', filestoanalyze, ');']);
    end
end


% --- Executes on button press in BehaviorTimecourse_PushButton.
function BehaviorTimecourse_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorTimecourse_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

listpos = get(handles.AnimalName_ListBox, 'Value');
list = get(handles.AnimalName_ListBox, 'String');
h1 = waitbar(0, 'Initializing...');
user = get(handles.figure1, 'UserData');

if length(listpos) ==1
    folder = dir('D:\Sara\All Behavioral Data');
    cd('D:\Sara\All Behavioral Data');
    filestoanalyze = [];
    animal = list{listpos};
    session = get(handles.Session_EditableText, 'String');
    sessioncounter = 0;
%     actual_sessions = [];
    for i = 1:length(folder)
        rightanimal = strfind(folder(i).name, [animal, '_']);
        wrongfile = strfind(folder(i).name, 'Summarized');
        if ~isempty(rightanimal) && isempty(wrongfile)
            load(folder(i).name)
            nameerror = regexp(folder(i).name, '\d'); % if the initials of the user are more than two letters, the save name excludes one; account for this by finding when the first digit character occurs (e.g. position 3 vs. 2 for GLB001 vs. NH001)
            if nameerror(1) > 2
%                 filestoanalyze = [filestoanalyze, ',', folder(i).name(2:end-4)];
                date = regexp(folder(i).name, '\d{6}\w*_Behavior', 'match');
                wrkspc = who;
                filestoanalyze = [filestoanalyze, ',', wrkspc{~cellfun(@isempty, regexp(who, date))}];
            else
                filestoanalyze = [filestoanalyze, ',', folder(i).name(1:end-4)];
            end
            sessioncounter = sessioncounter +1;
%             actual_sessions = [actual_sessions; folder(i).Session];
        end
        waitbar(i/length(folder), h1, 'Looking for files...')
    end
    close(h1)
    filestoanalyze = filestoanalyze(2:end);
    seshcheck = str2num(session);
    if sessioncounter ~= length(seshcheck)
        disp('Make sure to enter the correct session numbers!');
        session = num2str(1:sessioncounter);
    end

%     if length(session)~=length(filestoanalyze)
%         session = 1:length(filestoanalyze)
%     end
    eval(['NHanalySummarizeBehavior(', filestoanalyze, ',[', session, ']);']);
else
    if strcmpi(user, 'Nathan')
        folder = dir('C:\Users\Komiyama\Desktop\Output Data');
        cd('C:\Users\Komiyama\Desktop\Output Data');
    elseif strcmpi(user, 'Giulia')
        folder = dir('C:\Users\komiyama\Desktop\Giulia\All Behavioral Data');
        cd('C:\Users\komiyama\Desktop\Giulia\All Behavioral Data');
    end
    filestoanalyze = [];
    animal = list(listpos);
    for i = 1:length(folder);
        rightanimal = regexp(folder(i).name, animal);
        rightfile = strfind(folder(i).name, 'Summarized');
        if sum(cell2mat(rightanimal)) ~= 0 && ~isempty(rightfile)
            load(folder(i).name)
%             filestoanalyze = [filestoanalyze, ',', folder(i).name(1:end-4)];
            nameerror = regexp(folder(i).name, '\d'); % if the initials of the user are more than two letters, the save name excludes one; account for this by finding when the first digit character occurs (e.g. position 3 vs. 2 for GLB001 vs. NH001)
            if nameerror(1) > 2
                user = regexp(folder(i).name, '[A-Z]{2,3}', 'match'); user = user{1};
                animalnum = regexp(folder(i).name, '0{1,3}[A-Z,0-9]*', 'match'); animalnum = animalnum{1};
                wrkspc = who;
                combos = nchoosek(user,2);  %%% Cycle through all combinations of the intitials to find the closest match
                for n = 1:size(combos,1)
                    successfullyloaded = [];
                    options{n} = [combos(n,1:2), animalnum, '_SummarizedBehavior'];
                    if sum(~cellfun(@isempty, regexp(who,options{n})))
                        filestoanalyze = [filestoanalyze, ',', wrkspc{~cellfun(@isempty, regexp(who,options{n}))}];
                        successfullyloaded = wrkspc{~cellfun(@isempty, regexp(who,options{n}))};
                    else
                    end
                end
%                 if isempty(successfullyloaded)
%                 end
            else
                filestoanalyze = [filestoanalyze, ',', folder(i).name(1:end-4)];
            end

        end
        waitbar(i/length(folder), h1, 'Looking for files')
    end
    close(h1)
    startpoint = regexp(filestoanalyze, '\w'); startpoint = startpoint(1);
    filestoanalyze = filestoanalyze(startpoint:end);
    eval(['AverageBehavior(', filestoanalyze, ')']);
end


% --- Executes on button press in AlignActivty_PushButton.
function AlignActivty_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to AlignActivty_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global LeverTracePlots

listpos = get(handles.AnimalName_ListBox, 'Value');
list = get(handles.AnimalName_ListBox, 'String');
h1 = waitbar(0, 'Initializing...');

Behavior = cell(14,1);
Activity = cell(14,1);

Beh_folder = dir('C:\Users\Komiyama\Desktop\Behavioral Data\All Summarized Behavior Files list');
% storeddata = getappdata(KomiyamaLabHub);
% datafolder = storeddata.Folder;
datafolder = 'C:\Users\Komiyama\Desktop\ActivitySummary_UsingRawData';
Act_folder = dir(datafolder);
Output_folder = dir('C:\Users\Komiyama\Desktop\Output Data');

scsz = get(0, 'ScreenSize');

LeverTracePlots.figure = figure('Position', scsz);

if length(listpos) ==1
    animal = list{listpos};
    for i = 1:length(Act_folder)
        cd(datafolder);
        ismatch = strfind(Act_folder(i).name, animal);
        wrongfile = strfind(Act_folder(i).name, 'Poly');

        if ~isempty(ismatch) && isempty(wrongfile) 
            load(Act_folder(i).name)
            eval(['currentsession = ', Act_folder(i).name(1:end-4), '.Session;'])
            Activity{currentsession} = Act_folder(i).name(1:end-4);
            cd('C:\Users\Komiyama\Desktop\Behavioral Data\All Summarized Behavior Files list');
            if currentsession == 10
                a = 1;
            end
            for j = 1:length(Beh_folder)
                areboth = strncmp(Beh_folder(j).name, Activity(currentsession), 12);
                if areboth
                    load(Beh_folder(j).name)
                    Behavior{currentsession} = Beh_folder(j).name(1:end-4);
                end
            end
            if ~isempty(Activity{currentsession}) && ~isempty(Behavior{currentsession})
                eval(['[', animal, '_Correlations{', num2str(currentsession), '},', ...
                    animal, '_StatClassified{', num2str(currentsession), '},', ...
                    animal, '_TrialInformation{', num2str(currentsession),...
                    '}] = NHanalyAlignBehavior(', Activity{currentsession},',', Behavior{currentsession}, ');'])
                clear(Activity{currentsession})
                clear(Behavior{currentsession})
            elseif ~isempty(Activity{currentsession}) && isempty(Behavior{currentsession})
                clear(Activity{currentsession})
            elseif isempty(Activity{currentsession}) && ~isempty(Behavior{currentsession})
                clear(Behavior{currentsession})
            else
            end
        end
        waitbar(i/(length(Act_folder)), h1, 'Finding and aligning files...')
    end
    close(h1)
    fnameCorrelations = [animal, '_Correlations'];
    fnameStatClass = [animal, '_StatClassified'];
    fnameTrial = [animal, '_TrialInformation'];
    
    cd('C:\Users\Komiyama\Desktop\Output Data');
    
    save(fnameCorrelations, fnameCorrelations);
    save(fnameStatClass, fnameStatClass);
    save(fnameTrial, fnameTrial);
    toclear = who(['*', animal, '*']);
    for c = 1:length(toclear)
        clear(toclear{c})
    end
    disp(['Animal ', animal, ' alignment complete'])
else
    animals = list(listpos);
    count = 0;
    for f = 1:length(listpos)
        for i = 1:length(Act_folder)
            animal = animals{f};
            cd('C:\Users\Komiyama\Desktop\ActivitySummary');
            ismatch = strfind(Act_folder(i).name, animal);
            wrongfile = strfind(Act_folder(i).name, 'Poly');
            
            if ~isempty(ismatch) && isempty(wrongfile) 
                load(Act_folder(i).name)
                eval(['currentsession = ', Act_folder(i).name(1:end-4), '.Session;'])
                Activity{currentsession} = Act_folder(i).name(1:end-4);
                cd('C:\Users\Komiyama\Desktop\Behavioral Data\All Summarized Behavior Files list');
                for j = 1:length(Beh_folder)
                    areboth = strncmp(Beh_folder(j).name, Activity(currentsession), 12);
                    if areboth
                        load(Beh_folder(j).name)
                        Behavior{currentsession} = Beh_folder(j).name(1:end-4);
                    end
                end
                if ~isempty(Activity{currentsession}) && ~isempty(Behavior{currentsession})
                    eval(['[', animal, '_StatClassified{', num2str(currentsession), '},', ...
                        animal, '_TrialInformation{', num2str(currentsession),...
                        '}] = NHanalyAlignBehavior(', Activity{currentsession},',', Behavior{currentsession}, ');'])
                    clear(Activity{currentsession})
                    clear(Behavior{currentsession})
                elseif ~isempty(Activity{currentsession}) && isempty(Behavior{currentsession})
                    clear(Activity{currentsession})
                elseif isempty(Activity{currentsession}) && ~isempty(Behavior{currentsession})
                    clear(Behavior{currentsession})
                else
                end
            end
            count = count+1;
            waitbar(count/(length(Act_folder)*length(listpos)), h1, 'Finding and aligning files...')
        end
        fnameStatClass = [animal, '_StatClassified'];
        fnameTrial = [animal, '_TrialInformation'];

        cd('C:\Users\Komiyama\Desktop\Output Data');

        save(fnameStatClass, fnameStatClass);
        save(fnameTrial, fnameTrial);
        clear(fnameStatClass)
        clear(fnameTrial)
        toclear = who(['*', animal, '*']);
    %     fileclear = who(['*', fname, '*']);
        for c = 1:length(toclear)
            clear(toclear{c})
        end
    %     for c = 1:length(fileclear)
    %         clear(fileclear{c})
    %     end
        disp(['Animal ', animal, ' alignment complete'])
    end
end

% for i = 1:14
%     if ~isempty(Activity{i}) && ~ isempty(Behavior{i})
%         eval([animal, '_ActBehCorrelations{', num2str(i), '} = NHanalyAlignBehavior(', Activity{i},',', Behavior{i}, ');'])
%     end
%     waitbar((length(Act_folder)+length(Beh_folder)+i)/(length(Act_folder)+length(Beh_folder)+14), h1, ['Aligning activity and behavior of session', num2str(i)])
% end


% --------------------------------------------------------------------
function OpenCode_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to OpenCode_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ActivityCode_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to ActivityCode_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function BehaviorCode_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorCode_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ActivityBehavior_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to ActivityBehavior_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Alignment_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to Alignment_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pname = 'Z:\People\Nathan\Matlab\';
fname = 'KomiyamaLabHub.m';

matlab.desktop.editor.openAndGoToFunction([pname, fname], 'AlignActivty_PushButton_Callback')

edit NHanalyAlignBehavior

% --------------------------------------------------------------------
function ActivityBehaviorClustering_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to ActivityBehaviorClustering_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pname = 'Z:\People\Nathan\Matlab\';
fname = 'KomiyamaLabHub.m';

matlab.desktop.editor.openAndGoToFunction([pname, fname], 'Clustering_PushButton_Callback')

edit ClusterBehaviorCorrelations


% --------------------------------------------------------------------
function BehaviorTimecourse_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorTimecourse_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pname = 'Z:\People\Nathan\Matlab\';
fname = 'KomiyamaLabHub.m';

matlab.desktop.editor.openAndGoToFunction([pname, fname], 'Behavior_PushButton_Callback')

edit NHanalySummarizeBehavior

% --------------------------------------------------------------------
function ActivityTimecourse_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to ActivityTimecourse_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pname = 'Z:\People\Nathan\Matlab\';
fname = 'KomiyamaLabHub.m';

matlab.desktop.editor.openAndGoToFunction([pname, fname], 'Timecourse_PushButton_Callback')

edit SummarizeActivity

% --------------------------------------------------------------------
function ActivityClustering_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to ActivityClustering_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pname = 'Z:\People\Nathan\Matlab\';
fname = 'KomiyamaLabHub.m';

matlab.desktop.editor.openAndGoToFunction([pname, fname], 'Clustering_PushButton_Callback')

edit NHanalyClusteringAnalysis.m


% --------------------------------------------------------------------
function Diagnostics_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to Diagnostics_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SpineDendOverlap_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to SpineDendOverlap_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CalcSpineDendOverlap;


% --------------------------------------------------------------------
function SpineDendFitting_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to SpineDendFitting_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotRandomSpineDendFits;


% --------------------------------------------------------------------
function InspectRandomTraces_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to InspectRandomTraces_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plotspineanddenddata;


% --- Executes on button press in DendSubtracted_CheckBox.
function DendSubtracted_CheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to DendSubtracted_CheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DendSubtracted_CheckBox

set(handles.DendExcluded_CheckBox, 'value', 0);


% --- Executes on button press in DendExcluded_CheckBox.
function DendExcluded_CheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to DendExcluded_CheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DendExcluded_CheckBox

set(handles.DendSubtracted_CheckBox, 'value', 0);


% --------------------------------------------------------------------
function Redo_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to Redo_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function AllAnalysis_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to AllAnalysis_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

RedoAnalysis

% --------------------------------------------------------------------
function DendSub_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to DendSub_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

RedoSubtraction

% --------------------------------------------------------------------
function CorrShuff_DropDown_Callback(hObject, eventdata, handles)
% hObject    handle to CorrShuff_DropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ShuffleSpines;
