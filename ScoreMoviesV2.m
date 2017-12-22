







function ScoreMovies()

clc;

global CurrentFrame;
global Cs;
global ElapsedTime;
global Frame;
global FrameInspector;
global FrameMask;
global HandIsOff;
global HandOffVec;
global HandToggle;
global Initials;
global NumFrames;
global obj;
global ResetToggle;
global ResetVecs;
global SaveFileName;
global T;
global MovieImage;


[fname p] = uigetfile('*.avi*','SELECT MOVIE TO SCORE','SELECT MOVIE TO SCORE');
if fname ==0
    disp('Program exited.');
    return;
end
MovieFileName = [p fname];

% [fname p] = uigetfile('*.mat*','SELECT DATA FILE','SELECT DATA FILE');

disp('Loading movie...')
disp('May take up to 5 minutes.');
obj=VideoReader(MovieFileName);
disp('Video loaded.')
disp(' ');
disp('Commands:');
disp('H - mark hand off');
disp('R - reset frames');
disp('Right arrow - move forward 1 frame');
disp('Left arrow - move backward 1 frame');
T = tic;



NumFrames=obj.NumberOfFrames;
CurrentFrame=1;

Cs = [    0.8941    0.1020    0.1098;0.2157    0.4941    0.7216;0.3020    0.6863    0.2902;0.5961    0.3059    0.6392;1.0000    0.4980         0];
FrameInspector = zeros(1,NumFrames);
HandOffVec = zeros(1,NumFrames);
HandToggle = 0;
ResetToggle = 0;

MatFileName = [p fname(1:end-3) 'mat'];
SaveFileName = MatFileName;

load(MatFileName,'FrameInspector','FrameMask','HandOffVec','Initials','ElapsedTime','LastFrame');
Initials = input('Initials: ','s');

FrameMask = FrameMask(:,:,1);

if LastFrame ~= NumFrames
    disp('DATA FILE DOES NOT MATCH MOVIE FILE. EXITING.');
    return;
end

fig = figure;
set(fig,'WindowScrollWheelFcn' , @scrollMove);
set(fig,'KeyPressFcn',@key_move);
set(fig,'units','normalized');

Frame = uicontrol('style','edit',...
    'units','normalized',...
    'position',[0.1300    0.1100    0.1566    0.2157],...
    'string','Frame',...
    'callback',{@setFrame},...
    'fontsize',40);


ResetVecs = uicontrol('style','pushbutton',...
    'units','normalized',...
    'position',[0.3361    0.1100    0.1566    0.2157],...
    'string','Reset',...
    'callback',{@mark_reset_vecs},...
    'fontsize',20);

HandIsOff = uicontrol('style','pushbutton',...
    'units','normalized',...
    'position',[0.5422    0.1100    0.1566    0.2157],...
    'string','Hand Off',...
    'callback',{@mark_hand_off},...
    'fontsize',20);



uicontrol('style','pushbutton',...
    'units','normalized',...
    'position',[0.7484    0.1100    0.1566    0.2157],...
    'string','Save',...
    'callback',{@saveVecs},...
    'fontsize',20);

displayFrame(CurrentFrame,1)


function scrollMove(varargin)
global NumFrames;
global CurrentFrame;

frame_move = varargin{2}.VerticalScrollCount;
if CurrentFrame+frame_move<1 || CurrentFrame+frame_move>NumFrames
    frame_move = 0;
end

CurrentFrame = CurrentFrame+frame_move;
displayFrame(CurrentFrame)

function key_move(varargin)
global NumFrames;
global CurrentFrame;
global StimVec

tmp = varargin{2};

if strcmp('rightarrow',tmp.Key)
    frame_move = 1;
elseif strcmp('leftarrow',tmp.Key)
    frame_move = -1;
elseif strcmp('h',tmp.Key)
    mark_hand_off();
    return;
elseif strcmp('r',tmp.Key)
    mark_reset_vecs();
    return;
elseif strcmp('end',tmp.Key)
    CurrentFrame = find(StimVec,1,'last');
    frame_move = 0;
    return;
elseif strcmp('space',tmp.Key)
    frame_move = 10;
elseif strcmp('home',tmp.Key);
    frame_move = -CurrentFrame+1;
else
    return;
end

if CurrentFrame+frame_move<1 || CurrentFrame+frame_move>NumFrames
    frame_move = 0;
end

CurrentFrame = CurrentFrame+frame_move;

displayFrame(CurrentFrame)

function setFrame(source,event)
global Frame;
global CurrentFrame
CurrentFrame = str2double(get(Frame,'string'));
displayFrame(str2double(get(Frame,'string')));

function mark_hand_off(source,event)
global HandIsOff
global CurrentFrame
global HandOffVec
global HandToggle;

if HandToggle == 0
    HandToggle = CurrentFrame;
    set(HandIsOff,'string','Hand On','backgroundcolor','r');
else
    HandOffVec(min([HandToggle CurrentFrame]):max([HandToggle CurrentFrame])) = 1;
    HandToggle = 0;
    set(HandIsOff,'string','Hand On','backgroundcolor','w');
    saveVecs();
    displayFrame(CurrentFrame);
    
end

function mark_reset_vecs(source,event)
global ResetVecs
global CurrentFrame
global ResetToggle
global HandOffVec

if ResetToggle == 0
    ResetToggle = CurrentFrame;
    set(ResetVecs,'string','RESETTING','backgroundcolor','r');
else
    HandOffVec(min([ResetToggle CurrentFrame]):max([ResetToggle CurrentFrame])) = 0;
    ResetToggle = 0;
    set(ResetVecs,'string','Reset','backgroundcolor','w');
    saveVecs();
    displayFrame(CurrentFrame);
    
end

function displayFrame(InputFrame,FirstFlag)
global obj;
global Frame;
global HandOffVec;
global Cs;
global FrameMask;
global FrameInspector;
global CurrentFrame;
global MovieImage;
global SmallHandOffPlot;
global LargeHandOffPlot;
global SmallHandOffPlotVerticalLine;
global LargeHandOffPlotVerticalLine;
global S4;
global S8;

img=read(obj,InputFrame);

if nargin == 1
    FirstFlag = 0;
end

if FirstFlag
    subplot(3,4,[1 2 3 5 6 7]);cla;
    MovieImage = imshow(img(:,:,1).*FrameMask,[]);
    set(gca,'xlim',[find(any(FrameMask(:,:,1)),1,'first'),find(any(FrameMask(:,:,1)),1,'last')]);
    set(gca,'ylim',[find(any(FrameMask(:,:,1)'),1,'first'),find(any(FrameMask(:,:,1)'),1,'last')]);
    
    S4 = subplot(3,4,4);cla;hold(S4,'on');
    SmallHandOffPlot = plot(HandOffVec,'color',Cs(1,:),'linewidth',2);
    SmallHandOffPlotVerticalLine = plot([InputFrame InputFrame],[-0.1 1.5],'k:','linewidth',2);
    set(gca,'xtick',[InputFrame-300 InputFrame InputFrame+300]);
    set(gca,'xticklabel',{num2str(InputFrame-300),'Current Frame',num2str(InputFrame+300)});
    set(gca,'ytick',[]);
    set(gca,'xlim',[InputFrame-300 InputFrame+300]);
    set(gca,'ylim',[-0.1 1.5]);
    
    S8 = subplot(3,4,8);cla;hold(S8,'on');
    LargeHandOffPlot = plot(HandOffVec,'color',Cs(1,:),'linewidth',2);
    LargeHandOffPlotVerticalLine = plot([InputFrame InputFrame],[-0.1 1.5],'k:','linewidth',2);
    set(gca,'xlim',[1 length(HandOffVec)]);
    set(gca,'ylim',[-0.1 1.5]);
    set(gca,'xtick',[1 length(HandOffVec)]);
    set(gca,'xticklabel',{'Beginning of Movie','End of Movie'});
    set(gca,'ytick',[]);
    
else
    set(MovieImage,'cdata',img(:,:,1).*FrameMask);
    set(SmallHandOffPlot,'xdata',1:length(HandOffVec),'ydata',HandOffVec)
    set(LargeHandOffPlot,'xdata',1:length(HandOffVec),'ydata',HandOffVec)
    set(SmallHandOffPlotVerticalLine,'xdata',[InputFrame InputFrame]);
    set(LargeHandOffPlotVerticalLine,'xdata',[InputFrame InputFrame]);
    set(S4,'xlim',[InputFrame-300 InputFrame+300]);
    set(S4,'xtick',[InputFrame-300 InputFrame InputFrame+300]);
    set(S4,'xticklabel',{num2str(InputFrame-300),'Current Frame',num2str(InputFrame+300)});
end

set(Frame,'string',num2str(InputFrame));

FrameInspector(CurrentFrame) = 1;

drawnow;


function saveVecs(varargin)

global SaveFileName
global HandOffVec;
global FrameInspector;
global Initials;
global CurrentFrame;
global ElapsedTime;
global T;
ElapsedTime = ElapsedTime+toc(T);
save(SaveFileName,'HandOffVec','FrameInspector','Initials','CurrentFrame','ElapsedTime','-append');


