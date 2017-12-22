function tephys

%dh5.16.07 changed from mirrorControl_testEphys.m to  dephys...
%dh101607 adapted dephys for new DOM PC

%Clear the workspace (prompt first).
%clearWorkspace;

%Set up a wait bar to show the progress.
wb = waitbarWithCancel(0, 'Starting ephys...', 'Name', 'Loading Physiology software...');
pos = get(wb, 'Position');
pos(2) = pos(2) - pos(4);
set(wb, 'Position', pos);

%--------------------------------------------------
%Set up the triggering.   <<<<<<<<<<---------- CONFIG
acqJob = daqjob('acquisition');
scopeJob = daqjob('scope');
%NEW TRIGGERING CONFIGURATION
setTriggerOrigin(acqJob, '/dev3/port0/line0');
setTriggerOrigin(scopeJob, '/dev3/port0/line0');
setTriggerDestination(acqJob, 'PFI1');
setTriggerDestination(scopeJob, 'PFI1');
%ORIGINAL TRIGGERING CONFIGURATION
% setTriggerOrigin(acqJob, '/dev1/port0/line0');
% setTriggerOrigin(scopeJob, '/dev1/port0/line0');
% setTriggerDestination(acqJob, 'PFI0');
% setTriggerDestination(scopeJob, 'PFI0');



%----- Include following from startup.m to reset digital lines
% % to low after MATLAB bug sets them high upon creation of AI/AO objects: 
% % 7/27/06: DHO,TO
% disp('dephys: Setting digital IO lines low ...')
% %TO072406B - Modified to cover all lines on all boards, and be compatible across systems. -- Tim O'Connor 7/24/06
% info = daqhwinfo('nidaq');
% for i = 1 : size(info.ObjectConstructorName, 1)
%     if ~isempty(info.ObjectConstructorName{i, 3})
%         dio = eval(info.ObjectConstructorName{i, 3});
%         %     dio = digitalio('nidaq', i);
%         for j = 0 : 7
%             line = addLine(dio, j, 'out');
%             putvalue(line, 0);
%         end
%         delete(dio);
%     end
% end
% 

%--------------------------------------------------
%Start ephys
waitbar(0.05, wb, 'Opening ephys...');
if isWaitbarCancelled(wb)
    delete(wb);
    return;
end

ep = program('ephys', 'ephys', 'ephys');

openprogram(progmanager, ep);

%--------------------------------------------------
%Set up amplifier(s).   <<<<<<<<<<---------- CONFIG
waitbar(0.10, wb, 'Creating amplifiers...');
if isWaitbarCancelled(wb)
    delete(wb);
    return;
end

%The scaledOutput is the "output" of the amplifier, as stated on the BNC
%connector of the Axon hardware, even though it is seen as an input from
%the software point of view.
%The control signal, vCom, is an output from the software's point of view.
patch{1} = axopatch_200B('scaledOutputBoardID', 3, 'scaledOutputChannelID', 4, 'gain_daq_board_id', 3, 'gain_channel', 1, 'mode_daq_board_id', 3, 'mode_channel', 2,  ...
    'v_hold_daq_board_id', 3,  'v_hold_channel', 3, 'vComBoardID', 3, 'vComChannelID', 0);%<<<<<<<<<<---------- CONFIG
%Sample @dumbamp config. With an input and an output.
% patch{2} = dumbamp('name', 'EEG/LED', 'input_gain', 1, 'input_offset', 0, 'inputBoardID', 2, 'inputChannelID', 5, 'inputName', 'EEG', ...
%     'output_gain', 1, 'output_offset', 0, 'outputBoardID', 3, 'outputChannelID', 3, 'outputName', 'LED'); 
%See TO012808A for changes related to implementing @dumbamp.
%patch{2} = dumbamp('name', 'EEG-In', 'input_gain', 1, 'input_offset', 0, 'inputBoardID', 3, 'inputChannelID', 5, 'inputName', 'EEG-In'); 
% patch{3} = dumbamp('name', 'LED-In', 'input_gain', 1, 'input_offset', 0, 'inputBoardID', 3, 'inputChannelID', 6, 'inputName', 'LED-In'); 

for i = 1 : length(patch)
    bindToDaqJob(patch{i}, acqJob);
end

%--------------------------------------------------
%Add the amplifier(s) to ephys.
waitbar(0.12, wb, 'Mounting amplifiers in ephys...');
if isWaitbarCancelled(wb)
    delete(wb);
    return;
end

% %Create shared aomux/aimux objects.
% aim = aimux(getDaqmanager);
% aom = aomux(getDaqmanager);

% %----- Include following from startup.m to reset digital lines
% % to low after MATLAB bug sets them high upon creation of AI/AO objects: 
% % 7/27/06: DHO,TO
% disp('dephys: Setting digital IO lines low ...')
% %TO072406B - Modified to cover all lines on all boards, and be compatible across systems. -- Tim O'Connor 7/24/06
% info = daqhwinfo('nidaq');
% for i = 1 : size(info.ObjectConstructorName, 1)
%     if ~isempty(info.ObjectConstructorName{i, 3})
%         dio = eval(info.ObjectConstructorName{i, 3});
%         %     dio = digitalio('nidaq', i);
%         for j = 0 : 7
%             line = addLine(dio, j, 'out');
%             putvalue(line, 0);
%         end
%         delete(dio);
%     end
% end

% setGlobal(progmanager, 'aimux', 'ephys', 'ephys', aim);
% setGlobal(progmanager, 'aomux', 'ephys', 'ephys', aom);

ep = getGlobal(progmanager, 'hObject', 'ephys', 'ephys');

ephys_setAmplifiers(ep, patch);

% %----- Include following from startup.m to reset digital lines
% % to low after MATLAB bug sets them high upon creation of AI/AO objects: 
% % 7/27/06: DHO,TO
% disp('tephys: Setting digital IO lines low ...')
% %TO072406B - Modified to cover all lines on all boards, and be compatible across systems. -- Tim O'Connor 7/24/06
% info = daqhwinfo('nidaq');
% for i = 1 : size(info.ObjectConstructorName, 1)
%     if ~isempty(info.ObjectConstructorName{i, 3})
%         dio = eval(info.ObjectConstructorName{i, 3});
%         %     dio = digitalio('nidaq', i);
%         for j = 0 : 7
%             line = addLine(dio, j, 'out');
%             putvalue(line, 0);
%         end
%         delete(dio);
%     end

%---------------
% end-----------------------------------
%Open the stimulator gui.
waitbar(0.20, wb, 'Opening stimulator...');
if isWaitbarCancelled(wb)
    delete(wb);
    return;
end

stim = program('stimulator', 'stimulator', 'stimulator');
openprogram(progmanager, stim);




%--------------------------------------------------
%Configure stimulator channels
waitbar(0.25, wb, 'Mounting channels in stimulator...');
if isWaitbarCancelled(wb)
    delete(wb);
    return;
end

%Here's the channel(s) for the stimulator.
stimChannels(1).channelName = 'LED1';%<<<<<<<<<<---------- CONFIG
stimChannels(1).boardID = 4;%<<<<<<<<<<---------- CONFIG
stimChannels(1).channelID = 1;%<<<<<<<<<<---------- CONFIG
% 
% %Here's the channel(s) for the stimulator.
% stimChannels(2).channelName = 'Airpuff';%<<<<<<<<<<---------- CONFIG
% stimChannels(2).boardID = 4;%<<<<<<<<<<---------- CONFIG
% stimChannels(2).channelID = 2;%<<<<<<<<<<---------- CONFIG
% 
% %Here's the channel(s) for the stimulator.
% stimChannels(3).channelName = 'reserve';%<<<<<<<<<<---------- CONFIG
% stimChannels(3).boardID = 4;%<<<<<<<<<<---------- CONFIG
% stimChannels(3).channelID = 3;%<<<<<<<<<<---------- CONFIG

% setGlobal(progmanager, 'aomux', 'stimulator', 'stimulator', aom);
stim_setChannels(stim, stimChannels);
% 
% %----- Include following from startup.m to reset digital lines
% % to low after MATLAB bug sets them high upon creation of AI/AO objects: 
% % 7/27/06: DHO,TO
% disp('tephys: Setting digital IO lines low ...')
% %TO072406B - Modified to cover all lines on all boards, and be compatible across systems. -- Tim O'Connor 7/24/06
% info = daqhwinfo('nidaq');
% for i = 1 : size(info.ObjectConstructorName, 1)
%     if ~isempty(info.ObjectConstructorName{i, 3})
%         dio = eval(info.ObjectConstructorName{i, 3});
%         %     dio = digitalio('nidaq', i);
%         for j = 0 : 7
%             line = addLine(dio, j, 'out');
%             putvalue(line, 0);
%         end
%         delete(dio);
%     end
% end


%--------------------------------------------------
%Open the acquirer gui.
waitbar(0.80, wb, 'Opening acquirer...');
if isWaitbarCancelled(wb)
    delete(wb);
    return;
end

acq = program('acquirer', 'acquirer', 'acquirer');
openprogram(progmanager, acq);
%cycler_registerProgram(cyclerObj, acq);
%registerLoopable(lm, {@acq_loopListener, acq}, 'acquirer')

%--------------------------------------------------
%Configure acquirer channels
waitbar(.85, wb, 'Mounting channels in acquirer...');
if isWaitbarCancelled(wb)
    delete(wb);
    return;
end

%Here's the channel(s) for the acquirer.
acqChannels(1).channelName = 'Trial_number';%<<<<<<<<<<---------- CONFIG
acqChannels(1).boardID = 2;%<<<<<<<<<<---------- CONFIG
acqChannels(1).channelID = 5;%<<<<<<<<<<---------- CONFIG

% acqChannels(2).channelName = 'input_1';%<<<<<<<<<<---------- CONFIG
% acqChannels(2).boardID = 2;%<<<<<<<<<<---------- CONFIG
% acqChannels(2).channelID = 2;%<<<<<<<<<<---------- CONFIG
% 
% acqChannels(3).channelName = 'Behavior_Trial_Number';%<<<<<<<<<<---------- CONFIG
% acqChannels(3).boardID = 2;%<<<<<<<<<<---------- CONFIG
% acqChannels(3).channelID = 3;%<<<<<<<<<<---------- CONFIG

% acqChannels(4).channelName = 'input_3';%<<<<<<<<<<---------- CONFIG
% acqChannels(4).boardID = 2;%<<<<<<<<<<---------- CONFIG
% acqChannels(4).channelID = 3;%<<<<<<<<<<---------- CONFIG
% 

%setGlobal(progmanager, 'aimux', 'acquirer', 'acquirer', aim);
acq_setChannels(acq, acqChannels);


%--------------------------------------------------
%Open the xsg.
waitbar(0.30, wb, 'Opening experimentSavingGui...');
if isWaitbarCancelled(wb)
    delete(wb);
    return;
end

xsg_getFilename;

%--------------------------------------------------
%Open the headerGUI.
waitbar(0.33, wb, 'Starting headerGui...');
if isWaitbarCancelled(wb)
    delete(wb);
    return;
end

hGui = program('headerGUI', 'headerGUI', 'headerGUI');
openprogram(progmanager, hGui);


%--------------------------------------------------
%Open the scope.
waitbar(0.35, wb, 'Opening scopeGui...');
if isWaitbarCancelled(wb)
    delete(wb);
    return;
end

%Start the ephys-enabled scope.
scg = program('scopeGui', 'scopeGui', 'scopeGui', 'ephysScopeAccessory', 'ephysScopeAccessory');
openprogram(progmanager, scg);
ephysAcc = getGlobal(progmanager, 'hObject', 'ephysScopeAccessory', 'ScopeGui');

%set(startmanager('ephysAcc'), 'triggerBoardID', get(sm, 'triggerBoardID'), 'triggerChannelID', get(sm, 'triggerChannelID'));%Pull these parameters from the 'acquisition' instance.
% %----- Include following from startup.m to reset digital lines
% % to low after MATLAB bug sets them high upon creation of AI/AO objects: 
% % 7/27/06: DHO,TO
% disp('tephys: Setting digital IO lines low ...')
% %TO072406B - Modified to cover all lines on all boards, and be compatible across systems. -- Tim O'Connor 7/24/06
% info = daqhwinfo('nidaq');
% for i = 1 : size(info.ObjectConstructorName, 1)
%     if ~isempty(info.ObjectConstructorName{i, 3})
%         dio = eval(info.ObjectConstructorName{i, 3});
%         %     dio = digitalio('nidaq', i);
%         for j = 0 : 7
%             line = addLine(dio, j, 'out');
%             putvalue(line, 0);
%         end
%         delete(dio);
%     end
% end


%--------------------------------------------------
%Add the amplifier(s) to the scope.
waitbar(0.45, wb, 'Mounting amplifiers in scopeGui...');
if isWaitbarCancelled(wb)
    delete(wb);
    return;
end

for i = 1 : length(patch)
    bindToDaqJob(patch{i}, scopeJob);
end


ephysAcc_setAmplifiers(ephysAcc, patch);

%----- Include following from startup.m to reset digital lines
% % to low after MATLAB bug sets them high upon creation of AI/AO objects: 
% % 7/27/06: DHO,TO
% disp('tephys: Setting digital IO lines low ...')
% %TO072406B - Modified to cover all lines on all boards, and be compatible across systems. -- Tim O'Connor 7/24/06
% info = daqhwinfo('nidaq');
% for i = 1 : size(info.ObjectConstructorName, 1)
%     if ~isempty(info.ObjectConstructorName{i, 3})
%         dio = eval(info.ObjectConstructorName{i, 3});
%         %     dio = digitalio('nidaq', i);
%         for j = 0 : 7
%             line = addLine(dio, j, 'out');
%             putvalue(line, 0);
%         end
%         delete(dio);
%     end
% end

% %--------------------------------------------------
% %Open the mapper.
% waitbar(0.50, wb, 'Opening mapper...');
% if isWaitbarCancelled(wb)
%     delete(wb);
%     return;
% end
% 
% dm = getDaqmanager;
% % if ~hasChannel(dm, 'pockelsCell')
% %     nameOutputChannel(dm, 2, 2, 'pockelsCell');%<<<<<<<<<<---------- CONFIG
% %     enableChannel(dm, 'pockelsCell');
% % end
% % if ~hasChannel(dm, 'shutter0')
% %     nameOutputChannel(dm, 3, 1, 'shutter0');
% %     enableChannel(dm, 'shutter0');
% % end
% if ~hasChannel(dm, 'photodiode1')
%     nameInputChannel(dm, 3, 0, 'photodiode1');%<<<<<<<<<<---------- CONFIG
%     enableChannel(dm, 'photodiode1');
% end
% % if ~hasChannel(dm, 'xMirror')
% %     nameOutputChannel(dm, 1, 0, 'xMirror');
% %     enableChannel(dm, 'xMirror');
% % end
% % if ~hasChannel(dm, 'yMirror')
% %     nameOutputChannel(dm, 1, 1, 'yMirror');
% %     enableChannel(dm, 'yMirror');
% % end
% 
% mapperObj = program('mapper', 'mapper', 'mapper');
% openprogram(progmanager, mapperObj);
% 
% setPreprocessor(getLocal(progmanager, stim, 'aomux'), 'pockelsCell', ...
%     {@mapper_pockelsCellPreprocessor, mapperObj}, 'mapper_pockelsCellPreprocessor');
% 
% %TEST -- TO091406TEST1
% setPreprocessor(getLocal(progmanager, stim, 'aomux'), 'xMirror', ...
%     {@mapper_mirrorChannelPreprocessor, mapperObj, 'X'}, 'mapper_mirrorChannelPreprocessorX');
% setPreprocessor(getLocal(progmanager, stim, 'aomux'), 'yMirror', ...
%     {@mapper_mirrorChannelPreprocessor, mapperObj, 'Y'}, 'mapper_mirrorChannelPreprocessorY');
% 
% moap = program('mapperOnlineAnalysisParameters', 'mapperOnlineAnalysisParameters', 'mapperOnlineAnalysisParameters');
% openprogram(progmanager, moap);

%--------------------------------------------------
%Open the pulseEditor.
waitbar(0.55, wb, 'Opening pulseEditor...');
if isWaitbarCancelled(wb)
    delete(wb);
    return;
end

pe = program('pulseEditor', 'pulseEditor', 'pulseEditor');
openprogram(progmanager, pe);

%TO022406D
try
    peCbm = getLocal(progmanager, pe, 'callbackManager');
    addCallback(peCbm, 'pulseCreation', {@ephys_pulseCreation, ep}, 'ephys_pulseCreation');
    addCallback(peCbm, 'pulseDeletion', {@ephys_pulseCreation, ep}, 'ephys_pulseDeletion');
    addCallback(peCbm, 'pulseSetCreation', {@ephys_pulseSetCreation, ep}, 'ephys_pulseSetCreation');
    addCallback(peCbm, 'pulseSetDeletion', {@ephys_pulseSetDeletionn, ep}, 'ephys_pulseSetDeletion');
    addCallback(peCbm, 'pulseUpdate', {@ephys_pulseCreation, ep}, 'ephys_pulseUpdate');
    addCallback(peCbm, 'pulseCreation', {@stim_pulseCreation, stim}, 'stim_pulseCreation');
    addCallback(peCbm, 'pulseDeletion', {@stim_pulseDeletion, stim}, 'stim_pulseDeletion');
    addCallback(peCbm, 'pulseSetCreation', {@stim_pulseSetCreation, stim}, 'stim_pulseSetCreation');
    addCallback(peCbm, 'pulseSetDeletion', {@stim_pulseSetDeletionn, stim}, 'stim_pulseSetDeletion');
    addCallback(peCbm, 'pulseUpdate', {@stim_pulseUpdate, stim}, 'stim_pulseUpdate');
catch
    warning('Error registering callbacks for pulseEditor events.');
end

%--------------------------------------------------
% %Open the cycler.
% waitbar(0.60, wb, 'Opening cycler...');
% if isWaitbarCancelled(wb)
%     delete(wb);
%     return;
% end
% 
% cyclerObj = program('cycler', 'cycler', 'cycler');
% openprogram(progmanager, cyclerObj);
% %Register programs that may be affected by the cycler.
% cycler_registerProgram(cyclerObj, ep);
% cycler_registerProgram(cyclerObj, getLocal(progmanager, stim, 'hObject'));
% bindCompletionListener(startmanager('acquisition'), {@cycler_Iterate, cyclerObj}, 'cyclerIterate');

%--------------------------------------------------
%Open the pulseEditor
waitbar(0.65, wb, 'Registering loopable components...');
if isWaitbarCancelled(wb)
    delete(wb);
    return;
end

[lg lm] = lg_factory;%This opens the loop gui.
registerLoopable(lm, {@ephys_loopListener, ep}, 'ephys');

[lg lm] = lg_factory;%This opens the loop gui.
registerLoopable(lm, {@shared_loopListener, ep}, 'ephys');
registerLoopable(lm, {@shared_loopListener, acq}, 'acquirer');
registerLoopable(lm, {@shared_loopListener, stim}, 'stimulator');

% % registerLoopable(lm, {@testLoopTriggering, startmanager('acquisition')}, 'trigger');
% registerLoopable(lm, {@cycler_loopListener, cyclerObj}, 'cycler');
% registerLoopable(lm, {@stim_loopListener, stim}, 'stimulator');

%--------------------------------------------------
%Open the userFcns gui.
waitbar(0.70, wb, 'Opening userFcns...');
if isWaitbarCancelled(wb)
    delete(wb);
    return;
end

userFcns = program('userFcns', 'userFcns', 'userFcns');
openprogram(progmanager, userFcns);

%--------------------------------------------------
%Open the autonotes gui.
waitbar(0.75, wb, 'Opening autonotes...');
if isWaitbarCancelled(wb)
    delete(wb);
    return;
end

userFcns = program('autonotes', 'autonotes', 'autonotes');
openprogram(progmanager, userFcns);



% %----- Include following from startup.m to reset digital lines
% % to low after MATLAB bug sets them high upon creation of AI/AO objects: 
% % 7/27/06: DHO,TO
% disp('dephys: Setting digital IO lines low ...')
% %TO072406B - Modified to cover all lines on all boards, and be compatible across systems. -- Tim O'Connor 7/24/06
% info = daqhwinfo('nidaq');
% for i = 1 : size(info.ObjectConstructorName, 1)
%     if ~isempty(info.ObjectConstructorName{i, 3})
%         dio = eval(info.ObjectConstructorName{i, 3});
%         %     dio = digitalio('nidaq', i);
%         for j = 0 : 7
%             line = addLine(dio, j, 'out');
%             putvalue(line, 0);
%         end
%         delete(dio);
%     end
% end

% %--------------------------------------------------
% %Open the mirror configuration gui.
% waitbar(.90, wb, 'Initializing imaging system...');
% if isWaitbarCancelled(wb)
%     delete(wb);
%     return;
% end
% 
% %TEST - TO091406TEST1
% % xMirrorBoardID = 1;%<<<<<<<<<<---------- CONFIG
% % xMirrorChannelID = 0;%<<<<<<<<<<---------- CONFIG
% % yMirrorBoardID = 1;%<<<<<<<<<<---------- CONFIG
% % yMirrorChannelID = 1;%<<<<<<<<<<---------- CONFIG
% % 
% % nameOutputChannel(dm, xMirrorBoardID, xMirrorChannelID, 'xMirror');
% % enableChannel(dm, 'xMirror');
% % nameOutputChannel(dm, yMirrorBoardID, yMirrorChannelID, 'yMirror');
% % enableChannel(dm, 'yMirror');
% 
% imagingSys = program('imagingSys', 'imagingSys', 'imagingSystemConfiguration');
% openprogram(progmanager, imagingSys);
% 
% setLocalBatch(progmanager, imagingSys, 'xBoardID', xMirrorBoardID, 'xChannelID', xMirrorChannelID, ...
%     'yBoardID', yMirrorBoardID, 'yChannelID', yMirrorChannelID, 'name', 'Mapper-Scanner01');
% 
% %Set the um/pixel for the video image here.
% setLocalBatch(progmanager, mapper, 'xVideoScaleFactor', 1820, 'yVideoScaleFactor', 1385);%<<<<<<<<<<---------- CONFIG

% %--------------------------------------------------
% %Open the mirror configuration gui.
% waitbar(.93, wb, 'Initializing photodiode configuration...');
% if isWaitbarCancelled(wb)
%     delete(wb);
%     return;
% end
% 
% pdiodeConfig = program('photodiode', 'photodiode', 'photodiodeConfiguration');
% openprogram(progmanager, pdiodeConfig);
% 
% pdiode = photodiode;
% set(pdiode, 'boardID', acqChannels(1).boardID, 'channelID', acqChannels(1).channelID, 'name', acqChannels(1).channelName);
% setLocalBatch(progmanager, pdiodeConfig, 'boardID', acqChannels(1).boardID, 'channelID', acqChannels(1).channelID, 'photodiodeName', acqChannels(1).channelName, 'photodiodeObject', pdiode);
% 
% setLocal(progmanager, mapper, 'photodiodeObject', pdiode);

%--------------------------------------------------
%Open the pulse hijacking gui.

waitbar(0.98, wb, 'Opening pulseJacker...');

pj = program('pulseJacker', 'pulseJacker', 'pulseJacker');
openprogram(progmanager, pj);

% addCallback(peCbm, 'pulseCreation', {@pj_pulseCreation, pj}, 'pj_pulseCreation');
% addCallback(peCbm, 'pulseDeletion', {@pj_pulseDeletion, pj}, 'pj_pulseDeletion');
% addCallback(peCbm, 'pulseSetCreation', {@pj_pulseSetCreation, pj}, 'pj_pulseSetCreation');
% addCallback(peCbm, 'pulseSetDeletion', {@pj_pulseSetDeletionn, pj}, 'pj_pulseSetDeletion');

pj_setPrograms(pj, {ep, stim, acq});  
registerLoopable(lm, {@pj_loopListener, pj}, 'pulseJacker');

hs = program('hotswitch', 'hotswitch', 'hotswitch', 'hs_config', 'hs_config');
openprogram(progmanager, hs);

%--------------------------------------------------
%Load a configuration (if requested).
waitbar(0.99, wb, 'Loading configuration...');
loadConfigurations(progmanager);

%----- Include following from startup.m to reset digital lines
% to low after MATLAB bug sets them high upon creation of AI/AO objects: 
% % 7/27/06: DHO,TO
% disp('dephys: Setting digital IO lines low ...')
% %TO072406B - Modified to cover all lines on all boards, and be compatible across systems. -- Tim O'Connor 7/24/06
% info = daqhwinfo('nidaq');
% for i = 1 : size(info.ObjectConstructorName, 1)
%     if ~isempty(info.ObjectConstructorName{i, 3})
%         dio = eval(info.ObjectConstructorName{i, 3});
%         %     dio = digitalio('nidaq', i);
%         for j = 0 : 7
%             line = addLine(dio, j, 'out');
%             putvalue(line, 0);
%         end
%         delete(dio);
%     end
% end

fprintf(1, '\nLoading Completed.\n\n');

%Kill the waitbar.
delete(wb);


%dh 5.30.07 ======================hack from tim

%See TO030207A

global triggerRepeatHack;

triggerRepeatHack = 1;

fprintf(1, 'Executing:\n `global triggerRepeatHack;`\n `triggerRepeatHack = 1;`\nSee TO030207A for details.\n');


