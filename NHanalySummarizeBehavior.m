function a = NHanalySummarizeBehavior(varargin)

%%%% This function accepts all Behavior files that are generated from the
%%%% function 'NHanalyLeverPressBehavior', which simply summarizes the
%%%% behavioral results for a specific day. To summary the behavior over
%%%% multiple days, use all of the corresponding files as input arguments
%%%% to this function. **** If you want to assign session numbers to each
%%%% of the files used, make the last input argument an array of said
%%%% numbers. Otherwise, the program will assume sequential order.


if ~isstruct(varargin{end})
    animalname = regexp(inputname(1), '[A-Z]{2,3}0{1,3}[1-9,A-Z]{0,2}', 'match');
    used_sessions = varargin{end};
    ns = length(varargin)-1;
    if length(varargin)-1 > 14
        explength = length(varargin)-1;
    else
        explength = 14;
    end
else
    animalname = regexp(inputname(1), '[A-Z]{2,3}0{1,3}\w{1,2}', 'match');
    used_sessions = 1:length(varargin);
    ns = length(varargin);
end


for i = 1:ns
    File{i} = varargin{i};
    try
        beh_sample_freq(i,1) = length(File{i}.lever_force_smooth)/File{i}.Frame_Times(end);
        actualTime(i,1) = File{i}.Frame_Times(end);
    catch
        beh_sample_freq(i,1) = 1000;
        actualTime(i,1) = length(File{i}.lever_force_smooth)/beh_sample_freq(i,1);
    end
    actualTime(i,1) = actualTime(i,1)/60;   %%% Convert to minutes 
end

%%% Performance %%%

scrsz = get(0, 'ScreenSize');

figure('Position', scrsz); 

rewards = zeros(ns,1);
rxnTime = cell(1,length(used_sessions));
cuestart2reward = cell(1,length(used_sessions));
trial_length = cell(1,length(used_sessions));


for session = 1:ns
    ch = find(strcmp(File{session}.xsg_data.channel_names,'Trial_number'));
    bitcode = parse_behavior_bitcode(File{session}.xsg_data.channels(:,ch));
    if ~isempty(File{session}.Behavior_Frames)
        current_session = used_sessions(session);
        movements_only = File{session}.lever_force_smooth.*File{session}.lever_active;
        boundary_frames = find(diff([Inf; File{session}.lever_active;Inf])~=0);
        for trialnumber = 1:length(File{session}.Behavior_Frames)
            if ~isempty(File{session}.Behavior_Frames{trialnumber}.states.reward)
                rewards(session,1) = rewards(session,1)+1;
                reward_time = round(File{session}.Frame_Times(round(File{session}.Behavior_Frames{trialnumber}.states.reward(1)))*1000);
                cue_start = round(File{session}.Frame_Times(round(File{session}.Behavior_Frames{trialnumber}.states.cue(1)))*1000);
                if trialnumber ~= length(File{session}.Behavior_Frames)         %%% The last behavioral trial must be treated differently, since there is no future cue to use as a reference
                    next_cue = round(File{session}.Frame_Times(round(File{session}.Behavior_Frames{trialnumber+1}.states.cue(1)))*1000);
                else
                    next_cue = round(File{session}.Frame_Times(end))*1000;
                end
                end_trial = next_cue;
                File{session}.movement{rewards(session,1)} = File{session}.lever_force_smooth(cue_start:next_cue);
                File{session}.PastThreshRewTrials{rewards(session,1)} = File{session}.lever_force_smooth(cue_start:next_cue).*File{session}.lever_active(cue_start:next_cue);  %%% Binarizes lever motion for a particular cue period
                
                %%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%
                [File,tlength, cs2r, rxntime, fault] = ProfileRewardedMovements(File, boundary_frames,session, trialnumber, rewards,cue_start, reward_time, end_trial);
                if fault
                    continue
                end
                %%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%
                
                trial_length{session}(rewards(session,1),1) = tlength;
                rxnTime{session}(rewards(session,1),1) = rxntime;
                CuetoRew{session}(rewards(session,1),1) = cs2r;
            else
            end
        end
        AveRxnTime(session,1) = nanmean(rxnTime{session});
        AveCueToRew(session,1) = nanmean(CuetoRew{session});
        trial_length{session}(trial_length{session} == 0) = NaN;
        if ~isempty(trial_length{session})
            range(session,1) = min(trial_length{session});
        else
            range(session,1) = nan;
        end
        subplot(2,ns,ns+session); hold on;
        for rewardedtrial = 1:rewards(session,1) 
            try
                if ~isempty(File{session}.SuccessfulMovements{rewardedtrial})
                    MovementMat{current_session}(rewardedtrial,1:range(session,1)) = File{session}.SuccessfulMovements{rewardedtrial}(1:range(session,1));
                    MovementMat{current_session}(MovementMat{current_session}==0) = nan;
                    plot(MovementMat{current_session}(rewardedtrial,1:3000), 'k');
                    drawnow;
                else
                    MovementMat{current_session}(rewardedtrial,1:3001) = nan(1,3001);
                end
            catch
                MovementMat{current_session}(rewardedtrial,1:3001) = nan(1,3001);
                disp(['Movement was not tracked for trial ', num2str(rewardedtrial), ' from session ', num2str(session)])
            end
        end
        if rewards(session,1) ~= 0
            MovementAve(current_session,:) = nanmean(MovementMat{current_session}(:,1:3000),1);
            plot(MovementAve(current_session,1:3000), 'r', 'Linewidth', 2); drawnow;
            ylim([-2.5 0])
            title(['Session', num2str(current_session)])
        else
            MovementAve(current_session,:) = nan(1,3000);
        end
        trials(session,1) = length(File{session}.Behavior_Frames);
    else  %%% This section is for using data that is not aligned to imaging frames
        current_session = used_sessions(session);
        movements_only = File{session}.lever_force_smooth.*File{session}.lever_active;
        boundary_frames = find(diff([Inf; File{session}.lever_active;Inf])~=0);
        if boundary_frames(1) == 1;
            boundary_frames = boundary_frames(2:end);
        end
        
        bitcode_offset = mode([bitcode.behavior_trial_num]-(1:length(bitcode)));
        for trialnumber = 2:length(File{session}.DispatcherData.saved_history.ProtocolsSection_parsed_events) %%% For every trial(the first trial often has incorrect information)
            if trialnumber > length(bitcode)
                continue
            end
%             i_bitcode=find([bitcode.behavior_trial_num]==trialnumber,1);
%             if(isempty(i_bitcode)) continue; end
            i_bitcode = trialnumber-bitcode_offset;
            if i_bitcode<=0
                continue
            end
            start_trial = round(bitcode(i_bitcode).xsg_sec*1000);
            t0 = varargin{session}.DispatcherData.saved_history.ProtocolsSection_parsed_events{trialnumber}.states.bitcode(1);
            end_trial = start_trial+round((File{session}.DispatcherData.saved_history.ProtocolsSection_parsed_events{trialnumber}.states.state_0(2,1)-t0)*1000); %%% Subtract the starting point of each trial  (t0) on Dispatcher's terms so as to align it with the start point according to the .xsg data

            if ~isempty(File{session}.DispatcherData.saved_history.ProtocolsSection_parsed_events{trialnumber}.states.reward) %%% If the trial was rewarded
                rewards(session,1) = rewards(session,1)+1;
                reward_time = round(start_trial+(File{session}.DispatcherData.saved_history.ProtocolsSection_parsed_events{trialnumber}.states.reward(1)-t0)*1000);
                cue_start = round(start_trial+(File{session}.DispatcherData.saved_history.ProtocolsSection_parsed_events{trialnumber}.states.cue(1)-t0)*1000);
                if cue_start >= length(File{session}.lever_force_smooth) || reward_time >= length(File{session}.lever_force_smooth)
                    continue
                end
                
                if end_trial>length(File{session}.lever_force_smooth)
                    end_trial = length(File{session}.lever_force_smooth);
                end
                
                File{session}.movement{rewards(session,1)} = File{session}.lever_force_smooth(cue_start:end_trial);
                PastThreshRewTrials{rewards(session,1)} = File{session}.lever_force_smooth(cue_start:end_trial).*File{session}.lever_active(cue_start:end_trial);  %%% Binarizes lever motion for a particular cue period
                
                %%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%
                [File,tlength, cs2r, rxntime, fault] = ProfileRewardedMovements(File, boundary_frames,session, trialnumber, rewards,cue_start, reward_time, end_trial);
                if fault
                    continue
                end
                %%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%
                
                trial_length{session}(rewards(session,1),1) = tlength;
                rxnTime{session}(rewards(session,1),1) = rxntime;
                CuetoRew{session}(rewards(session,1),1) = cs2r;
            else
            end
        end
        AveRxnTime(session,1) = nanmean(rxnTime{session});
        AveCueToRew(session,1) = nanmean(CuetoRew{session});
        trial_length{session}(trial_length{session} == 0) = NaN;
        try
            range(session,1) = min(trial_length{session});
        catch
            range(session,1) = nan;
        end
        subplot(2,ns,ns+session); hold on;
        for trialnumber = 1:rewards(session,1) 
            try
                if ~isempty(RecordedSuccessfulMovement{trialnumber})
                    MovementMat{current_session}(trialnumber,1:range(session,1)) = RecordedSuccessfulMovement{trialnumber}(1:range(session,1));
                    MovementMat{current_session}(MovementMat{current_session}==0) = nan;
                    plot(MovementMat{current_session}(trialnumber,1:3000), 'k');
                else
                    MovementMat{current_session}(trialnumber,1:3001) = nan(1,3001);
                    disp(['Movement was not tracked for trial ', num2str(trialnumber), ' from session ', num2str(session)])
                    drawnow;
                end
            catch
                MovementMat{current_session}(trialnumber,1:3001) = nan(1,3001);
                disp(['Movement was not tracked for trial ', num2str(trialnumber), ' from session ', num2str(session)])
            end
        end
        if rewards(session,1) ~= 0
            MovementAve(current_session,:) = nanmean(MovementMat{current_session}(:,1:3000),1);
            plot(MovementAve(current_session,1:3000), 'r', 'Linewidth', 2); drawnow;
            ylim([-2.5 0])
            title(['Session', num2str(current_session)])
        else
            MovementAve(current_session,:) = nan(1,3000);
        end
        trials(session,1) = length(File{session}.DispatcherData.saved_history.ProtocolsSection_parsed_events);
        if rewards(session,1) == 0
            File{session}.SuccessfulMovements = [];
            File{session}.PastThreshRewTrials = [];
        else
            File{session}.SuccessfulMovements = RecordedSuccessfulMovement;
            File{session}.PastThreshRewTrials = PastThreshRewTrials;
        end
    end
end

temp = nan(explength,1); 
temp(used_sessions) = (rewards./trials).*100;
temp(temp == 0) = nan;
rewards = temp;

subplot(2,ns,1:round(ns/4))
plot(1:explength, rewards(1:end,1))
rewards;
title('Correct Trials')
xlabel('Session')
ylabel('Rewards')

temp = nan(explength,1); 
temp(used_sessions) = AveRxnTime;
temp(temp == 0) = nan;
AveRxnTime = temp;

temp = nan(explength,1); 
temp(used_sessions) = AveCueToRew;
temp(temp == 0) = nan;
AveCueToRew = temp;

subplot(2,ns,round(ns/4)+1:round(ns/2))
plot(1:explength, AveRxnTime, 'k'); hold on;
plot(1:explength, AveCueToRew, 'r');
title('Reaction Time')
xlabel('Session')
legend({'Cue to movement', 'Cue to reward'})

subplot(2,ns, round(ns/2)+2:ns);

%%% Concatenate all the movement traces

unused_days = explength-ns;

total = 0;
for session = 1:ns
    total = total+size(MovementMat{session},1);
end

cat_data = nan(3001,total+unused_days);

counter = used_sessions(1); %%% Start from the first day that was actually used, leave preceding days blank
for session = 1:ns
    current_session = used_sessions(session);
    if size(MovementMat{current_session},1)
        cat_data(:,counter:counter+size(MovementMat{current_session},1)-1) = MovementMat{current_session}';
        counter = counter + size(MovementMat{current_session},1);
    else
    end
%     if session < ns
%         if isempty(MovementMat{current_session+1})
%             addon = 1;
%             for trialnumber = current_session+2:ns
%                 if isempty(MovementMat{trialnumber})
%                     addon = addon+1;
%                 end
%             end
%             counter = counter + addon;
%         end
%     end
end

%%% Find the correlation within and between individual movements
[r, p] = corrcoef(cat_data, 'rows', 'pairwise'); 

% [r_lever p_lever] = corrcoef(MovementAve');


%%% Find the median of each block of data correlations

r_lever = nan(explength,explength);

counter1 = 1;
for session = 1:ns
    current_session_row = used_sessions(session);
    temp1 = counter1:counter1+size(MovementMat{current_session_row},1)-1;
    counter2 = counter1; %%% to step down the diagonal, make counter2 start where counter 1 does!
        for trialnumber = session:ns
            current_session_column = used_sessions(trialnumber);
            temp2 = counter2:counter2+size(MovementMat{current_session_column},1)-1;
            r_lever(current_session_row,current_session_column) = nanmean(nanmean(r(temp1,temp2))); 
            r_lever(current_session_column,current_session_row) = nanmean(nanmean(r(temp1,temp2))); %%% Accounts for the symmetry of heatmaps (only half needs to be calculated, the rest can just be filled in, as done here)
            counter2 = counter2+size(MovementMat{current_session_column},1);
        end
    counter1 = counter1 + size(MovementMat{current_session_row},1);
end

imagesc(r_lever);
set(gcf, 'ColorMap', hot)
set(gca, 'CLim', [min(min(r_lever)), max(max(r_lever))])
colorbar
ylabel('Session')
xlabel('Session')
title('Movement correlation over sessions')


a.rewards = rewards;
a.ReactionTime = AveRxnTime;
a.MovementAverages = MovementAve;
a.MovementCorrelation = r_lever;
a.UsedSessions = used_sessions;
a.CuetoReward = AveCueToRew;

try
    cd('C:\Users\Komiyama\Desktop\Output Data');
catch
    cd('C:\Users\komiyama\Desktop\Giulia\All Behavioral Data')
end
eval([animalname{1}, '_SummarizedBehavior = a']);
save([animalname{1}, '_SummarizedBehavior'], [animalname{1}, '_SummarizedBehavior']);

figure; plot(diag(r_lever), 'k');
hold on;
plot(diag(r_lever,1),'Color', [0.6 0.6 0.6])
ylabel('Correlations')
xlabel('Session')
legend({'Within sessions', 'Across sessions'})

