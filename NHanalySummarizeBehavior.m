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
cs2r = cell(1,length(used_sessions));
trial_length = cell(1,length(used_sessions));


for i = 1:ns
    ch = find(strcmp(File{i}.xsg_data.channel_names,'Trial_number'));
    bitcode = parse_behavior_bitcode(File{i}.xsg_data.channels(:,ch));
    if ~isempty(File{i}.Behavior_Frames)
        current_session = used_sessions(i);
        movements_only = File{i}.lever_force_smooth.*File{i}.lever_active;
        boundary_frames = find(diff([Inf; File{i}.lever_active;Inf])~=0);
        for j = 1:length(File{i}.Behavior_Frames)
            if ~isempty(File{i}.Behavior_Frames{j}.states.reward)
                rewards(i,1) = rewards(i,1)+1;
                reward_time = round(File{i}.Frame_Times(round(File{i}.Behavior_Frames{j}.states.reward(1)))*1000);
                cue_start = round(File{i}.Frame_Times(round(File{i}.Behavior_Frames{j}.states.cue(1)))*1000);
                if j ~= length(File{i}.Behavior_Frames)         %%% The last behavioral trial must be treated differently, since there is no future cue to use as a reference
                    next_cue = round(File{i}.Frame_Times(round(File{i}.Behavior_Frames{j+1}.states.cue(1)))*1000);
                else
                    next_cue = round(File{i}.Frame_Times(end))*1000;
                end
                File{i}.movement{rewards(i,1)} = File{i}.lever_force_smooth(cue_start:next_cue);
                File{i}.PastThreshRewTrials{rewards(i,1)} = File{i}.lever_force_smooth(cue_start:next_cue).*File{i}.lever_active(cue_start:next_cue);  %%% Binarizes lever motion for a particular cue period
                
                %%% Discard any trials for which the animal is already moving
                try    
                    if File{i}.PastThreshRewTrials{rewards(i,1)}(1) ~= 0        
                        disp('Animal was moving at the beginning of trial ', num2str(i));
                        continue
                    end
                catch
                    continue
                end
                if length(File{i}.movement{rewards(i,1)}) < 1000 || cue_start == next_cue
                    File{i}.MovementInitiation{rewards(i,1)} = [];
                    continue
                end
                
                File{i}.CueStarttoReward{rewards(i,1)} = File{i}.lever_force_smooth(cue_start:reward_time);
                cs2r{i}(rewards(i,1)) = length(File{i}.CueStarttoReward{rewards(i,1)})/1000;
                %%% The following are all ways to determine the beginning of
                %%% a movement... select based on what works!
                    temp = find(boundary_frames < reward_time);
                    if boundary_frames(temp(end))<400
                        baseLine_start = length(1:boundary_frames(temp(end)));
                    else
                        baseLine_start = 400;
                    end
                    successful_mvmt_start = boundary_frames(temp(end))-baseLine_start; %%% Defines the entire contiguous movement that leads to rewards, including failed attempts!
                    if successful_mvmt_start == 0
                        successful_mvmt_start = 1;
                    end
                    rise = find(File{i}.CueStarttoReward{rewards(i,1)}> median(File{i}.CueStarttoReward{rewards(i,1)}));        %%% Only includes successful attempts! Find the last value that is below the threshold (in this case, since lever position is negative, find the last value that is greater that the median)
                if isempty(rise)
                    rise = 1;
                end
    %             File{i}.MovementInitiation{rewards(i,1)} = File{i}.movement{rewards(i,1)}(rise(end):end);
                File{i}.MovementInitiation{rewards(i,1)} = File{i}.lever_force_smooth(successful_mvmt_start:successful_mvmt_start+3000);
                trial_length{i}(rewards(i,1),1) = length(File{i}.MovementInitiation{rewards(i,1)});
                if trial_length{i}(rewards(i,1),1) == 0
                    dbstop
                end
                try
                    rxnTime{i}(rewards(i,1),1) = find(File{i}.PastThreshRewTrials{rewards(i,1)}, 1, 'first')/1000; %%% Reaction time in seconds (Starts at cue and ends at motion start) (note that this data is downsampled from ephus' 10,000Hz to 1,000Hz by Andy's code)
                catch
                    rxnTime{i}(rewards(i,1),1) = NaN;
                end
            else
            end
        end
        AveRxnTime(i,1) = nanmean(rxnTime{i});
        CueToRew(i,1) = nanmean(cs2r{i});
        trial_length{i}(trial_length{i} == 0) = NaN;
        try
            range(i,1) = min(trial_length{i});
        catch
            range(i,1) = nan;
        end
        subplot(2,ns,ns+i); hold on;
        for j = 1:rewards(i,1) 
            try
                if ~isempty(File{i}.MovementInitiation{j})
                    MovementMat{current_session}(j,1:range(i,1)) = File{i}.MovementInitiation{j}(1:range(i,1));
                    MovementMat{current_session}(MovementMat{current_session}==0) = nan;
                    plot(MovementMat{current_session}(j,1:3000), 'k');
                else
                    MovementMat{current_session}(j,1:3001) = nan(1,3001);
                end
            catch
                MovementMat{current_session}(j,1:3001) = nan(1,3001);
                disp(['Movement was not tracked for trial ', num2str(j), ' from session ', num2str(i)])
            end
        end
        if rewards(i,1) ~= 0
            MovementAve(current_session,:) = nanmean(MovementMat{current_session}(:,1:3000),1);
            plot(MovementAve(current_session,1:3000), 'r', 'Linewidth', 2);
            ylim([-2.5 0])
            title(['Session', num2str(current_session)])
        else
            MovementAve(current_session,:) = nan(1,3000);
        end
        trials(i,1) = length(File{i}.Behavior_Frames);
    else  %%% This section is for using data that is not aligned to imaging frames
        current_session = used_sessions(i);
        movements_only = File{i}.lever_force_smooth.*File{i}.lever_active;
        boundary_frames = find(diff([Inf; File{i}.lever_active;Inf])~=0);
        if boundary_frames(1) == 1;
            boundary_frames = boundary_frames(2:end);
        end
        
        bitcode_offset = mode([bitcode.behavior_trial_num]-(1:length(bitcode)));
        for j = 2:length(File{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events) %%% For every trial(the first trial often has incorrect information)
            if j > length(bitcode)
                continue
            end
%             i_bitcode=find([bitcode.behavior_trial_num]==j,1);
%             if(isempty(i_bitcode)) continue; end
            i_bitcode = j-bitcode_offset;
            if i_bitcode<=0
                continue
            end
            start_trial = round(bitcode(i_bitcode).xsg_sec*1000);
            t0 = varargin{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events{j}.states.bitcode(1);
            end_trial = start_trial+round((File{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events{j}.states.state_0(2,1)-t0)*1000); %%% Subtract the starting point of each trial  (t0) on Dispatcher's terms so as to align it with the start point according to the .xsg data

            if ~isempty(File{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events{j}.states.reward) %%% If the trial was rewarded
                rewards(i,1) = rewards(i,1)+1;
                reward_time = round(start_trial+(File{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events{j}.states.reward(1)-t0)*1000);
                cue_start = round(start_trial+(File{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events{j}.states.cue(1)-t0)*1000);
                if cue_start >= length(File{i}.lever_force_smooth) || reward_time >= length(File{i}.lever_force_smooth)
                    continue
                end
                
                if end_trial>length(File{i}.lever_force_smooth)
                    end_trial = length(File{i}.lever_force_smooth);
                end
                
                File{i}.movement{rewards(i,1)} = File{i}.lever_force_smooth(cue_start:end_trial);
                PastThreshRewTrials{rewards(i,1)} = File{i}.lever_force_smooth(cue_start:end_trial).*File{i}.lever_active(cue_start:end_trial);  %%% Binarizes lever motion for a particular cue period
                
                %% Discard any trials for which the animal is already moving
%                 try    
%                     if PastThreshRewTrials{rewards(i,1)}(1) ~= 0        
%                         disp(['Animal was moving at the beginning of trial ', num2str(i)]);
%                         continue
%                     end
%                 catch
%                     continue
%                 end
                %%
                
                if length(File{i}.movement{rewards(i,1)}) < 1000 
                    RecordedSuccessfulMovement{rewards(i,1)} = [];
                    continue
                end
                
                File{i}.CueStarttoReward{rewards(i,1)} = File{i}.lever_force_smooth(cue_start:reward_time);
                cs2r{i}(rewards(i,1)) = length(File{i}.CueStarttoReward{rewards(i,1)})/1000;
                
                %%% Find the beginning of the rewarded movement
                %%%
                    temp = find(boundary_frames < reward_time); %%% Finds the movement periods that come before the reward
                    if isempty(temp)
                        continue
                    end
                    if boundary_frames(temp(end))<400
                        baseLine_start = length(1:boundary_frames(temp(end)));  %%% finds the index of the end of the last movement period before the reward
                    else
                        baseLine_start = 500;
                    end
                    successful_mvmt_start = boundary_frames(temp(end))-baseLine_start; %%% Defines the entire contiguous movement that leads to rewards, including failed attempts!
                    if successful_mvmt_start == 0
                        successful_mvmt_start = 1;
                    end
                    rise = find(File{i}.CueStarttoReward{rewards(i,1)}> median(File{i}.CueStarttoReward{rewards(i,1)}));        %%% Only includes successful attempts! Find the last value that is below the threshold (in this case, since lever position is negative, find the last value that is greater that the median)
                if isempty(rise)
                    rise = 1;
                end
    %             MovementInitiation{rewards(i,1)} = File{i}.movement{rewards(i,1)}(rise(end):end);
                if baseLine_start<400
                    RecordedSuccessfulMovement{rewards(i,1)} = [ones(1,400-baseLine_start), File{i}.lever_force_smooth(successful_mvmt_start:successful_mvmt_start+3000)];
                else
                    RecordedSuccessfulMovement{rewards(i,1)} = File{i}.lever_force_smooth(successful_mvmt_start:successful_mvmt_start+3000);
                end
                trial_length{i}(rewards(i,1),1) = length(RecordedSuccessfulMovement{rewards(i,1)});
                if trial_length{i}(rewards(i,1),1) == 0
                    dbstop
                end
                try
                    rxnTime{i}(rewards(i,1),1) = find(PastThreshRewTrials{rewards(i,1)}, 1, 'first')/1000; %%% Reaction time in seconds (Starts at cue and ends at motion start) (note that this data is downsampled from ephus' 10,000Hz to 1,000Hz by Andy's code)
                catch
                    rxnTime{i}(rewards(i,1),1) = NaN;
                end
            else
            end
        end
        AveRxnTime(i,1) = nanmean(rxnTime{i});
        CueToRew(i,1) = nanmean(cs2r{i});
        trial_length{i}(trial_length{i} == 0) = NaN;
        try
            range(i,1) = min(trial_length{i});
        catch
            range(i,1) = nan;
        end
        subplot(2,ns,ns+i); hold on;
        for j = 1:rewards(i,1) 
            try
                if ~isempty(RecordedSuccessfulMovement{j})
                    MovementMat{current_session}(j,1:range(i,1)) = RecordedSuccessfulMovement{j}(1:range(i,1));
                    MovementMat{current_session}(MovementMat{current_session}==0) = nan;
                    plot(MovementMat{current_session}(j,1:3000), 'k');
                else
                    MovementMat{current_session}(j,1:3001) = nan(1,3001);
                    disp(['Movement was not tracked for trial ', num2str(j), ' from session ', num2str(i)])
                end
            catch
                MovementMat{current_session}(j,1:3001) = nan(1,3001);
                disp(['Movement was not tracked for trial ', num2str(j), ' from session ', num2str(i)])
            end
        end
        if rewards(i,1) ~= 0
            MovementAve(current_session,:) = nanmean(MovementMat{current_session}(:,1:3000),1);
            plot(MovementAve(current_session,1:3000), 'r', 'Linewidth', 2);
            ylim([-2.5 0])
            title(['Session', num2str(current_session)])
        else
            MovementAve(current_session,:) = nan(1,3000);
        end
        trials(i,1) = length(File{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events);
        if rewards(i,1) == 0
            File{i}.MovementInitiation = [];
            File{i}.PastThreshRewTrials = [];
        else
            File{i}.MovementInitiation = RecordedSuccessfulMovement;
            File{i}.PastThreshRewTrials = PastThreshRewTrials;
        end
    end
end

temp = nan(ns,1); 
temp(used_sessions) = (rewards./trials).*100;
temp(temp == 0) = nan;
rewards = temp;

subplot(2,ns,1:round(ns/4))
plot(1:ns, rewards(1:end,1))
rewards;
title('Correct Trials')
xlabel('Session')
ylabel('Rewards')

temp = nan(ns,1); 
temp(used_sessions) = AveRxnTime;
temp(temp == 0) = nan;
AveRxnTime = temp;

temp = nan(ns,1); 
temp(used_sessions) = CueToRew;
temp(temp == 0) = nan;
CueToRew = temp;

subplot(2,ns,round(ns/4)+1:round(ns/2))
plot(1:ns, AveRxnTime, 'k'); hold on;
plot(1:ns, CueToRew, 'r');
title('Reaction Time')
xlabel('Session')
legend({'Cue to movement', 'Cue to reward'})

subplot(2,ns, round(ns/2)+2:ns);

%%% Concatenate all the movement traces

unused_days = 14-ns;

total = 0;
for i = 1:ns
    total = total+size(MovementMat{i},1);
end

cat_data = nan(3001,total+unused_days);

counter = used_sessions(1); %%% Start from the first day that was actually used, leave preceding days blank
for i = 1:ns
    current_session = used_sessions(i);
    if size(MovementMat{current_session},1)
        cat_data(:,counter:counter+size(MovementMat{current_session},1)-1) = MovementMat{current_session}';
        counter = counter + size(MovementMat{current_session},1);
    else
    end
%     if i < ns
%         if isempty(MovementMat{current_session+1})
%             addon = 1;
%             for j = current_session+2:ns
%                 if isempty(MovementMat{j})
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

r_lever = nan(14,14);

counter1 = 1;
for i = 1:ns
    current_session_row = used_sessions(i);
    temp1 = counter1:counter1+size(MovementMat{current_session_row},1)-1;
    counter2 = counter1; %%% to step down the diagonal, make counter2 start where counter 1 does!
        for j = i:ns
            current_session_column = used_sessions(j);
            temp2 = counter2:counter2+size(MovementMat{current_session_column},1)-1;
            r_lever(current_session_row,current_session_column) = nanmedian(nanmedian(r(temp1,temp2))); 
            r_lever(current_session_column,current_session_row) = nanmedian(nanmedian(r(temp1,temp2))); %%% Accounts for the symmetry of heatmaps (only half needs to be calculated, the rest can just be filled in, as done here)
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
a.CuetoReward = CueToRew;

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

