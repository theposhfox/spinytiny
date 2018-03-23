function [File,trial_length,cs2r, rxnTime, fault] = ProfileRewardedMovements(File, boundary_frames,session, trialnumber, rewards,cue_start, reward_time, end_trial)


%%% Discard any trials for which the animal is already moving 
% if File{session}.PastThreshRewTrials{rewards(session,1)}(1) ~= 0        
%     disp(['Animal was moving at the beginning of trial ', num2str(trialnumber), ' from session ', num2str(session)]);
%     File{session}.SuccessfulMovements{rewards(session,1)} = [];
%     trial_length = [];
%     rxnTime = [];
%     fault = 1;
%     return
% end
% if length(File{session}.movement{rewards(session,1)}) < 1000 || cue_start == end_trial
%     File{session}.SuccessfulMovements{rewards(session,1)} = [];
%     trial_length = [];
%     rxnTime = [];
% 	return
% end

File{session}.CueStarttoReward{rewards(session,1)} = File{session}.lever_force_smooth(cue_start:reward_time);
cs2r = length(File{session}.CueStarttoReward{rewards(session,1)})/1000;
%%% The following are all ways to determine the beginning of
%%% a movement... select based on what works!
    temp = find(boundary_frames < reward_time);
    if boundary_frames(temp(end))<400
        baseLine_start = length(1:boundary_frames(temp(end)));
    else
        baseLine_start = 400;
    end
    successful_mvmt_start = boundary_frames(temp(end))-baseLine_start; %%% Defines the entire contiguous movement that leads to rewards, ****including failed attempts!**** (this means that
    if successful_mvmt_start == 0
        successful_mvmt_start = 1;
    end
    rise = find(File{session}.CueStarttoReward{rewards(session,1)}> median(File{session}.CueStarttoReward{rewards(session,1)}));        %%% Only includes successful attempts! Find the last value that is below the threshold (in this case, since lever position is negative, find the last value that is greater that the median)

if isempty(rise)
    rise = 1;
end
%           
if baseLine_start<400
    shift = abs(baseLine_start-400);
else
    shift = 0;
end
File{session}.SuccessfulMovements{rewards(session,1)} = [nan(shift,1);File{session}.lever_force_smooth(successful_mvmt_start:successful_mvmt_start+(3000-shift))];
trial_length = length(File{session}.SuccessfulMovements{rewards(session,1)});
if trial_length == 0
    dbstop
end
try
    rxnTime = find(File{session}.PastThreshRewTrials{rewards(session,1)}, 1, 'first')/1000; %%% Reaction time in seconds (Starts at cue and ends at motion start) (note that this data is downsampled from ephus' 10,000Hz to 1,000Hz by Andy's code)
catch
    rxnTime = NaN;
end

if isempty(rxnTime)
    rxnTime = 0;
end

fault = 0;
