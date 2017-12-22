function AverageBehavior(varargin)

rewards = nan(length(varargin), 14);
ReactionTime = nan(length(varargin),14);
CuetoReward = nan(length(varargin),14);
MovementCorrelation = nan(14,14,length(varargin));

%%%%%%%%%%%%%%%%%%%%%%
%%Color Information%%%
%%%%%%%%%%%%%%%%%%%%%%

    lgray = [0.50 0.51 0.52];   brown = [0.28 0.22 0.14];
    gray = [0.50 0.51 0.52];    lbrown = [0.59 0.45 0.28];
    yellow = [1.00 0.76 0.05];  orange = [0.95 0.40 0.13];
    lgreen = [0.55 0.78 0.25];  green = [0.00 0.43 0.23];
    lblue = [0.00 0.68 0.94];   blue = [0.00 0.33 0.65];
    magenta = [0.93 0.22 0.55]; purple = [0.57 0.15 0.56];
    red = [0.93 0.11 0.14];     black = [0 0 0];
    colorj = {red,lblue,green,lgreen,gray,brown,yellow,blue,purple,magenta,orange,brown,lbrown};
    
%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:length(varargin)
    rewards(i,varargin{i}.UsedSessions) = varargin{i}.rewards(varargin{i}.UsedSessions);
    ReactionTime(i,varargin{i}.UsedSessions) = varargin{i}.ReactionTime(varargin{i}.UsedSessions);
    CuetoReward(i,varargin{i}.UsedSessions) = varargin{i}.CuetoReward(varargin{i}.UsedSessions);
%     MovementAverages(i,varargin{i}.UsedSessions) = 
    MovementCorrelation(:,:,i) = varargin{i}.MovementCorrelation(1:14,1:14);
end

for i = 1:14
    rewardsSEM(1,i) = nanstd(rewards(:,i),0,1)/sqrt(length(~isnan(rewards(:,i))));
    RTSEM(1,i) = nanstd(ReactionTime(:,i),0,1)/sqrt(length(~isnan(ReactionTime(:,i))));
    CtRSEM(1,i) = nanstd(CuetoReward(:,i),0,1)/sqrt(length(~isnan(CuetoReward(:,i))));
end

scrsz = get(0, 'ScreenSize');
figure('Position', scrsz);

subnum = ceil(sqrt(length(varargin)));

for i = 1:length(varargin)
    ax(i) = subplot(subnum,subnum,i); hold on;
    plot(varargin{i}.UsedSessions,varargin{i}.rewards(varargin{i}.UsedSessions)/100, 'Color', black, 'Linewidth', 2)
    plot(varargin{i}.UsedSessions,varargin{i}.ReactionTime(varargin{i}.UsedSessions)./nanmax(varargin{i}.ReactionTime(varargin{i}.UsedSessions)), 'Color',red, 'Linewidth', 2)
    plot(varargin{i}.UsedSessions,varargin{i}.CuetoReward(varargin{i}.UsedSessions)./nanmax(varargin{i}.CuetoReward(varargin{i}.UsedSessions)), 'Color', lblue, 'Linewidth', 2)
    withinsessions = diag(MovementCorrelation(1:14,1:14,i));
    plot(1:14,withinsessions./nanmax(withinsessions), 'Color', green, 'Linewidth', 2)
    acrosssessions = diag(MovementCorrelation(1:14,1:14, i),+1);
    plot(2:14,acrosssessions./nanmax(acrosssessions), 'Color', lgreen, 'Linewidth', 2)
    set(gca, 'XTick', 1:14)
%     ylim([-0.05 1.05])
    ylabel('Percent Max')
    xlabel('Session')
    file = regexp(inputname(i), '[ABCDEFGHIJKLMNOPQRSTUVWXYZ]{2,3}0+[A-Z,0-9]*', 'match');
    title(file{1});
    plot(1:14, ones(1,14),'--k')
end

linkaxes(ax);

legend({'Percent Rewarded', 'Reaction Time', 'Cue to Reward', 'Mov Corr within', 'Mov Corr Across'},'Location','BestOutside')


scrsz = get(0, 'ScreenSize');
figure('Position', scrsz);

subplot(6,6,[1:3, 7:9]); plot(nanmean(rewards,1), 'b'); xlim([0 15]); ylabel('% Rewarded'); xlabel('Session');
r_errorbar(1:14, nanmean(rewards, 1), rewardsSEM, 'b');
set(gca, 'XTick', 1:14)
subplot(6,6,[4:6, 10:12]); reactplot = plot(nanmean(ReactionTime,1), 'k'); xlim([0 15]); ylabel('Time (s)'); xlabel('Session');
hold on; c2rplot = plot(nanmean(CuetoReward,1), 'r');
r_errorbar(1:14, nanmean(ReactionTime,1), RTSEM, 'k'); 
r_errorbar(1:14, nanmean(CuetoReward,1), CtRSEM, 'r');
legend([reactplot, c2rplot], {'Reaction Time', 'Cue to Reward'});
set(gca, 'XTick', 1:14, 'YTick', 1:14)
subplot(6,6,[13:15, 19:21, 25:27, 31:33]); imagesc(nanmean(MovementCorrelation, 3))
set(gcf, 'ColorMap', hot); colorbar;
ylabel('Session', 'Fontsize', 14)
xlabel('Session', 'Fontsize', 14)
set(gca, 'XTick', 1:14, 'YTick', 1:14)
title('Movement correlation over sessions')

for i = 1:size(MovementCorrelation,3)
    within(i,1:14) = diag(MovementCorrelation(:,:,i));
end
for i = 1:size(MovementCorrelation,3)
    across(i,1:13) = diag(MovementCorrelation(:,:,i),1);
end

subplot(6,6,[16:18, 22:24, 28:30, 34:36]);
withinplot = plot(1:14,nanmean(within,1),'k', 'Linewidth', 2); hold on;
acrossplot = plot(2:14, nanmean(across,1), 'Color', [0.5 0.5 0.5], 'Linewidth', 2);
for i = 1:14
    within_SEM(1,i) = nanstd(within(:,i))/sqrt(sum(~isnan(within(i))));
end
for i = 1:13
    across_SEM(1,i) = nanstd(across(:,i))/sqrt(sum(~isnan(across(i))));
end
r_errorbar(1:14, nanmean(within,1), within_SEM, 'k')
r_errorbar(2:14, nanmean(across,1), across_SEM, [0.5 0.5 0.5])
legend([withinplot, acrossplot], {'Within sessions', 'Across sessions'})
ylabel('Correlation')
xlabel('Session')
xlim([0 15])
set(gca, 'XTick', 1:14)
