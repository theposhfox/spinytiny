t0=varargin{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events{1}.states.bitcode(1);
figure;hold on;
for z = 1:99
plot((varargin{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events{z}.states.bitcode-t0)*1000, -1.5, 'ok', 'MarkerFaceColor', 'k')
plot((varargin{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events{z}.states.cue-t0)*1000, -1.5, 'vk', 'MarkerFaceColor', 'g')
try
    plot((varargin{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events{z}.states.reward-t0)*1000, -1.5, 'ok', 'MarkerFaceColor', 'c')
catch
    plot((varargin{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events{z}.states.cue-t0)*1000, -1.5, 'ok', 'MarkerFaceColor', 'r')
end
end;
plot(varargin{i}.xsg_data.channels(1:10:end,3))
plot(varargin{i}.xsg_data.channels(1:10:end,2))
title('Dispatcher Data')

bitcode = parse_behavior_bitcode(varargin{i}.xsg_data.channels(:,3));

for z=1:100
    dispatcher_bitcode(z)=varargin{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events{z}.states.bitcode(1);
end
% dispatcher_bitcode(100)=NaN
figure;plot([bitcode(1:100).xsg_sec],dispatcher_bitcode(1:100),'.')
figure;plot([bitcode(1:100).xsg_sec]-dispatcher_bitcode(1:100),'.')

figure; hold on; 

for z = 1:100
start = bitcode(z).xsg_sec;
t0=varargin{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events{z}.states.bitcode(1);
plot((start+varargin{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events{z}.states.bitcode-t0)*1000, -1.5, 'ok', 'MarkerFaceColor', 'k')
plot((start+varargin{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events{z}.states.cue-t0)*1000, -1.5, 'vk', 'MarkerFaceColor', 'g')
try
    plot((start+varargin{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events{z}.states.reward-t0)*1000, -1.5, 'ok', 'MarkerFaceColor', 'c')
catch
    plot((start+varargin{i}.DispatcherData.saved_history.ProtocolsSection_parsed_events{z}.states.cue-t0)*1000, -1.5, 'ok', 'MarkerFaceColor', 'r')
end
end;

plot(varargin{i}.xsg_data.channels(1:10:end,3))
plot(varargin{i}.xsg_data.channels(1:10:end,2))
title('Using bitcode')
