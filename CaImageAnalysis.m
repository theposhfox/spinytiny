function CaImageAnalysis(file);

%%% pixelperMicronat20xzoom = 10.666;

Total_ROIs = size(file.raw)-3;

nearby_spine_ROIs = file.spines(2:end);
Dendrite_ROIs = file.dendrite;

if nearby_spine_ROIs ~= 4
    for i = 1:numel(nearby_spine_ROIs)
        nearby_spine(:,i) = (file.raw(:,nearby_spine_ROIs(i)+3)/mean(file.raw(1:32,nearby_spine_ROIs(i)+3)))./file.red(:,nearby_spine_ROIs(i)+3);
    end
end

for i = 1:numel(Dendrite_ROIs)
    dendrite(:,i) = (file.raw(:,Dendrite_ROIs(i)+3)/mean(file.raw(1:32,Dendrite_ROIs(i)+3)))./file.red(:,Dendrite_ROIs(i)+3);
end

Stim_Spine_Delta = (file.raw(:,4)/mean(file.raw(1:32,4)))./file.red(:,4);

Time = [0:0.125:70];
Time = Time-4;
Time = Time(1:560);

figure(1); plot(Time,Stim_Spine_Delta, 'g'); hold on;

for i = 1:numel(nearby_spine_ROIs);
    plot(Time,nearby_spine(:,i), 'k'); hold on;
end

subplot1 = 1
subplot2 = numel(Dendrite_ROIs)

figure(2);
for i= 1:numel(Dendrite_ROIs)
    h(i) = subplot(subplot1, subplot2, i);
        plot(Time,dendrite(:,i), 'k');
        xlim([-35 75]);
        y_scale(i,:) = get(gca,'YLim');
end

y_scale = max(y_scale);
for i = 1:numel(Dendrite_ROIs);
    set(h(i), 'YLim', [0, y_scale(2)]);
end


