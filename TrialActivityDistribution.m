function [ output_args ] = TrialActivityDistribution(TrialInformation)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if length(TrialInformation)<14
    TrialInformation{14} = [];
end

for i = 1:14
    if ~isempty(TrialInformation{i})
        availablesessiondata = TrialInformation{5}(~cell2mat(cellfun(@(x) isempty(x), TrialInformation{5}, 'UniformOutput', false)));
        cuelength(i) = max(cell2mat(cellfun(@(x) length(x.trialactivity(:,1:x.CueEnd)), availablesessiondata, 'UniformOutput', false)));
        sessioncuedistribution = zeros(1,cuelength(i));
        for j = 1:length(availablesessiondata)
            trialcuelength = availablesessiondata{j}.CueEnd;
            sessioncuedistribution(1:trialcuelength) = sessioncuedistribution(1:trialcuelength) + fliplr(availablesessiondata{j}.trialbinaryactivity(1:trialcuelength));
        end
    end
end

figure; hist(sessioncuedistribution)



