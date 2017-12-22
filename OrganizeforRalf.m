function [Distances, AllActCorrelations, MovementCorrelations, QuiescentCorrelations, AllRates, MovementRates, QuietRates] = OrganizeforRalf(data, correlations, session)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Distances = data.DistanceHeatMap;

numspines = size(Distances,1);

samedendcorrection = Distances; samedendcorrection(~isnan(samedendcorrection)) = 1;

AllActCorrelations = correlations{session}.DendSubtractedSpineCorrelations(10:10+numspines-1, 10:10+numspines-1).*samedendcorrection;

MovementCorrelations = correlations{session}.SpineDuringMovePeriods.*samedendcorrection;
QuiescentCorrelations = correlations{session}.SpineDuringStillPeriods.*samedendcorrection;

expTime = length(data.SynapseOnlyBinarized_DendriteSubtracted)/30.5/60;

AllRates = diff(data.SynapseOnlyBinarized_DendriteSubtracted,1,2);
AllRates(AllRates<0) = 0;
AllRates = sum(AllRates,2)/expTime;

movementMat = correlations{session}.BinarizedBehavior;

if length(movementMat)~=length(data.SynapseOnlyBinarized_DendriteSubtracted)
    if length(movementMat)>length(data.SynapseOnlyBinarized_DendriteSubtracted)
        movementMat = movementMat(1:length(data.SynapseOnlyBinarized_DendriteSubtracted));
    else
        movementMat = [movementMat; zeros(abs(length(data.SynapseOnlyBinarized_DendriteSubtracted)-length(movementMat)),1)];
    end
end

movementMat = repmat(movementMat,1,numspines);

MovementAct = movementMat'.*data.SynapseOnlyBinarized_DendriteSubtracted;

MovementRates = diff(MovementAct,1,2);
MovementRates(MovementRates<0) = 0;
MovementRates = sum(MovementRates,2)/expTime;

QuietAct = ~movementMat'.*data.SynapseOnlyBinarized_DendriteSubtracted;

QuietRates = diff(QuietAct,1,2);
QuietRates(QuietRates<0) = 0;
QuietRates = sum(QuietRates,2)/expTime;

figure; 
sub1 = 5;
sub2 = 3;
subplot(sub1,sub2,1:3)
imagesc(Distances)
subplot(sub1,sub2,4:6)
imagesc(AllActCorrelations)
subplot(sub1,sub2,7:9)
imagesc(MovementCorrelations)
subplot(sub1,sub2,10:12)
imagesc(QuiescentCorrelations)
subplot(sub1,sub2,13)
imagesc(AllRates)
subplot(sub1,sub2,14)
imagesc(MovementRates)
subplot(sub1,sub2,15)
imagesc(QuietRates)


end

