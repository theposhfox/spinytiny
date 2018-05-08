function NewSpineAnalysis(experimentnames)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


if strcmpi(getenv('computername'), 'Nathan-Lab-PC')
    cd(['C:\Users\Komiyama\Desktop\Output Data', filesep, experimentnames, ' New Spine Analysis'])
end

fieldsource = fastdir(cd, 'Field');

for i = 1:length(fieldsource)
    load(fieldsource{i})
    fieldnumber = regexp(fieldsource{i}, '\d+.Spine');
    eval(['FieldData{', fieldsource{i}(fieldnumber), '} = SpineRegistry']);
    clear SpineRegistry
end

for i = 1:length(FieldData)
    FieldChanges{i} = diff(FieldData{i}.Data,1,2);
end

if strcmpi(getenv('computername'), 'Nathan-Lab-PC')
    cd(['C:\Users\Komiyama\Desktop\ActivitySummary_UsingRawData'])
end

activitydata = fastdir(cd, [experimentnames, '.+_Summary']);

for i = 1:length(activitydata)
    load(activitydata{i})
end

wrkspc = whos;
for i = 1:length(FieldData)
    for j = 1:length(FieldData{i}.DatesAcquired)
        locate =(regexp(who, FieldData{i}.DatesAcquired{j}));
        FieldData{i}.CalciumData{j} = eval(wrkspc(~cell2mat(cellfun(@isempty, locate, 'uni',false))).name);
    end
end



end

