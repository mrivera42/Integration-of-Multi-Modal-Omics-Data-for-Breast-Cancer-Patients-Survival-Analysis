function [survivalFinal, os, ostime, patientid, sampleid] = preprocessing_survival()
% open file 
survival = tdfread('survival.tsv', '\t')

% loop through rows, delete any repeated samples

% separate data into 4 separate vectors 
sampleid = cellstr(survival(:).sample);
os = num2cell(survival(:).OS); 
patientid = cellstr(survival(:).x_PATIENT); 
ostime = num2cell(survival(:).OS0x2Etime); 

% use unique() on sample ID 
[~, ind] = unique(sampleid,'stable');

% use unique mask to other columns 
os = os(ind);
patientid = patientid(ind);
ostime = ostime(ind);

% concatenate back together 
survivalFinal = [sampleid os patientid ostime]; 
end
