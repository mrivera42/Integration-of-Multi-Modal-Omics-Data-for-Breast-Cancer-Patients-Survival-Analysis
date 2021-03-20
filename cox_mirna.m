
function [B, survivalFinal, riskscore, mirnaTest] = cox_mirna(mirnaTopData, mirnaTopFeatures, mirnaTest, fieldNamestrain, cellid1)
[~, os, ostime, patientid, sampleid] = preprocessing_survival();
%implement coxphfit model for mirna

%X and Y need same number of rows
%transpose mirnaTop so columns (number of samples) become rows 
    %1202x3 array
    
X = mirnaTopData';
%survival data where rows and number of samples
    %create array of sampleid and os where sampleid is row and os is column
    %1202x1 array
    %we have to make sure it is same order as sample IDs in X
    
samplenames = fieldNamestrain;

%loop through each sample id and if it is in mirna then keep it in
%survival, otherwise remove it 
    %compare samplenames (mirna) to sampleid (survival)
    %reorder samplenames and os so that it is in same order as mirna
    index = ismember(sampleid,samplenames); 
    mirnasampleid = sampleid(index);
    ind1 = ismember(samplenames, mirnasampleid);
    mirnaTopData =  mirnaTopData(:,ind1);
    mirnaTopFeatures = cellstr(mirnaTopFeatures);
    samplenames = samplenames(ind1);
    patientmirna = patientid(index);
    osmirna =os(index);
    ostimemirna = ostime(index);
    %output final array for data 1179x4
   mirnaTop = [samplenames; num2cell(mirnaTopData)];
   mirnaFeatures = mirnaTopFeatures;
   mirnaTop = [mirnaFeatures mirnaTop];
   survivalFinal = [samplenames' osmirna patientmirna ostimemirna];

%Cox Proportional Hazard Model
    %X is nxp where n is number of samples and p is number of features
    %T survival time of nx1 where n is number of samples is ostime
    %B should give coefficient for each feature 
X = mirnaTopData';
T = cell2mat(ostimemirna)
osmirna = logical(cell2mat(osmirna));
B = coxphfit(X, T,'Censoring',osmirna);

%Use coefficients in Cox function to generate risk score per patient per
%modality
%then average modalities to get single risk score per patient
%need to pull in test data, then need to pull features that were chosen 

%edit test data for only features we want
mirnaTest = [cellid1 mirnaTest];
mirnaTestFeatures = mirnaTest(2:end,1);
mirnaTopFeatures = mirnaTopFeatures(2:end);
[~, indices] = ismember(mirnaTopFeatures,mirnaTestFeatures);
%indices should output the indices in mirnaTestFeatures we want to keep

%apply feature mask to test data
mirnaTest = mirnaTest(indices,:);
mirnaTest1 = mirnaTest(:,2:end);
[~, col] = size(mirnaTest1);

%get risk score
riskscore = [];
mirnaTest1 = cell2mat(mirnaTest1);

for x = 1:col
    sample = mirnaTest1(:,x);
    %multiple features and B value
    XB = sample' * B;
    XB = exp(XB);
    riskscore = [riskscore XB];
end
end