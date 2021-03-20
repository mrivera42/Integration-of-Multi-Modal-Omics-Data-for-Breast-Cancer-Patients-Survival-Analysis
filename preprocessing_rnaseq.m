function [mirnafinal, cellid, cellid1, mirnaNorm, fieldNamescell] = preprocessing_rnaseq()

% open file 
mirna = tdfread('rnaseq-aaa.tsv', '\t');

% get an array of the fieldnames 
fieldNames = char(fieldnames(mirna));
featuresFieldName = fieldNames(1,:);
featuresFieldName(featuresFieldName == ' ') = [];
id = mirna(:).(featuresFieldName) ; 
cellid = cellstr(id);


% use fieldnames to loop through struture and take out double values
fieldNames_length = length(fieldNames);
mirnaNew = [];
for x = 2:fieldNames_length % 1203 = length of fieldNames from workspace
    
    currField = fieldNames(x,:);
  
    currSample = [mirna(:).(currField)];
    
    mirnaNew = [mirnaNew currSample];
    
    
end

% looping through rows, remove rows containing all zeros
 %remove those row labels from feature names also
[rowSize, col] = size(mirnaNew);

ind = 1;
while ind <= rowSize
    mask = mirnaNew(ind, :) == 0;
    
    if sum(mask) == col
        mirnaNew(ind,:) = [];
        cellid(ind) = [];
    else
        ind = ind+1;
    end
    [rowSize,~] = size(mirnaNew);
end

 %remove rows with any empty values
[rowSize, ~] = size(mirnaNew);

ind = 1;
while ind <= rowSize
    if any(isnan(mirnaNew(ind,:))) 
        mirnaNew(ind,:) = [];
        cellid(ind) = [];
    else
        ind = ind+1;
    end
    [rowSize,~] = size(mirnaNew);
end

% remove outliers    
[mirnaOut,outlierindex] = rmoutliers(mirnaNew, 'mean');  
cellid(outlierindex) = [];

% change range to [0,1]

mirnaNorm = rescale(mirnaOut);    
   
%remove repeated samples



%fix labeling of sample ID
fieldNamescell = cellstr(fieldNames);
for x = 2:fieldNames_length
    field_str = fieldNamescell{x};
    field_str = replace(field_str,'0x2D','-');
    fieldNamescell(x) = cellstr(field_str);
end
[~, ~] = unique(fieldNamescell,'stable');


%make the double array a cell array and concatenate accordingly
mirnaNormcell = num2cell(mirnaNorm);
cellid1 = [fieldNamescell(1);cellid];
fieldNamescell(1) = [];
mirnafinal = [fieldNamescell'; mirnaNormcell];
mirnafinal = [cellid1 mirnafinal];
end

