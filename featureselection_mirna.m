function [mirnaTopData, mirnaTopFinal, mirnaTopFeatures, fieldNamestrain, cellid1] = featureselection_mirna(cellid, cellid1, mirnaNorm, fieldNamestrain)
% load the preprocessed mirna data
% [mirnafinal, cellid, cellid1, mirnaNorm, fieldNamescell] = preprocessing_mirna();

% looping through rows, removing rows (genes) that aren't significant 
[rowSize, ~] = size(mirnaNorm);
% cellid = cellstr(cellid); 
% cellid(1) = []; 
pvals = []; 
ind = 1;
while ind <= rowSize
    
    % get the current row, find median 
    currRow = mirnaNorm(ind,:); 
    med = median(currRow); 
    
    % group 1, above the median 
    mask1 = currRow > med; 
    above = currRow(mask1);
    
    % group 2, below the median
    mask2 = currRow < med; 
    below = currRow(mask2);
    
    % compare groups, generate p-value 
    [~,p,~,~] = ttest2(above,below,'Vartype','unequal');
    
    % create p val vector 
    pvals = [pvals p]; 
    
    % delete row if insignificant 
    if p > 0.05 || isnan(p)
        
        mirnaNorm(ind,:) = []; 
        cellid(ind,:) = []; 
        
    else
        ind = ind + 1;
    end
  
    % update row size 
    [rowSize,~] = size(mirnaNorm);
    
end

mirnaTopData = mirnaNorm;
mirnaTopcell = num2cell(mirnaTopData);
mirnaTopFinal = [fieldNamestrain; mirnaTopcell];
label = cellid1(1);
cellid = [label; cellid];
mirnaTopFinal = [cellid mirnaTopFinal];
mirnaTopFeatures = cellid;

end 