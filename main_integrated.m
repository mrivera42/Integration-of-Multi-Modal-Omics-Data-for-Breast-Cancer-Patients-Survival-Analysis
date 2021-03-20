

% pull sample ids from survival data 
[~, ~, ~, ~, sampleid] = preprocessing_survival();

% get indices for training and testing 
cvIndices = crossvalind('kfold',sampleid,5);
chosenTestActual = sampleid(cvIndices == 1);

% initialize vectors for cross validation and external validation
% performance
pValues = []; 
ciCV = []; 
ciEV = []; 
cv = [1 2 3 4]; 
% create a crossval loop 
for fold = 2:5

    
    % current Testing and Training indices for this loop 
    chosenTest = sampleid(cvIndices == fold);
    chosenTrain = sampleid(cvIndices ~= fold & cvIndices ~= 1);
    
    mirnaTrain = []; 
    mirnaTest = [];
    mirnaTestActual = []; 

% run mirna preprocessing with no inputs 
[mirnafinal, cellid, cellid1, mirnaNorm, fieldNamescell] = preprocessing_mirna();   



% filter through mirna samples and separate into train and test 
for x = 1:length(fieldNamescell)
    
    name = fieldNamescell(x); 
    name2 = chosenTrain(1);
    
    if ismember(name, chosenTrain) == 1
        
        
        chosenCol = num2cell(mirnaNorm(:,x));
        chosenCol = [fieldNamescell(x); chosenCol];
        mirnaTrain = [mirnaTrain chosenCol];
       
        
    elseif ismember(name, chosenTest) == 1
        
        chosenCol = num2cell(mirnaNorm(:,x));
        chosenCol = [fieldNamescell(x); chosenCol];
        mirnaTest = [mirnaTest chosenCol];
        
    elseif ismember(name, chosenTestActual) == 1
        chosenCol = num2cell(mirnaNorm(:,x)); 
        chosenCol = [fieldNamescell(x); chosenCol]; 
        mirnaTestActual = [mirnaTestActual chosenCol]; 
        
    end
   
end

% split mirna into test, train, and validation
fieldNamestrain = mirnaTrain(1,:);
fieldNamestest = mirnaTest(1,:);
fieldNamestestActual = mirnaTestActual(1,:); 
mirnaTrainData = cell2mat(mirnaTrain(2:end,:));
mirnaTestData = cell2mat(mirnaTest(2:end,:)); 
mirnaTestActual = cell2mat(mirnaTestActual(2:end,:));

% call FS, with inputs training 
[mirnaTopData, mirnaTopFinal, mirnaTopFeatures, fieldNamestrain, cellid1] = featureselection_mirna(cellid, cellid1, mirnaTrainData, fieldNamestrain);

% call cox with inputs as training and testing
[B, survivalFinal, riskscore, mirnaTest] = cox_mirna(mirnaTopData, mirnaTopFeatures, mirnaTest, fieldNamestrain, cellid1);
B_mirna = B; 
survivalFinal_mirna = survivalFinal; 
riskscore_mirna = riskscore; 
mirnaTest_mirna = mirnaTest; 
fieldNamestest_mirna = fieldNamestest;

% RNASEQ 
 % preprocessing 
[mirnafinal, cellid, cellid1, mirnaNorm, fieldNamescell] = preprocessing_rnaseq();

mirnaTrain = []; 
mirnaTest = [];
% filter through mirna samples and separate into train and test 
for x = 1:length(fieldNamescell)
    
    name = fieldNamescell(x); 
    name2 = chosenTrain(1);
    
    if ismember(name, chosenTrain) == 1
        
        
        chosenCol = num2cell(mirnaNorm(:,x));
        chosenCol = [fieldNamescell(x); chosenCol];
        mirnaTrain = [mirnaTrain chosenCol];
       
        
    elseif ismember(name, chosenTest) == 1
        
        chosenCol = num2cell(mirnaNorm(:,x));
        chosenCol = [fieldNamescell(x); chosenCol];
        mirnaTest = [mirnaTest chosenCol];
        
    end
   
end
% split mirna into test and train 
fieldNamestrain = mirnaTrain(1,:);
fieldNamestest = mirnaTest(1,:);
mirnaTrainData = cell2mat(mirnaTrain(2:end,:));
mirnaTestData = cell2mat(mirnaTest(2:end,:)); 
% fs
[mirnaTopData, mirnaTopFinal, mirnaTopFeatures, fieldNamestrain, cellid1] = featureselection_mirna(cellid, cellid1, mirnaTrainData, fieldNamestrain);
% cox
[B, survivalFinal, riskscore, mirnaTest] = cox_mirna(mirnaTopData, mirnaTopFeatures, mirnaTest, fieldNamestrain, cellid1);
B_rnaseq = B; 
survivalFinal_rnaseq = survivalFinal; 
riskscore_rnaseq = riskscore;
mirnaTest_rnaseq = mirnaTest;
fieldNamestest_rnaseq = fieldNamestest;


% % only use samples common in both modalities
% if length(riskscore_mirna) > length(riskscore_rnaseq)
%     ind = ismember(fieldNamestest_mirna,fieldNamestest_rnaseq);
%     riskscore_mirna = riskscore_mirna(ind);
% elseif length(riskscore_rnaseq) > length(riskscore_mirna) 
%     ind = ismember(fieldNamestest_rnaseq,fieldNamestest_mirna); 
%     riskscore_rnaseq = riskscore_rnaseq(ind); 
% end

% updated code
if length(riskscore_mirna) > length(riskscore_rnaseq)
    riskscore_mirna = riskscore_mirna(1:length(riskscore_rnaseq));
elseif length(riskscore_rnaseq) > length(riskscore_mirna) 
    riskscore_rnaseq = riskscore_rnaseq(1:length(riskscore_mirna));
end

% % old code
%     ind = ismember(fieldNamestest_mirna,fieldNamestest_rnaseq);
%     riskscore_mirna = riskscore_mirna(ind);
% elseif length(riskscore_rnaseq) > length(riskscore_mirna) 
%     ind = ismember(fieldNamestest_rnaseq,fieldNamestest_mirna); 
%     riskscore_rnaseq = riskscore_rnaseq(ind); 
% end    
    
% average risk scores 
riskscore_integrated = (riskscore_mirna + riskscore_rnaseq)./2; 


% split risk scores into high and low risk groups 
[testSurvivalFinal, survivalTimeHigh, survivalTimeLow, censoredFlagHigh, censoredFlagLow, p] = kaplanmeier(riskscore_integrated, fieldNamestest, fold);

% plot the survival functions 
subplot(1,5,fold-1)
ecdf(survivalTimeHigh, 'Censoring', censoredFlagHigh, 'Function', 'survivor'); % high risk group
hold on
ecdf(survivalTimeLow, 'Censoring', censoredFlagLow,'Function', 'survivor'); % low risk group
str = sprintf('Fold %d : miRNA + mRNA', fold-1);
title(str)
xlabel('Survival Time (days)')
ylabel('Probability of Survival')
legend('High Risk Group', 'Low Risk Group')
txt = ['p = ',num2str(p)];
text(300,0.8,txt)
pValues = [pValues p]; 




% c index
survivalTimes = survivalFinal(:,4);
survivalTimes = survivalTimes(1:length(riskscore_integrated));
survivalTimes = cell2mat(survivalTimes)';
ci = concordanceIndex(survivalTimes,riskscore_integrated);

% ev vs cv 
 
ciCV = [ciCV ci]; 

% do the same thing for the actual test data 
[testSurvivalFinal, survivalTimeHigh, survivalTimeLow, censoredFlagHigh, censoredFlagLow, p] = kaplanmeier(riskscore_integrated, fieldNamestestActual, fold);

% c-index for external validation
survivalTimes = survivalFinal(:,4);
survivalTimes = survivalTimes(1:length(riskscore_integrated));
survivalTimes = cell2mat(survivalTimes)';
ciTest = concordanceIndex(survivalTimes,riskscore_integrated);
ciEV = [ciEV ciTest];

% plot the survival functions 
subplot(2,5,fold + 4)
ecdf(survivalTimeHigh, 'Censoring', censoredFlagHigh, 'Function', 'survivor'); % high risk group
hold on
ecdf(survivalTimeLow, 'Censoring', censoredFlagLow,'Function', 'survivor'); % low risk group
str = sprintf('Fold %d : miRNA + mRNA', fold-1);
title(str)
xlabel('Survival Time (days)')
ylabel('Probability of Survival')
legend('High Risk Group', 'Low Risk Group')
txt = ['p = ',num2str(p)];
text(300,0.8,txt)


end

% plot concordance interval for cross validation  
subplot(2,5,5)
scatter(cv, ciCV)
title('CV vs CI')
xlabel('Cross Validation')
ylabel('c-index')
axis([1 4 0 1])


% plot cross validation vs external validation
subplot(2,5, 10)
scatter(ciCV, ciEV, 'filled')
title('Cross Validation vs External Validation')
xlabel('Training Performance')
ylabel('Testing Performance')
axis([0.45 0.55 0.45 0.55])
hold on 
p = polyfit(ciCV,ciEV,1);
f = polyval(p,ciCV);
plot(ciCV,f)
hold off





