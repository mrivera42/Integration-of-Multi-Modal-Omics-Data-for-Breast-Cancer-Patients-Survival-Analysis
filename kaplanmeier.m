function [testSurvivalFinal, survivalTimeHigh, survivalTimeLow, censoredFlagHigh, censoredFlagLow, p] = kaplanmeier(riskscore, fieldNamestest, fold)

[survivalFinal, ~, ~, ~, sampleid] = preprocessing_survival();

%calculate median
medi = median(riskscore);

% match up the survival times to risk score
%ismember and save the indices of where the test names are in survival
testSampleIDs = fieldNamestest(2:end);
testSampleIDs = testSampleIDs';
survivalSampleIDs = sampleid;

[~, indices] = ismember(testSampleIDs, survivalSampleIDs); %(sampleids from test, sampleids from survival)
%apply mask to survival
testSurvivalFinal = survivalFinal(indices);


%make table
risk_survival = []; %FILL IN

highrisk = [];
lowrisk = [];
highriskscore = [];
lowriskscore = [];
for x = 1:length(riskscore)
    if riskscore(x) > medi
       %move row to high risk group
       risk_survival = survivalFinal(x,:);
       highrisk = [highrisk; risk_survival];
       highriskscore = [highriskscore; riskscore(x)];
    else
        %move row to low risk group
        risk_survival = survivalFinal(x,:);
        lowrisk = [lowrisk; risk_survival];
        lowriskscore = [lowriskscore; riskscore(x)];
    end
end


% pull ostime and os for each risk group
survivalTimeHigh = cell2mat(highrisk(:,4));
censoredFlagHigh = cell2mat(highrisk(:,2));
survivalTimeLow = cell2mat(lowrisk(:,4));
censoredFlagLow = cell2mat(lowrisk(:,2));


% compute pvalues
[~,p,~,~] = ttest2(survivalTimeHigh,survivalTimeLow,'Vartype','unequal');


% % plot the survival functions 
% subplot(1,2,1)
% ecdf(survivalTimeHigh, 'Censoring', censoredFlagHigh, 'Function', 'survivor'); % high risk group
% hold on
% ecdf(survivalTimeLow, 'Censoring', censoredFlagLow,'Function', 'survivor'); % low risk group
% str = sprintf('KM Survival Probability Fold %d : miRNA', fold);
% title(str)
% xlabel('Survival Time (days)')
% ylabel('Probability of Survival')
% legend('High Risk Group', 'Low Risk Group')
% txt = ['p = ',num2str(p)];
% text(300,0.8,txt)
% 
% hold off

% % plot the cumulative hazard functions 
% subplot(1,2,2)
% ecdf(survivalTimeHigh, 'Censoring', censoredFlagHigh,'Function', 'cumulative hazard'); % high risk group
% hold on
% ecdf(survivalTimeLow, 'Censoring', censoredFlagLow,'Function', 'cumulative hazard'); % low risk group
% str = sprintf('KM Cumulative Hazard Fold %d : miRNA', fold);
% title(str)
% xlabel('Survival Time (days)')
% ylabel('Cumulative Hazard')
% legend('High Risk Group', 'Low Risk Group')
% txt = ['p = ',num2str(p)];
% text(4.5,0.8,txt)
% hold off

end
