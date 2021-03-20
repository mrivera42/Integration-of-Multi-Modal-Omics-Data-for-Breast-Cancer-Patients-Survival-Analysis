% create and display heatmap of model performances across 4-fold 
cdata = [0.4999 0.4999 0.5053; 0.4502 0.4514 0.4798 ; 0.5322 0.4794 0.5028; 0.5050 0.5052 0.4748];
xvalues = {'miRNA','mRNA','miRNA + mRNA'};
yvalues = {'Fold 1','Fold 2','Fold 3','Fold 4'};
h = heatmap(xvalues,yvalues,cdata);

h.Title = '4-fold C-Index Performance';
h.XLabel = 'Modality';
h.YLabel = 'Cross Validation';

% create and display table of model performances and parameters 
h={'Fold' 'miRNA' 'mRNA' 'miRNA + mRNA'};
data=[1    0.4999     0.4999        0.5053
  2        0.4502     0.4514        0.4798
  3        0.5322     0.4794        0.5028
  4        0.5050     0.5052        0.4748 ];
f=figure;
t1=uitable(f,'data',data,'columnname',h);

% create and display table of p-values across folds 
h={'Fold' 'miRNA' 'mRNA' 'miRNA + mRNA'};
data=[1    0.61961     0.79023        0.71143
  2        0.45737     0.0076219        0.28302
  3        0.53458     0.53735        0.51501
  4        0.29056     0.69345        0.72231 ];
f=figure;
t2=uitable(f,'data',data,'columnname',h);