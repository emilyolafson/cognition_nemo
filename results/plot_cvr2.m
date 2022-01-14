% read cv results.

sc=readtable('GIT/cognition_nemo/results/100perm_ridge_spearmancorr_featureselect_r2_SC.txt')
fc=readtable('GIT/cognition_nemo/results/100perm_ridge_spearmancorr_featureselect_r2_FC.txt')

sc=table2array(sc)
fc=table2array(fc)

violinplot([sc, fc])
xticklabels(["SC", "FC"])
xlabel('Input')
ylabel('average test R^2 over 10 folds')

%writematrix(array, 'GIT/cognition_nemo/results/100perm_ridge_spearmancorr_featureselect_r2_SC.txt')
title('WAIS prediction')
text(0.6, 0.25,sprintf('mean R^2 = %0.4f', mean(sc)))
text(1.5, 0.25,sprintf('mean R^2 = %0.4f', mean(fc)))

violinplot(array)
mean(array)