% load cleaned dataframe
%subs with no lesion but cog scores

nolesion={'ca089', 'ca069', 'ca124', 'ca094'};
nolesionindex=[64, 82, 87, 115];
index=zeros(149, 1)
index(nolesionindex)=1
df = readtable('dataframe_noduplicates_NaNsincluded.csv');
df=df(:,2:end); % remove index variable
colnamesall=df.Properties.VariableNames
tmt=colnamesall(1:8)
wasi=colnamesall(9:end)

clear ratio
for i=1:90
    df = readtable('dataframe_noduplicates_NaNsincluded.csv');
    ids = readtable('IDlist_noduplicates_NaNsincluded.csv');
    ids=ids(~index,:);

    ids=ids(2:end,:); % remove index variable
    ids=table2array(ids(:,2));
    df=df(:,2:end); % remove index variable
    df=df(~index,:);
    nans = readtable('nans.txt');
    nans=table2array(nans);
    nans=nans(~index,:);

    % find measures with NaNs in at least 40 subjects, remove them.
    nans_cols=sum(nans,1);
    nans_cols_thresh=nans_cols>i;
    
    nans=nans(:,~nans_cols_thresh);

    % find subjects with missing data in the remaining measures.
    nans_rows=sum(nans,2);
    nans_rows_thresh=nans_rows>=1;
    nans=nans(~nans_rows_thresh,:);

    % remaining subjects & measures:
    ratio{i}=size(nans)
    df=df(:, ~nans_cols_thresh);
    df=df(~nans_rows_thresh,:);
    colnames=df.Properties.VariableNames;
    tmtsum(i)=sum(ismember(colnames, tmt));
    wasisum(i)=sum(ismember(colnames, wasi));
end


ratio=cell2mat(ratio')
% look at trade off between n and p 
figure;
plot(1:90, ratio)
set(gca, 'FontSize', 14)
xlabel('Cutoff')
ylabel('Count')
hold on;
bar(1:90, [tmtsum;wasisum])

legend('Number of subjects', 'Number of variables', 'TMT', 'WASI')

set(gca, 'FontSize', 14)
xlabel('Cutoff')
ylabel('Count')
legend('TMT', 'WASI')


% implement best tradeoff 
besti=31
df = readtable('dataframe_noduplicates_NaNsincluded.csv');
ids = readtable('IDlist_noduplicates_NaNsincluded.csv');
ids=ids(2:end,:); % remove index variable
ids=table2array(ids(:,2));
ids=ids(~index,:);
df=df(:,2:end); % remove index variable
df=df(~index,:);

nans = readtable('nans.txt');
nans=table2array(nans);
nans=nans(~index,:);

% find measures with NaNs in at least 31 subjects, remove them.
nans_cols=sum(nans,1);
nans_cols_thresh=nans_cols>besti;
nans=nans(:,~nans_cols_thresh);

% find subjects with missing data in the remaining measures.
nans_rows=sum(nans,2);
nans_rows_thresh=nans_rows>=1;
nans=nans(~nans_rows_thresh,:);

% PCA
df=df(:, ~nans_cols_thresh);
df=df(~nans_rows_thresh,:);


colnames=df.Properties.VariableNames
df=table2array(df);


[coeff,score,latent,tsquared,explained,mu]=pca(zscore(df));
figure;
bar(explained);title('Variance explained by PCA component')
figure
biplot(coeff(:,1:2),'scores',score(:,1:2),'varlabels',colnames');




% Output new list of subject IDs
ids=ids(~nans_rows_thresh,:);
writecell(ids, '/Users/emilyolafson/GIT/cognition_nemo/final_subIDs.csv')

% Load demographic data
demogs = readtable('demogs_noduplicates_NaNsincluded.csv');
demogs=demogs(:,3:end)
demogs=removevars(demogs,{'Dates_of_neuropsychological_testing'})

demogs=demogs(~nans_rows_thresh,:); % remove rows of subjects excluded 
colnames=demogs.Properties.VariableNames

writecell(colnames, '/Users/emilyolafson/GIT/cognition_nemo/colnames.csv')

writetable(demogs, '/Users/emilyolafson/GIT/cognition_nemo/demographic_final.csv')

writematrix(table2array(df), '/Users/emilyolafson/GIT/cognition_nemo/df_minimalfeatures.csv')

writematrix(df.T_MT_A_TIME, '/Users/emilyolafson/GIT/cognition_nemo/T_MT_A_TIME.csv')
writematrix(df.WAIS_IV_DS_TOTAL_RAW, '/Users/emilyolafson/GIT/cognition_nemo/WAIS_IV_DS_TOTAL_RAW.csv')

% save matrix with cognitive scores for only subset of subjects 
df = readtable('dataframe_noduplicates_NaNsincluded.csv');
df=df(:,2:end); % remove index variable

test=df(~nans_rows_thresh,:)
final_dataframe=table2array(test)
writematrix(final_dataframe, '/Users/emilyolafson/GIT/cognition_nemo/dataframe_allscores_101subs.csv')


nancols=sum(isnan(final_dataframe),2)./40
