% load cleaned dataframe
%subs with no lesion but cog scores

nolesion={'ca089', 'ca069', 'ca124', 'ca094'};

nolesionindex=[117];
index=zeros(149, 1)
index(nolesionindex)=1
df = readtable('df_all.csv');
df=df(:,2:end); % remove index variable
colnamesall=df.Properties.VariableNames
wasi=colnamesall(30:72)

clear ratio
for i=1:90
    df=readtable("~/GIT/cognition_nemo/df_all.csv");
    ids = readtable('ids_all.csv');

    ids=ids(2:end,:); % remove index variable
    ids=table2array(ids(:,2));
    df=df(:,2:end); % remove index variable
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
legend('TMT', 'WAIS')

% implement best tradeoff 
besti=31
df = readtable('dataframe_noduplicates_NaNsincluded_chronic.csv');

ids = readtable('IDlist_noduplicates_NaNsincluded_chronic.csv');
ids=ids(2:end,:); % remove index variable
ids=table2array(ids(:,2));
ids=ids(~index,:);
df=df(:,2:end); % remove index variable
df=df(~index,:);

nans = readtable('nans_chronic.txt');
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


%% only chronic scores
df = readtable('dataframe_noduplicates_NaNsincluded.csv');
df=df(:,3:end)
nans_rows_thresh=isnan(df.("WAIS_IV_DS_TOTAL_RAW"))
sum(nans_rows_thresh)

% Output new list of subject IDs
ids = readtable('IDlist_noduplicates_NaNsincluded_chronic.csv');
ids=ids(2:end,:); % remove index variable
ids=ids(~nans_rows_thresh,"Var3");
writetable(ids, '/Users/emilyolafson/GIT/cognition_nemo/final_subIDs.csv')


colnames=demogs.Properties.VariableNames
writecell(colnames, '/Users/emilyolafson/GIT/cognition_nemo/colnames.csv')

demogs=demogs(~nans_rows_thresh,:); % remove rows of subjects excluded 
writetable(demogs, '/Users/emilyolafson/GIT/cognition_nemo/demographic_final.csv')

% save matrix with cognitive scores for only subset of subjects 
df = readtable('dataframe_noduplicates_NaNsincluded_chronic.csv');
df=df(:,3:end); % remove index variable

test=df(~nans_rows_thresh,:)
colnames=test.Properties.VariableNames

x=table2array(test)
zscore = (x - nanmean(x))./nanstd(x);

[coeff,score,latent,tsquared,explained,mu] = pca(zscore,'algorithm','als')
pc1=score(:,1)
corr(pc1,table2array(test2))

bar(explained);title('Variance explained by PCA component')
figure
biplot(coeff(:,1:2),'scores',score(:,1:2), 'varlabels', colnames)


test2=df(~nans_rows_thresh,"WAIS_IV_DS_TOTAL_RAW")
final_dataframe=table2array(test)
writematrix(final_dataframe, '/Users/emilyolafson/GIT/cognition_nemo/WAIS_IV_DS_TOTAL_RAW_79.csv')

pc1=zscore(pc1)
writematrix(pc1, '/Users/emilyolafson/GIT/cognition_nemo/pc1_79.csv')



%%
df_chronic_wais=readtable('df_chronic_wais.csv')
df_acute_wais=readtable('df_acute_wais.csv')

ids_chronic=readtable('ids_chronic.csv')
ids_acute=readtable('ids_acute.csv')

ids_chronic=ids_chronic(2:end,:);
ids_acute=ids_acute(2:end,:);

ids_all= [ids_chronic; ids_acute]

ids_all=sortrows(ids_all)
ids_all_a=table2array(ids_all(:,2))

ids_chronic_a=table2array(ids_chronic(:,2))
ids_acute_a=table2array(ids_acute(:,2))

writecell(ids_chronic_a, '/Users/emilyolafson/GIT/cognition_nemo/ids_chronic_a.csv')

writecell(ids_acute_a, '/Users/emilyolafson/GIT/cognition_nemo/ids_acute_a.csv')

writecell(table2array(ids_all(:,2)), '/Users/emilyolafson/GIT/cognition_nemo/ids_all.csv')



%% check to see if distribution of the two batches of chaco are different?

before=load('/Users/emilyolafson/GIT/cognition_nemo/SC/txtfiles/ca065_nemo_output_chacoconn_fs86subj_nemoSC_volnorm.txt')
after=load('/Users/emilyolafson/GIT/cognition_nemo/SC/txtfiles_rem/ca144_nemo_output_chacoconn_fs86subj_nemoSC_volnorm.txt')

corr(reshape(before, 1, [])', reshape(after, 1, [])')
plot(reshape(before, 1, [])', reshape(after, 1, [])')




%%

df=readtable("~/GIT/cognition_nemo/data/spreadsheets/df_all.csv");
ids = readtable('ids_all.csv');



IDs=df.RedID
colnames=df.Properties.VariableNames
colnames=colnames(5:end)
df1=table2array(df(:,5:end))


isnan(df1)

zscore=normalize(df1)
[coeff,score,latent,tsquared,explained,mu] = pca(zscore,'algorithm','als')

bar(explained);title('Variance explained by PCA component')
set(gca, 'FontSize', 14)

figure
biplot(coeff(:,1:2),'scores',score(:,1:2), 'varlabels', colnames)


set(gca, 'FontSize', 14)
after=load('/Users/emilyolafson/GIT/cognition_nemo/SC/txtfiles_rem/ca144_nemo_output_chacoconn_fs86subj_nemoSC_volnorm.txt')

