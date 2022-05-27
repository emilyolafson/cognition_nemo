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

df=readtable("~/GIT/cognition_nemo/df_all.csv");
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

d2=readtable("~/GIT/cognition_nemo/df_all.csv");

%% Aaron email: get dominant variable from correlation
df=readtable("~/GIT/cognition_nemo/data/spreadsheets/cog_stroke_alldata_waisbinary.csv");
ids=readtable("~/GIT/cognition_nemo/ids_all.csv");

testname=df(:,38:80).Properties.VariableNames
subid=df.RedID

wais=df(:,38:80);

wais(:,40)=[]
wais(:,40)=[]

wais=table2array(wais)
%
newwais=wais(:,~(sum(isnan(wais),1)>125))
newnames=testname(:,~(sum(isnan(wais),1)>125))
newwais2=newwais(sum(isnan(newwais),2)==0,:)
zscore=normalize(newwais2)

tiledlayout(1,2,'padding','none')
nexttile;
imagesc(corr(zscore))
xticks(1:1:26)
xticklabels(newnames)
xtickangle(90)
yticks(1:1:26)
yticklabels(newnames)
set(gca,'TickLabelInterpreter','none')
colorbar
caxis([-1,1])
title('Correlation')
set(gca, 'FontSize', 14)

nexttile;
answer=corr(zscore,'rows', 'pairwise')
imagesc(inv(answer))
xticks(1:1:26)
xticklabels(newnames)
xtickangle(90)
yticks(1:1:26)
yticklabels(newnames)
set(gca,'TickLabelInterpreter','none')
colormap(flipud(brewermap([], 'RdBu')))
colorbar
%caxis([-80,80])
title('Inverse Correlation (Precision)')
set(gca, 'FontSize', 14)


newwais=wais(:,~(sum(isnan(wais),1)>100))
newnames=testname(:,~(sum(isnan(wais),1)>100))
newwais2=newwais(sum(isnan(newwais),2)==0,:)
zscore=normalize(newwais2)

[coeff,score,latent,tsquared,explained,mu]=pca(zscore)

imagesc(coeff)
xticks(1:1:size(coeff,1))
%ytickangle(90)
yticks(1:1:size(coeff,1))
yticklabels(newnames)
set(gca,'TickLabelInterpreter','none')
colormap(flipud(brewermap([], 'RdBu')))

answer=corr(zscore)
[~, D, W]=eig(answer)

[coeff2,score2,latent,tsquared,explained2,mu]=pca(sqrt(D)*W')
imagesc(coeff)
xticks(1:1:size(coeff,1))
xtickangle(90)
yticks(1:1:size(coeff,1))
yticklabels(newnames)
set(gca,'TickLabelInterpreter','none')
colormap(flipud(brewermap([], 'RdBu')))

diag(D./sum(sum(D)))

imagesc(inv(answer))

xticks(1:1:size(answer,1))
xticklabels(newnames)
xtickangle(90)
yticks(1:1:size(answer,1))
yticklabels(newnames)
set(gca,'TickLabelInterpreter','none')
colormap(flipud(brewermap([], 'RdBu')))


[coeff,score,latent,tsquared,explained,mu] = pca(answer, 'algorithm', 'als')

imagesc(score)
xtickangle(90)
yticks(1:1:41)
yticklabels(testname)
set(gca,'TickLabelInterpreter','none')
colormap(flipud(brewermap([], 'RdBu')))




%% Mark Bowren email: get dominant variable from WAIS correlation matrix
% core tests of WAIS
%Verbal Comprehension Index has three core subtests, which are: Similarities,
%Vocabulary, and Information.

t=readtable('~/GIT/cognition_nemo/data/WAIS_corr_matrix_full.csv')

tcorrs=t(:, 2:end)
tcorrs=table2array(tcorrs)
tcorrs=[ones(1,24)*NaN;tcorrs]

tcorrs(logical(eye(24)))=1 % fill diagonal with 1
A=tcorrs
B = triu(A.',1) + tril(A)

imagesc(B)
xticks(1:1:24)
xticklabels(colnames)
yticks(1:1:24)
yticklabels(colnames)

colnames=t(:,1)
colnames=table2array(colnames)
colnames=[{'BD'};colnames]

% PCA only of the working memory scores 

ntest=3

Bsub=[B(3, 3),B(3, 6),B(3, 11); B(6, 3),B(6, 6),B(6, 11); B(11, 3),B(11, 6),B(11, 11)];

% WM only: digit span, arithmetic, letter number 
% 3, 6, 11
% VC only: similarity, vocab, information
% 2, 5, 9
% PR only: block design, matrix reasoning, visual puzzles
% 1, 4, 8
% PS only: symbol search, coding, cancellation
% 7, 10, 14

colnamessub=[colnames_full(3),colnames_full(6),colnames_full(11)]
imagesc(Bsub)
yticks(1:1:3)
yticklabels(colnamessub)
xticks(1:1:3)
xticklabels(colnamessub)
xtickangle(90)

set(gca,'TickLabelInterpreter','none')
colormap(brewermap([], 'Reds'))

bar(mean(Bsub))
xticks(1:1:ntest)
xticklabels(colnamessub)
xtickangle(90)

% get eigenvectors and eigenvalues
[~, D, W]=eig(Bsub)
[coeff,~,~,~,explained,~]=pca(sqrt(D)*W')

imagesc(coeff)
yticks(1:1:ntest)
yticklabels(colnamessub)

figure
subplot(1,2,1)
bar(W(:,3))
xticks(1:1:ntest)
xticklabels(colnamessub)
xtickangle(90)
set(gca, 'FontSize', 15)

% variance explained between working memory subtests

clear sigsq
for x=1:3
    Bsub=B([3, 6, 11, 2, 5, 9, 1, 4, 8, 7, 10, 14], [3, 6, 11, 2, 5, 9, 1, 4, 8, 7, 10, 14]);
    ntest=3
    for y=1:3
        if x==y
            disp('test')
            continue
        end
        X=Bsub(x,x);
        Y=Bsub(x,y);
        XTX=X.'*X;
        XTY=X.'*Y;
        YTX=Y.*X;
        sigsq(x,y)=1-YTX*inv(XTX)*XTY;
    end
    if x==y
        continue
    end
end

rsq=1-sigsq
rsq(logical(eye(3)))=NaN
bar(mean(rsq, 'omitnan'))
xticks(1:1:3)
colnamessub=colnames_full([3, 6, 11])
%colnamessub=colnames_full(1:ntest)
xticklabels(colnamessub)
xtickangle(90)
ylabel('R^2')
set(gca, 'FontSize', 15)



%% Plot cognitive scores (TMT and WAIS), PCA with missing data.

df=readtable("~/GIT/cognition_nemo/data/spreadsheets/cog_stroke_alldata_waisbinary.csv");
ids=readtable("~/GIT/cognition_nemo/ids_all.csv");

testname=df(:,31:81).Properties.VariableNames
subid=df.RedID

% sort columns by number of nans
[c,idx]=sort(sum(isnan(wais)))
out=wais(:,idx)
names=testname(:,idx)

colormap(flipud(hot))
out2=(out - nanmean(out))./nanstd(out)
imagesc(out2)
xticks(1:1:51)
xticklabels(names)
xtickangle(90)
set(gca,'TickLabelInterpreter','none')
ylabel('Subject', 'FontSize', 15)
xlabel('Cognitive test', 'FontSize', 15)
ax = gca;
ax.FontSize = 16; 

% find at which test there's > 30% missing data
sum(~isnan(out2),1)/size(out2,1)

% select data up to that test
out2_reduced=out2(:,1:8)
names_reduced=names(:,1:8)

zscore=normalize(out2_reduced)
[coeff,score,latent,tsquared,explained,mu] = pca(zscore,'algorithm','als')

bar(explained);ylabel('% variance explained')
set(gca, 'FontSize', 14)

figure
biplot(coeff(:,1:2),'scores',score(:,1:2), 'varlabels', names_reduced)

set(gca,'FontSize', 15)

bar(coeff(:,1));xticks(1:1:11)
xticklabels(names_reduced)

set(gca,'TickLabelInterpreter','none')
ylabel('PC score', 'fontsize', 15)

xtickangle(90)
ax = gca;
ax.FontSize = 16; 


corr(one,two, 'rows','complete')


imagesc(corr(zscore,'rows','complete'))
xticks(1:1:11)
xticklabels(names_reduced)
xtickangle(90)
yticks(1:1:11)

yticklabels(names_reduced)
set(gca,'TickLabelInterpreter','none')
colormap(flipud(brewermap([], 'RdBu')))

sum(corr(out2_reduced,'rows','complete'))


subid=df.RedID
final=tmtwais(:,23)
nans=isnan(final)
final2=final(~nans,:)
subid2=subid(~nans,:)




[coeff,score,latent,tsquared,explained,mu] = pca(zscore,'algorithm','als')

bar(explained);ylabel('% variance explained')
set(gca, 'FontSize', 14)

figure
biplot(coeff(:,1:2),'scores',score(:,1:2), 'varlabels', names_reduced)

set(gca,'FontSize', 15)

bar(coeff(:,1));xticks(1:1:11)
xticklabels(names_reduced)

set(gca,'TickLabelInterpreter','none')
ylabel('PC score', 'fontsize', 15)

xtickangle(90)
ax = gca;
ax.FontSize = 16; 



%% Get final list (unbiased)
df=readtable("~/GIT/cognition_nemo/df_all.csv");
chronicity=df.('acute_chronic')
testname=df(:,31:81).Properties.VariableNames
subid=df.RedID

tmtwais=df(:,31:81);
tmtwais=table2array(tmtwais)

subid=df.RedID
final=tmtwais(:,23)
nans=isnan(final)
final2=final(~nans,:)
subid2=subid(~nans,:)
chronicity2=chronicity(~nans,:)
writecell(subid2, '~/GIT/cognition_nemo/test.csv')

writecell(chronicity2, '~/GIT/cognition_nemo/test.csv')
