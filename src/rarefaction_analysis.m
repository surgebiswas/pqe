% Rarefaction analysis
clear;
rng('default');
cd('~/GitHub/pqe/data/rarefaction');

vals = fliplr([1, 5, 10, 30, 50, 80, 90]);
refDC = load('DC_100.txt');

R = zeros(length(vals),3);
for i = 1 : length(vals)
    dc = load(['DC_', num2str(vals(i)), '.txt']);
    R(i,:) = abs(diag(corr(dc, refDC)));
end


semilogx(vals,R, '-', 'LineWidth', 3);
legend('DC1', 'DC2', 'DC3', 'Location', 'SouthEast');
axis square
xlabel('Percent of reads kept for each cell');
ylabel('Pearson correlation of DC to full data DC')
plotSave('~/GitHub/pqe/figures/rarefaction_analysis_DC_correlation.png');
close

md = dataset('file', '~/GitHub/pqe/data/expression_filtered_and_DE_genes_design_mat.txt',  ...
    'ReadVarNames', true, 'ReadObsNames', true);
GFP = md.EGFP;
E_stage = md.EStage;
E_num = str2double(strrep(E_stage, 'E', ''));


cm = redgreencmap(100);
cmg = flipud(cm(1:50,:)).^3;
hm = colormap('hot');
cmaps = {cmg; hm(1:50,:)};

colorby = {standardize(GFP), E_num};
plot_PCA_summary(dc, [], colorby, cmaps);
plotSave('~/GitHub/pqe/figures/rarefaction_analysis_DC_plot_of_1percent_rarefied.png');
close