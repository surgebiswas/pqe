%% Load the data
clear;
rng('default');
cd('~/GitHub/pqe/src');

d = dataset('file', '~/GitHub/pqe/data/expression_filtered_and_DE_genes_expression_mat.txt', ... 
    'ReadVarNames', true, 'ReadObsNames', true);
md = dataset('file', '~/GitHub/pqe/data/expression_filtered_and_DE_genes_design_mat.txt',  ...
    'ReadVarNames', true, 'ReadObsNames', true);

y = double(d)';
sy = standardize(y); % standardized log(RPKM + 0.001)
genes = get(d, 'ObsNames');
cells = get(d, 'VarNames');
GFP = md.EGFP;
E_stage = md.EStage;
E_num = str2double(strrep(E_stage, 'E', ''));

% batch ?
c = char(cells'); b = double(c(:,1));
batch = b - min(b);


%% PCA
% Set up plotting parameters
fexp = @(g) fetch_expression(sy, genes, g);

% cm = redgreencmap(100);
% cmg = flipud(cm(1:50,:)).^3;
cmg = prgn;
hm = colormap('hot');
cmaps = {cmg; hm(1:50,:)};


% Compute PCA on noisy expression matrix.
% PCA is computed on standardized expression matrix so PCs aren't dominated
% by high expression genes. 
[~,s,~,~,pexp] = pca(sy, 'NumComponents', 3);
lmx1a = fexp('Lmx1a');
colorby = {standardize(lmx1a), E_num};
plot_PCA_summary(s, pexp, colorby, cmaps);
plotSave('../figures/PCA_summary_original_matrix.png');
close


% Use MAGIC to perform imputation 
path(genpath('~/GitHub/magic'), path)
npca = 20; % ususally between 10 and 200
ka = 3; % can be smaller, eg 3
k = 9; % can be smaller, eg 9
t = 6; % usually between 6 and 12, smaller ka/k requires bigger t
y_with_gfp = [GFP, y];
y_imputed_with_gfp = run_magic(y_with_gfp, t, 'npca', npca, 'ka', ka, 'k', k, 'rescale_to', 0);
sy_imputed_with_gfp = standardize(y_imputed_with_gfp);
sy_imputed = sy_imputed_with_gfp(:,2:end);
GFP_imputed = sy_imputed_with_gfp(:,1);


% Compute PCA on the standardized imputed expression matrix. 
fexp2 = @(g) fetch_expression(sy_imputed, genes, g);
[~,s_magic,~,~,pexp_magic] = pca(sy_imputed, 'NumComponents', 3);
lmx1a_imputed = fexp2('Lmx1a');
colorby = {standardize(lmx1a_imputed), E_num};
plot_PCA_summary(s_magic, pexp_magic, colorby, cmaps);
plotSave('../figures/PCA_summary_imputed_matrix.png');
close


%% Hierarchical clustering
[ri,ci] = hclust(sy);

figure;
imagesc(sy(ri,ci), [-2,2]); 
set(gca, 'XTick', []);
set(gca, 'YTick', []);
xlabel('Genes');
ylabel('Cells');
plotSave('../figures/clustered_original.png');
close

figure;
imagesc(sy_imputed(ri,ci), [-2,2]);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
xlabel('Genes');
ylabel('Cells');
plotSave('../figures/clustered_imputed.png');
close

figure
imagesc(lmx1a_imputed(ri)); colormap(cmg);
colorbar;
set(gca, 'XTick', []);
set(gca, 'YTick', []);
plotSave('../figures/lmx1a_annotation_sorted.png');
close

figure
imagesc(E_num(ri)); colormap(hm(1:50,:));
colorbar;
set(gca, 'XTick', []);
set(gca, 'YTick', []);
plotSave('../figures/Estage_annotation_sorted.png');
close








%% Let's examine some genes in "Axis 2". 

axis2_list = {'Barhl1', 'Barhl2', 'Dbx1', 'Epha3', 'Wnt8b', 'Rspo2', ...
     'Pitx2', 'Lmx1a'};
 
gene_plots(axis2_list, fexp2, E_num, lmx1a_imputed, cmg);
plotSave('../figures/axis_2_gene_plots_imputed.png');
close


% Other interesting genes.
% Cell cycle.
glist = {'Camk2a', 'Myb', 'Brca2', 'Cdk2', 'Casp3', 'Brca1', 'Itgb1', 'Skp2'};
gene_plots(glist, fexp2, E_num, lmx1a_imputed, cmg);
plotSave('../figures/cell_cycle_gene_plots.png');
close


% Those that are predictive of neurons being DA
glist = {'Barhl1', 'Ebf2', 'Ebf3', 'Foxb1', 'Myt1', 'Neurod1', 'Neurod2'};
gene_plots(glist, fexp2, E_num, lmx1a_imputed, cmg);
plotSave('../figures/DA_gene_plots.png')


glist = {'Barhl1', 'Ebf2', 'Ebf3', 'Foxb1', 'Myt1', 'Neurod1', 'Neurod2'};
gene_plots(glist, fexp, E_num, lmx1a, cmg);
plotSave('../figures/DA_gene_plots_original.png')

%% Export the imputed expression matrix for downstream analysis.
md.EGFP_imputed = GFP_imputed;
dout = mat2dataset(y_imputed_with_gfp(:,2:end)', 'VarNames', cells, 'ObsNames', genes);
export(dout, 'file', '../data/expression_filtered_and_DE_genes_imputed_expression_mat.txt');
export(md, 'file', '../data/expression_filtered_and_DE_genes_design_mat.txt');








