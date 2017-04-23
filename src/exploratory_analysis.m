%% Load the data
clear;
rng('default');
cd('~/GitHub/src/');

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


%% PCA
% Set up plotting parameters
fexp = @(g) fetch_expression(y, genes, g);

cm = redgreencmap(100);
cmg = flipud(cm(1:50,:)).^3;
cmaps = {cmg; jet};


% Compute PCA on noisy expression matrix.
% PCA is computed on standardized expression matrix so PCs aren't dominated
% by high expression genes. 
[~,s,~,~,pexp] = pca(sy, 'NumComponents', 3);
colorby = {standardize(GFP), E_num};
plot_PCA_summary(s, pexp, colorby, cmaps);
plotSave('../figures/PCA_summary_original_matrix.png');
close


% Use MAGIC to perform imputation 
path(genpath('~/GitHub/magic'), path)
npca = 20; % ususally between 10 and 200
ka = 3; % can be smaller, eg 3
k = 9; % can be smaller, eg 9
t = 10; % usually between 6 and 12, smaller ka/k requitres bigger t
sy_with_gfp = [standardize(GFP), sy];
sy_imputed_with_gfp = run_magic(sy_with_gfp, t, 'npca', npca, 'ka', ka, 'k', k, 'rescale_to', 0);
sy_imputed = standardize(sy_imputed_with_gfp(:,2:end));
GFP_imputed = sy_imputed_with_gfp(:,1);


% Compute PCA on the standardized imputed expression matrix. 
[~,s_magic,~,~,pexp_magic] = pca(sy_imputed, 'NumComponents', 3);
colorby = {GFP_imputed, E_num};
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
imagesc(GFP_imputed(ri)); colormap(cmg);
colorbar;
set(gca, 'XTick', []);
set(gca, 'YTick', []);
plotSave('../figures/gfp_annotation_sorted.png');
close

figure
imagesc(E_num(ri)); colormap(jet);
colorbar;
set(gca, 'XTick', []);
set(gca, 'YTick', []);
plotSave('../figures/Estage_annotation_sorted.png');
close








%% Let's examine some genes in "Axis 2". 
fexp2 = @(g) fetch_expression(sy_imputed, genes, g);
axis2_list = {'Barhl1', 'Barhl2', 'Dbx1', 'Epha3', 'Wnt8b', 'Rspo2', ...
     'Pitx2'};
gx = {};
for i = 1 : length(axis2_list)
    gx{i} = fexp2(axis2_list{i});
end
gx{end+1} = gfp_imputed;
axis2_list{end+1} = 'EGFP';
 
 
idx = 1;
for i = 1 : length(axis2_list)
    subplot(ceil(length(axis2_list)/2), 2, idx);
    
    if ~isempty(gx{i})
        plot(E_num + 0.3*rand(size(E_num)) - 0.15, gx{i}, '.k')
        
        hold on
        ue = unique(E_num);
        day_means = zeros(1,length(ue));
        for j = 1 : length(day_means)
            day_means(j) = mean(gx{i}(E_num == ue(j)));
        end
        plot(ue, day_means, '-or', 'MarkerFaceColor', 'r');
        set(gca, 'XTick', ue);
        
        
        
        idx = idx + 1;
        title(axis2_list{i});
    end
end
plotSave('../figures/axis_2_gene_plots_imputed.png');
close












