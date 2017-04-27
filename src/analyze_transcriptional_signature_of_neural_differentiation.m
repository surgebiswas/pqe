%% Load the data
clear;
rng('default');
cd('~/GitHub/pqe/src');

d = dataset('file', '~/GitHub/pqe/data/expression_filtered_and_DE_genes_imputed_expression_mat.txt', ... 
    'ReadVarNames', true, 'ReadObsNames', true);
md = dataset('file', '~/GitHub/pqe/data/expression_filtered_and_DE_genes_design_mat.txt',  ...
    'ReadVarNames', true, 'ReadObsNames', true);

y = double(d)';
sy = standardize(y); % standardized log(RPKM + 0.001)
genes = get(d, 'ObsNames');
cells = get(d, 'VarNames');
GFP = md.EGFP_imputed;
E_stage = md.EStage;
E_num = str2double(strrep(E_stage, 'E', ''));


% Transcription factors
tf = dataset('file', '~/GitHub/pqe/data/tf_list.txt', 'ReadVarNames', false, 'ReadObsNames', true);
tfs = tf.Var2;
[int_tfs, ia] = intersect(genes, tfs);
sy_tfs = sy(:,ia);


% Load diffusion component data
dc = load('~/GitHub/pqe/data/rarefaction/DC_100.txt');




% Plot with lmx1a
lmx1a = sy(:,find(strcmpi(genes, 'Lmx1a')));
scatter3(dc(:,1), dc(:,2), dc(:,3), 40, lmx1a, 'filled');
title('Diffusion plot');
colormap(prgn);
colorbar
xlabel('DC1');
ylabel('DC2');
zlabel('DC3');
c = colorbar;
c.Label.String = 'Lmx1a expression';
axis tight;
buffer_axis
savefig('~/GitHub/pqe/figures/Diffusion_plot_imputed_expression_matrix_lmx1a.fig')


% plot with E_num
hm = hot;
scatter3(dc(:,1), dc(:,2), dc(:,3), 40, E_num, 'filled');
title('Diffusion plot');
colormap(hm(1:50,:));
colorbar
xlabel('DC1');
ylabel('DC2');
zlabel('DC3');
c = colorbar;
c.Label.String = 'Embryonic stage';
axis tight;
buffer_axis
savefig('~/GitHub/pqe/figures/Diffusion_plot_imputed_expression_matrix_Enum.fig')


% 2D diffusion plot
figure
scatter(dc(:,1), dc(:,2), 40, lmx1a, 'filled');
colormap(prgn);
colorbar
xlabel('DC1');
ylabel('DC2');
c = colorbar;
c.Label.String = 'Lmx1a expression';
axis square
axis tight;
buffer_axis
plotSave('~/GitHub/pqe/figures/Diffusion_plot_2D_with_lmx1a_imputed.png');
close

% Axis 2 gene signature overlaid onto diffusion 
a2genes = {'Barhl1', 'Barhl2', 'Dbx1', 'Epha3', 'Wnt8b', 'Rspo2', 'Pitx2'};
[~,a2s] = pca(sy(:,find(steq(genes, a2genes))));
a2sig = a2s(:,1);
scatter(dc(:,1), dc(:,2), 40, a2sig, 'filled');
colormap(prgn);
colorbar
xlabel('DC1');
ylabel('DC2');
c = colorbar;
c.Label.String = 'Axis 2 signature expression';
axis square
axis tight;
buffer_axis
plotSave('~/GitHub/pqe/figures/Diffusion_plot_2D_with_axis2_signature_imputed.png');
close



% HOw does pca look again?
[~,s] = pca(sy);
scatter(s(:,1), s(:,2), 40, lmx1a, 'filled');
colormap(prgn);
colorbar
xlabel('PC1');
ylabel('PC2');
c = colorbar;
c.Label.String = 'Lmx1a expression';
axis square;
axis tight;
buffer_axis
plotSave('~/GitHub/pqe/figures/PCA_plot_2D_with_lmx1a_imputed.png');
close


% PCA of TFs
[~,stf] = pca(sy_tfs);
scatter(stf(:,1), stf(:,2), 40, lmx1a, 'filled');
colormap(prgn);
colorbar
xlabel('PC1');
ylabel('PC2');
c = colorbar;
c.Label.String = 'Lmx1a expression';
axis square;
axis tight;
buffer_axis
plotSave('~/GitHub/pqe/figures/PCA_plot_2D_of_TFs_with_lmx1a_imputed.png');
close


% TF heatmap
[ri,ci] = hclust(sy_tfs);
[~,sidx] = sort(dc(:,1));

row_order = sidx; order_name = 'DC';
figure;
imagesc(sy_tfs(row_order,ci), [-2,2]);
colormap(prgn);
colorbar
axis equal
axis tight
set(gca, 'XTick', []);
set(gca, 'YTick', []);
plotSave('~/GitHub/pqe/figures/TF_heatmap.png');
close


figure
annot = [standardize(E_num), lmx1a, standardize(dc(:,1:3))];
imagesc(annot(row_order,:));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
plotSave('~/GitHub/pqe/figures/TF_heatmap_annotation.png');
close




% Find TFs that are DE wrt time AND lmx1a expression:
interaction = lmx1a.*dc(:,1);
R = corr(interaction, sy_tfs, 'type', 'spearman');
imagesc(R(ci),[-1,1]);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
plotSave('~/GitHub/pqe/figures/TF_correlation_to_interaction.png');
close

% Overlay onto Tradict paper's data.
load('Tradict_sub_data.mat');
mask = abs(R) > 0.5;
tfs_hot = int_tfs(mask);

[tfs_hot, num2cell(R(mask))']

idx = find(steq(tids_pqe, tfs_hot));
lY_pqe_sub = lY_pqe(:,idx);
[csig,sig] = pca(standardize(lY_pqe_sub), 'NumComponents', 1);

scatter(s_pqe(:,1), s_pqe(:,2), 4, sig, 'filled');
hm = hot;
colormap(hm(1:50,:));
caxis([0 8.2]);
box on
xlabel('PC1');
ylabel('PC2');
set(gca, 'XTick', []);
set(gca, 'YTick', []);
axis square
axis tight
buffer_axis
colorbar
plotSave('~/GitHub/pqe/figures/tradict_extension_interaction_TFs.png');


% What if we relax to include all TFs?
mask = abs(R) > 0;
tfs_hot = int_tfs(mask);

idx = find(steq(tids_pqe, tfs_hot));
lY_pqe_sub = lY_pqe(:,idx);
[~,sig] = pca(standardize(lY_pqe_sub), 'NumComponents', 1);

scatter(s_pqe(:,1), s_pqe(:,2), 4, sig, 'filled');
hm = hot;
colormap(hm(1:50,:));
caxis([0 11]);
box on
xlabel('PC1');
ylabel('PC2');
set(gca, 'XTick', []);
set(gca, 'YTick', []);
axis square
axis tight
buffer_axis
colorbar
plotSave('~/GitHub/pqe/figures/tradict_extension_all_TFs.png');
close






