imputed_word = 'original';
d = dataset('file', ['~/GitHub/pqe/data/diffusion_components_', imputed_word, '.txt'], 'ReadVarNames', false, 'ReadObsNames', true);
md = dataset('file', '~/GitHub/pqe/data/expression_filtered_and_DE_genes_design_mat.txt',  ...
    'ReadVarNames', true, 'ReadObsNames', true);

EV = double(d);

GFP = md.EGFP;
E_stage = md.EStage;
E_num = str2double(strrep(E_stage, 'E', ''));


colorby = {standardize(GFP), E_num};
cm = redgreencmap(100);
cmg = flipud(cm(1:50,:)).^3;
hm = colormap('hot');
cmaps = {cmg; hm(1:50,:)};


plot_PCA_summary(EV, [], colorby, cmaps);
plotSave(['../figures/diffusion_plot_', imputed_word, '.png'])
close