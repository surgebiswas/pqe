%% Load the data
clear;
rng('default');
cd('~/GitHub/pqe/src');

d = dataset('file', '~/GitHub/pqe/data/expression_filtered_and_DE_genes_expression_mat.txt', ... 
    'ReadVarNames', true, 'ReadObsNames', true);
md = dataset('file', '~/GitHub/pqe/data/expression_filtered_and_DE_genes_design_mat.txt',  ...
    'ReadVarNames', true, 'ReadObsNames', true);

y = double(d)';

y_rpkm = 2.^y - 0.001; y_rpkm(y_rpkm < 1e-6) = 0;

% assume each gene is approximately 1 kb .
% the following is an approximate count table.
map_rate = 0.64; %average
y_rpk = round((y_rpkm*1e6)./repmat(map_rate*md.ApproxDepth, 1, size(y_rpkm,2)))';


rarefy_amts = fliplr([0.01 0.05 0.1 0.3 0.5 0.8 0.9 1]); % 
for i = 1 : length(rarefy_amts)
    if rarefy_amts(i) == 1
        rar_mat = y_rpk;
    else
        rar_mat = y_rpk;
        for j = 1 : size(rar_mat,2)
            rar_mat(:,j) = rarefy(rar_mat(:,j), round(sum(y_rpk(:,j))*rarefy_amts(i)) ) ;
        end
    end
    
    % Convert back to rpkm
    rar_rpkm = rar_mat.*(map_rate*repmat(md.ApproxDepth', size(rar_mat,1), 1))/(1e6*rarefy_amts(i));
    dout = mat2dataset(rar_rpkm, 'VarNames', get(d, 'VarNames'), 'ObsNames', get(d, 'ObsNames'));
    export(dout, 'file', ['~/GitHub/pqe/data/rarefaction/rarefied_', num2str(100*rarefy_amts(i)), '.txt']);
    disp(i);
end





