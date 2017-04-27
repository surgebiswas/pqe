function gene_plots( gene_list, fexp, E_num, colorby, cmap )

fexp2 = fexp;
axis2_list = gene_list;
Lmx1a = colorby;
cmg = cmap;

gx = {};
for i = 1 : length(axis2_list)
    gx{i} = fexp2(axis2_list{i});
end
gx{end+1} = Lmx1a;
axis2_list{end+1} = 'Lmx1a';
 
 
idx = 1;
for i = 1 : length(axis2_list)
    subplot(ceil(length(axis2_list)/2), 2, idx);
    
    if ~isempty(gx{i})
        scatter(E_num + 0.3*rand(size(E_num)) - 0.15, gx{i}, 4, Lmx1a)
        colormap(cmg)
        
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


end

