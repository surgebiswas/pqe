function plot_PCA_summary( s, pexp, colorby, cmaps )

idx = 1;
figure
for c = 1 : length(colorby)
    for i = 1 : 3
        for j = i+1 : 3
            subplot(length(colorby), 3, idx);

            scatter(s(:,i), s(:,j), 7, colorby{c}, 'filled');
            axis tight
            buffer_axis;
            xlabel(sprintf('PC%0.0f (%0.2f%%)', i, pexp(i)));
            ylabel(sprintf('PC%0.0f (%0.2f%%)', j, pexp(j)));
            set(gca, 'XTick', []);
            set(gca, 'YTick', []);
            colormap(gca, cmaps{c});
            %colorbar

            idx = idx + 1;
        end
    end
end




end

