function [ gy ] = fetch_expression( y, genes, gene )

gy = y(:, strcmpi(genes, gene) );

end

