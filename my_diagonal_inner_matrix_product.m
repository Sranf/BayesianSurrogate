function [diagonal] = my_diagonal_inner_matrix_product(M, H);

[mi, mj] = size(M);
diagonal = zeros([mi,1]);
for i=1:mi
    diagonal(i) =  M(i,:)*H(:,i);
end

