function [M, I] = permn(V, N)
% Returns all permutations with repetition of choose N out of V

nV = numel(V) ;
if nV == 0 || N == 0
    M = zeros(nV, N) ;
    I = zeros(nV, N) ;
elseif N == 1
    % return column vectors
    M = V(:) ;
    I = (1:nV).' ;
else
    % NOTE: this is fast but may cause memory issues
    [Y{N:-1:1}] = ndgrid(1:nV) ;
    I = reshape(cat(N+1, Y{:}), [], N) ;
    M = V(I) ;
end


 