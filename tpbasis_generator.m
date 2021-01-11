function [multivariate_polynomials] = tpbasis_generator(dimension, order, x, setsize, multi_idx)
% Creates a basis of multivariate orthonormal polynomials with a tensor
% product
if nargin == 3, setsize = size(x,1); end
if setsize ~= size(x,1)
    error(('tpbasisgenerator:dimension_mismatch(samplesize+stochastic_dimension):transpose_inputvector'));
end

univariate_polynomials = zeros(([setsize, order+1, dimension]));

for n = 1:dimension
    for k=0:order
        L1 = legendre(k, x(:,n)  );
        univariate_polynomials(:,k+1,n) = L1(1,:);
    end
end

cardinality = factorial(order+dimension)/(factorial(order)*factorial(dimension)); %includes k=0-term!!
multivariate_polynomials = zeros(([setsize, cardinality]));

for m = 1:cardinality
    Q = 1;
    for n = 1:dimension
        Q = Q.*univariate_polynomials(:,  multi_idx(m,n)+1 , n);
    end
    multivariate_polynomials(:,m) = Q;
    
end

 