% Author: Sascha Ranftl, 11.01.2021

 
warning('on', 'MATLAB:singularMatrix')
warning('on', 'MATLAB:nearlySingularMatrix')

seed = 13;
rng(seed);

dimension = 4;
order = 3;
nugget = 10e-10; 
legeps = 10^-6;

N = 1:100;
Ndata = numel(N);

flag_check_surrogateuncertainty = 0;

N_pivots_array{1} = []; % Needs be defined by user


limits_a = limits_a_array{1};
for i = 1:dimension
    atrain(:,i) = []; % Training parameter set samples need be specified as list (matrix with Ndata x dimension)
end


% Get multi-index of basis functions
multi_idx = permn(0:order, dimension);    % Get all the combinations, equals n-choose-k with repetition
multi_idx = multi_idx(sum(multi_idx, 2) <= order,:); % Drop all the combinations for which the combined order exceeds the total order
cardinality = factorial(order+dimension)/(factorial(order)*factorial(dimension)); %includes k=0-term!!
Nq = cardinality;
if size(multi_idx,1) ~= cardinality, warning('Multi-Index does not match cardinality');end % Check wheter everything is alrigh


[Mtrain] = tpbasis_generator(dimension, order, atrain, Ndata, multi_idx);
inv_MtM_train = inv(Mtrain'*Mtrain);
Htrain = inv_MtM_train*Mtrain';


pivots = cell(1,dimension);
for i = 1:dimension
    N_pivots_i = N_pivots_array{idx_select_patient+1}(i);
    pivots{i} = linspace(limits_a(i,1), limits_a(i,2), N_pivots_i);
end

pivots_rescaled =cell(dimension,1); % Pivots generally need be rescaled to range of Legendre polynomials or whatever polynomials are used
for i = 1:dimension
    pivots_rescaled{i} =  calc_rescale_forward(pivots{i}, limits_a(i,:))*(1-legeps);
end


logpriorpdf = calc_prior() ; %  Prior in p(a/d_exp) needs be defined by user
auxlnp = max(logpriorpdf(:));
priorpdf = exp(logpriorpdf-auxlnp);
priorpdf = priorpdf./sum(priorpdf(:));

loglikelihood = calc_likelihood(); %  likelihood  in p(a/d_exp) needs be defined by user
 
logposterior = loglikelihood + logpriorpdf;
posteriorpdf = exp(logposterior -max(logposterior(:)));
NormAll = sum(posteriorpdf(:));
posteriorpdf = posteriorpdf./sum(posteriorpdf(:));

atest = zeros([prod(N_pivots)] ,4);
posterior_reordered =zeros([prod(N_pivots)] ,1);
cnt = 1;

% atest and posterior need be reordered in into single index for efficient
% vectorisation
for i1=1:N_pivots(1)
    for i2=1:N_pivots(2)
        for i3=1:N_pivots(3)
            for i4=1:N_pivots(4)
                atest(cnt, :) = [pivots_rescaled{1}(i1),   pivots_rescaled{2}(i2),  pivots_rescaled{3}(i3),  pivots_rescaled{4}(i4)];
                posterior_reordered(cnt) = posteriorpdf(i1,i2,i3,i4);
                cnt = cnt+1;
            end
        end
    end
end



[Mtest] = tpbasis_generator(dimension, order, atest, N_pivots(1)^4, multi_idx);
MtM_test = Mtest'*Mtest;
inv_MtM_test = inv(MtM_test+nugget*eye(size(MtM_test)));
Htest = inv_MtM_test *Mtest';

my_diagonal_inner_matrix_product_Mtest_Htest = my_diagonal_inner_matrix_product(Mtest*inv_MtM_train, Mtest');

est_z = zeros(numel(T),dimension);
std_z = zeros(numel(T),dimension);

z_data=zeros(Nd,1); % Define z_data here




est_c =  Htrain*z;
Phi = z'*(eye(size(Mtrain*Htrain))-Mtrain*Htrain)*z;
est_DeltaSQ = Phi/(Ndata-Nq-2) ;
est_z_condA = Mtest*est_c;
est_z_condA(est_z_condA<0) = 0;

std_z_condA =  (est_DeltaSQ * my_diagonal_inner_matrix_product_Mtest_Htest + est_z_condA.^2);
if flag_check_surrogateuncertainty
    aux_check = [ est_DeltaSQ * my_diagonal_inner_matrix_product_Mtest_Htest , est_z_condA.^2  ]' * posterior_reordered;
    if aux_check(1)*100>aux_check(2)
        warning('surrogate uncertainties non-negligible')
        sprintf(['Ratio of surrogate uncertainty to rheological uncertainty: ', num2str(aux_check(1)/aux_check(2))])
    end
end
est_zabs = est_z_condA' * posterior_reordered;
std_zabs = sqrt(std_z_condA' * posterior_reordered - est_zabs.^2);


