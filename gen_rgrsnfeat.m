function [X_train, X_cv] = gen_rgrsnfeat(feature_mtx,train_ratio,randomness)

% generates regression features for training and cross validation.
% Needs input of a m x n matrix of the raw features (feature_mtx), m = length of each
% matrix, n = the number of features. *** this matrix cannot be transposed ***
% do not need to include an array of 1s for bias. it will be included
% automatically.
% Also needs input of the ratio of original feature set that one wants use
% as training (train_ratio), the rest goes into CV, or shared by CV and testing.
% set randomness to 1 picks the training set randomly to avoid bias,
% default is 1.

% Outputs two structures, each containing the processed feature matrices,
% normalized feature matrices, mean, and standard deviation.

if ~exist('randomness', 'var')
    randomness = 1;
end

samp_size = size(feature_mtx,1);
samp_idx = 1:samp_size;
if randomness == 1
    X_train.idx = randperm(samp_size, ceil(samp_size*train_ratio))';
else
    X_train.idx = [1:ceil(samp_size*train_ratio)]';
end

X_cv.idx = setdiff(samp_idx, X_train.idx)';

X_train.rawmtx = [ones(length(X_train.idx),1) feature_mtx(X_train.idx,:)];
[X_train.normmtx, X_train.mu, X_train.sigma] = featureNormalize(X_train.rawmtx);

if size(X_train.rawmtx,1) < size(feature_mtx,1)
    X_cv.rawmtx = [ones(length(X_cv.idx),1) feature_mtx(X_cv.idx,:)];
    [X_cv.normmtx, X_cv.mu, X_cv.sigma] = featureNormalize(X_cv.rawmtx);
end


end