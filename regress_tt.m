function [coeff, train_rsq, test_rsq] = regress_tt(y, x, train_ratio)

% the function looks like:
% [coeff, train_rsq, test_rsq] = regress_tt(y, x, train_ratio). 
% 
% coeff = the regression coefficients
% train_rsq = the rsq of the training data set if coeff is applied to the predictors in
% the training dataset
% test_rsq = the rsq of the test data set if coeff is applied to the predictors in
% the test dataset
% y = depedent variable
% x = independent variable
% train_ratio = the fraction of x-y pairs that will be used to determine the regression
% coefficients

if size(x,2) == size(y,1)
    x = x';
end


dat_sz = length(y); % the size of the entire data set
train_set_sz = ceil(dat_sz*train_ratio); % the size of the training set
train_set_idx = randperm(dat_sz,train_set_sz)'; % randomly pick the indices

test_set_idx = setdiff(1:dat_sz, train_set_idx)'; % indices not picked before will be tested

coeff = regress(y(train_set_idx),x(train_set_idx,:));

y_test_hat = x(test_set_idx,:)*coeff;
y_test = y(test_set_idx);

y_train_hat = x(train_set_idx,:)*coeff;
y_train = y(train_set_idx);

notnan_idx_test = ~isnan(y_test_hat);
notnan_idx_train = ~isnan(y_train_hat);

test_rsq = 1 - nansum((y_test_hat-y_test).^2)/nansum((y_test(notnan_idx_test)-nanmean(y_test(notnan_idx_test))).^2);
train_rsq = 1 - nansum((y_train_hat-y_train).^2)/nansum((y_train(notnan_idx_train)-nanmean(y_train(notnan_idx_train))).^2);

end