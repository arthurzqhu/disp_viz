function [coeff, train_rsq, test_rsq]=regress_tt(y, x, train_ratio)

% the function looks like:
% [coeff, train_rsq, test_rsq]=regress_tt(y, x, train_ratio). 
% 
% coeff=the regression coefficients
% train_rsq=the rsq of the training data set if coeff is applied to the predictors in
% the training dataset
% test_rsq=the rsq of the test data set if coeff is applied to the predictors in
% the test dataset
% y=depedent variable
% x=independent variable
% train_ratio=the fraction of x-y pairs that will be used to determine the regression
% coefficients

if size(x,2) == size(y,1)
    x=x';
end

dat_sz=length(y); % the size of the entire data set
train_set_sz=ceil(dat_sz*train_ratio); % the size of the training set
train_set_idx=randperm(dat_sz,train_set_sz)'; % randomly pick the indices

test_set_idx=setdiff(1:dat_sz, train_set_idx)'; % indices not picked before will be tested

coeff=regress(y(train_set_idx),x(train_set_idx,:));

y_test_hat=x(test_set_idx,:)*coeff;
y_test=y(test_set_idx);

y_train_hat=x(train_set_idx,:)*coeff;
y_train=y(train_set_idx);

vidx_ts=~isnan(y_test_hat+y_test);
vidx_tr=~isnan(y_train_hat+y_train);

test_rsq=1-sum((y_test_hat(vidx_ts)-y_test(vidx_ts)).^2)/...
   sum((y_test(vidx_ts)-mean(y_test(vidx_ts))).^2);
train_rsq=1-sum((y_train_hat(vidx_tr)-y_train(vidx_tr)).^2)/...
   sum((y_train(vidx_tr)-mean(y_train(vidx_tr))).^2);

end