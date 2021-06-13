function y = wmean(A,w)
% y = wmean(A, w). A = random variable. w = weights.
if any(w<0)
    w(w<0)=0;
    warning('negative weights are set to 0')
end

if sum(w) ~= 1
    w = w/sum(w);
end

A(isnan(A)) = 0;
w(isnan(A)) = 0;

if size(A,1) > size(A,2)
    A = A';
end

if size(A,2) == size(w,1)
    y = A*w;
elseif size(A,2) == size(w,2)
    y = A*w';
else
    error('weights must have the same dimension as the data')
end




end