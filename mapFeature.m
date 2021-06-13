function out = mapFeature(X1, X2, deg)
% MAPFEATURE Feature mapping function to polynomial features
%
%   MAPFEATURE(X1, X2) maps the two input features
%   to quadratic features used in the regularization exercise.
%
%   Returns a new feature array with more features, comprising of 
%   X1, X2, X1.^2, X2.^2, X1*X2, X1*X2.^2, etc..
%
%   Inputs X1, X2 must be the same size
%
out=[];
% out = ones(size(X1(:,1)));
for i = 0:deg
    for j = -i:i
        out(:, end+1) = (X1.^(i-j)).*(X2.^j);
    end
end

end