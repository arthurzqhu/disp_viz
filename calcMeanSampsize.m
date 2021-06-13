function [mean, sampsize] = calcMeanSampsize(raw_arr, logic_arr)

if nargin == 1
    logic_arr = true(length(raw_arr),1);
end

filtered_arr = raw_arr(logic_arr);
mean = nanmean(filtered_arr);
sampsize = length(find(~isnan(filtered_arr)));

end