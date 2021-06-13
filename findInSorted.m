function [a,b] = findInSorted(arr,min,max)
    % input: arr, min, (max)
    % output: a, b
    % output two boundary indices (a, b) of an interval (min,max] in a sorted array 'arr' using binary search
    % function will return the two indices adjacent to min if max is not specified
    
    if nargin < 2
        error('Please specify the sorted array and/or the minimum and/or maximum value')
    end
    
    L = 1; % lower bound
    R = length(arr); % upper bound
    m = floor((L+R)/2); % mid value
    
    while R>L+1 % while R and L is more than 1 apart
        if arr(m)<=min
            L = m;
            m = ceil((L+R)/2); % ceil() is used here and floor() below to avoid skipping indices
        elseif arr(m)>min 
            R = m;
            m = floor((L+R)/2);
        else
            break
        end
    end
    
    a=R; % set a to the upper bound
    

    if nargin == 2
        b = a;
        a = a-1;
    elseif nargin == 3
        L = a;
        R = length(arr);
        m = floor((L+R)/2);
    
        while R>L+1
            if arr(m)<=max
                L = m;
                m = ceil((L+R)/2);
            elseif arr(m)>max
                R = m;
                m = floor((L+R)/2);
            else
                break
            end
        end
    
    b=L;
    end
    
    if  isnan(min) || any(isnan(arr))
        a=nan;
        b=nan;
    end
   
end