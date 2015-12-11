function [ smallestNumberIndex ] = findSmallestNumbers( numberArray, percentage )
%function [ smallestNumberIndex ] = findSmallestNumbers( numberArray,
%percentage )
%   Return the indices of the 'percentage' (0,1) smallest numbers of the
%   array 'numberArray'.


N = length(numberArray);
Ns = floor(percentage*N);

count = 0;
smallestNumberIndex = zeros(Ns,1);
while count < Ns
    [~,idx] = min(numberArray);
    smallestNumberIndex(count+1) = idx;
    numberArray(idx) = realmax;
    count = count +1;
end


end

