function [ gratestNumberIndex ] = findGratestNumbers( numberArray, percentage )
%function [ gratestNumberIndex ] = findGratestNumbers( numberArray,
%percentage )
%   Return the indices of the 'percentage' (0,1) gratest numbers of the
%   array 'numberArray'.


N = length(numberArray);
Ns = floor(percentage*N);

count = 0;
gratestNumberIndex = zeros(Ns,1);
while count < Ns
    [~,idx] = max(numberArray);
    gratestNumberIndex(count+1) = idx;
    numberArray(idx) = -realmax;
    count = count +1;
end


end

