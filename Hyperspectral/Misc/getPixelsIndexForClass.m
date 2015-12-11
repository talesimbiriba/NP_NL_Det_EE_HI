function [ idx ] = getPixelsIndexForClass(gt, class)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

idx = zeros(1,2);

[n,p] = size(gt);

count = 1;
for i=1:n,
    for j=1:p,
        if (gt(i,j) == class)
            idx(count,:) = [i j];
            count = count + 1;
        end
    end
end



end

