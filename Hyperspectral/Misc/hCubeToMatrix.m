function [ hypMatrix ] = hCubeToMatrix(hCube)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


[n1,n2,L] = size(hCube);


hypMatrix = zeros(L,n1*n2);

count = 1;

for i=1:n1,
    for j=1:n2
        hypMatrix(:,count) = squeeze(hCube(i,j,:));
        count = count + 1;
    end
end

end

