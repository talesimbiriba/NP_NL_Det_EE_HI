function [ HCube ] = matrixToHCube( Y, nRow, nCol)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


[L,N] = size(Y);

if (N ~= nRow*nCol)
    error('Size missmatch! Check number of columns and rows.')
end



HCube = zeros(nRow, nCol, L);
count =1;

for i=1:nRow, 
    for j=1:nCol,
        HCube(i,j,:) = Y(:,count);
        count = count +1;
    end
end


end

