function [ hcube_error ] = computeReconstructionErrorIMG( HCube,HCube_est )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[lin,col, L] = size(HCube);

hcube_error = zeros(lin,col);

% tt = HCube - HCube_est;
for i=1:lin
    for j=1:col
        hcube_error(i,j) = RMSEAndSTDForMatrix(squeeze(HCube(i,j,:)),squeeze(HCube_est(i,j,:)));
    end
end

end

