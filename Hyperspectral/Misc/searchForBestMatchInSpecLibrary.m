function [ libraryIDX ] = searchForBestMatchInSpecLibrary(M, spec_library, method, oneForEach )
%[ libraryIDX ] = searchForBestMatchInSpecLibrary(M, spec_library, method )
%
%   Search the best match for the columns of M in the spectral library spec_library
%
%   Output: libraryIDX is a vector with the indices of spec_library that
%   best match the vectors in M.
%   Inputs:
%       M is a LxR matrix with R column interest vectors 
%       spec_library is a LxD matrix with D spectra
%       method = '1' to use the angle and '0' for squared error. (default =
%       1).
%       oneForEach = '1' to find spectral vectors without replacement (default)
%       oneForEach = '0' to let it repeat the same endmember for different
%       vectors in M.
%
%   Author: Tales Imbiriba
%   November 2015.

if nargin < 4
    oneForEach = 1;
    if nargin < 3
        method=1;
        if nargin <2 
            error('Number of arguments cannot be smaller than 2!');
        end
    end
end

R = size(M,2);
D=size(spec_library,2);

if oneForEach ==1 
    
    if method==1
        theta = zeros(D,R);
        for i=1:R,
            for j=1:D
                theta(j,i) = acos(M(:,i)'*spec_library(:,j)/(norm(M(:,i),2)*norm(spec_library(:,j),2))); 
            end
        end
        for i=1:R
            [val,idx]= min(theta);
            [~,ii] = min(val);
            libraryIDX(ii)=idx(ii);
            theta(idx(ii),:) = ones(1,R)*2;
            theta(:,ii) = ones(D,1)*2;
        end
        
        %[~,libraryIDX] = min(theta);

    else
        sErr = zeros(D,R);
        for i=1:R,
            for j=1:D
                diff = M(:,i)-spec_library(:,j);
                sErr(j,i) = diff'*diff; 
            end
        end
        %[~,libraryIDX] = min(sErr);
        for i=1:R
            [val,idx]= min(sErr);
            [~,ii] = min(val);
            libraryIDX(ii)=idx(ii);
            sErr(idx(ii),:) = ones(1,R)*1e10;
            sErr(:,ii) = ones(D,1)*1e10;
        end
    end

    
    
    
    
    
    
    
else
    if method==1
        theta = zeros(D,R);
        for i=1:R,
            for j=1:D
                theta(j,i) = acos(M(:,i)'*spec_library(:,j)/(norm(M(:,i),2)*norm(spec_library(:,j),2))); 
            end
        end
        [~,libraryIDX] = min(theta);

    else
        sErr = zeros(D,R);
        for i=1:R,
            for j=1:D
                diff = M(:,i)-spec_library(:,j);
                sErr(j,i) = diff'*diff; 
            end
        end
        [~,libraryIDX] = min(sErr);
    end

end





end

