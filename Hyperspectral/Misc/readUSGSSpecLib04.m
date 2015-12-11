function [ M, labels] = readUSGSSpecLib04(fileName)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
 
disp('Warning: The band 206 is missing in this spectral library!')


% important strings 

dummy = 'ZZZZZZZZZZZZZZZZZZZZ';
maxZZZCount = 5;

nameDelimiter1 = 'ititl:';       % ititl: Mineral name =? AVIRIS...
nameDelimiter2 = '=';

lineBeforeData = 'tempd:';               %line before vector
lineAfterData = 'stat:complete';

fid = fopen(fileName,'r');
s='';
M = zeros(223,1);
labelCount=1;
dataCount=1;
count =0;
countZZZ = 0;
while ~feof(fid)
    line = fgetl(fid);
    if isempty(line), 
        continue, 
    end
    
    [tk,remain] = strtok(line, ' ');
    
    switch tk 
        case dummy
            countZZZ=countZZZ+1;   
            line = fgetl(fid);
            continue;
        case nameDelimiter1
            if countZZZ >=maxZZZCount
                [~,labels{labelCount}] = strtok(line,' ');
                labelCount = labelCount +1;
            end
        case lineBeforeData
            cc=1;
            %ttvec = 0; 
            ttvec=zeros(223,1); % if we now the vector sizes.
            while 1
                line=fgetl(fid);
                [tk,remain]=strtok(line,' ');
                if strcmp(tk,lineAfterData), break, end
                if strcmp(tk,'206'), continue, end
               ttvec(cc) =  str2num(remain);
               cc = cc+1;
            end
            
            M(:,dataCount) = ttvec;
            dataCount = dataCount +1;
        otherwise
            continue;
    end
    
%     if strcmp(line,dummy), 
%         countZZZ=countZZZ+1;   
%         line = fgetl(fid);
%         continue;
%     end
    
%     if strcmp(line,nameDelimiter1)&& countZZZ >=maxZZZCount, 
%         [~,labels{labelCount}] = strtok(line,' ');        
%     end
%     
%     if strcmp(line, lineAfterData)
%         continue 
%     end
    
%     if strcmp(line,'MAMAE'), 
%         break, 
%     end
%     
%     if ~mod(count,500) 
%         line
%     end

%     s = strvcat(s,line);
%     count = count +1;
end

fclose(fid);
