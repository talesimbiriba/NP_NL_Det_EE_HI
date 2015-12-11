clear all
close all

nCores = feature('numCores');
isOpen = matlabpool('size') > 0;
if (isOpen)
    %matlabpool close;    
end


dF = 4;


useJavaProgressMonitor = 1;


% loading indian pines image and ground truth (from 0 (background) to 16).
load Indian_pines_corrected.mat
Yo = indian_pines_corrected(:,:,1:dF:end)/(2^15); %(max(max(max(indian_pines_corrected))));
clear indian_pines_corrected;
load Indian_pines_gt

[n1, n2, L] = size(Yo);


nlMap = zeros(size(indian_pines_gt));

% indian pines has 16 classes. 
% I will divide in 5 sub-images. 
% the first 4 images contain 3 classes each. The last one has 4 classes.
% 
% subimage 1 % numOfEndmembers = 3;
% n9 = 20       % Oats (looks like just one signature)
% n7 = 28       % Grass-pasture-mowed ( looks like 2 sig)
% 
% subimage 2
% n1 = 46       % alfafa (maybe one maybe 2 maybe it has nonlinear
%                           interactions)
% n13 = 205     % Wheat (looks like it has just one sig.)
% n4 = 237      % Corn  (i think it has 3 to 4 different signatures. -
%                           Looking we can see at least 2 different tones)
%
% subimage 3 % numOfEndmembers = 3;
% n16 = 93      % Stone-steel-towers
%
% subimage 4  % numOfEndmembers = 4;
% n15 = 386     % Buildings-Grass-Trees-Drives
%
% subimage 5  % numOfEndmembers = 3;
% n5 = 483      % Grass-Pasture  (Looks like 2 endmembers maybe becouse 2 diff areas)
%
%subimage 6   % numOfEndmembers = 3;
% n8 = 478      % Hay-windrowed  (two)
% n12 = 593     % Soybean-clean (three)
% 
%subimage 7  % numOfEndmembers = 4;
% n3 = 830      % Corn-mintill (2 or 3)
% n6 = 730      % Grass-trees  (looks like 2)
% n10 = 972     % Soybean-notill (2)
% 
% subimage 8  % numOfEndmembers = 4;
% n14 = 1265    % Woods     (2)
% n2 = 1428     % Corn-notill  (three sigs, maybe due to different areas)
% n11 = 2455    % Soybean-mintill  (three sigs, e.g. there is a path in the middle
%                                     of one area.)

clim = [0, 18];

PFA = 0.001;
maxNumOfRuns = 10;

%% processing subimage 1 with classes 9 and 7.

[id1,id2] = find(indian_pines_gt==9 | indian_pines_gt==7);
numOfEndmembers = 3;
figure(1);
imagesc(indian_pines_gt, clim)
hold on;
plot(id2,id1,'ok')
title('Processing black resgions')

% building a smaller image
count =1;
Y = zeros(L,length(id1));
for i=1:length(id1),
    Y(:,count) = squeeze(Yo(id1(i),id2(i),:));
    count = count +1;
end

% iterative joint endmember estimation and nonlinearity detection

[endmembers, nlPixIdx] = iterativeEndmemberEstAndNonlinDetect(Y,numOfEndmembers, PFA, maxNumOfRuns);

for i=1:length(id1),
nlMap(id1(i),id2(i)) = nlPixIdx(i);
end
figure(2)
imagesc(nlMap)
title('Nonlinear map (nonlinearly mixed pixels in red)')

fullM = endmembers;


%return;
%% processing subimage 2 with classes 1, 13 and 4.

[id1,id2] = find(indian_pines_gt==1 | indian_pines_gt==13 | indian_pines_gt ==4);
numOfEndmembers = 3;
figure(1);
clf;
imagesc(indian_pines_gt, clim)
hold on;
plot(id2,id1,'ok')
title('Processing black resgions')

% building a smaller image
count =1;
Y = zeros(L,length(id1));
for i=1:length(id1),
    Y(:,count) = squeeze(Yo(id1(i),id2(i),:));
    count = count +1;
end

% iterative joint endmember estimation and nonlinearity detection

[endmembers, nlPixIdx] = iterativeEndmemberEstAndNonlinDetect(Y,numOfEndmembers, PFA, maxNumOfRuns);

for i=1:length(id1),
nlMap(id1(i),id2(i)) = nlPixIdx(i);
end
figure(2)
imagesc(nlMap)
title('Nonlinear map (nonlinearly mixed pixels in red)')

fullM = endmembers;



%% processing subimage 3 with classes 16 and 13.

[id1,id2] = find(indian_pines_gt==16);
numOfEndmembers = 3;
figure(1);
clf;
imagesc(indian_pines_gt, clim)
hold on;
plot(id2,id1,'ok')
title('Processing black resgions')

% building a smaller image
count =1;
Y = zeros(L,length(id1));
for i=1:length(id1),
    Y(:,count) = squeeze(Yo(id1(i),id2(i),:));
    count = count +1;
end

% iterative joint endmember estimation and nonlinearity detection

[endmembers, nlPixIdx] = iterativeEndmemberEstAndNonlinDetect(Y,numOfEndmembers, PFA, maxNumOfRuns);

for i=1:length(id1),
nlMap(id1(i),id2(i)) = nlPixIdx(i);
end
figure(2)
imagesc(nlMap)
title('Nonlinear map (nonlinearly mixed pixels in red)')

fullM = [fullM, endmembers];


%% processing subimage 4 with classes 15.

[id1,id2] = find(indian_pines_gt==15);
numOfEndmembers = 4;
figure(1);
clf;
imagesc(indian_pines_gt, clim)
hold on;
plot(id2,id1,'ok')
title('Processing black resgions')


% building a smaller image
count =1;
Y = zeros(L,length(id1));
for i=1:length(id1),
    Y(:,count) = squeeze(Yo(id1(i),id2(i),:));
    count = count +1;
end

% iterative joint endmember estimation and nonlinearity detection

[endmembers, nlPixIdx] = iterativeEndmemberEstAndNonlinDetect(Y,numOfEndmembers, PFA, maxNumOfRuns);

for i=1:length(id1),
nlMap(id1(i),id2(i)) = nlPixIdx(i);
end
figure(2)
imagesc(nlMap)
title('Nonlinear map (nonlinearly mixed pixels in red)')

fullM = [fullM, endmembers];


%% processing subimage 5 with classe  5.

[id1,id2] = find(indian_pines_gt==5);
numOfEndmembers = 3;
figure(1);
clf;
imagesc(indian_pines_gt, clim)
hold on;
plot(id2,id1,'ok')
title('Processing black resgions')


% building a smaller image
count =1;
Y = zeros(L,length(id1));
for i=1:length(id1),
    Y(:,count) = squeeze(Yo(id1(i),id2(i),:));
    count = count +1;
end

% iterative joint endmember estimation and nonlinearity detection

[endmembers, nlPixIdx] = iterativeEndmemberEstAndNonlinDetect(Y,numOfEndmembers, PFA, maxNumOfRuns);

for i=1:length(id1),
nlMap(id1(i),id2(i)) = nlPixIdx(i);
end
figure(2)
imagesc(nlMap)
title('Nonlinear map (nonlinearly mixed pixels in red)')

fullM = [fullM, endmembers];


%% processing subimage 6 with classes  8 and 12.

[id1,id2] = find(indian_pines_gt==8 | indian_pines_gt==12);
numOfEndmembers = 3;
figure(1);
clf;
imagesc(indian_pines_gt, clim)
hold on;
plot(id2,id1,'ok')
title('Processing black resgions')


% building a smaller image
count =1;
Y = zeros(L,length(id1));
for i=1:length(id1),
    Y(:,count) = squeeze(Yo(id1(i),id2(i),:));
    count = count +1;
end

% iterative joint endmember estimation and nonlinearity detection

[endmembers, nlPixIdx] = iterativeEndmemberEstAndNonlinDetect(Y,numOfEndmembers, PFA, maxNumOfRuns);

for i=1:length(id1),
nlMap(id1(i),id2(i)) = nlPixIdx(i);
end
figure(2)
imagesc(nlMap)
title('Nonlinear map (nonlinearly mixed pixels in red)')

fullM = [fullM, endmembers];




%% processing subimage 7 with classes 3, 6 and 10.

[id1,id2] = find(indian_pines_gt==3 |indian_pines_gt==6 | indian_pines_gt==10);
numOfEndmembers = 3;
figure(1);
clf;
imagesc(indian_pines_gt, clim)
hold on;
plot(id2,id1,'ok')
title('Processing black resgions')


% building a smaller image
count =1;
Y = zeros(L,length(id1));
for i=1:length(id1),
    Y(:,count) = squeeze(Yo(id1(i),id2(i),:));
    count = count +1;
end

% iterative joint endmember estimation and nonlinearity detection

[endmembers, nlPixIdx] = iterativeEndmemberEstAndNonlinDetect(Y,numOfEndmembers, PFA, maxNumOfRuns);

for i=1:length(id1),
nlMap(id1(i),id2(i)) = nlPixIdx(i);
end
figure(2)
imagesc(nlMap)
title('Nonlinear map (nonlinearly mixed pixels in red)')

fullM = [fullM, endmembers];



%% processing subimage 8 with classes 14, 2 and 11.

[id1,id2] = find(indian_pines_gt==14 | indian_pines_gt==2 | indian_pines_gt ==11);
numOfEndmembers = 3;
figure(1);
clf;
imagesc(indian_pines_gt, clim)
hold on;
plot(id2,id1,'ok')
title('Processing black resgions')

% building a smaller image
count =1;
Y = zeros(L,length(id1));
for i=1:length(id1),
    Y(:,count) = squeeze(Yo(id1(i),id2(i),:));
    count = count +1;
end

% iterative joint endmember estimation and nonlinearity detection

[endmembers, nlPixIdx] = iterativeEndmemberEstAndNonlinDetect(Y,numOfEndmembers, PFA, maxNumOfRuns);

for i=1:length(id1),
nlMap(id1(i),id2(i)) = nlPixIdx(i);
end
figure(2)
imagesc(nlMap)
title('Nonlinear map (nonlinearly mixed pixels in red)')

fullM = [fullM, endmembers];


%%
% closing matlab pool
matlabpool close;    

% plots
% example to plot specific pixel class 
% [id1,id2] =find(indian_pines_gt==3);
% plot(id2,id1,'ok')
% note that the plot inverted the id1 and id2 order! 

% 1) the whole image 
figure;
imagesc(squeeze(Yo(1:end,1:end,[35,15,5]))/max(max(max(Yo))))

% 2) the labels  
[tid1,tid2] = find(nlMap ~= 0);
h = figure;
%imagesc(indian_pines_gt)
%imagesc(indian_pines_gt,[0 18])
imm = indian_pines_gt;
imm(find(nlMap~=0)) = 17;
cmap = jet(18);
cmap(end,:)=zeros(1,3);
imagesc(imm,[0, 18])
colormap(cmap)
% 
% hold on
% for i=1:length(tid1),
%     plot(tid2(i),tid1(i),'k.','linewidth',1)
% end

%filename = '/home/tales/Dropbox/simulations/indianPinesDetectionMapPFA0p001_5RunsOOO.eps';
%print(h, '-depsc', filename);



figure;
imshow(indian_pines_gt)
hold on
h =imshow(nlMap)
set( h, 'AlphaData', .5 )
