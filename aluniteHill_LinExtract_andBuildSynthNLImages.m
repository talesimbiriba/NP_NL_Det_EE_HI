clear all;
close all;


load cuprite_ROI_tales.mat
load wavelength.txt

[nRow, nCol, L] = size(cuprite_hsimg);

% 
% cd '~/Work/matlab/gpml-matlab-v3.3-2013-10-19/'
% startup;
% 
% cd /home/tales/Work/matlab/sedumi-master;
% install_sedumi;





SNR=30;
nlDegree = 0.3;
R=3;
df = 1;
scaleFactor = 1.1;
%startBand = 139;
%plotBands = round([1,20,30]);
plotBands = round([30,70,100]/df);
%plotBands = round([10,20,30]/df);
numOfNlPixels = floor(0.5*nRow*nCol);
nlmodel=1; % 1 for GBM and 2 for PNMM



% convert HCube to matrix 
%cuprite_hsimg = cuprite_hsimg(:,:,startBand:end);
cuprite_hsimg = cuprite_hsimg(:,:,1:df:end);
Y = hCubeToMatrix(cuprite_hsimg);
% normalize
Yn = Y./(scaleFactor*max(max(Y)));

[L,N] = size(Yn);

% extract R=3 endmembers using MVES
M = MVES(Yn,R,1);
%M = hyperVca(Yn,R);

% estimate abundances using fully constrained linear model
a = constrainedAbundanceEstimation(Yn, M, 0.001);

% convert to cube 
acube = matrixToHCube(a,nRow,nCol);

% plot abundances 
figure
%subplot(3,1,1)
%imagesc(acube(:,:,1))
%imshow(flipdim(imrotate(acube(:,:,2),-90),2))
imagesc(flipdim(imrotate(acube(:,:,1),-90),2))
colormap 'gray'

figure
%subplot(3,1,2)
%imagesc(acube(:,:,2))
imagesc(flipdim(imrotate(acube(:,:,2),-90),2))
colormap 'gray'

figure
%subplot(3,1,3)
%imagesc(acube(:,:,3))
imagesc(flipdim(imrotate(acube(:,:,3),-90),2))
colormap 'gray'

% create sinthetic linear mixed image
Ysin = M*a;

% convert to HCube
Ycube = matrixToHCube(Ysin,nRow,nCol);

% plot original and LMM synth image
figure
%imagesc(squeeze(cuprite_hsimg(:,:,[10 30 60]))/(max(max(max(cuprite_hsimg(:,:,[10 30 60]))))))
%imagesc(squeeze(flipdim(imrotate(cuprite_hsimg(:,:,[95 105 125])/(1.1*max(max(max(cuprite_hsimg)))),-90),2)));
imagesc(squeeze(flipdim(imrotate(cuprite_hsimg(:,:,plotBands)/(scaleFactor *max(max(max(cuprite_hsimg)))),-90),2)));
figure;
%imagesc(squeeze(flipdim(imrotate(Ycube(:,:,[10,30,60]),-90),2)));
imagesc(squeeze(flipdim(imrotate(Ycube(:,:,plotBands),-90),2)));


% select half of the pixels randomly
%idx = datasample(1:448,224,'Replace',false);
idx = datasample(1:N,numOfNlPixels,'Replace',false);
%idx = sort(idx);
% idx=1:224;
Ar = (1-nlDegree)/nlDegree;
k = sqrt(Ar/(1 + Ar));


linPart = M*a(:,idx);
   
    
%linPart = M*a(:,idx(i));
       
nonlinPart = zeros(L,length(idx));

if nlmodel==1
%compute the nonlinear contributions  
for i=1:length(idx)
    count = 1;
    for t=1:R-1,
        for j=1+t:R;
            M_bilin(:,count) = M(:,t).*M(:,j);
            a_bilin(count,1) = a(t,idx(i))*a(j,idx(i));
            count = count+1;
        end
    end
    nonlinPart(:,i) = (M_bilin*a_bilin);   
end


% Sample mean of the Energies
El = trace(linPart'*linPart)/N;
Enl = trace(nonlinPart'*nonlinPart)/N;
Elnl = trace(linPart'*nonlinPart)/N;
% El = diag(linPart'*linPart)/N;
% Enl = diag(nonlinPart'*nonlinPart)/N;
% Elnl = diag(linPart'*nonlinPart)/N;

% compute the mean K
%K = ( -2*gamma*Elnl + sqrt( 4*(gamma^2)*(Elnl^2) -4*El*( (gamma^2)*Enl-El) ) )/ (2*El);
B = (1-k^2)*El;

gamma = (-2*k*Elnl + sqrt(4*(k^2)*Elnl^2 + 4*B*Enl))/(2*Enl);
%gamma = (-2*k*Elnl + sqrt(4*(k^2)*Elnl.^2 + 4*B.*Enl))./(2*Enl);

% generate the noisy signal.
%Y = K*linPart + gamma*nonlinPart + randn(numOfBands,numOfSamples)*sqrt(noiseVar);    
Ysin(:,idx) = k*linPart + gamma*nonlinPart;    
%Ysin(:,idx) = k*linPart + (diag(gamma)*nonlinPart')';    
% compute the nonlinearity degree.
%nlDegree =  (2*K*gamma*Elnl + (gamma^2)*Enl)/( (K^2)*El + 2*K*gamma*Elnl + (gamma^2)*Enl );    

else


disp('PNMM')
p = 0.3;    
linPart = M*a;
nlPart = (M*a).^p;
El = trace(linPart'*linPart)/N;
Enl = trace(nlPart'*nlPart)/N;
Elnl = trace(linPart'*nlPart)/N;


B = (1-k^2)*El;
gamma = (-2*k*Elnl + sqrt(4*(k^2)*Elnl^2 + 4*B*Enl))/(2*Enl);

%Y = k*linPart + gamma*nlPart + randn(L,N)*sqrt(noiseVar);
Y = k*linPart + gamma*nlPart;


%nld = trace(sconst*2*k*linPart'*nonlinPart + (sconst^2)*nonlinPart'*nonlinPart)/trace((k^2)*linPart'*linPart + sconst*2*k*linPart'*nonlinPart + (sconst^2)*nonlinPart'*nonlinPart);
nld =  (2*k*gamma*Elnl + (gamma^2)*Enl)/( (k^2)*El + 2*k*gamma*Elnl + (gamma^2)*Enl );
end


[L,N]= size(Ysin);

%noisePw = (norm(Ysin,'fro')^2/(N*L))/(10^(SNR/10));
noisePw = noisePowerSNR(M,SNR);

Ysin = Ysin + randn(L,N)*sqrt(noisePw);

Ycube = matrixToHCube(Ysin,nRow,nCol);
figure;
imagesc(squeeze(flipdim(imrotate(Ycube(:,:,plotBands),-90),2)));


simpHand = figure;
draw3dSimplex(M,simpHand)
hold on;
ymnl = M'*Ysin;
scatter3(ymnl(1,:),ymnl(2,:),ymnl(3,:),'.k')

Mvca = hyperVca(Ysin, R);
Mmves = MVES(Ysin,R,0);
mm1 = Mvca'*M;
scatter3(mm1(:,1),mm1(:,2),mm1(:,3),'og','fill')
mm1 = Mmves'*M;
scatter3(mm1(:,1),mm1(:,2),mm1(:,3),'or','fill')


[endmembers, nlPixIdx] = iterativeEndmemberEstAndNonlinDetect(Ysin,R, 0.1, 10, 0.9 ,1, M);


figHandsimp2 = figure;
draw3dSimplex(M,figHandsimp2)
hold on
draw3dSimplex(endmembers,figHandsimp2,M,'b')
scatter3(ymnl(1,:),ymnl(2,:),ymnl(3,:),'.k')

tt =M - repmat(endmembers(:,1),1,3);
[~,i1] = min(diag(tt'*tt));
tt =M - repmat(endmembers(:,2),1,3);
[~,i2] = min(diag(tt'*tt));
tt =M - repmat(endmembers(:,3),1,3);
[~,i3] = min(diag(tt'*tt));
endmembers = [endmembers(:,i1),endmembers(:,i2),endmembers(:,i3)];

tt =M - repmat(Mmves(:,1),1,3);
[~,i1] = min(diag(tt'*tt));
tt =M - repmat(Mmves(:,2),1,3);
[~,i2] = min(diag(tt'*tt));
tt =M - repmat(Mmves(:,3),1,3);
[~,i3] = min(diag(tt'*tt));
Mmves = [Mmves(:,i1),Mmves(:,i2),Mmves(:,i3)];


tt =M - repmat(Mvca(:,1),1,3);
[~,i1] = min(diag(tt'*tt));
tt =M - repmat(Mvca(:,2),1,3);
[~,i2] = min(diag(tt'*tt));
tt =M - repmat(Mvca(:,3),1,3);
[~,i3] = min(diag(tt'*tt));
Mvca = [Mvca(:,i1),Mvca(:,i2),Mvca(:,i3)];


for i=1:R
    figure;
%     plot(wavelength(174:end,1), M(:,i),'linewidth',2)
plot(M(:,i),'linewidth',2)
    hold on
%     plot(wavelength(174:end,1), endmembers(:,i),'r','linewidth',2)
%     plot(wavelength(174:end,1), Mmves(:,i),'g','linewidth',2)
%     plot(wavelength(174:end,1), Mvca(:,i),'k','linewidth',2)
plot(endmembers(:,i),'r','linewidth',2)
    plot(Mmves(:,i),'g','linewidth',2)
    plot(Mvca(:,i),'k','linewidth',2)
    legend('Endmember','Estimated proposed','Estimated MVES','Estimated VCA')    
end


[Y_est_lin, Aest_lin] = estimateLinearModel(Ysin, endmembers, 1);
[Aest_skp,~,Y_est_skp] = tskHype(Ysin,endmembers,[],1000,10);

Aest_detect = zeros(R,N);
for i=1:N,
    if nlPixIdx(i) == 1
        Aest_detect(:,i) = Aest_skp(:,i);
    else
        Aest_detect(:,i) = Aest_lin(:,i);
    end
end

norm(M-endmembers,'fro')^2
norm(M-Mmves,'fro')^2
norm(M-Mvca,'fro')^2

[rmse_d,std_d] = RMSEAndSTDForMatrix(a, Aest_detect)
[rmse_lin, std_lin] = RMSEAndSTDForMatrix(a, Aest_lin)
[rmse_skp, std_skp] = RMSEAndSTDForMatrix(a, Aest_skp)