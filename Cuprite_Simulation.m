


clear all;
close all;

% load cuprite scene
load hcube_tales.mat

BANDS = BANDS; % to avoid matlab complain 

df = 2; % decimation factor 

% % taking a smaller portion of the image.
% smallIMG = hcube(160:250,110:190,:);

% virtual dimension 
R = HFC(smallIMG,0.0001);

Y = hCubeToMatrix(smallIMG);
Y = Y(1:df:end,:);
Mvca = hyperVca(Y,R);
Mmvs = MVES(Y,R,0);

% estimate endmembers and detect nonlinearly mixed pixels usign the
% proposed method (Algorithm 1)
PFA = 0.0001; rf = 0.7; Nmax = 10;
[M_ours,nlPix]= iterativeEndmemberEstAndNonlinDetect(Y,R,PFA,10,rf,1);

% estimate endmembers and detect nonlinearly mixed pixels usign the
% Algorithm 1 with the Robust LS based detector.
rf2 = 1.2;
[M_rls,nlPix_rls]= iterativeEndmemberEstAndNonlinDetect_RLS(Y,R,PFA,Nmax,rf2,1);

% unmix using skyphe and FCLS
[Askp,~,Yskp] = tskHype(Y,M_ours);
[Ylin, Alin] = estimateLinearModel(Y, M_ours, 1);




% Load USGS Spectral Library
load splib04.mat
% selecting 17 minerals known to exist in the cuprite mining field.
% labels
cupriteLabels = mineral_labels(sort([374 311 321 235 35 295 303 296 18 67 71 234 300 425 135 20 239]));
% spectra
cupriteLib = spec_library(:,sort([374 311 321 235 35 295 303 296 18 67 71 234 300 425 135 20 239]));

% find the best match between the USGS spectra and the estimated 
libIDX = searchForBestMatchInSpecLibrary(M_ours,cupriteLib(BANDS(1:df:end),:),1);
Mtrue = cupriteLib(BANDS(1:df:end),libIDX);
Mlabels = cupriteLabels(libIDX);


% reordering colunms of estimated endmember matrices 
ord =searchForBestMatchInSpecLibrary(Mvca,Mtrue,1);
Mvca_reord = Mvca(:,ord);
ord =searchForBestMatchInSpecLibrary(Mmvs,Mtrue,1);
Mmves_reord = Mmvs(:,ord);

ord =searchForBestMatchInSpecLibrary(M_rls,Mtrue,1);
M_rls_reord = M_rls(:,ord);

BANDS2 = BANDS(1:df:end);

figFontSize=20;


%%%%% PLOT ENDMEMBERS %%%%%%
% the scale may be needed to change. It depends on the order of endmembers extracted by MVES and VCA 
figure;
set(gca,'FontSize',figFontSize)
plot(BANDS2,Mtrue(:,1),'linewidth',2)
hold on
plot(BANDS2,M_ours(:,1)*2.5,'r','linewidth',2)
%plot(BANDS,Mvca(:,4)*1.2,'g','linewidth',2)
plot(BANDS2,M_rls_reord(:,1)*8,'g','linewidth',2)
%plot(BANDS2,Mmves_reord(:,1)*1.7,'g','linewidth',1)
ylim([0.05,0.9])
xlim([0 224])
%legend('USGS', 'Estimated Proposed','VCA','MVES','location','SouthEast')
legend('USGS', 'Estimated Proposed','Estimated LS','location','SouthEast');
xlabel('Bands')
ylabel('Reflectance')
title(Mlabels{1})
% savefig('em1.fig')
% print('em1','-depsc')

figure;
set(gca,'FontSize',figFontSize)
plot(BANDS2,Mtrue(:,2),'linewidth',2)
hold on
plot(BANDS2,M_ours(:,2)*1.05,'r','linewidth',2)
plot(BANDS2,M_rls_reord(:,2)*0.9,'g','linewidth',2)
%plot(BANDS,Mvca(:,3)*1.55,'g','linewidth',2)
%plot(BANDS2,Mmves(:,2)*1,'g','linewidth',1)
xlim([0 224])
%legend('USGS', 'Estimated Proposed','VCA','MVES','location','SouthEast')
legend('USGS', 'Estimated Proposed','Estimated LS','location','South');
xlabel('Bands')
ylabel('Reflectance')
title(Mlabels{2})
% savefig('em2.fig')
% print('em2','-depsc')

figure;
set(gca,'FontSize',figFontSize)
plot(BANDS2,Mtrue(:,3),'linewidth',2)
hold on
plot(BANDS2,M_ours(:,3)*1.8,'r','linewidth',2)
plot(BANDS2,M_rls_reord(:,3)*1.55,'g','linewidth',2)
%plot(BANDS,Mvca(:,2)*1.7,'g','linewidth',2)
%plot(BANDS2,Mmves(:,3)*1,'g','linewidth',1)
ylim([0 1.2])
xlim([0 224])
%legend('USGS', 'Estimated Proposed','VCA','MVES','location','SouthEast')
legend('USGS', 'Estimated Proposed','Estimated LS','location','South');
xlabel('Bands')
ylabel('Reflectance')
title(Mlabels{3})
% savefig('em3.fig')
% print('em3','-depsc')


figure;
set(gca,'FontSize',figFontSize)
plot(BANDS2,Mtrue(:,4),'linewidth',2)
hold on
plot(BANDS2,M_ours(:,4)*0.9,'r','linewidth',2)
plot(BANDS2,M_rls_reord(:,4)*0.9,'g','linewidth',2)
%plot(Mvca(:,1)*1.25,'g','linewidth',2) 
%plot(Mmves(:,2)*0.86,'r','linewidth',2)
ylim([0 0.6])
xlim([0 224])
% legend('USGS', 'Estimated Proposed','VCA','MVES','location','SouthEast')
legend('USGS', 'Estimated Proposed','Estimated LS','location','South');
xlabel('Bands')
ylabel('Reflectance')
title(Mlabels{4})
% savefig('em4.fig')
% print('em4','-depsc')


figure;
set(gca,'FontSize',figFontSize)
plot(BANDS2,Mtrue(:,5),'linewidth',2)
hold on
plot(BANDS2,M_ours(:,5)*1.9,'r','linewidth',2)
plot(BANDS2,M_rls_reord(:,5)*1.6,'g','linewidth',2)
% plot(Mvca(:,5)*2.2,'g','linewidth',2)
% plot(Mmves(:,5)*1.44,'r','linewidth',2)
xlim([0 224])
ylim([0 0.9])
%legend('USGS', 'Estimated Proposed','VCA','MVES','location','SouthEast')
legend('USGS', 'Estimated Proposed','Estimated LS','location','SouthEast');
xlabel('Bands')
ylabel('Reflectance')
title(Mlabels{5})
% savefig('em5.fig')
% print('em5','-depsc')

%ord = searchForBestMatchInSpecLibrary(Mmves,Mtrue,1)
%Mvca_reord = Mvca(:,[4 3 2 5]);
% reorder:  ord = searchForBestMatchInSpecLibrary(Mmves,Mtrue,1)
%Mmves_reord = Mmves(:,ord);

err_our = zeros(R,1);
err_rls = zeros(R,1);
err_vca = zeros(R,1);
err_mves = zeros(R,1);

[L,N] = size(Y);

for i=1:R
    err_our(i) = acos((M_ours(:,i)'*Mtrue(:,i))/(norm(M_ours(:,i),2)*norm(Mtrue(:,i),2)));
    err_rls(i) = acos((M_rls_reord(:,i)'*Mtrue(:,i))/(norm(M_rls_reord(:,i),2)*norm(Mtrue(:,i),2)));
    err_vca(i) = acos((Mvca_reord(:,i)'*Mtrue(:,i))/(norm(Mvca_reord(:,i),2)*norm(Mtrue(:,i),2)));
    err_mves(i) = acos((Mmves_reord(:,i)'*Mtrue(:,i))/(norm(Mmves_reord(:,i),2)*norm(Mtrue(:,i),2)));
       
end


fprintf('\n')
for i=1:R
    fprintf('%s & %1.4f & %1.4f & %1.4f & %1.4f  \\\\\n', Mlabels{i}(1:25), err_our(i), err_rls(i), err_vca(i), err_mves(i));
end


Ymix = Ylin;
Amix = Alin;
for i=1:length(nlPix)
    if nlPix(i)==1;
        Amix(:,i) = Askp(:,i);
        Ymix(:,i) = Yskp(:,i);
    end
end

acube= matrixToHCube(Askp,Lines,Columns);
%acube= matrixToHCube(Amix,91,81);
for i=1:R
    figure;
    imagesc(acube(:,:,i))
    caxis([0 1])
    set(gca,'XtickLabel',[],'YtickLabel',[]);
    set(gca,'visible','off')
    %figName = ['abundaces_em',num2str(i)];
    %print(figName,'-depsc');
end

[~,~,Yvca] = tskHype(Y,Mvca);
fprintf('\n RMSE_{prop} = %.4f\n',RMSEAndSTDForMatrix(Ymix,Y))
fprintf('\n RMSE_{VCA} = %.4f\n',RMSEAndSTDForMatrix(Yvca,Y))



for i=1:N,
    tt(i) = sqrt((Y(:,i)-Ymix(:,i))'*(Y(:,i)-Ymix(:,i))/L);
end
%size(tt)
figure;
tcube= matrixToHCube(tt,Lines,Columns);
imagesc(tcube)
caxis([0 0.03])
title('RMSE Proposed')

for i=1:N,
    tt(i) = sqrt((Y(:,i)-Yvca(:,i))'*(Y(:,i)-Yvca(:,i))/L);
end
%size(tt)
figure;
tcube= matrixToHCube(tt,Lines,Columns);
imagesc(tcube)
caxis([0 0.03])
title('RMSE VCA')
