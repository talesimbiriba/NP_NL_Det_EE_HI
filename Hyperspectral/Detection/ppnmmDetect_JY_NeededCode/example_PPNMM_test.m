clear all
close all
clc


%% PPNMM test: linearly mixed pixels
display('Detection procedure for a linearly mixed image')
load image_synth_LMM_2500_bruite.mat
MPlus=MPlus(:,1:3);
[N L]=size(y);
R=3;

% Unmixing procedurels

% for i=1:N
%     i
% [al1 b compt]=subgradient_poly2b_u2(y(i,:)',MPlus);
% al(:,i)=al1;
% B(i)=b;
% sigma2_est_PPNMM(i)=1/L*norm(y(i,:)'-gene_poly2(MPlus,al(:,i),B(i)))^2;
% end
% 
% save result_unmix_PPNMM_LMM.mat al B sigma2_est_PPNMM

% Detection
load result_unmix_PPNMM_LMM.mat al B sigma2_est_PPNMM
PFA=0.05;
for n=1:N
    d(n)=nonlinearity_detection_PPNMM(MPlus,al(:,n),B(n),sigma2_est_PPNMM(n),PFA);
end
PFA_emp=sum(d)/N


%%%% PPNMM test: nonlinearly mixed pixels
display('Detection procedure for a nonlinearly mixed image')
load image_synth_FM_2500_bruite.mat
MPlus=MPlus(:,1:3);
[N L]=size(y);
R=3;

% Unmixing procedure

% for i=1:N
%     i
% [al1 b compt]=subgradient_poly2b_u2(y(i,:)',MPlus);
% al(:,i)=al1;
% B(i)=b;
% sigma2_est_PPNMM(i)=1/L*norm(y(i,:)'-gene_poly2(MPlus,al(:,i),B(i)))^2;
% end
%
% save result_unmix_PPNMM_FM.mat al B sigma2_est_PPNMM

% Detection

load result_unmix_PPNMM_FM.mat al B sigma2_est_PPNMM
PFA=0.05;
for n=1:N
    d(n)=nonlinearity_detection_PPNMM(MPlus,al(:,n),B(n),sigma2_est_PPNMM(n),PFA);
end
PD_emp=sum(d)/N
