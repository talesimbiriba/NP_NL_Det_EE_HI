function CCRB= CCRB_theta(MPlus,alpha,b,sigma2)

[L R]=size(MPlus);
c=[ones(R,1)' 0 0]';
CRB_u= CRB_theta(MPlus,alpha,b,sigma2);

Q=eye(R+2)-CRB_u*c*pinv(c'*CRB_u*c)*c';

CCRB=Q*CRB_u;