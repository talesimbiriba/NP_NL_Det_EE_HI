function CRB= CRB_theta(MPlus,alpha,b,sigma2)

[L R]=size(MPlus);
a=alpha;
Ma=MPlus*a;
h=(MPlus*a).^2;

J=zeros(R+2,R+2);

d=zeros(L,R+1);
for r=1:R
    d(:,r)=MPlus(:,r) + 2*b * (MPlus*a).*MPlus(:,r);
end
d(:,R+1)=h;

for i=1:R+1
    for j=1:R+1
        J(i,j)=d(:,i)'*d(:,j);
    end
end
J(R+2,R+2)=L/(2*sigma2);
J=J/sigma2;
% J(R+2,R+2)=L/(sigma2);
% J=J/(2*sigma2);
CRB=pinv(J);

