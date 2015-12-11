function [alpha,b,compt] = subgradient_poly2b_u2(y,MPlus)
%alpha0= FCLS_single_pixel(MPlus,y);
alpha0 = constrainedAbundanceEstimation(y,MPlus);  % Tales: modify since the code was missing.
alpha0(find(alpha0<0))=10^-2;
alpha0(find(alpha0>1))=1-10^-2;
alpha0=alpha0/sum(alpha0);
b0=0;
Niter=3000;
y02=gene_poly2(MPlus,alpha0,b0);
ERR1=norm(y-y02)^2;
ERR1=exp(ERR1);
ERR2=0;
dERR=1;
compt=0;
compt2=0;
alpha=alpha0;
b=b0;
R=size(MPlus,2);
while (dERR>10^-8 && compt<Niter) || (dERR==0 && compt2<20)
    
    k=unidrnd(R);
    for r=1:R
    Mp(:,r)=MPlus(:,r)-MPlus(:,k);
    end
    Mp(:,k)=MPlus(:,k);
    
    I=setdiff(1:R,k) ;
    compt=compt+1;
    e(compt)=ERR1;
    
    %Gradient computation
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %update alpha_1
    
    
%current estimation of the pixel
yp=gene_poly2(MPlus,alpha,b);

    for r=I;
        z=y-MPlus*alpha;
        h=(MPlus*alpha).^2;
        b=(z'*h)/(h'*h);
        
        hp=2*Mp(:,r).* (Mp(:,k)+ Mp*alpha -Mp(:,k)*alpha(k));
        zp=-Mp(:,r);
        bp=( (zp'*h + z'*hp)*(h'*h) - (z'*h)*(2*hp'*h) ) / ( (h'*h)^2 );
        dy= -zp + bp*h + b*hp;   
    df=(yp-y)'*dy;
    df=df/norm(df);

% Computation of the maximal move to respect the constraints
lambda_max=0;
if df~=0
    if df>0
    lambda_max=(alpha(r))/df;
    else
    lambda_max=-alpha(k)/df;
    end
end

 
lambda_M=max(lambda_max,0);
lambda_m=0;

if lambda_M~=0
alpha1=alpha;
alpha2=alpha;

while ((lambda_M-lambda_m) >10^-6)
        lambda1=lambda_m+0.4*(lambda_M-lambda_m);
        lambda2=lambda_m+0.6*(lambda_M-lambda_m);
        alpha1(r)=alpha(r)-lambda1*df;
        alpha1(k)=1-sum(alpha1)+alpha1(k);
        b1=(y-MPlus*alpha1)'*((MPlus*alpha1).^2)/(norm((MPlus*alpha1).^2)^2);
        y1=gene_poly2(MPlus,alpha1,b1);
        alpha2(r)=alpha(r)-lambda2*df;
        alpha2(k)=1-sum(alpha2)+alpha2(k);
        b2=(y-MPlus*alpha2)'*((MPlus*alpha2).^2)/(norm((MPlus*alpha2).^2)^2);

        y2=gene_poly2(MPlus,alpha2,b2);
        err1=norm(y-y1)^2;
        err2=norm(y-y2)^2;
        if err2<err1
            lambda_m=lambda1;
        else
            lambda_M=lambda2;
        end
end  

alpha(r)=alpha(r)-(lambda_M+lambda_m)/2*df;
alpha(k)=1-sum(alpha)+alpha(k);

end
    end

b=(y-MPlus*alpha)'*((MPlus*alpha).^2)/(norm((MPlus*alpha).^2)^2);
%b=max(b,-0.5);
%b=min(b,1); 

ERR2=norm(y-gene_poly2(MPlus,alpha,b))^2;
ERR2=exp(ERR2);
dERR=ERR1-ERR2;
if dERR==0
    compt2=compt2+1;
end
ERR1=ERR2;
end

% figure
 %plot(e)