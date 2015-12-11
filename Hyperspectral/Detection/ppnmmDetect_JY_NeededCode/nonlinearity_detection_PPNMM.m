function d = nonlinearity_detection_PPNMM(MPlus,alpha,b,sigma2,PFA)
[L R]=size(MPlus);
CCRB=CCRB_theta(MPlus,alpha,b,sigma2);
sigma2b=CCRB(R+1,R+1);
eta=(-norminv(PFA/2)).^2;
d=0;
if (b^2/sigma2b > eta)
    d=1;
end



