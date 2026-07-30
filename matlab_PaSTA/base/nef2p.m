function p = nef2p(rX, nef)
%Compute p value that taking account of autocorrelation based on effective
%sample size nef
%   rX: pearson correlation coefficient
%   nef: effective sample size
df=max(0,nef-2); %see Bretherton for using nef-1 rather than nef-2
t=rX*sqrt(df)/sqrt(1-rX^2);
p=2*tcdf(abs(t),df,'upper');
end