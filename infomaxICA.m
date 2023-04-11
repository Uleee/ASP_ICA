function [w,wz,n]=infomaxICA(x)

N=size(x,1); P=size(x,2); M=N;	 %M is dimension of the ICA output
mx=mean(x');%subtracting mean
x=x-(ones(P,1)*mx)';
c=cov(x');%'calculating whitening
wz=2*inv(sqrtm(c));
x=wz*x;%whitening

xx=inv(wz)*x;                     %mean extracted.


w=eye(N); 
count=0; 
perm=randperm(P); 
sweep=0; 
Id=eye(M);
oldw=w; 
olddelta=ones(1,N*M); 
angle=1000; 
change=1000;

 
B=50; L=0.001; F=5000;  %data manipulations parameter 
for I=1:200%I means iteration
    x=x(:,perm);
    sweep=sweep+1; t=1;
    noblocks=fix(P/B);
    BI=B*Id;
    for t=t:B:t-1+noblocks*B
        count=count+B;
        u=w*x(:,t:t+B-1);
        w=w+L*(BI+(1-2*(1./(1+exp(-u))))*u')*w;
    end
    n(1,I)=norm(w-wz, 'fro');
end
 