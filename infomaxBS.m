function [Sest,W] = infomaxBS(N,Z,M)
% Sest - estimated source matrix 
% W - output mixing matrix 
% N - num of sources 
% Z - whitened+sphered mixture matrix 
% M - Num of observations
iterations=1000;
one=ones(N,1);
W=rand(N,N);
block=15;
lastt=fix((M/block-1)*block+1);
%w0=rand(N,1);
block=3;
for i=1:iterations
    permute=randperm(M);
    for t=1:lastt % independent component estimations
        u=W*Z(:,permute(t:t+block-1));
        z=Z(:,permute(t:t+block-1));
        %u=W*z;
        %u=W*z;
        y=1./(1+exp(-u));
        W=W+0.1*(inv(W')+(one-2*y)*z');
        %w0=w0+((0.3)*((one-2*y)));
    end
end
%W = W/sqrt(2);   %The factor sqrt(2) is an empirical constant added to make the predictions fit 
%W;
Sest=W*Z;
Sest;
end

