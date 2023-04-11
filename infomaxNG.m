function [Sest,W] = infomaxNG(N,Z,M)
iterations=500;
W=ones(N,N)+rand(N,N);
% w0=rand(N,1);
gamma=zeros(N,1);
lr=0.00001;
for ind=1:M % independent component estimations
    z=Z(:,ind);
    for iter=1:iterations
        y=W*z;
        for i=1:N
                gamma(i)=(1-0.2)*gamma(i)+0.2*mean(-tanh(y(i))*y(i)+(1-(tanh(y(i).^2))));
        end
        if gamma<0
            g=tanh(y)-y;
        end
        if gamma>0
            g=-2*tanh(y);
        end
        W=W+(lr*(eye(N)+g*y')*W);
    end
end
Sest=W*Z;
end

