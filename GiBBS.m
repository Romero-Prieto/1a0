function [sample] = GiBBS(y,X,seed)

rng(seed);
[n,k]  = size(X);
yX     = y'*X;
XX     = X'*X;
 
priorB = [zeros(k,1),ones(k,1)*10^-4]; 
priorT = [1,1]*10^-2;
B      = zeros(k,1);
T      = 1*10^-4;
sample = [B',T'];
 
for i = 1:12500
    for j = 1:k
        zx      = yX(:,j) - B'*XX(:,j) + XX(j,j)*B(j);
        S       = pinv(priorB(j,2) + T*XX(j,j));
        R       = priorB(j,2)*priorB(j,1) + T*zx;
        B(j)    = normrnd(S*R,sqrt(S));
        pM(i,j) = S*R;
        pS(i,j) = S;
        %clear S R
    end
    e               = y - X*B;
    T               = gamrnd(priorT(1) + n/2,1/(priorT(2) + e'*e/2));
    sample(i + 1,:) = [B',T'];
end
sample = sample(2502:end,:);