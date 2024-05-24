function [maTCh,q,La,E,MSE] = Model(B,xo,x,SET,match,lambda)
ReG     = B{2};
B       = B{1};
tables  = size(xo,1);
Z       = table2array(ReG(1,5));
lAB     = x{2};
x       = x{1};
X       = [xo(:,1).^(0:Z),xo(:,2:end)];
order   = Z + size(xo,2);
B       = B(:,1:order);

if ~iscell(B)
    ages = size(B,1);
    B    = reshape(B',[numel(B) 1]);
    X    = kron(eye(ages),X);
    fx   = reshape(min(exp(X*B),1),[tables ages])';
    clear X B
else
    y  = 0;
    for i = 1:size(X,2)
        y = y + B{i}.*X(:,i)';
    end
    fx   = exp(y);
    clear y X B
end

if isequal(char(table2array(ReG(1,4))),'q')
    l = 1 - [zeros(1,tables);fx];
    d = [l(1:end-1,:) - l(2:end,:);l(end,:)];
    l = cumsum(d,'reverse');
    q = 1 - l./l(1,:);
elseif isequal(char(table2array(ReG(1,4))),'m')
    m = fx;
    n = x(2:end) - x(1:end-1);
    q = 1 - [ones(1,tables);cumprod(exp(-m.*n),1)];
end
 
[maTCh] = Set(q,{x,lAB},SET,xo);
w       = SET{1};
w       = table2array(w(:,4));
w       = w/sum(w);

E       = recode(log(match{1}./maTCh{1}),[-inf,inf],[0;0])*sqrt(w);
MSE     = (recode(log(match{1}./maTCh{1}),[-inf,inf],[0;0]).^2)*w;
La      = (MSE) - lambda.*recode(log(maTCh{2}./match{2}),[-inf,inf],[0;0])*ones(size(lambda,2),1);