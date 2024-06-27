function [maTCh,q,La,E,MSE] = Model(B,xo,x,SET,match,lambda)
tables  = size(xo,1);
Z       = table2array(B{2}(1,5));
X       = [xo(:,1).^(0:Z),xo(:,2:end)];
order   = Z + size(xo,2);
fx      = exp(B{1}(:,1:order)*X');

if isequal(char(table2array(B{2}(1,4))),'q')
    q = [zeros(1,tables);fx]
elseif isequal(char(table2array(B{2}(1,4))),'m')
    m = fx;
    n = x{1}(2:end) - x{1}(1:end-1);
    q = 1 - [ones(1,tables);cumprod(exp(-m.*n),1)];
end
 
[maTCh] = Set(q,x,SET,xo);
w       = table2array(SET{1}(:,4));
w       = w/sum(w);

E       = recode(log(match{1}./maTCh{1}),[-inf,inf],[0;0])*sqrt(w);
MSE     = (recode(log(match{1}./maTCh{1}),[-inf,inf],[0;0]).^2)*w;
La      = (MSE) - lambda.*recode(log(maTCh{2}./match{2}),[-inf,inf],[0;0])*ones(size(lambda,2),1);
