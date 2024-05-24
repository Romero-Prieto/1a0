function [beta,Q,M,x,info] = Coef(ReG)

load('data.mat','WS');
ages     = max(WS.index);
tables   = size(WS,1)/ages;
x        = [(0:7:28)/365.25,[(2:1:12),(15:3:24),(36:12:60)]/12]';
LAB      = [string(0);string([(7:7:28)';(2:1:11)';(12:3:24)';(36:12:60)']) + char([kron('d',ones(4,1));kron('m',ones(18,1))])];

n        = x(2:end) - x(1:end-1);
x        = {x,LAB};
set      = [WS.mx_f,WS.mx_m,WS.mx_t];
for j = 1:size(set,2)
    M{1,j}  = reshape(set(:,j),ages,tables)';
    Q{1,j}  = 1 - [ones(1,tables);exp(-cumsum(n.*M{j}',1))]';
    maTCh   = Set(Q{1,j}',x,{ReG(1,1:3),[]},[]);
    maTCh   = maTCh{1};
    xx{1,j} = maTCh;
    clear maTCh 
end

country  = reshape(WS.country,ages,tables);
country  = country(1,:)';
year     = reshape(WS.year,ages,tables);
year     = year(1,:)';
info     = table(country,year);
clear WS country year

Z        = table2array(ReG(1,5));
for k = 1:size(Q,2)
    if isequal(char(table2array(ReG(1,4))),'q')
        q    = Q{1,k};
        q    = q(:,2:end);
        y    = log(reshape(q,[tables*ages,1]));
        xo   = log(xx{1,k});
        clear q
    elseif isequal(char(table2array(ReG(1,4))),'m')
        m    = M{1,k};
        y    = log(reshape(m,[tables*ages,1]));
        xo   = log(xx{1,k});
        clear m
    end
    I       = diag(sparse(1 - ismember(xo,-inf)));
    xo      = recode(xo,-inf,0);
    for z = 0:Z
        X(:,z+1) = xo.^z;
    end    
    X       = kron(sparse(eye(ages)),X);
    I       = min(kron(sparse(eye(ages)),I),diag(sparse(1 - ismember(y,-inf))));
    y       = recode(y,-inf,0);
    B       = pinv(full(X'*I*X))*full(X'*I*y);
    E       = I*y - full(I*X*B);
    B       = reshape(B,[Z + 1 ages])';
    E       = reshape(E,tables,ages);
    tI      = reshape(diag(I),tables,ages);
    [u,s,~] = svd(E'*E./(tI'*tI - Z - 1));
    B       = [B,u];
    beta{k} = {B,ReG(1,:)};
    spl{3}  = B;
    spl{4}  = E;
    s       = diag(s);
    s(1:5)/sum(s)
    clear y X I B E u s z xo
end