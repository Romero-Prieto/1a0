function [XO,q,MSE,maTCh,dLa,lAMBda] = Match(match,B,x,SET,xo,R)
% Minimizes the MSE for a given set of non-linear equality constraints
tables            = max(size(match{1},1),size(match{1},1)); % Number of LTs to be matched
order             = max(size(match{2},2),size(xo,2));
ReG               = B{2};
Z                 = table2array(ReG(1,5));
B{1}              = B{1}(:,1:order + Z);                    % Set the coefficients to be used
lambda            = ones(size(match{2}));
iNdeX             = (1:tables)';
mIn               = zeros(tables,1);                        % Define an indicator of convergent solution
d                 = 0.0001;                                 % Define the proximity value to calculate a numerical derivative
XA                = 1.25*xo;
XO                = NaN(size(xo));
lAMBda            = lambda;
mATCh             = match;
tABLes            = tables;
Betas             = B;

for r = 1:R
    dLa          = ModelD(B,xo,x,SET,match,lambda);
    M            = reshape(dLa',[numel(dLa),1]);
    both         = [xo,lambda];
    for h = 1:size(both,2) % numerical derivative
        bothU      = both;
        bothU(:,h) = bothU(:,h) + d;
        dLU        = ModelD(B,bothU(:,1:size(xo,2)),x,SET,match,bothU(:,size(xo,2) + 1:end));
        bothL      = both;
        bothL(:,h) = bothL(:,h) - d;
        dLL        = ModelD(B,bothL(:,1:size(xo,2)),x,SET,match,bothL(:,size(xo,2) + 1:end));        
        M(:,h+1)   = reshape((dLU - dLL)'/(2*d),[numel(dLa),1]);
        clear dLU dLL bothU bothL
    end
    M            = mat2cell(M,ones(tables,1)*(size(M,2) - 1),[1 size(M,2) - 1]);    
    for k = 1:tables % solving the set of equations for each life table
        test           = sum(sum(cell2mat(M(k,:))));
        test2          = max(abs(XA(iNdeX(k),:) - xo(k,:)));
        XA(iNdeX(k),:) = xo(k,:);
        if max(abs(M{k,1})) > 10^-10 && ~isnan(test) && ~isinf(test) && test2 > 10^-10 % tolerance of errors
            [u,W,v]            = svd(M{k,2});
            change(k,:)        = -real(v*pinv(W)*u'*M{k,1})';
        else
            mIn(k)             = 1;
            change(k,:)        = zeros(1,size(both,2));
            XO(iNdeX(k),:)     = xo(k,:);
            lAMBda(iNdeX(k),:) = lambda(k,:);
        end 
    end
    clear k test
    
    bothK              = both + change;
    xo                 = bothK(:,1:size(xo,2));
    lambda             = bothK(:,size(xo,2) + 1:end);
    clear both bothK change
    
    if max(mIn) == 1
        xo     = xo(mIn ~= 1,:);
        lambda = lambda(mIn ~= 1,:);
        match  = {match{1}(mIn ~= 1,:),match{2}(mIn ~= 1,:)};
        if iscell(B{1})
            for h = 1:numel(B{1})
                B{1}{h} = B{1}{h}(:,mIn ~= 1);
            end
        end
                
        iNdeX  = iNdeX(mIn ~= 1);
        tables = size(xo,1);
        mIn    = zeros(tables,1);
    end
    clc;
    table([{'iteration';'Life Tables';'Convergence';'%'},{r;tABLes;tABLes-tables;floor((1 - tables/tABLes)*100)}],'VariableNames',{'Information'})
    if numel(mIn) == 0
        char('all done!')
        break
    end
end

XO(iNdeX,:)       = xo;
lAMBda(iNdeX,:)   = lambda;
[maTCh,q,~,~,MSE] = Model(Betas,XO,x,SET,mATCh,lAMBda);
dLa               = ModelD(Betas,XO,x,SET,mATCh,lAMBda);
end


function dLa = ModelD(B,xo,x,SET,match,lambda)
maTCh = Model(B,xo,x,SET,match,lambda);
d     = 0.0001;
for h = 1:size(xo,2)
    xU       = xo;
    xU(:,h)  = xU(:,h) + d;
    [~,~,LU] = Model(B,xU,x,SET,match,lambda);
    xL       = xo;
    xL(:,h)  = xL(:,h) - d;
    [~,~,LL] = Model(B,xL,x,SET,match,lambda);
    dLa(:,h) = (LU - LL)/(2*d);
    clear xU xL LU LL
end
dLa   = [dLa,recode(-log(maTCh{2}./match{2}),-inf,0)];
end
