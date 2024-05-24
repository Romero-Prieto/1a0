function [maTCh,lAB] = Set(q,x,SET,K)

n     = x{1}(2:end) - x{1}(1:end-1);
l     = 1 - q;
d     = [l(1:end-1,:) - l(2:end,:);l(end,:)];
m     = (1./n).*log(l(1:end-1,:)./l(2:end,:));
L     = [d(1:end-1,:)./m;zeros(1,size(q,2))];
L     = min(L,[n;0].*l);

y     = (1:numel(x{1}))';
maTCh = [];

for j = 1:1:numel(SET)
    set        = SET{j};
    match      = zeros(size(q,2),0);
    lab        = {};
    for k = 1:1:size(set,1)
        if iscell(set{k,3})
            BOx = x{2};
        else
            BOx = x{1};
        end
        s              = y'*[ismember(BOx,table2array(set(k,2))),ismember(BOx,table2array(set(k,3)))];
        clear BOx
        if isequal(char(table2array(set(k,1))),'q')
            match(:,end + 1) = 1 - ((1 - q(s(2),:))./(1 - q(s(1),:)))';
        elseif isequal(char(table2array(set(k,1))),'m')
            match(:,end + 1) = (sum(d(s(1):s(2)-1,:),1)./sum(L(s(1):s(2)-1,:),1))';
        elseif isequal(char(table2array(set(k,1))),'a')
            match(:,end + 1) = 1./(sum(d(s(1):s(2)-1,:),1)./sum(L(s(1):s(2)-1,:),1))' - (x{1}(s(2)) - x{1}(s(1)))*((1 - q(s(2),:))./(1 - q(s(1),:)))'./(1 - ((1 - q(s(2),:))./(1 - q(s(1),:)))');
        elseif isequal(char(table2array(set(k,1))),'z')
            match(:,end + 1) = (log((1 - q(s(1),:))./(1 - q(1,:)))./log((1 - q(s(2),:))./(1 - q(1,:))))';
        elseif isequal(char(table2array(set(k,1))),'b')
            match(:,end + 1) = (q(s(1),:)./q(s(2),:))';
        else
            temp             = char(table2array(set(k,1)));
            if isequal(temp(1),'k')
                match(:,end + 1) = exp(K(:,str2num(temp(2))));
            end
        end
        lab{end + 1,1} = char("$\mathit{" + table2array(set(k,1)) + "}\mathrm{[" + sscanf(x{2}{s(1)},'%f') + "}\mathit{" + x{2}{s(1)}(~ismember(x{2}{s(1)},char(string(sscanf(x{2}{s(1)},'%f'))))) + "}\mathrm{, " + sscanf(x{2}{s(2)},'%f') + "}\mathit{" + x{2}{s(2)}(~ismember(x{2}{s(2)},char(string(sscanf(x{2}{s(2)},'%f'))))) + "}\mathrm{)}$");          
    end
    maTCh{1,j} = match;
    lAB{1,j}   = lab;
    clear set match lab temp k j
end