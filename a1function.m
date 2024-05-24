function a = a1function(x,sex,model,input)

load('data.mat','models_a1');
models = models_a1;
t      = ones(size(x,1),1)*find(models.model == model & models.sex == sex);
S      = [models.two(t),models.one(t)];
R      = cell2mat(models.range(t));
R      = [zeros(size(x,1),1),R];

if ismember(model,"RT")
    a = 0.5;
    if isequal(input,"m")
        q = x(:,1:4)./(1 + (1 - a).*x(:,1:4));
        l = [ones(size(x,1),1),cumprod(1 - q,2)];
        d = l(:,1:end - 1) - l(:,2:end);
        L = l(:,2:end) + a*d;
        a = (sum(L,2) - 4*(1 - sum(d,2)))./sum(d,2);
        clear q l d L
    else
        l = [ones(size(x,1),1),cumprod(1 - x(:,1:4),2)];
        d = l(:,1:end - 1) - l(:,2:end);
        L = l(:,2:end) + a*d;
        a = (sum(L,2) - 4*(1 - sum(d,2)))./sum(d,2);
        clear q l d L
    end
elseif ismember(model,"PHG")
    if isequal(input,"m")
        t = sum(x(:,1) >= R,2);
        for j = 1:1:numel(t)
            A(j,:) = S{j,t(j)};
        end
        a = A(:,1) + A(:,2).*x(:,1);
        clear A j t
    else
        a = aOfunction(x,sex,model,input);
        m = x(:,1)./(1 - x(:,1).*(1 - a));
        t = sum(m >= R,2);
        for j = 1:1:numel(t)
            A(j,:) = S{j,t(j)};
        end
        a = A(:,1) + A(:,2).*m;
        clear A j t
    end
else
    if isequal(input,"m")
        a = aOfunction(x,sex,model,input);
        q = x(:,1)./(1 + (1 - a).*x(:,1));
        t = sum(q >= R,2);
        for j = 1:1:numel(t)
            A(j,:) = S{j,t(j)};
        end
        a = A(:,1) + A(:,2).*q;
        clear A j t
    else
        t = sum(x(:,1) >= R,2);
        for j = 1:1:numel(t)
            A(j,:) = S{j,t(j)};
        end
        a = A(:,1) + A(:,2).*x(:,1);
        clear A j t
    end
end