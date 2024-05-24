function a = aOfunction(x,sex,model,input)

load('data.mat','models_a0');
models = models_a0;

for i = 4:1:size(models,2) - 1
    for j = 1:1:size(models,1)
        if numel(cell2mat(table2array(models(j,i)))) == 0
            models(j,i) = models(j,i - 1);
            if numel(models.range{j}) == 1
                models.range{j} = ones(1,2)*models.range{j};
            end
        end
    end
end

t      = ones(size(x,1),1)*find(models.model == model & models.sex == sex);
S      = [models.three(t),models.two(t),models.one(t)];
R      = cell2mat(models.range(t));
R      = [zeros(size(x,1),1),R(:,2),R(:,1)];

if ismember(model,["K" "PHG" "AKm"])
    if isequal(input,"m")
        t = sum(x(:,1) >= R,2);
        for j = 1:1:numel(t)
            A(j,:) = S{j,t(j)};
        end
        a = A(:,1) + A(:,2).*x(:,1);
        clear A j t
    else
        a = 0;
        for i = 1:1:25
            m = x(:,1)./(1 - x(:,1).*(1 - a));            
            t = sum(m >= R,2);
            for j = 1:1:numel(t)
                A(j,:) = S{j,t(j)};
            end
            a = A(:,1) + A(:,2).*m;
            clear A j t
        end
    end
elseif isequal(model,"AR")
    if isequal(input,"q")
        t = sum(x(:,1) >= R,2);
        for j = 1:1:numel(t)
            A(j,:) = S{j,t(j)};
        end
        a = A(:,1) + A(:,2).*x(:,1) + A(:,3).*(x(:,1)./x(:,2));
        clear A j t 
    else
        c = 4*x(:,2)./(1 + (4 - 1.5)*x(:,2));
        a = 0;
        for i = 1:1:25
            q(:,1) = x(:,1)./(1 + (1 - a).*x(:,1));
            q(:,2) = 1 - (1 - q(:,1)).*(1 - c);
            t      = sum(q(:,1) >= R,2);
            for j = 1:1:numel(t)
                A(j,:) = S{j,t(j)};
            end
            a      = A(:,1) + A(:,2).*q(:,1) + A(:,3).*(q(:,1)./q(:,2));
            clear A j t
        end
    end
else
    if isequal(input,"q")
        t = sum(x(:,1) >= R,2);
        for j = 1:1:numel(t)
            A(j,:) = S{j,t(j)};
        end
        a = A(:,1) + A(:,2).*x(:,1);
        clear A j t
    else
        a = 0;
        for i = 1:1:25
            q = x(:,1)./(1 + (1 - a).*x(:,1));
            t = sum(q >= R,2);
            for j = 1:1:numel(t)
                A(j,:) = S{j,t(j)};
            end            
            a = A(:,1) + A(:,2).*q;
            clear A j t
        end
    end
end