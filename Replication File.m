%% Table 1 %%
clear;
pATh              = ""; %Ajust path accordingly%

ReG               = table({'q'},{'0'},{'60m'},{'q'},2);
[beta,Q,M,x,info] = Coef(ReG);
F                 = table([],[],[],ones(0,1));
A                 = table({'a';'a'},{'0';'12m'},{'12m';'60m'});
G                 = table({'q';'q';'q';'q';'q';'m';'m';'m';'z';'z'},{'0';'0';'0';'12m';'0';'0';'12m';'0';'28d';'3m'},{'28d';'3m';'12m';'60m';'60m';'12m';'60m';'60m';'12m';'12m'});
Z                 = {[6 10],[6 9],[3 5],[6 7],[],[],[],3,4,5,6,7,8,[],[],[],[],[],[]};
Y                 = {[3 5],[1 3],[2 3]};
sex               = {'female','male','total'};

country = string(info.country(1));
year    = "";
temp    = info.year(1);
z       = temp;
tables  = 1;
for i = 2:size(info,1)
    if isequal(info.country(i),country(end))
        if isequal(info.year(i),info.year(i - 1) + 1)
            if isequal(floor(str2num(string(temp))/100),floor(info.year(i)/100))
                y = extractAfter(string(info.year(i)),0);
            else
                y = string(info.year(i));
            end                
            z = temp + "-" + y;
        else
            year(end,1) = year(end,1) + z;
            temp        = ", " + string(info.year(i));
            z           = temp;
        end
        tables(end,1)      = tables(end,1) + 1;
    else
        year(end,1)        = year(end,1) + z;
        temp               = string(info.year(i));
        z                  = temp;
        tables(end + 1,1)  = 1;
        country(end + 1,1) = string(info.country(i));
        if ~isequal(i,size(info,1))
            year(end + 1,1) = "";
        else
            year(end + 1,1) = temp;
        end
    end
end
year(end) = year(end) + z;
clear temp z info y

load('data.mat','Name')
for i = 1:size(country,1)
    country(i) = Name.Name(country(i) == Name.S);
end

country(end + 1)  = "Total";
year(end + 1)     = "";
sum(tables)
tables            = string(tables);
tables(end + 1)   = "1,219";

lABs              = {num2cell(1:25) {26}};
tABleSuMm({},{{'Years','Life Tables'}},{'left','right';10,11},lABs,{'Country','$\mathit{Under}$-$\mathit{5\ Mortality\ Database}$, https://web.sas.upenn.edu/global-age-patterns-under-five-mortality/data/'},country,[year tables],[.125 .340 .0])
saveas(gcf,char(pATh + "Results/Table 1.png"));
clear year country tables


%% Figure 1 and Figure 2 %% 
clear;
pATh              = ""; %Ajust path accordingly%
ReG               = table({'q'},{'0'},{'60m'},{'q'},2);
[beta,Q,~,x]      = Coef(ReG);
B                 = beta{3};
Q                 = Q{3}';
F                 = table([],[],[],ones(0,1));
G                 = table({'q';'q';'a';'a'},{'0';'0';'0';'12m'},{'12m';'60m';'12m';'60m'});
G2                = table({'k2'},{'0'},{'0'});
match             = Set(Q,x,{F,G});
match             = match{2};
sex               = {'female','male','total'};
k                 = [" 0.0";"+0.5";"+1.0";"-0.5";"-1.0"];
K                 = str2num(char(k));
width             = {1.2,0.9,0.9,0.9,0.9};
sTyle             = {'-','-.',':','-.',':'};
TiTLe             = {{'$\textbf{a. The log-quadratic model vs.}$','$\textbf{equations depending on one parameter}$'},{'$\textbf{b. The log-quadratic model vs.}$','$\textbf{equations depending on two parameters}$'}};

init              = log(.2./(2.^[7 0]));
n                 = 250;
h                 = (min(init):(max(init)-min(init))/(n - 1):max(init))';
H                 = exp(h);
init              = [kron(ones(size(K)),h),kron(K,ones(size(h)))];
models            = {zeros(size(init,1),0),exp(init)};
[~,q]             = Match(models,B,x,{F,[G(1,:);G2]},init,50);
models            = Set(q,x,{G([3 4],:)});
models            = models{1};
a0                = reshape(models(:,1),n,numel(K));
a1                = reshape(models(:,2),n,numel(K));
models            = Set(q,x,{G([1 2],:)});
models            = models{1};

set               = ["AK";"PHG";"K";"CD_W";"CD_E"];
for j = 1:numel(set)
    a0S(:,j) = aOfunction(H,sex{3},set{j},'q');
end
a0D               = aOfunction(models,sex{3},"AR",'q');
a0D               = reshape(a0D,n,numel(K));

set               = ["PHG";"CD_N";"CD_S";"CD_W";"CD_E"];
for j = 1:numel(set)
    a1S(:,j) = a1function(H,sex{3},set{j},'q');
end
clear set models

pix                      = 1/37.7952755906;
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0.25 0.25 17 10]/pix;
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

for i = 1:2
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 8.5;
    ax{i}.YAxis.TickLabelFormat = '%.2f';
    ax{i}.XScale                = 'log';
    ax{i}.YMinorGrid            = 'off';
    ax{i}.XMinorGrid            = 'off';
    
    ylabel('$\mathbf{_1}\textbf{\emph{a}}\mathbf{_0}$ $\textbf{(years)}$','Interpreter','latex','FontName','Times New Roman','FontSize',12);
    xlabel('$\mathbf{_1}\textbf{\emph{q}}\mathbf{_0}$ $\textbf{(log scale)}$','Interpreter','latex','FontName','Times New Roman','FontSize',12);
    ax{i}.XTick                 = .2./(2.^(7:-1:0));
    ax{i}.XAxis.MinorTickValues = .2./(2.^(7:-1:0));
    xlim(.2./(2.^[8 0]))
    
    t                           = ax{i}.XAxis.TickValues;
    T                           = ax{i}.XAxis.TickLabels;
    for j = 1:numel(t)
        T{j} = char("$\mathbf{" + string(T{j}) + "}$");
    end
    ax{i}.XAxis.TickLabels           = T;
    ax{i}.XAxis.TickLabelInterpreter = 'latex';
    
    grid on;
    box on;
    hold on;
end

coloR             = {[1.0 0.0 0.0],[0.0 0.0 1.0],[0.0 0.6 0.0],[1.0 0.0 1.0],[0.5 0.0 0.5]};
sTyle2            = {'-','-','-','--','-'};
for i = 1:2
    nexttile(i)
    ylim([0 0.4])
    ax{i}                       = gca;
    t                           = ax{i}.YAxis.TickValues;
    T                           = ax{i}.YAxis.TickLabels;    
    for j = 1:numel(t)
        T{j} = char("$\mathbf{" + string(T{j}) + "}$");
    end
    ax{i}.YAxis.TickLabels           = T;
    ax{i}.YAxis.TickLabelInterpreter = 'latex';

    U{i}{1}       = scatter(match(:,1),match(:,3),2.5,'filled','MarkerFaceColor',[0.0 0.6 0.6],'MarkerFaceAlpha',.30,'MarkerEdgeColor',[0.0 0.0 0.0],'MarkerEdgeAlpha',.30);
    U{i}{end + 1} = title(TiTLe{i},'Interpreter','latex','FontSize',10);
    if i == 1
        for j = 1:size(a0S,2)
            U{i}{end + 1} = plot(H,a0S(:,j),'color',coloR{j},'lineWidth',1.0,'LineStyle',sTyle2{j});
        end
        
        for j = 1:size(a0,2)
            U{i}{end + 1} = plot(H,a0(:,j),'color','k','lineWidth',width{j},'LineStyle',sTyle{j});
            U{i}{end + 1} = text(H(1)*17/32,a0(1,j),char("$\textbf{\emph{k }}\mathbf{=" + string(k{j}) + "}$"),'Interpreter','latex','FontName','Times New Roman','FontAngle','oblique','FontSize',6.5);
        end
        leGend        = {'$\textbf{U5MD}$';'$\textbf{Andreev and Kingkade (2015)}$';'$\textbf{Preston et al. (2001)}$';'$\textbf{Keyfitz (1970)}$';'$\textbf{Coale and Demeny (1966)}$';'$\textbf{Coale and Demeny (1966)\textendash East}$';'$\textbf{Log-quadratic model}$'};
        U{i}{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','Box','off');
    elseif i == 2
        for j = 1:3
            U{i}{end + 1} = plot(H,a0D(:,j),'color','r','lineWidth',width{j},'LineStyle',sTyle{j});
        end
        for j = 1:3
            U{i}{end + 1} = plot(H,a0(:,j),'color','k','lineWidth',width{j},'LineStyle',sTyle{j});
            U{i}{end + 1} = text(H(1)*17/32,a0(1,j),char("$\textbf{\emph{k }}\mathbf{=" + string(k{j}) + "}$"),'Interpreter','latex','FontName','Times New Roman','FontAngle','oblique','FontSize',6.5);
        end
        for j = 4:size(a0D,2)
            U{i}{end + 1} = plot(H,a0D(:,j),'color','r','lineWidth',width{j},'LineStyle',sTyle{j});
        end
        for j = 4:size(a0,2)
            U{i}{end + 1} = plot(H,a0(:,j),'color','k','lineWidth',width{j},'LineStyle',sTyle{j});
            U{i}{end + 1} = text(H(1)*17/32,a0(1,j),char("$\textbf{\emph{k }}\mathbf{=" + string(k{j}) + "}$"),'Interpreter','latex','FontName','Times New Roman','FontAngle','oblique','FontSize',6.5);
        end        
        leGend        = {'$\textbf{U5MD}$';'$\textbf{Alexander and Root (2022)}$';'$\textbf{\emph{k }}\mathbf{= \pm 0.5}$';'$\textbf{\emph{k }}\mathbf{= \pm 1.0}$';'$\textbf{Log-quadratic model, \emph{k }= 0}$';'$\textbf{\emph{k }}\mathbf{= \pm 0.5}$';'$\textbf{\emph{k }}\mathbf{= \pm 1.0}$'};
    end
    U{i}{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','Box','off');
end
saveas(gcf,char(pATh + "Results/Figure 1.png"));

for i = 1:numel(U)
    for j = 1:numel(U{i})
        delete(U{i}{j});
    end
end
clear U t T

TiTLe             = {{'$\textbf{a. The log-quadratic model vs.}$','$\textbf{equations depending on one parameter}$'},{'$\textbf{b. Predicted values using}$','$\textbf{the log-quadratic model of under-five mortality}$'}};
for i = 1:2
    nexttile(i)
    ylim([1 2])
    ax{i}                       = gca;
    t                           = ax{i}.YAxis.TickValues;    
    for j = 1:numel(t)
        T{j} = char("$\mathbf{" + string(sprintf('%.2f',t(j))) + "}$");
    end
    ax{i}.YAxis.TickLabels           = T;
    ax{i}.YAxis.TickLabelInterpreter = 'latex';    
    
    ylabel('$\mathbf{_4}\textbf{\emph{a}}\mathbf{_1}$ $\textbf{(years)}$','Interpreter','latex','FontName','Times New Roman','FontSize',12);    
    U{i}{1}       = scatter(match(:,1),match(:,4),2.5,'filled','MarkerFaceColor',[0.0 0.6 0.6],'MarkerFaceAlpha',.30,'MarkerEdgeColor',[0.0 0.0 0.0],'MarkerEdgeAlpha',.30);
    U{i}{end + 1} = title(TiTLe{i},'Interpreter','latex','FontSize',10);
    if i == 1
        coloR  = {[0.0 0.0 1.0],[0.0 0.6 0.0],[1.0 0.0 0.0],[1.0 0.0 1.0],[0.5 0.0 0.5]};
        for j = 1:size(a1S,2)
            U{i}{end + 1} = plot(H,a1S(:,j),'color',coloR{j},'lineWidth',1,'LineStyle',sTyle2{j});
        end
        for j = 1:size(a1,2)
            U{i}{end + 1} = plot(H,a1(:,j),'color','k','lineWidth',width{j},'LineStyle',sTyle{j});
        end
        leGend = {'$\textbf{U5MD}$';'$\textbf{Preston et al. (2001)}$';'$\textbf{Coale and Demeny (1966)\textendash North}$';'$\textbf{Coale and Demeny (1966)\textendash South}$';'$\textbf{Coale and Demeny (1966)\textendash West}$';'$\textbf{Coale and Demeny (1966)\textendash East}$';'$\textbf{Log-quadratic model}$'};
    elseif i == 2
        coloR  = {[0.0 0.0 0.0];[1.0 0.0 1.0];[1.0 0.0 1.0];[0.0 0.0 1.0];[0.0 0.0 1.0]};
        for j = 1:size(a1,2)
            U{i}{end + 1} = plot(H,a1(:,j),'color',coloR{j},'lineWidth',width{j},'LineStyle',sTyle{j});
        end
        leGend = {'$\textbf{U5MD}$';'$\textbf{Log-quadratic model, \emph{k }= 0}$';char("$\textbf{\emph{k }}\mathbf{=" + string(k{2}) + "}$");char("$\textbf{\emph{k }}\mathbf{=" + string(k{3}) + "}$");char("$\textbf{\emph{k }}\mathbf{=" + string(k{4}) + "}$");char("$\textbf{\emph{k }}\mathbf{=" + string(k{5}) + "}$")};
    end
    U{i}{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','Box','off');
end
saveas(gcf,char(pATh + "Results/Figure 2.png"));

%% Figure 3 %%
clear;
pATh              = ""; %Ajust path accordingly%
ReG               = table({'q'},{'0'},{'60m'},{'q'},2);
[beta,Q,~,x,info] = Coef(ReG);
B                 = beta{3};
F                 = table([],[],[],ones(0,1));
G                 = table({'m';'q';'q';'q';'z';'z'},{'0';'0';'0';'0';'28d';'3m'},{'12m';'12m';'28d';'3m';'12m';'12m'});
leGend            = {'$\mathbf{_1}\textbf{\emph{m}}\mathbf{_0}\mathbf{(_5}\textbf{\emph{q}}\mathbf{_0,}\textbf{\emph{k}}\mathbf{)}$';'$\mathbf{_1}\textbf{\emph{q}}\mathbf{_0}\mathbf{(_5}\textbf{\emph{q}}\mathbf{_0,}\textbf{\emph{k}}\mathbf{)}$';'$\textbf{\emph{q}}\mathbf{(28}\textbf{\emph{d};}\mathbf{_5}\textbf{\emph{q}}\mathbf{_0,}\textbf{\emph{k}}\mathbf{)}$';'$\textbf{\emph{q}}\mathbf{(3}\textbf{\emph{m};}\mathbf{_5}\textbf{\emph{q}}\mathbf{_0,}\textbf{\emph{k}}\mathbf{)}$';'$\textbf{\emph{z}}\mathbf{(28}\textbf{\emph{d};}\mathbf{_5}\textbf{\emph{q}}\mathbf{_0,}\textbf{\emph{k}}\mathbf{)}$';'$\textbf{\emph{z}}\mathbf{(3}\textbf{\emph{m};}\mathbf{_5}\textbf{\emph{q}}\mathbf{_0,}\textbf{\emph{k}}\mathbf{)}$'};

W                 = {[3 2],[1 5],[3 2 1 5],[4 2 1 6]};
coloR             = {[0.00 0.00 0.75],[0.00 0.75 0.75],[0.05 0.05 0.05],[0.95 0.00 0.95]};
TiTLe             = {{'$\textbf{a. Neonatal mortality rate}$ $\textbf{\emph{q}}\mathbf{(28}\textbf{\emph{d}}\mathbf{)}$','$\textbf{and infant mortality rate}$ $\mathbf{_1}\textbf{\emph{q}}\mathbf{_0}$'},{'$\textbf{b. Death rate from 0 to 1}$ $\mathbf{_1}\textbf{\emph{m}}\mathbf{_0}$ $\textbf{and the observed proportion}$','$\textbf{of infant deaths during the first month of life}$ $\textbf{\emph{z}}\mathbf{(28}\textbf{\emph{d}}\mathbf{)}$'},{' ',' ','$\textbf{c. Matching points, Norway 1975 (I):}$','$\textbf{\emph{q}}\mathbf{(28}\textbf{\emph{d}}\mathbf{)}$ $\textbf{vs.}$ $\textbf{\emph{z}}\mathbf{(28}\textbf{\emph{d}}\mathbf{)}$'},{' ',' ','$\textbf{d. Matching points, Norway 1975 (II):}$','$\textbf{\emph{q}}\mathbf{(3}\textbf{\emph{m}}\mathbf{)}$ $\textbf{vs.}$ $\textbf{\emph{z}}\mathbf{(3}\textbf{\emph{m}}\mathbf{)}$'}};

init              = log(.2./(2.^[8 -1]));
n                 = 100;
h                 = (min(init):(max(init)-min(init))/(n - 1):max(init))';
H                 = exp(h);
K                 = 1.25;
K                 = (-K:2*K/(n - 1):K)';
init              = [kron(h,ones(size(K))),kron(ones(size(h)),K)];
lambda            = {zeros(size(init,1),size(F,1)),zeros(size(init,1),size(G,1))};
[co,qts]          = Model(B,init,x,{F,G},lambda,lambda{2});

co                = co{2};
qts               = (1 - qts(2:end,:))./(1 - qts(1:end - 1,:));
test1             = max(qts >= 1);
test2             = max(qts == 0);
test3             = max(isnan(qts));
ok                = reshape(1 - max([test1;test2;test3])',numel(K),numel(H));

for i = 1:1:size(co,2)
    CO{i} = reshape(co(:,i),numel(K),numel(H));
end
clear co qts test1 test2 test3 n h

n                 = 20;
h                 = (min(init):(max(init)-min(init))/(n - 1):max(init))';
init              = [h,zeros(size(h))];
lambda            = {zeros(size(init,1),size(F,1)),zeros(size(init,1),size(G,1))};
co                = Model(B,init,x,{F,G},lambda,lambda{2});
co                = co{2};
co(:,end - 1:end) = (0:.05:.95)'*ones(1,2);

country           = Q{3}';
country           = Set(country,x,{F,G});
country           = country{2};

pix                      = 1/37.7952755906;
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0.0 0.0 17 16]/pix;
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

for i = 1:4
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 8.5;
    ax{i}.YAxis.TickLabelFormat = '%.2f';
    ax{i}.XScale                = 'log';
    ax{i}.YMinorGrid            = 'off';
    ax{i}.XMinorGrid            = 'off';
    xlim(.2./(2.^[8 0]))
    ylim([min(K) max(K)])
    
    ylabel('$\textbf{\emph{k}}$','Interpreter','latex','FontName','Times New Roman','FontSize',12);
    xlabel('$\mathbf{_5}\textbf{\emph{q}}\mathbf{_0}$ $\textbf{(log\,scale)}$','Interpreter','latex','FontName','Times New Roman','FontSize',12);
    ax{i}.XTick                 = .2./(2.^(8:-1:0));
    ax{i}.XAxis.MinorTickValues = .2./(2.^(8:-1:0));
    ax{i}.XAxis.TickLabelInterpreter = 'latex';
    ax{i}.YTick                 = min(K):.25:max(K);
    ax{i}.YAxis.MinorTickValues = min(K):.25:max(K);
    t                           = ax{i}.XAxis.TickValues;
    T                           = ax{i}.XAxis.TickLabels;
    for j = 1:numel(t)
        T{j} = char("$\mathbf{" + string(T{j}) + "}$");
    end
    ax{i}.XAxis.TickLabels           = T;
    ax{i}.XAxis.TickLabelInterpreter = 'latex';

    T                           = ax{i}.YAxis.TickLabels;
    t                           = ax{i}.YAxis.TickValues;    
    for j = 1:numel(t)
        T{j} = char("$\mathbf{" + string(T{j}) + "}$");
    end
    ax{i}.YAxis.TickLabels           = T;
    ax{i}.YAxis.TickLabelInterpreter = 'latex';
    
    grid on;
    box on;
    hold on;
end

for i = 1:numel(W)
    nexttile(i)
    w             = W{i};
    if i < 3
        [~,U{i}{1}] = contour(H,K,CO{w(1)},co(:,w(1)),'color',coloR{2*i - 1});
        if i ~= 2
            [~,U{i}{2}]       = contour(H,K,CO{w(2)},co(:,w(2)),'color',coloR{2*i});
        else
            [U{i}{2},U{i}{3}] = contour(H,K,CO{w(2)},co(:,w(2)),'color',coloR{2*i},'ShowText','on','LabelSpacing',1200);
            clabel(U{i}{2},U{i}{3},'FontSize',8,'FontName','Times New Roman','color',coloR{2*i},'Interpreter','latex');
        end
    else
        J           = 952;
        [~,U{i}{1}] = contour(H,K,CO{w(1)},country(J,w(1))*ones(2,1),'-','color',coloR{1},'LineWidth',.9);
        [~,U{i}{2}] = contour(H,K,CO{w(2)},country(J,w(2))*ones(2,1),'-','color',coloR{2},'LineWidth',1.5);
        [~,U{i}{3}] = contour(H,K,CO{w(3)},country(J,w(3))*ones(2,1),'--','color',coloR{3},'LineWidth',1.5);
        [~,U{i}{4}] = contour(H,K,CO{w(4)},country(J,w(4))*ones(2,1),'-','color',coloR{4},'LineWidth',.9);
        info(J,:)
    end
    U{i}{end + 1} = legend(leGend{w},'Interpreter','latex','FontName','Times New Roman','FontSize',11,'FontAngle','oblique','Location','southoutside','NumColumns',2,'Box','off');
    U{i}{end + 1} = title(TiTLe{i},'Interpreter','latex','FontSize',10);
end
saveas(gcf,char(pATh + "Results/Figure 3.png"));






%% Figure A1 %%
clear;
pATh              = ""; %Ajust path accordingly%
ReG               = table({'q'},{'0'},{'60m'},{'q'},2);
[~,Q,~,x,info]    = Coef(ReG);
F                 = table([],[],[],ones(0,1));
G                 = table({'a';'z';'z'},{'0';'3m';'28d'},{'12m';'12m';'12m'});
for i = 1:3
    match         = Set(Q{i}',x,{F,G});
    match         = match{2};
    y             = match(:,1);
    X             = [ones(size(match,1),1),match(:,2)];    
    sample28d{i}  = GiBBS(y,X,0);
    X             = [ones(size(match,1),1),match(:,3)];    
    sample3m{i}   = GiBBS(y,X,0);
end
match             = Set(Q{3}',x,{F,G});
match             = match{2};
TiTLe             = {'a. The first trimester of life','b. The first month of life'};

pix                      = 1/37.7952755906;
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0.25 0.25 15 7.5]/pix;
axes1                    = axes('Parent',fi,'Position',[0.05 0.05 0.95 0.95]);
hold(axes1,'on');
TL                       = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

for i = 1:2
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 9;
    ax{i}.YAxis.TickLabelFormat = '%.2f';
    ax{i}.XAxis.TickLabelFormat = '%.2f';
    ax{i}.YMinorGrid            = 'off';
    ax{i}.XMinorGrid            = 'off';
    ylim([0 0.4])

    ylabel('$\mathrm{_1}\mathit{a}\mathrm{_0}$ (in years)','Interpreter','latex','FontName','Times New Roman','FontSize',12);    
    grid on;
    box on;
    hold on;
end

N                 = 100;
for i = 1:2
    nexttile(i)
    if i == 1
        xlim([0.5 1.0])
        X             = [ones(N,1),(0.3:0.7/(N - 1):1.0)'];
        y             = X*sample28d{3}(:,1:2)' + normrnd(0,1,N,size(sample28d{3},1)).*sqrt(1./sample28d{3}(:,3))';
        xlabel('$\mathit{z}\mathrm{(3}\mathit{m}\mathrm{)}$','Interpreter','latex','FontName','Times New Roman','FontSize',12);
        U{i}{1}       = scatter(match(:,2),match(:,1),2.5,'filled','MarkerFaceColor',[0.0 0.0 0.6],'MarkerFaceAlpha',.3,'MarkerEdgeColor',[0.0 0.0 0.0],'MarkerEdgeAlpha',.3);
        U{i}{end + 1} = plot(X(:,2),prctile(y',50),'color',[.90 0.00 .90],'lineWidth',1.0,'LineStyle','-');
    elseif i == 2
        xlim([0.3 1.0])        
        X             = [ones(N,1),(0.3:0.7/(N - 1):1.0)'];
        y             = X*sample3m{3}(:,1:2)' + normrnd(0,1,N,size(sample3m{3},1)).*sqrt(1./sample3m{3}(:,3))';
        xlabel('$\mathit{z}\mathrm{(28}\mathit{d}\mathrm{)}$','Interpreter','latex','FontName','Times New Roman','FontSize',12);
        U{i}{1}       = scatter(match(:,3),match(:,1),2.5,'filled','MarkerFaceColor',[0.0 0.0 0.6],'MarkerFaceAlpha',.3,'MarkerEdgeColor',[0.0 0.0 0.0],'MarkerEdgeAlpha',.3);
        U{i}{end + 1} = plot(X(:,2),prctile(y',50),'color',[.90 0.00 .90],'lineWidth',1.0,'LineStyle','-');
    end
    
    U{i}{end + 1} = legend('U5MD','linear fitting','Interpreter','latex','FontName','Times New Roman','FontSize',10,'FontAngle','oblique','Location','southoutside','NumColumns',2,'Box','off');
    U{i}{end + 1} = title(TiTLe{i},'Interpreter','latex','FontSize',10);
end
saveas(gcf,char(pATh + "Appendix/Figure A1.png"));











%% Table A3, Table 2, Table 3, and Figure A3 - A8 %%
clear;
pATh              = ""; %Ajust path accordingly%
ReG               = table({'q'},{'0'},{'60m'},{'q'},2);
[beta,Q,M,x,info] = Coef(ReG);
F                 = table([],[],[],ones(0,1));
A                 = table({'a';'a'},{'0';'12m'},{'12m';'60m'});
G                 = table({'q';'q';'q';'q';'q';'m';'m';'m';'z';'z'},{'0';'0';'0';'12m';'0';'0';'12m';'0';'28d';'3m'},{'28d';'3m';'12m';'60m';'60m';'12m';'60m';'60m';'12m';'12m'});
Z                 = {[6 10],[6 9],[3 5],[6 7],[],[],[],3,4,5,6,7,8,[],[],[],[],[],[]};
Y                 = {[3 5],[1 3],[2 3]};
sex               = {'female','male','total'};

for i = 1:numel(sex)
    B             = beta{i};
    q             = Q{i}';
    init          = [log(q(end,:))',zeros(size(q,2),1)];
    sET           = Set(q,x,{A,G});
    a             = sET{1};
    sET           = sET{2};
    for j = 1:numel(Z)
        if numel(Z{j}) > 0
            match     = {zeros(size(sET,1),0),sET(:,Z{j})};
            [~,q]     = Match(match,B,x,{F,G(Z{j},:)},init(:,1:numel(Z{j})),50);
            aH        = Set(q,x,{A});
            a0(:,j)   = aH{1}(:,1);
            a1(:,j)   = aH{1}(:,2);
            clear match q aH
        else
            a0(:,j)   = NaN;
            a1(:,j)   = NaN;            
        end
    end
    
    a0(:,5)       = aOfunction(sET(:,[3 5]),sex{i},'AR','q');
    
    X             = [ones(size(sET,1),1),sET(:,10)];
    gibbs         = GiBBS(a(:,1),X,0);
    a0(:,6)       = X*prctile(gibbs(:,1:end - 1),50)';
    X             = [gibbs(:,1),gibbs(:,2),NaN(size(gibbs,1),1)];    
    coEF{1,i}     = prctile(X,[50 2.5 97.5])';
    X             = [ones(size(sET,1),1),sET(:,9)];
    gibbs         = GiBBS(a(:,1),X,0);
    a0(:,7)       = X*prctile(gibbs(:,1:end - 1),50)';
    X             = [gibbs(:,1),NaN(size(gibbs,1),1),gibbs(:,2)];
    coEF{2,i}     = prctile(X,[50 2.5 97.5])';
    
    X             = [ones(size(sET,1),1),sET(:,[10 6])];
    gibbs         = GiBBS(a(:,1),X,0);
    X             = [gibbs(:,1),gibbs(:,2),NaN(size(gibbs,1),1),gibbs(:,3)];
    coEF{3,i}     = prctile(X,[50 2.5 97.5])';
    X             = [ones(size(sET,1),1),sET(:,[9 6])];
    gibbs         = GiBBS(a(:,1),X,0);
    X             = [gibbs(:,1),NaN(size(gibbs,1),1),gibbs(:,2),gibbs(:,3)];
    coEF{4,i}     = prctile(X,[50 2.5 97.5])';    
    X             = [ones(size(sET,1),1),sET(:,3),sET(:,3)./sET(:,5)];
    gibbs         = GiBBS(a(:,1),X,0);
    coEF{5,i}     = prctile(gibbs(:,1:end - 1),[50 2.5 97.5])';
        
    a0(:,14)      = aOfunction(sET(:,3),sex{i},'AK','q');
    a0(:,15)      = aOfunction(sET(:,6),sex{i},'AKm','m');
    a0(:,16)      = aOfunction(sET(:,6),sex{i},'PHG','m');
    a0(:,17)      = aOfunction(sET(:,6),sex{i},'K','m');
    a0(:,18)      = aOfunction(sET(:,3),sex{i},'CD_W','q');
    a0(:,19)      = aOfunction(sET(:,3),sex{i},'CD_E','q');
    
    a1(:,14)      = a1function(sET(:,6),sex{i},'PHG','m');    
    a1(:,15)      = 1.5;    
    a1(:,16)      = a1function(sET(:,3),sex{i},'CD_N','q');
    a1(:,17)      = a1function(sET(:,3),sex{i},'CD_S','q');
    a1(:,18)      = a1function(sET(:,3),sex{i},'CD_W','q');
    a1(:,19)      = a1function(sET(:,3),sex{i},'CD_E','q');    

    OuT(:,i)      = {[a(:,1),a0];[a(:,2),a1]};
    RElaTiVe(:,i) = {[];[]};
    RElaTiVS(:,i) = {[];[]};
    clear init q B sET a a0 a1 sample X
end

sEt               = {'Direct estimation, using $\mathit{z}\mathrm{(3}\mathit{m}\mathrm{)}$','Direct estimation, using $\mathit{z}\mathrm{(28}\mathit{d}\mathrm{)}$'};
vARs              = {'Female','Male','Both sexes'};
vARs              = {vARs,vARs};
foRMaT            = {'%+0.4f','%+0.4f','%+0.4f'};
lABs              = {{1} {2} {3}};
nOTe              = {'','$\mathrm{p50/}\mathit{[p2.5,p97.5]}$'};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,{{'$\delta\mathrm{_0}$';''},{'$\delta\mathrm{_1}^{|\mathit{z}\mathrm{(3}\mathit{m}\mathrm{)}}$';''},{'$\delta\mathrm{_1}^{|\mathit{z}\mathrm{(28}\mathit{d}\mathrm{)}}$';''}},[cell2mat(coEF(1,:)) cell2mat(coEF(2,:))],0.05,0.075,[])
saveas(gcf,char(pATh + "Appendix\Table A3.png"));


for i = 1:numel(OuT)
    set         = OuT{i};
    E           = set(:,2:end) - set(:,1);
    e           = log(set(:,2:end)./set(:,1));
    O           = ones(size(E,1),1);
    for j = 1:size(E,2)
        if isnan(E(1,j))
            indexA(j,:) = NaN(1,9);
            indexR(j,:) = NaN(1,9);
        else            
            s           = GiBBS(e(:,j),O,j);
            s           = [s,sqrt(s(:,1).^2 + 1./s(:,2))];
            indexR(j,:) = [prctile(s(:,1),[50 2.5 97.5]),prctile(s(:,2),[50 2.5 97.5]),prctile(s(:,3),[50 2.5 97.5])];
            SR{j}       = s;
        end
        clc;
        [i,j]
    end
    RElaTiVe{i} = indexR;
    RElaTiVS{i} = SR;
    clear set E e O s indexR SR
end

coloR             = {[0.05 0.05 0.05],[0.95 0.00 0.95],[0.00 0.75 0.75],[0.00 0.00 0.75]};
TiTLe             = {{'a. Direct (linear regression) vs. indirect (model)','estimation of $\mathrm{_1}\mathit{a}\mathrm{_0}$, using $\mathit{z}\mathrm{(3}\mathit{m}\mathrm{)}$ as input'},{'b. Direct (linear regression) vs. indirect (model)','estimation of $\mathrm{_1}\mathit{a}\mathrm{_0}$, using $\mathit{z}\mathrm{(28}\mathit{d}\mathrm{)}$ as input'},{'c. New approaches using traditional inputs','to estimate $\mathrm{_1}\mathit{a}\mathrm{_0}$'},{'d. The systematic bias of','classic approaches estimating $\mathrm{_1}\mathit{a}\mathrm{_0}$'}};
W                 = {[6 1],[7 2],[5 15 11 13],[16 17 18];3,3,3,1};
laBEl             = {'$\mathit{\mu}$ (bias - $\mathrm{_1}\mathit{a}\mathrm{_0}$)','$\mathit{\varphi}$ (precision - $\mathrm{_1}\mathit{a}\mathrm{_0}$)','RMSE (predicting - $\mathrm{_1}\mathit{a}\mathrm{_0}$)'};

leGend1a0         = {
    {'Log-quadratic model, $\mathrm{_1}\mathit{m}\mathrm{_0}$ and $\mathit{z}\mathrm{(3}\mathit{m}\mathrm{)}$';''};
    {'Log-quadratic model, $\mathrm{_1}\mathit{m}\mathrm{_0}$ and $\mathit{z}\mathrm{(28}\mathit{d}\mathrm{)}$';''};
    {'Log-quadratic model, $\mathrm{_1}\mathit{q}\mathrm{_0}$ and $\mathrm{_5}\mathit{q}\mathrm{_0}$';''};
    {'Log-quadratic model, $\mathrm{_1}\mathit{m}\mathrm{_0}$ and $\mathrm{_4}\mathit{m}\mathrm{_1}$';''};
    {'Alexander and Root (2022), $\mathrm{_1}\mathit{q}\mathrm{_0}$ and $\mathrm{_1}\mathit{q}\mathrm{_0}$/$\mathrm{_5}\mathit{q}\mathrm{_0}$';''};
    {'Linear regression, $\mathit{z}\mathrm{(3}\mathit{m}\mathrm{)}$';''};
    {'Linear regression, $\mathit{z}\mathrm{(28}\mathit{d}\mathrm{)}$';''};
    {'Log-quadratic model, $\mathrm{_1}\mathit{q}\mathrm{_0}$ and $\mathit{k}$ = 0';''};
    {'Log-quadratic model, $\mathrm{_4}\mathit{q}\mathrm{_1}$ and $\mathit{k}$ = 0';''};
    {'Log-quadratic model, $\mathrm{_5}\mathit{q}\mathrm{_0}$ and $\mathit{k}$ = 0';''};
    {'Log-quadratic model, $\mathrm{_1}\mathit{m}\mathrm{_0}$ and $\mathit{k}$ = 0';''};
    {'Log-quadratic model, $\mathrm{_4}\mathit{m}\mathrm{_1}$ and $\mathit{k}$ = 0';''};
    {'Log-quadratic model, $\mathrm{_5}\mathit{m}\mathrm{_0}$ and $\mathit{k}$ = 0';''};
    {'Andreev and Kingkade (2015), $\mathrm{_1}\mathit{q}\mathrm{_0}$';''};
    {'Andreev and Kingkade (2015), $\mathrm{_1}\mathit{m}\mathrm{_0}$';''};
    {'Preston et al. (2001), $\mathrm{_1}\mathit{m}\mathrm{_0}$';''};
    {'Keyfitz (1970), $\mathrm{_1}\mathit{m}\mathrm{_0}$';''};
    {'Coale and Demeny (1966), $\mathrm{_1}\mathit{q}\mathrm{_0}$';''};
    {'Coale and Demeny (1966) - East, $\mathrm{_1}\mathit{q}\mathrm{_0}$';''}};

leGend4a1                = {
    {'Log-quadratic model, $\mathrm{_1}\mathit{m}\mathrm{_0}$ and $\mathit{z}\mathrm{(3}\mathit{m}\mathrm{)}$';''};
    {'Log-quadratic model, $\mathrm{_1}\mathit{m}\mathrm{_0}$ and $\mathit{z}\mathrm{(28}\mathit{d}\mathrm{)}$';''};
    {'Log-quadratic model, $\mathrm{_1}\mathit{q}\mathrm{_0}$ and $\mathrm{_5}\mathit{q}\mathrm{_0}$';''};
    {'Log-quadratic model, $\mathrm{_1}\mathit{m}\mathrm{_0}$ and $\mathrm{_4}\mathit{m}\mathrm{_1}$';''};
    {'';''};
    {'';''};
    {'';''};
    {'Log-quadratic model, $\mathrm{_1}\mathit{q}\mathrm{_0}$ and $\mathit{k}$ = 0';''};
    {'Log-quadratic model, $\mathrm{_4}\mathit{q}\mathrm{_1}$ and $\mathit{k}$ = 0';''};
    {'Log-quadratic model, $\mathrm{_5}\mathit{q}\mathrm{_0}$ and $\mathit{k}$ = 0';''};
    {'Log-quadratic model, $\mathrm{_1}\mathit{m}\mathrm{_0}$ and $\mathit{k}$ = 0';''};
    {'Log-quadratic model, $\mathrm{_4}\mathit{m}\mathrm{_1}$ and $\mathit{k}$ = 0';''};
    {'Log-quadratic model, $\mathrm{_5}\mathit{m}\mathrm{_0}$ and $\mathit{k}$ = 0';''};
    {'Prestonet al. (2001), $\mathrm{_1}\mathit{m}\mathrm{_0}$';''};
    {'Keyfitz and Flieger(1971), $\mathrm{_4}\mathit{a}\mathrm{_1}$ = 1.5';''};
    {'Coale and Demeny (1966) - North, $\mathrm{_1}\mathit{q}\mathrm{_0}$';''};
    {'Coale and Demeny (1966) - South, $\mathrm{_1}\mathit{q}\mathrm{_0}$';''};
    {'Coale and Demeny (1966) - West, $\mathrm{_1}\mathit{q}\mathrm{_0}$';''};
    {'Coale and Demeny (1966) - East, $\mathrm{_1}\mathit{q}\mathrm{_0}$';''}};

pix                      = 1/37.7952755906;
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0.25 0.25 15 15]/pix;
axes1                    = axes('Parent',fi,'Position',[0.05 0.05 0.95 0.95]);
hold(axes1,'on');
TL                       = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
u                        = 3;
for i = 1:4
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 9;
    ax{i}.YAxis.TickLabelFormat = '%.0f';
    ax{i}.XAxis.TickLabelFormat = '%.3f';
    
    ylabel('$\mathit{kernel\ density}$','Interpreter','latex','FontName','Times New Roman','FontSize',12);
    xlabel(laBEl{W{2,i}},'Interpreter','latex','FontName','Times New Roman','FontSize',12);
    grid on;
    box on;
    hold on;
end

for h = 1:size(RElaTiVS,2)
    for i = 1:4
        nexttile(i)
        w = W{1,i};
        for j = 1:numel(w)
            [f{j},xi{j}]  = ksdensity(RElaTiVS{1,h}{w(j)}(:,W{2,i}));
            U{i}{j}       = plot(xi{j},f{j},'color',coloR{j},'LineWidth',1.00);
        end        
        for j = 1:numel(w)
            U{i}{end + 1} = fill(xi{j},f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
        end
        clear f xi
        for j = 1:numel(w)
            [f,xi(:,j)]   = ksdensity(RElaTiVS{1,h}{w(j)}(:,W{2,i}),prctile(RElaTiVS{1,h}{w(j)}(:,W{2,i}),[50 2.5 97.5 0.001 99.999]));
            U{i}{end + 1} = plot(xi(1,j)*ones(2,1),[f(1) 0],'color',coloR{j},'LineWidth',1.0,'LineStyle',':');
            U{i}{end + 1} = plot(xi(2,j)*ones(2,1),[f(2) 0],'color',coloR{j},'LineWidth',0.5,'LineStyle','-');
            U{i}{end + 1} = plot(xi(3,j)*ones(2,1),[f(3) 0],'color',coloR{j},'LineWidth',0.5,'LineStyle','-');
            leGend{j}     = leGend1a0{w(j)}{1};
        end
        xlim([min(xi(4,:)) max(xi(5,:))])
        clear f xi
        U{i}{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','Box','off');
        U{i}{end + 1} = title(TiTLe{i},'Interpreter','latex','FontSize',10);
        coloR         = coloR([3 4 1 2]);
        clear leGend
    end
    
    saveas(gcf,char(pATh + "Appendix/Figure A" + string(u) + ", " + string(sex{h}) + ".png"));
    u = u + 1;
    for i = 1:numel(U)
        for j = 1:numel(U{i})
            delete(U{i}{j});
        end
    end
    clear U
end

for i = 3:4
    nexttile(i)
    ax{i}.XAxis.TickLabelFormat = '%.0f';
    ax{i}.YAxis.TickLabelFormat = '%.2f';
end

TiTLe             = {{'a. Indirect (model) estimation of $\mathrm{_4}\mathit{a}\mathrm{_1}$,','using one or two inputs'},{'b. Indirect (model) estimation of $\mathrm{_4}\mathit{a}\mathrm{_1}$,','using one or two inputs (Continuation)'},{'c. Direct (classic approaches) vs. indirect (model)','estimation of $\mathrm{_4}\mathit{a}\mathrm{_1}$, using $\mathit{z}\mathrm{(3}\mathit{m}\mathrm{)}$ as input'},{'d. Direct (classic approaches) vs. indirect (model)','estimation of $\mathrm{_4}\mathit{a}\mathrm{_1}$, using $\mathit{z}\mathrm{(3}\mathit{m}\mathrm{)}$ as input (Continuation)'}};
W                 = {[1 11 2 4],[1 11 13 8],[1 11 14 15],[1 16 17 19];3,3,2,2};
laBEl             = {'$\mathit{\mu}$ (bias - $\mathrm{_4}\mathit{a}\mathrm{_1}$)','$\mathit{\varphi}$ (precision - $\mathrm{_4}\mathit{a}\mathrm{_1}$)','RMSE (predicting - $\mathrm{_4}\mathit{a}\mathrm{_1}$)'};
for h = 1:size(RElaTiVS,2)
    for i = 1:4
        nexttile(i)
        w = W{1,i};
        for j = 1:numel(w)
            [f{j},xi{j}]  = ksdensity(RElaTiVS{2,h}{w(j)}(:,W{2,i}));
            U{i}{j}       = plot(xi{j},f{j},'color',coloR{j},'LineWidth',1.00);
        end        
        for j = 1:numel(w)
            U{i}{end + 1} = fill(xi{j},f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
        end
        clear f xi
        for j = 1:numel(w)
            [f,xi(:,j)]   = ksdensity(RElaTiVS{2,h}{w(j)}(:,W{2,i}),prctile(RElaTiVS{2,h}{w(j)}(:,W{2,i}),[50 2.5 97.5 0.001 99.999]));
            U{i}{end + 1} = plot(xi(1,j)*ones(2,1),[f(1) 0],'color',coloR{j},'LineWidth',1.0,'LineStyle',':');
            U{i}{end + 1} = plot(xi(2,j)*ones(2,1),[f(2) 0],'color',coloR{j},'LineWidth',0.5,'LineStyle','-');
            U{i}{end + 1} = plot(xi(3,j)*ones(2,1),[f(3) 0],'color',coloR{j},'LineWidth',0.5,'LineStyle','-');
            leGend{j}     = leGend1a0{w(j)}{1};
        end
        xlabel(laBEl{W{2,i}},'Interpreter','latex','FontName','Times New Roman','FontSize',12);
        xlim([min(xi(4,:)) max(xi(5,:))])
        clear f xi
        U{i}{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','Box','off');
        U{i}{end + 1} = title(TiTLe{i},'Interpreter','latex','FontSize',10);
        clear leGend
    end
    saveas(gcf,char(pATh + "Appendix/Figure A" + string(u) + ", " + string(sex{h}) + ".png"));
    u = u + 1;
    for i = 1:numel(U)
        for j = 1:numel(U{i})
            delete(U{i}{j});
        end
    end
    clear U
end

sEt               = {'Female','Male','Both sexes combined'};
vARs              = {'$\mathit{\mu}$ (bias)','$\mathit{\varphi}$ (precision)','RMSE'};
foRMaT            = {'%+0.4f','%0.2f','%0.4f'};
lABs              = {{1 2 3 4} {5} {6 7} {8 11 13} {15} {16} {17} {18 19}};
nOTe              = {'$\mathrm{method/}\mathit{input}$','$\mathrm{p50/}\mathit{[p2.5,p97.5]}$'};
S                 = ones(15,3);
S(13,[1 2])       = 0;
S([11 12],3)      = 0;
S([14 15],3)      = 0;
tABleBAyEs(sEt,{vARs,vARs,vARs},foRMaT,lABs,nOTe,leGend1a0,[RElaTiVe{1,1},RElaTiVe{1,2},RElaTiVe{1,3}],0.225,0.0625,S)
saveas(gcf,char(pATh + "Results/Table 2.png"));


lABs              = {{1 2 3 4} {8 11 13} {14} {15} {16 17 18 19}};
S                 = ones(13,3);
S(9,[1 2])        = 0;
S(8,3)            = 0;
S(10:13,3)        = 0;
tABleBAyEs(sEt,{vARs,vARs,vARs},foRMaT,lABs,nOTe,leGend4a1,[RElaTiVe{2,1},RElaTiVe{2,2},RElaTiVe{2,3}],0.225,0.0625,S)
saveas(gcf,char(pATh + "Results/Table 3.png"));


%% Figure A2 %%
clear;
pATh              = ""; %Ajust path accordingly%
ReG               = table({'q'},{'0'},{'60m'},{'q'},2);
[beta,Q,~,x]      = Coef(ReG);
B                 = beta{3};
Q                 = Q{3}';
F                 = table([],[],[],ones(0,1));
G                 = table({'m';'m';'q';'q';'a';'a'},{'0';'12m';'0';'12m';'0';'12m'},{'12m';'60m';'12m';'60m';'12m';'60m'});
G2                = table({'k2'},{'0'},{'0'});
match             = Set(Q,x,{F,G});
match             = match{2};
sex               = {'female','male','total'};
k                 = [" 0.0";"+1.00";"+0.75";"+0.50";"+0.25"];
K                 = str2num(char(k));


init              = log(.2./(2.^[7 0]));
n                 = 250;
h                 = (min(init):(max(init)-min(init))/(n - 1):max(init))';
H                 = exp(h);
init              = [kron(ones(size(K)),h),kron(K,ones(size(h)))];
models            = {zeros(size(init,1),0),exp(init)};
[~,q]             = Match(models,B,x,{F,[G(1,:);G2]},init,50);
models            = Set(q,x,{G});

mx{1}             = reshape(models{1}(:,1),n,numel(K));
mx{2}             = reshape(models{1}(:,2),n,numel(K));
qx{1}             = reshape(models{1}(:,3),n,numel(K));
qx{2}             = reshape(models{1}(:,4),n,numel(K));
ax{1}             = reshape(models{1}(:,5),n,numel(K));
ax{2}             = reshape(models{1}(:,6),n,numel(K));
sEL               = [2 3 4 5];
N                 = [1 4];

for i = 1:numel(mx)
    OuT{i}             = (ax{i}(:,sEL) - ax{i}(:,1))./ax{i}(:,1)*100;
    OuT{i + numel(mx)} = OuT{i}.*qx{i}(:,sEL).*ax{i}(:,sEL)/N(i);
end

xtick                    = {.2./(2.^(7:-1:0)),.2./(2.^(10:-1:2))};
ylaB                     = {'(\%)'};
xlaB                     = {'$\mathrm{_1}\mathit{m}\mathrm{_0}$ (log scale)','$\mathrm{_4}\mathit{m}\mathrm{_1}$ (log scale)'};
TiTLe                    = {'$\mathrm{a.\,Percentage\,correction\,in}$ $\mathrm{_1}\mathit{a}\mathrm{_0}$','$\mathrm{b.\,Percentage\,correction\,in}$ $\mathrm{_4}\mathit{a}\mathrm{_1}$','$\mathrm{c.\,Expected\,adjustment\,in}$ $\mathrm{_1}\mathit{q}\mathrm{_0}$','$\mathrm{d.\,Expected\,adjustment\,in}$ $\mathrm{_4}\mathit{q}\mathrm{_1}$'};
coloR                    = {[0.00 0.00 0.75],[0.95 0.00 0.95],[0.85 0.35 0.01],[0.45 0.65 0.20],[0.65 0.10 0.20],[0.00 0.55 0.65],[0.05 0.05 0.05]};

pix                      = 1/37.7952755906;
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0.25 0.25 15 15]/pix;
axes1                    = axes('Parent',fi,'Position',[0.05 0.05 0.95 0.95]);
hold(axes1,'on');
TL                       = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

for i = 1:4
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 9;
    ax{i}.YAxis.TickLabelFormat = '%.2f';
    ax{i}.XScale                = 'log';
    ax{i}.YMinorGrid            = 'off';
    ax{i}.XMinorGrid            = 'off';
    
    title(TiTLe{i},'Interpreter','latex','FontSize',10);
    ylabel(ylaB{1},'Interpreter','latex','FontName','Times New Roman','FontSize',12);
    xlabel(xlaB{2 - mod(i,2)},'Interpreter','latex','FontName','Times New Roman','FontSize',12);
    ax{i}.XTick                 = xtick{2 - mod(i,2)};
    ax{i}.XAxis.MinorTickValues = xtick{2 - mod(i,2)};
    %ax{i}.YTick                 = min(K):.25:max(K);
    %ax{i}.YAxis.MinorTickValues = min(K):.25:max(K);
    xlim([min(xtick{2 - mod(i,2)}) max(xtick{2 - mod(i,2)})]);
    %ylim([min(K) max(K)])
    grid on;
    box on;
    hold on;
end

for i = 1:4
    nexttile(i)
    for j = 1:numel(sEL)
        plot(mx{2 - mod(i,2)}(:,sEL(j)),OuT{i}(:,j),'color',coloR{j},'lineWidth',1.0,'LineStyle','-')
    end
    
end

nexttile(3)
for j = 1:numel(sEL)
    tEXt{j} = char("$\mathit{k:\,}\mathit{from}\mathrm{\,0.0\,}\mathit{to}\mathrm{\," + string(k{sEL(j)}) + "}$");
end
legend(tEXt,'Interpreter','latex','FontName','Times New Roman','FontSize',11,'FontAngle','oblique','Location','southoutside','NumColumns',2,'Box','off');
saveas(gcf,char(pATh + "Appendix/Figure A2.png"));


%% Worked Examples %%
clear;
pATh              = ""; %Ajust path accordingly%
ReG          = table({'q'},{'0'},{'60m'},{'q'},2);
[beta,Q,~,x] = Coef(ReG);
B            = beta{3};

OBS          = 0.05;
xo           = [log(0.10) 0.5]';
for i = 1:3
    [R,dR,MOD]    = dRx(xo(:,end),B,x,OBS);
    xo(:,end + 1) = xo(:,end) - [pinv(dR)*R;0];
    
    TaBLE{i,1}    = [exp(xo(1,end - 1)),NaN(1,2)];
    TaBLE{i,2}    = [xo(2,end - 1),NaN(1,2)];
    TaBLE{i,3}    = [MOD(1),NaN(1,2)];
    TaBLE{i,4}    = [R(1),NaN(1,2)];
    TaBLE{i,5}    = [dR(1),NaN(1,2)];
    lABs{1}{i}    = i;
    leGend{i,1}   = "$\mathit{" + string(i) + "}$";
end
sEt          = {'$\textit{Worked example:}$ $\mathrm{_1}\mathit{m}\mathrm{_0\,=\,0.05}$ $\textit{and}$ $\mathit{k}\mathrm{\,=\,0.5}$'};
vARs         = {{'$\mathrm{_5}\mathit{q}\mathrm{_0}$','$\mathit{k}$','$\mathrm{_1}\mathit{m}\mathrm{_0(}\mathrm{_5}\mathit{q}\mathrm{_0,}\mathit{k}\mathrm{)}$','$\mathrm{R_1(}\mathrm{_5}\mathit{q}\mathrm{_0,}\mathit{k}\mathrm{)}$','$\mathrm{\partial}\mathrm{R_1(\cdot)/}\mathrm{\partial}\mathrm{ln[}\mathrm{_5}\mathit{q}\mathrm{_0]}$'}};
foRMaT       = {'%0.4f','%0.4f','%0.4f','%0.5f','%0.4f'};
nOTe         = {'$\mathit{iteration}$',''};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,leGend,cell2mat(TaBLE),0.01,0.085,[])
saveas(gcf,char(pATh + "Appendix/Table A1.png"));
[R,MOD,a]    = Rx(xo(:,end),B,x,OBS)
clear TaBLE


OBS          = [0.025 0.75]';
xo           = [log(0.05) 0.00]';
for i = 1:3
    [R,dR,MOD]    = dRx(xo(:,end),B,x,OBS);
    xo(:,end + 1) = xo(:,end) - pinv(dR)*R;
    
    TaBLE{i,1}    = [exp(xo(1,end - 1)),NaN(1,2)];
    TaBLE{i,2}    = [xo(2,end - 1),NaN(1,2)];
    TaBLE{i,3}    = [MOD(1),NaN(1,2)];
    TaBLE{i,4}    = [R(1),NaN(1,2)];    
    TaBLE{i,5}    = [MOD(2),NaN(1,2)];
    TaBLE{i,6}    = [R(2),NaN(1,2)];
    TaBLE{i,7}    = [det(dR),NaN(1,2)];
    TaBLE{i,8}    = [dR(1,1),NaN(1,2)];
    TaBLE{i,9}    = [dR(1,2),NaN(1,2)];
    TaBLE{i,10}   = [dR(2,1),NaN(1,2)];
    TaBLE{i,11}   = [dR(2,2),NaN(1,2)];
    
    lABs{1}{i}    = i;
    leGend{i,1}   = "$\mathit{" + string(i) + "}$";
end
sEt          = {'$\mathit{Worked\,example:}$ $\mathrm{_1}\mathit{m}\mathrm{_0\,=\,0.025}$ $\mathit{and}$ $\mathit{z}\mathrm{(3}\mathit{m}\mathrm{)}\mathrm{\,=\,0.75}$'};
vARs         = {{'$\mathrm{_5}\mathit{q}\mathrm{_0}$','$\mathit{k}$','$\mathrm{_1}\mathit{m}\mathrm{_0(}\mathrm{_5}\mathit{q}\mathrm{_0,}\mathit{k}\mathrm{)}$','$\mathrm{R_1(}\mathrm{_5}\mathit{q}\mathrm{_0,}\mathit{k}\mathrm{)}$','$\mathit{z}\mathrm{(3}\mathit{m}\mathrm{;}\mathrm{_5}\mathit{q}\mathrm{_0,}\mathit{k}\mathrm{)}$','$\mathrm{R_2(}\mathrm{_5}\mathit{q}\mathrm{_0,}\mathit{k}\mathrm{)}$','$\mathrm{det[J]}$','$\mathrm{\partial}\mathrm{R_1(\cdot)/}\mathrm{\partial}\mathrm{ln[}\mathrm{_5}\mathit{q}\mathrm{_0]}$','$\mathrm{\partial}\mathrm{R_1(\cdot)/}\mathrm{\partial}\mathit{k}$','$\mathrm{\partial}\mathrm{R_2(\cdot)/}\mathrm{\partial}\mathrm{ln[}\mathrm{_5}\mathit{q}\mathrm{_0]}$','$\mathrm{\partial}\mathrm{R_2(\cdot)/}\mathrm{\partial}\mathit{k}$'}};
foRMaT       = {'%0.4f','%0.4f','%0.4f','%0.5f','%0.4f','%0.5f','%0.4f','%0.4f','%0.4f','%0.4f','%0.4f'};
nOTe         = {'$\mathit{iteration}$','$\textrm{det[J] is the determinant of the Jacobian matrix}\mathrm{\,[first\,derivatives\,of\,R_j(}\mathrm{_5}\mathit{q}\mathrm{_0,}\mathit{k}\mathrm{)].}$'};
tABleBAyEs(sEt,{vARs{1}(1:7)},foRMaT,lABs,nOTe,leGend,cell2mat(TaBLE(:,1:7)),0.01,0.07,[])
saveas(gcf,char(pATh + "Appendix/Table A2.png"));
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,leGend,cell2mat(TaBLE),0.01,0.085,[])
[R,MOD,a]    = Rx(xo(:,end),B,x,OBS)


function [R,MOD,a] = Rx(xo,B,x,OBS)
n       = x{1}(2:end) - x{1}(1:end - 1);
q       = [0;exp(B{1}(:,1:4)*[1 xo(1) xo(1)^2 xo(2)]')];
m       = (log(1 - q(1:end - 1)) - log(1 - q(2:end)))./n;
d       = q(2:end) - q(1:end - 1);
L       = d./m;
a(1)    = (sum(L(1:15)) - 1*(1 - sum(d(1:15))))./sum(d(1:15));
a(2)    = (sum(L(16:end)) - 4*(1 - sum(d)))./sum(d(16:end));

m       = sum(d(1:15))/sum(L(1:15));
z       = log(1 - q(7))/log(1 - q(16));
MOD     = [m;z];
MOD     = MOD(1:numel(OBS));
R       = log(OBS) - log(MOD);
end

function [R,dR,MOD] = dRx(xo,B,x,OBS)
delta   = 10^-5;
[R,MOD] = Rx(xo,B,x,OBS);
for i = 1:numel(R)
    xa      = xo;
    xa(i)   = xa(i) + delta;
    dR(:,i) = (Rx(xa,B,x,OBS) - R)/delta;
end
end