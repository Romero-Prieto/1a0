function tABleSuMm(sEt,vARs,sIZe,lABs,nOTe,leGend,dATa,d)

if min(numel(vARs{1}),numel(sEt)) == 0
    S  = [0.060 1/50];
    Sx = [0.025 0.025];
else
    S  = [0.085 1/50];
    Sx = [0.025 0.050];
end
if numel(nOTe{2}) == 0
    S(2) = 0;
end

s       = 3/200;
f       = 0.03;                
r       = 2.5/200;

n       = NaN(numel(vARs),1);
for i = 1:numel(vARs)
    n(i) = numel(vARs{i});
end
R       = d(1);
d       = d(2:end);
D       = R + 2*r*numel(n) + numel(n)*sum(d);
xT      = NaN(1,numel(n));
x       = NaN(1,sum(n));
xlab    = num2cell(x);
xTab    = NaN(numel(n),2);
for i = 1:numel(n)
    for j = 1:n(i)
        x(sum(n(1:i - 1)) + j)    = (R + 2*r*(i - 1) + (i - 1)*sum(d) + r + sum(d(1:j - 1)))/D;
        xlab{sum(n(1:i - 1)) + j} = vARs{i}{j};
    end
    xT(i)     = (R + i*sum(d) + 2*r*i - (r + sum(d)/2))/D;
    xTab(i,:) = (R + (i - 1)*sum(d) + 2*r*(i - 1) + r + [0 sum(d)])/D;
end
xTab(1) = r/D; 


n       = NaN(numel(lABs),1);
for i = 1:numel(lABs)
    n(i) = numel(lABs{i});
end
F       = sum(S) + s*(numel(n) - 1) + sum(n)*f;
y       = NaN(1,sum(n));
ylab    = num2cell(x);
data    = string(NaN(numel(y),numel(x)));
for i = 1:numel(n)
    for j = 1:n(i)
        y(sum(n(1:i - 1)) + j)      = 1 - (S(1) + s*(i - 1) + f*sum(n(1:i - 1)) + f*(j - 1))/F;
        ylab{sum(n(1:i - 1)) + j}   = leGend{lABs{i}{j}};
        data(sum(n(1:i - 1)) + j,:) = dATa(lABs{i}{j},:);
    end
end

pix                      = 1/37.7952755906;
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0 0 27*D 15*F]/pix;
axes1                    = axes('Parent',fi,'Position',[0.00 0.00 1.00 1.00]);
hold(axes1,'on');
axis off
xlim([0 1])
ylim([0 1])

plot([0 1],.999*ones(2,1),'color','k','Linewidth',1.25);
plot([0 1],(S(2) + s + 0.005)/F*ones(2,1),'color','k','Linewidth',0.75);
plot([0 1],(S(2) + s)/F*ones(2,1),'color','k','Linewidth',1.25);
text(r/D,1 - Sx(2)/F,nOTe{1},'Interpreter','latex','HorizontalAlignment','left','FontName','Times New Roman','FontSize',12);
text(xTab(1,1),(S(2) - s/2)/F,nOTe{2},'Interpreter','latex','HorizontalAlignment','left','FontName','Times New Roman','FontSize',10);

for i = 1:numel(xT)
    if numel(sEt) > 0
        text(xT(i),1 - Sx(1)/F,sEt{i},'Interpreter','latex','HorizontalAlignment','center','FontName','Times New Roman','FontSize',13);
    end
    plot([xTab(i,1) xTab(i,2)],1 - (S(1) - 0.02)/F*ones(2,1),'color','k','Linewidth',0.75);
end

for i = 1:numel(y)
    text(r/D,y(i),ylab{i},'Interpreter','latex','HorizontalAlignment','left','FontName','Times New Roman','FontSize',11);
    if numel(y) > 5 && numel(y) < 27 
        text(r/4/D,y(i),char(64 + i),'Interpreter','latex','HorizontalAlignment','center','FontName','Times New Roman','FontSize',5);
        text(1 - r/4/D,y(i),char(64 + i),'Interpreter','latex','HorizontalAlignment','center','FontName','Times New Roman','FontSize',5);
    end
end

for i = 1:numel(x)
    text(x(i),1 - Sx(2)/F,xlab{i},'Interpreter','latex','HorizontalAlignment',sIZe{1,1 + mod(numel(sIZe) - 1 + i,numel(sIZe))},'FontName','Times New Roman','FontSize',11);
    for j = 1:numel(y)
        text(x(i),y(j),data(j,i),'Interpreter','latex','HorizontalAlignment',sIZe{1,1 + mod(numel(sIZe) - 1 + i,numel(sIZe))},'FontName','Times New Roman','FontSize',sIZe{2,1 + mod(numel(sIZe) - 1 + i,numel(sIZe))});
    end
end