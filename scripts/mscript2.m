% Pavel Loskot, October 2022
% experiment 3: variances of reaction occurrence

clearvars -except m

Models= ['A' 'B' 'C' 'D' 'E'];
model= Models(m);
fprintf('Model: %s\n',model);

dat= readmatrix(['model' model '.csv']);
dat= dat(:,1:2);
dat1= dat(:,2);
[rxn,frq]= stat(dat1);
N= length(dat1);
Nr= length(rxn);
  
% variance
g= 0;
Grp= [];
Vars= [];
for W=10:10:50
  fprintf('W=%d\n',W);
  g= g+1;
  x= floor(N/W);
  dat2= reshape(dat1(1:W*x),W,x)';
  dat2= sort(dat2,2,'descend');
  vars= zeros(1,x);
  for i=1:x
    [~,fr]= stat(dat2(i,:));
    fr= sort(fr,'descend');
    %vars(i)= var(fr);
    vars(i)= length(find(fr));
  end
  Vars= [Vars vars];
  Grp= [Grp g*ones(1,x)];  
end

% plot
boxplot(Vars,Grp);grid
xticklabels({'10','20','30','40','50'})
%title(sprintf('tt%d',m));
xlabel('xx')
ylabel('yy')
limsy=get(gca,'YLim');
set(gca,'Ylim',[0 limsy(2)])
  
% export figure
fname= sprintf('figs/rxnvars%d.eps',m);
q= input(sprintf('\n Export as %s? [y/N] ',fname),'s');
if q=='y'
  expfig(fname);
  fprintf(' -> Exported as %s\n',fname);
end
disp(' ')


