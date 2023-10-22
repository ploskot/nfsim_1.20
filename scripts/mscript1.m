% Pavel Loskot, October 2022
% experiment 2: reaction frequency and clustering

clearvars -except m

Models= ['A' 'B' 'C' 'D' 'E'];
model= Models(m);
fprintf('Model: %s\n',model);

dat= readmatrix(['model' model '.csv']);
dat= dat(:,1:2);
dat1= dat(:,2);
[rxn,frq]= stat(dat1);
N= length(dat1);
  
%subplot(5,1,m)
frq1= sort(100*frq/sum(frq),'descend');
bar(frq1);grid
%title(sprintf('Model %s',model))
xlabel('xx')
ylabel('yy')

Nr= length(rxn);
% reaction clusters [see m1a]
switch m
  case 1 % model A
    ci1= [1 2];
    ci2= [3 4];
    ci3= [5 6];
  case 2 % model B
    ci1= [1];
    ci2= [2:7];
    ci3= [8:22];
  case 3 % model C
    ci1= [1:2];
    ci2= [3:6];
    ci3= [7:16];    
  case 4 % model D
    ci1= [1:18];
    ci2= [19:24];
    ci3= [];    
  case 5 % model E
    ci1= [1:12];
    ci2= [13:18];
    ci3= [];    
end

% plot
clf
h(1)= bar(frq1(ci1),'b'); % cl1
hold on
h(2)= bar([zeros(length(ci1),1);frq1(ci2)],'r'); % cl2
if ~isempty(ci3)
  h(3)= bar([zeros(length([ci1 ci2]),1);frq1(ci3)],'k'); % cl3
end
hold off
grid;xlabel('xx');ylabel('yy')
if length(h)==3
  legend(h,'cluster #1','cluster #2','cluster #3','Location','NE')
else
  legend(h,'cluster #1','cluster #2','Location','NE')
end
xticks(1:Nr)
tmp= cell(1,Nr);
for i=1:Nr
  if i==1
    tmp{i}= '1';
  elseif mod(i,5)==0
    tmp{i}= sprintf('%d',i);
  end
end
if m>1, xticklabels(tmp); end

% export figure
fname= sprintf('figs/rxnfrq%d.eps',m);
q= input(sprintf('\n Export as %s? [y/N] ',fname),'s');
if q=='y'
  expfig(fname);
  fprintf(' -> Exported as %s\n',fname);
end
disp(' ')

