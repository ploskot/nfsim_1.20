% Pavel Loskot, October 2022
% experiment 7: generate anti-causal statements

clear all

m= 1;
atleast= 1;
rxn_type= 3;
W2= 25; % query (effect)
W1= 25; % response (cause)

Models= ['A' 'B' 'C' 'D' 'E'];
model= Models(m);
report('Model',model,': ',atleast, rxn_type)

if (atleast+5)>W1 || (atleast+5)>W2
  warning('atleast may be too large');
end

dat= readmatrix(['model' model '.csv']);
dat= dat(:,1:2);
dat1= dat(:,2);
[rxn,frq]= stat(dat1);
N= length(dat1);
Nr= length(rxn);
% check
if isequal(sort(rxn),0:Nr-1)
  error('Missing reaction index(es).')
end

% reaction clusters [see m1a]
% model A
ci1= [1 2];
ci2= [3 4];
ci3= [5 6];

[~,ii]= sort(frq,'descend');
rxn1= rxn(ii(ci1)); 
rxn2= rxn(ii(ci2)); 
rxn3= rxn(ii(ci3));
% check
tmp= [rxn1;rxn2;rxn3];
if length(tmp)~=Nr || length(intersect(tmp,rxn))~=Nr
  error('Reaction clustering error.')
end

PP= zeros(1,3);
W= W1+W2;
ii1= 1:W1;
ii2= (W1+1):W;  
x= floor(N/W);
dat2= reshape(dat1(1:W*x),W,x)';
% query
Q2= zeros(x,1);
switch rxn_type
  case 1
    for i=1:x
      Q2(i)= length(find(ismember(dat2(i,ii2),[rxn1]))); % select
    end
  case 2
    for i=1:x
      Q2(i)= length(find(ismember(dat2(i,ii2),[rxn2]))); % select
    end
  case 3
    for i=1:x
      Q2(i)= length(find(ismember(dat2(i,ii2),[rxn3]))); % select
    end
end
jj1= find(Q2>=atleast); % at least
jj2= find(Q2<atleast);

n1= length(jj1);
%if n1==0, continue; end    
n2= length(jj2);
if n1+n2~=x
  error('Wrong jj1 or jj2 length.')
end

% CAUSE (response)
resp1= dat2(jj1,ii1); % cause
resp2= dat2(jj2,ii1); % not cause
allrxn= [rxn1;rxn2;rxn3];
NR1= zeros(n1,Nr);
NR2= zeros(n2,Nr);
for i=1:Nr
  NR1(:,i)= sum(resp1==allrxn(i),2);
  NR2(:,i)= sum(resp2==allrxn(i),2);      
end
% unique
xNR1= unique(NR1,'rows');
xNR2= unique(NR2,'rows');
xn1= size(xNR1,1);
xn2= size(xNR2,1);
yNR12= intersect(xNR1,xNR2,'rows');
yNR1= setdiff(xNR1,yNR12,'rows');
yNR2= setdiff(xNR2,yNR12,'rows');    

% effect/response stats
tmp1= yNR1(:,ci1);    
cmin1= min(tmp1,[],1);
cmax1= max(tmp1,[],1);

tmp2= yNR1(:,ci2);    
cmin2= min(tmp2,[],1);
cmax2= max(tmp2,[],1);

tmp3= yNR1(:,ci3);    
cmin3= min(tmp3,[],1);
cmax3= max(tmp3,[],1);

% percentages
yn1= size(yNR1,1);
yn2= size(yNR2,1);
yn12= size(yNR12,1);
if (xn1+xn2~=yn1+yn2+2*yn12), error('Check NR1 and NR2.'); end
xyn= yn1+yn2+yn12;
if xn1==0, p1= '0.0%'; PP(1)= 0; else
  PP(1)= 100*yn1/xyn; p1= sprintf('%.1f%%',PP(1)); 
end
if xn2==0, p2= '0.0%'; PP(2)= 0; else
  PP(2)= 100*yn2/xyn; p2= sprintf('%.1f%%',PP(2));
end
if yn12==0, p3= '0.0%'; PP(3)= 0; else
  PP(3)= 100*yn12/xyn; p3= sprintf('%.f%%',PP(3));
end

% prints
fprintf(' n1=%-6d n2=%-5d %3d(%s)/%d %3d(%s)/%d %d(%s)\n', ...
        n1,n2,yn1,p1,xn1,yn2,p2,xn2,yn12,p3);
fprintf('n1=%d, xn1=%d\n',n1,xn1);
tmp= [cmin1,cmin2,cmin3;cmax1,cmax2,cmax3];
for i=1:Nr
  fprintf('R%d: %2d-%d\n',i,tmp(1,i),tmp(2,i));
end



