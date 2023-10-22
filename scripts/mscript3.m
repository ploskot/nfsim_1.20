% Pavel Loskot, October 2022
% experiment 4: query cause, get effects

clearvars -except m atleast

atleast= 1; % minimum occurrence
rxn_type= 2; % cluster to query
WW1= [5 10 15 20]; % query (cause) sub-sequence length
WW2= 1:15; % response (effect) sub-sequence length

Models= ['A' 'B' 'C' 'D' 'E'];
model= Models(m);
fprintf('Model %s: %d, %d\n', model,atleast,rxn_type);

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
[~,ii]= sort(frq,'descend');
rxn1= rxn(ii(ci1)); 
rxn2= rxn(ii(ci2)); 
rxn3= rxn(ii(ci3));
% check
tmp= [rxn1;rxn2;rxn3];
if length(tmp)~=Nr || length(intersect(tmp,rxn))~=Nr
  error('Reaction clustering error.')
end
if isempty(rxn3) && rxn_type==3
  error('Rxn3 is empty.')
end

% calculation
nW1= length(WW1);
nW2= length(WW2);
for w1=1:nW1
  W1= WW1(w1);
  report(W1)
  PP= zeros(nW2,3);
  for w2=1:nW2
    W2= WW2(w2);
    W= W1+W2;
    ii1= 1:W1;
    ii2= (W1+1):W;  
    x= floor(N/W);
    dat2= reshape(dat1(1:W*x),W,x)';
    % query
    Q1= zeros(x,1);
    switch rxn_type
      case 1
        for i=1:x
          Q1(i)= length(find(ismember(dat2(i,ii1),[rxn1]))); % select
        end
      case 2
        for i=1:x
          Q1(i)= length(find(ismember(dat2(i,ii1),[rxn2]))); % select
        end
      case 3
        for i=1:x
          Q1(i)= length(find(ismember(dat2(i,ii1),[rxn3]))); % select
        end
    end
    jj1= find(Q1>=atleast); % at least
    jj2= find(Q1<atleast);
    n1= length(jj1);  
    %if n1==0, continue; end
    n2= length(jj2);
    if n1+n2~=x
      error('Wrong jj1 or jj2 length.')
    end
    resp1= dat2(jj1,ii2); % cause
    resp2= dat2(jj2,ii2); % not cause

    allrxn= [rxn1;rxn2;rxn3];
    NR1= zeros(n1,Nr);
    NR2= zeros(n2,Nr);
    for i=1:Nr
      NR1(:,i)= sum(resp1==allrxn(i),2);
      NR2(:,i)= sum(resp2==allrxn(i),2);      
    end
    xNR1= unique(NR1,'rows');
    xNR2= unique(NR2,'rows');
    xn1= size(xNR1,1);
    xn2= size(xNR2,1);
    yNR12= intersect(xNR1,xNR2,'rows');
    yNR1= setdiff(xNR1,yNR12,'rows');
    yNR2= setdiff(xNR2,yNR12,'rows');    
    yn1= size(yNR1,1);
    yn2= size(yNR2,1);
    yn12= size(yNR12,1);
    if (xn1+xn2~=yn1+yn2+2*yn12)
      error('Check NR1 and NR2.')
    end
    xyn= yn1+yn2+yn12;
    if xn1==0, p1= '-'; PP(w2,1)= 0; else
      PP(w2,1)= 100*yn1/xyn; p1= sprintf('%.1f%%',PP(w2,1)); 
    end
    if xn2==0, p2= '-'; PP(w2,2)= 0; else
      PP(w2,2)= 100*yn2/xyn; p2= sprintf('%.1f%%',PP(w2,2));
    end
    if yn12==0, p3= '-'; PP(w2,3)= 0; else
      PP(w2,3)= 100*yn12/xyn; p3= sprintf('%.f%%',PP(w2,3));
    end
    %fprintf(' W2=%-2d n1=%-6d n2=%-5d %3d(%s)/%d %3d(%s)/%d %d(%s)\n', ...
    %        W2,n1,n2,yn1,p1,xn1,yn2,p2,xn2,yn12,p3);
    fprintf('.');
  end 
  fprintf('\n');

  % bar plot
  subplot(2,2,w1)
  bar(PP,'stacked');grid;axis tight
  if w1==3 || w1==4, xlabel('xx'); end %xlabel('W_2')
  if w1==1 || w1==3, ylabel('yy'); end % ylabel('%')
  title(sprintf('tt%d',w1)); % title(sprintf('W_1=%d',W1))
  xticklabels({'1','','','','5','','','','','10','','','','','15'})
end

% export figure
fname= sprintf('figs/pict1%d-%d%d.eps',m,rxn_type,atleast);
q= input(sprintf('\n Export as %s? [y/N] ',fname),'s');
if q=='y'
  expfig(fname);
  fprintf(' -> Exported as %s\n',fname);
end
disp(' ')


