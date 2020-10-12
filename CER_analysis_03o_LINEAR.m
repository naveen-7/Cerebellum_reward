%% function CER_analysis_03o_LINEAR
% analyses behavior critically ------------------
% written by naveen at JLG on 2/8/19

outcometable(:,1) = Infos(:,5);  % symbol
outcometable(:,2) = Infos(:,9); % movt
outcometable(:,3) = Infos(:,10); % reward

% % % WHICHHAND(find(outcometable(:,2)==0 & outcometable(:,3)==1),1)=1;
% % % WHICHHAND(find(outcometable(:,2)==1 & outcometable(:,3)==1),1)=2;
% % % WHICHHAND(find(outcometable(:,2)==1 & outcometable(:,3)==2),1)=1;
% % % WHICHHAND(find(outcometable(:,2)==0 & outcometable(:,3)==2),1)=2;
% % % 
% % % outcometable(:,2)=WHICHHAND;




if strcmp(TASK_TYPE,'N')
    CHANGE = find(abs(diff(outcometable(:,1)))>112,1)+1;
end

% UN = (unique(outcometable(:,1)));

PARSED = outcometable(CHANGE:end,:);
UN = (unique(PARSED(:,1)));
for i=1:length(UN)
   PARSED((PARSED(:,1)==UN(i)))=i; 
   SYMOUT{i,1} = PARSED((PARSED(:,1)==i),3);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c = jet(length(UN));

F = figure();
subplot(2,2,1)

w_bin = 5; % bin width
s_bin = 3;  % shift in bin
for i=1:length(UN)
    [x_trial_LC, per_corr_LC] = LCparams(SYMOUT{i,1},w_bin,s_bin,0);
    X_LC{i,1}=x_trial_LC;
    Y_LC{i,1}=per_corr_LC;
    
    
    hold on;
    plot(X_LC{i,1},Y_LC{i,1},'-','color',c(i,:),'LineWidth',0.7);
    plot(X_LC{i,1},Y_LC{i,1},'o','MarkerSize',4,'color','k');
   
    
end


xlabel('number of occurences'); ylabel('% Correct');

yLim = ylim;
if yLim(1)>20
    ylim([20 yLim(2)]);
end
yLim = ylim;
ylim([yLim(1) 100]);


cd(Results_dir)
filename = strcat(NOME,'_LCparams');
print(F, '-dpdf', filename, '-r400')

LCPARAMS.X = X_LC;
LCPARAMS.Y = Y_LC;



%% SOME NUMBERS %%%%%%%%%%%%%%%%%%%%%%%


h1 = PARSED(find(PARSED(:,1)==1 & PARSED(:,3)==1,1),2); % hand for symbol 1
h2 = PARSED(find(PARSED(:,1)==2 & PARSED(:,3)==1,1),2); % hand for symbol 1




n_bin = 0;  % number of bins
w_bin = 5; % bin width
s_bin = 3;  % shift in bin

% w_bin = 15; % bin width
% s_bin = 8;  % shift in bin



SIGNAL = PARSED(:,[1 3]);
% Part 2: after CHANGE
n_bin = 0;
m_bin = round(length(SIGNAL)/(s_bin))-1; % max number of bins
% clear per_corr2 x_trial2
clear winstay loseswitch winstay_next loseswitch_next
for i=1:s_bin:length(SIGNAL)
     if(i<=((m_bin-2)*s_bin)+0)
        n_bin = n_bin+1;
        signal = SIGNAL(i:i+w_bin-1,:);
        
        %%%%%%% winstay for any trials in a bin
        clear ind
        ind = find(signal(:,2)==1)+1;
        winstay(n_bin) = 0;
        if length(ind)>0
            if ind(end)>signal(end) ind = ind(1:end-1); end
            winstay(n_bin) = length(find(signal(ind,2)==1))/length(ind);
        end
        
        %%%%%%% winstay and loseswitch for pairs of trials
        clear ind
        DIFF = diff(signal(:,1));
        winstay_next(n_bin) = 0;
        loseswitch_next(n_bin)= 0;
        ind = find(DIFF==0);
          % given a trial is correct, what is the prob that an immedietly next trial with the same stimulus is correct?
          if length(ind)>0
              for jj=1:length(ind)
                  
                  if signal(ind(jj),2)==1
                      if signal(ind(jj)+1,2)==1
                          winstay_next(n_bin)=winstay_next(n_bin)+1;
                      end
                  end
                  
                  if signal(ind(jj),2)==2
                      if signal(ind(jj)+1,2)==1
                          loseswitch_next(n_bin)=loseswitch_next(n_bin)+1;
                      end
                  end
                  
                  
              end
              
              winstay_next(n_bin) = winstay_next(n_bin)/length(ind);
              loseswitch_next(n_bin)=loseswitch_next(n_bin)/length(ind);
          end
        
       
        
        clear ind
        ind = find(signal(:,2)==2)+1;
        loseswitch(n_bin) = 0;
        if length(ind)>0
            if ind(end)>signal(end) ind = ind(1:end-1); end
            loseswitch(n_bin) = length(find(signal(ind,2)==1))/length(ind);
        end
        
        trialcount(n_bin)=(i+(i+w_bin-1))/2;
        
     end
    
    
    
    
end




F = figure();
hold on;
plot(trialcount,winstay,'-b');
plot(trialcount,winstay_next,':b');
plot(trialcount,loseswitch,'-r');
plot(trialcount,loseswitch_next,':r');

cd(Results_dir)
filename = strcat(NOME,'_WSLS');
print(F, '-dpdf', filename, '-r400')


LCPARAMS.WS = winstay;
LCPARAMS.WS_next = winstay_next;
LCPARAMS.LS = loseswitch;
LCPARAMS.LS_next = loseswitch_next;
LCPARAMS.trialcount = trialcount;











 
%% REACTION TIME

RT = Infos(:,14);

RT = RT(CHANGE:CHANGE+20);
RTcomp{1,:} = RT(find(RT<median(RT)));
RTcomp{2,:} = RT(find(RT>median(RT)));

[~,p] = ttest(RTcomp{1,:},RTcomp{2,:})





save(POP_file,'LCPARAMS','-append');
save(MERGE_file,'LCPARAMS','-append');
save(ALLCELLS_file,'LCPARAMS','-append');

