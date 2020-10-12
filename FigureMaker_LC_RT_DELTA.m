% script to make a figure for LC and RT and Delta and correlate them.
% follow up code to CER_analysis_03b_LINEAR

% written by naveen at cumc on 7/27/17

%% FigureMaker_LC_RT_DELTA


%% First make LC and RT ---------------

Type8 = NaN(length(Infos),1);

for i=2:length(Infos)
    
    if (Infos(i-1,10)==1)
        Type8(i,1)=1;
    end
    
end



Signal = Infos(:,14);
W_bin = 15; % bin width
S_bin = 5;  % shift in bin
TYPE = 'a';

[x_trial, SIG_trial, per_corr] = LEARNING_CURVE_n(Signal, CHANGE, Infos, W_bin, S_bin, TYPE);

F = figure();
axes('box','off','tickdir','out','Linewidth',1.25,'FontSize',12)

subplot(3,3,1)
hold on;
clear ylim yLim

plot(x_trial,per_corr,'-','color',[0.5 0.5 0.5],'LineWidth',1.5);
plot(x_trial,per_corr,'o','MarkerSize',3,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
set(gca,'XTick',[])
set(gca,'XColor','w')
ylabel('% Correct');
xlim([1 x_trial(length(x_trial))+W_bin])
ylim([min(min(per_corr),20) 100]);

plot([CHANGE CHANGE], ylim,'--','lineWidth',1,'color',[0.4 0.4 0.4]);
set(gca,'FontSize',10,'LineWidth',1)


subplot(3,3,4)
hold on;
plot(x_trial,SIG_trial(:,1),'-','color',[1 0.4 0.4],'LineWidth',1.5);
errorbar(x_trial,SIG_trial(:,1),SIG_trial(:,3),'.k');
plot(x_trial,SIG_trial(:,1),'o','MarkerSize',3,'MarkerFaceColor',[1 0.4 0.4],'color','k');
xlabel([]);
ylabel('RT');
set(gca,'XTick',[])
set(gca,'XColor','w')
xlim([1 x_trial(length(x_trial))+W_bin])
ylim([nanmin(SIG_trial(:,1)-SIG_trial(:,3)) nanmax(SIG_trial(:,1)+SIG_trial(:,3))]);
plot([CHANGE CHANGE], ylim,'--','lineWidth',1,'color',[0.4 0.4 0.4]);
set(gca,'FontSize',10,'LineWidth',1)







%% Then do delta
    



if strcmp(ALIGN,'T')
    P_USE = P_T;
    TIME = time_T;
elseif strcmp(ALIGN,'M')
    P_USE = P_M;
    TIME = time_M;
end


temp_ind = find(TIME==EPOCH_START) : find(TIME==EPOCH_END);
clear P_delta;
for i=1:size(P_USE,1)
    P_delta(i,1) = nanmean(P_USE(i,temp_ind));
end



clear x_trial1 P_bef_DELTA1 P_curr_DELTA1
n_bin=0;
for i=1:S_bin:CHANGE-1
    if(i <= CHANGE-1+1-W_bin )
        n_bin = n_bin+1;
        temp_trials = find(Type8(i:i+W_bin-1,1)>=1)+i-1;
        P_bef_DELTA1(n_bin,1) = nanmean(P_delta(temp_trials-1));
        P_bef_DELTA1(n_bin,2) = nanstd(P_delta(temp_trials-1))/(sqrt(length(temp_trials)));
        
        P_curr_DELTA1(n_bin,1) = nanmean(P_delta(temp_trials));
        P_curr_DELTA1(n_bin,2) = nanstd(P_delta(temp_trials))/(sqrt(length(temp_trials)));
    end
end


% Part 2: after CHANGE
n_bin = 0;
clear per_corr2 P_bef_DELTA2 P_curr_DELTA2

for i=CHANGE:S_bin:size(Infos,1)
    if(i <= size(Infos,1)+1-W_bin )
        n_bin = n_bin+1;
        temp_trials = find(Type8(i:i+W_bin-1,1)>=1)+i-1;
        P_bef_DELTA2(n_bin,1) = nanmean(P_delta(temp_trials-1));
        P_bef_DELTA2(n_bin,2) = nanstd(P_delta(temp_trials-1))/(sqrt(length(temp_trials)));
        
        P_curr_DELTA2(n_bin,1) = nanmean(P_delta(temp_trials));
        P_curr_DELTA2(n_bin,2) = nanstd(P_delta(temp_trials))/(sqrt(length(temp_trials)));
    end
end

clear P_DELTA
P_DELTA = [P_curr_DELTA1; P_curr_DELTA2]-[P_bef_DELTA1; P_bef_DELTA2];







subplot(3,3,7)
hold on;

plot(x_trial,P_DELTA(:,1),'-','color',[1 0.4 0.4],'LineWidth',1.5);
errorbar(x_trial,P_DELTA(:,1),P_DELTA(:,2),'.k');
plot(x_trial,P_DELTA(:,1),'o','MarkerSize',3,'MarkerFaceColor',[1 0.4 0.4],'color','k');
xlabel('Trial number');
ylabel('Delta');
xlim([1 x_trial(length(x_trial))+W_bin])
ylim([nanmin(P_DELTA(:,1)-abs(P_DELTA(:,2))) nanmax(P_DELTA(:,1)+P_DELTA(:,2))]);
plot([CHANGE CHANGE], ylim,'--','lineWidth',1,'color',[0.4 0.4 0.4]);

set(gca,'FontSize',10,'LineWidth',1)
text(CHANGE,nanmin(P_DELTA(:,1)-P_DELTA(:,2)),num2str(CHANGE));   
    



[CORR_DELTA.LC_RT p_DELTA.LC_RT] = (corr(per_corr(:,1),SIG_trial(:,1)));
[CORR_DELTA.LC_DL p_DELTA.LC_DL] = (corr(per_corr(:,1),P_DELTA(:,1)));
[CORR_DELTA.RT_DL p_DELTA.RT_DL] = (corr(SIG_trial(:,1),P_DELTA(:,1)));



subplot(3,3,5)
axis off;
text(0.5,1,strcat('CORR'))
text(0.5,0.8,strcat('LC-RT = ',num2str(CORR_DELTA.LC_RT)))
text(0.5,0.6,strcat('LC-DL = ',num2str(CORR_DELTA.LC_DL)))   
text(0.5,0.4,strcat('DL-RT = ',num2str(CORR_DELTA.RT_DL)))   
    
subplot(3,3,8)
axis off;
text(0.5,1,strcat('p'))
text(0.5,0.8,strcat('LC-RT = ',num2str(p_DELTA.LC_RT)))
text(0.5,0.6,strcat('LC-DL = ',num2str(p_DELTA.LC_DL)))   
text(0.5,0.4,strcat('DL-RT = ',num2str(p_DELTA.RT_DL)))     
    
    
    
   

cd(Results_dir)
filename = strcat(FileName1(6:17),'_LC_RT_DL');

print(F, '-dpdf', filename, '-r600')


DELTA.CORR_RHO(IND)=CORR_DELTA;
DELTA.CORR_P(IND)=p_DELTA;   


% save(MERGE_file,'CORR_DELTA','p_DELTA','-append');
% save(ALLCELLS_file,'CORR_DELTA','p_DELTA','-append');




