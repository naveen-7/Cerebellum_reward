%% Part of a batch process Required after CER_analysis_01
%% Analyses Behavior
%% by taking pre calculated CHANGE and LEARNT

% Created by NAVEEN ON 05/13/16 at CUMC


% function CER_analysis_02Br_AFTER

% clc;
% clear all;
% close all;


%% Setup directories--------------------------------------------------------
% codes_dir = 'C:\NAVEEN_Work\Cerebellum\Codes\CER_codes';
% data_dir  = 'C:\NAVEEN_Work\Cerebellum\Data';


cd(data_dir)
disp('!!! CER_analysis_02Br_AFTER has started running !!!');





%%


%Method 2: dis-continuous --------

n_bin = 0;  % number of bins

% % % w_bin = 20; % bin width
% % % s_bin = 8;  % shift in bin

w_bin = round(0.1*(size(Infos,1))); % bin width
s_bin = round(0.05*(size(Infos,1)));  % shift in bin

clear per_corr1 x_trial1 RT_bef
clear per_corr2 x_trial2 RT_aft
x_trial1 = [];
x_trial2 = []; 
per_corr1 = []; 
per_corr2 = []; 



% Part 1: before CHANGE

m_bin = round((CHANGE-1)/(s_bin))-1; % max number of bins


for i=1:s_bin:CHANGE-1
    if(i<=(m_bin-1)*s_bin)
        n_bin = n_bin+1;
        count = 0;
        RT_bef(n_bin) = nanmean(Infos(i:i+w_bin-1,14));
        RT_std_bef(n_bin) = nanstd(Infos(i:i+w_bin-1,14));
        for j=i:i+w_bin-1
            if Infos(j,10)==1
                count = count+1;
            end
        end
        per_corr1(n_bin) = (count/w_bin)*100;
        x_trial1(n_bin) = (i+(i+w_bin-1))/2;
    end
end


% Part 2: after CHANGE
n_bin = 0;
m_bin = round((size(Infos,1)-CHANGE)/(s_bin))-1; % max number of bins
% clear per_corr2 x_trial2

for i=CHANGE:s_bin:size(Infos,1)
    if(i<=((m_bin-2)*s_bin)+CHANGE)
        n_bin = n_bin+1;
        count = 0;
        RT_aft(n_bin) = nanmean(Infos(i:i+w_bin-1,14));
        RT_std_aft(n_bin) = nanstd(Infos(i:i+w_bin-1,14));
        for j=i:i+w_bin-1
            if Infos(j,10)==1
                count = count+1;
            end
        end
        per_corr2(n_bin) = (count/w_bin)*100;
        x_trial2(n_bin) = (i+(i+w_bin-1))/2;


    end
end



per_corr = NaN(length(per_corr1)+length(per_corr2),1);
per_corr(1:length(per_corr1)) = per_corr1;
per_corr(length(per_corr1)+1:length(per_corr1)+length(per_corr2)) = per_corr2;

x_trial = NaN(length(x_trial1)+length(x_trial2),1);
x_trial(1:length(x_trial1)) = x_trial1;
x_trial(length(x_trial1)+1:length(x_trial1)+length(x_trial2)) = x_trial2;

RT_tot = NaN(length(RT_bef)+length(RT_aft),1);
RT_tot(1:length(RT_bef)) = RT_bef;
RT_tot(length(RT_bef)+1:length(RT_bef)+length(RT_aft)) = RT_aft;

RT_std_tot = NaN(length(RT_std_bef)+length(RT_std_aft),1);
RT_std_tot(1:length(RT_std_bef)) = RT_std_bef;
RT_std_tot(length(RT_std_bef)+1:length(RT_std_bef)+length(RT_std_aft)) = RT_std_aft;



RT_tot_LC = RT_tot;
RT_std_tot_LC = RT_std_tot;
per_corr_LC = per_corr;
x_trial_LC = x_trial;


% Plotting learning curve

F = figure();
axes('box','off','tickdir','out','Linewidth',1.25,'FontSize',12)
hold on;
clear ylim yLim

plot(x_trial,per_corr,'-','color',[1	0.7255	0.0588],'LineWidth',3);
plot(x_trial,per_corr,'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
xlabel('Trial number','fontweight','bold');
ylabel('% Correct','fontweight','bold');

yLim = ylim;
if yLim(1)>20
    ylim([20 yLim(2)]);
end
yLim = ylim;
plot([CHANGE CHANGE], ylim,'--','lineWidth',2,'color',[0.4 0.4 0.4]);

text(CHANGE,yLim(1)+1,num2str(CHANGE));
title('Learning curve');
set(gca,'FontSize',12,'LineWidth',3)
set(gca,'fontweight','bold')


plot([LEARNT LEARNT], ylim,'--','lineWidth',2,'color',[0.7 0.7 0.7]);
text(LEARNT,yLim(1)+1,num2str(LEARNT));
% xlim([1 152.5]);
% ylim([10 100]);
% set(gca,'FontSize',14,'LineWidth',4)
% set(gca,'fontweight','bold')


for KK=1:length(ALL_CHANGES)
    plot([ALL_CHANGES(KK) ALL_CHANGES(KK)], ylim,'-','lineWidth',1,'color',[0.7 0.7 0.7]);
end

set(gcf, 'PaperUnits','inches','PaperSize',[8 8],'PaperPosition',[1 1 6.65 5])
cd(Results_dir)
filename = 'Learning_curve';
print(F, '-dpdf', filename, '-r400')
%delete(gcf)





% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %
save(MERGE_file,'CHANGE','per_corr_LC','x_trial_LC','LEARNT','w_bin','s_bin','RT_tot_LC','RT_std_tot_LC','-append');
save(ALLCELLS_file,'CHANGE','per_corr_LC','x_trial_LC','LEARNT','w_bin','s_bin','RT_tot_LC','RT_std_tot_LC','-append');
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %





% Calculating Normalized Learning Curve --------------------------


n_bin = 0;  % number of bins
w_bin = round(0.1*(size(Infos,1))); % bin width
s_bin = round(0.02*(size(Infos,1)));  % shift in bin


% Part 1: before CHANGE

m_bin = round((CHANGE-1)/(s_bin))-1; % max number of bins
clear per_corr1 x_trial1 RT_1

for i=1:s_bin:CHANGE-1
    if(i <= CHANGE-1+1-w_bin ) 
        n_bin = n_bin+1;
        count = 0;
        RT_1(n_bin) = nanmean(Infos(i:i+w_bin-1,14));
        for j=i:i+w_bin-1
            if Infos(j,10)==1
                count = count+1;
            end
        end
        per_corr1(n_bin) = (count/w_bin)*100;
        x_trial1(n_bin) = (i+(i+w_bin-1))/2;
    end
end


% Part 2: after CHANGE
n_bin = 0;
m_bin = round((size(Infos,1)-CHANGE)/(s_bin)); % max number of bins
clear per_corr2 x_trial2 RT_2

for i=CHANGE:s_bin:size(Infos,1)
    if(i <= size(Infos,1)+1-w_bin )                  %((m_bin-2)*s_bin+CHANGE))
%         i
        n_bin = n_bin+1;
        count = 0;
        RT_2(n_bin) = nanmean(Infos(i:i+w_bin-1,14));
        for j=i:i+w_bin-1
            if Infos(j,10)==1
                count = count+1;
            end
        end
        per_corr2(n_bin) = (count/w_bin)*100;
        x_trial2(n_bin) = (i+(i+w_bin-1))/2;

    end

end


RT_tot = NaN(length(RT_1)+length(RT_2),1);
RT_tot(1:length(RT_1)) = RT_1;
RT_tot(length(RT_1)+1:length(RT_1)+length(RT_2)) = RT_2;

per_corr = NaN(length(per_corr1)+length(per_corr2),1);
per_corr(1:length(per_corr1)) = per_corr1;
per_corr(length(per_corr1)+1:length(per_corr1)+length(per_corr2)) = per_corr2;

x_trial = NaN(length(x_trial1)+length(x_trial2),1);
x_trial(1:length(x_trial1)) = x_trial1;
x_trial(length(x_trial1)+1:length(x_trial1)+length(x_trial2)) = x_trial2;


RT_tot_NEWX = RT_tot;
per_corr_NEWX = per_corr;
x_trial_NEWX = x_trial;





x_trial_N = x_trial_NEWX-CHANGE;
x_trial_TEMP = x_trial_N((find(x_trial_NEWX>CHANGE,1)):(find(x_trial_NEWX>CHANGE,1) +round( 0.75*( length(x_trial_NEWX) - (find(x_trial_NEWX>CHANGE,1)-1) ) )));
x_trial_NEW = (x_trial_TEMP/(max(x_trial_TEMP)))*100;


per_corr_TEMP = per_corr_NEWX((find(x_trial_NEWX>CHANGE,1)):(find(x_trial_NEWX>CHANGE,1) +round( 0.75*( length(x_trial_NEWX) - (find(x_trial_NEWX>CHANGE,1)-1) ) )));
RT_tot_TEMP = RT_tot_NEWX((find(x_trial_NEWX>CHANGE,1)):(find(x_trial_NEWX>CHANGE,1) +round( 0.75*( length(x_trial_NEWX) - (find(x_trial_NEWX>CHANGE,1)-1) ) )));

clear x_trial_Nbin per_corr_Nbin RT_tot_Nbin

count = 0;
for i=1:12:96
    count = count+1;
ind = find(i<=x_trial_NEW & x_trial_NEW<=i+11);

if ~isempty(ind)
    per_corr_Nbin(count,1) = nanmean(per_corr_TEMP(ind));
    RT_tot_Nbin(count,1) = nanmean(RT_tot_TEMP(ind));
else
    per_corr_Nbin(count,1)=NaN;
    RT_tot_Nbin(count,1) = NaN;
end

x_trial_Nbin(count,1)=i+5;

end

% figure
% plot(x_trial_Nbin,per_corr_Nbin)
% 







w_bin = round(0.1*(size(Infos,1))); % bin width
s_bin = round(0.05*(size(Infos,1)));  % shift in bin



save(MERGE_file,'CHANGE','x_trial_Nbin','per_corr_Nbin','RT_tot_Nbin','LEARNT','-append');
save(ALLCELLS_file,'CHANGE','x_trial_Nbin','per_corr_Nbin','RT_tot_Nbin','LEARNT','-append');







%% simple RT analysis


F = figure; hold on;
for i=1:length(Infos)
    if Infos(i,10)==1
        plot(i,Infos(i,14),'og')
    elseif Infos(i,10)==2
        plot(i,Infos(i,14),'or')
    end
end

plot([CHANGE CHANGE],ylim,'-k')


IND_C = find(Infos(:,10)==1);
IND_W = find(Infos(:,10)==2);

plot(IND_C,Infos(IND_C,14),'-g');
plot(IND_W,Infos(IND_W,14),'-r');

xlabel('Trial number');
ylabel('RT');
box off;
title(strcat(NOME(1:8),'-',NOME(10:12),'-RT-CORR-WRNG'));
cd(Results_dir)
filename = strcat(NOME,'_RT_CORR_WRNG');

print(F, '-dpdf', filename, '-r600')
    













disp('!!! END OF CODE 02 !!!');











% end