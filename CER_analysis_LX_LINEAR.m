%% Code for analysing Licks for reward expectation


% Created by NAVEEN ON 08/31/17 at CUMC


% function CER_analysis_LX_LINEAR

% % % % % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % FOR REGULAR, RELEASE HERE -----------------------------------------

clc;
clear 
close all;


codes_dir = fullfile('C:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW','LINEAR');
data_dir  = fullfile('C:','NAVEEN_Work','Cerebellum','Data','ALL_CELLS');

cd(data_dir)
disp('!!! CER_analysis_L_LINEAR has started running !!!');



% LOADING FILE--------------------------------------------------------
disp('*************************************************')
disp('*************************************************')
disp('LOAD THE "M" DATA FILE')
[FileName1,PathName1] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
DATAfile_M = cat(2,PathName1,FileName1);
disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
load(DATAfile_M);
cd(PathName1);




% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % FOR BATCH, RELEASE HERE -----------------------------------------
% % % 
% % % close all;
% % % File = DATAfile;
% % % FileName = Tempy;
% % % 
% % % % % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%
%% ONE: SIMPLE SPIKES ----------------------------------------------------


disp('*************************************************');
disp('*************************************************');
disp('!!! Analysing simple spikes !!!');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1: PSTHs --------------------



F = figure();

START = -1250;
END = 650;
StartLIM = START+150;
EndLIM = END-150;
temp_time = START:END;


subplot(2,2,1)
hold on;
Signal = Licks; Align = Infos(:,4); Start = START; End = END;
Cured = time2logic_n(Signal,Align,Start,End);
Percent = smooth(nansum(Cured)/size(Cured,1));
plot(START:END,Percent,'color',[0.5 0.5 0.5])
xlim([START END])


subplot(2,2,2)
hold on;
Signal = Licks; Align = Infos(:,11); Start = START; End = END;
Cured = time2logic_n(Signal,Align,Start,End);
Percent = smooth(nansum(Cured)/size(Cured,1));
plot(START:END,Percent,'color',[0.5 0.5 0.5])
xlim([START END])








% % % % % % 
% % % % % % P1 = PSTHe_n(Licks,Infos(:,11),START,END,50,[1.0000    0.7255    0.0588],1,1); % Aligned to reward
% % % % % % temp_ind = find(StartLIM<=temp_time & temp_time<=EndLIM);
% % % % % % errorline_n(StartLIM:EndLIM,P1(1,temp_ind),P1(3,temp_ind),1,[0.5 0.5 0.5],0.3,0,1)
% % % % % % xlim([StartLIM EndLIM])
% % % % % % hold on;
% % % % % % 
% % % % % % MIN = nanmin([nanmin(P1(1,find(temp_time>=-1100 & temp_time<=500))-P1(3,find(temp_time>=-1100 & temp_time<=500)))])-2;
% % % % % % 
% % % % % % MAX = nanmax([nanmax(P1(1,find(temp_time>=-1100 & temp_time<=500))+P1(3,find(temp_time>=-1100 & temp_time<=500)))])+2;
% % % % % % 
% % % % % % ylim([MIN MAX])
% % % % % % % 
% % % % % % % plot([StartLIM StartLIM],[MIN MIN+20],'-k','linewidth',1)
% % % % % % % plot([StartLIM StartLIM+200],[MIN MIN],'-k','linewidth',1)
% % % % % % % axis off;
% % % % % % 
% % % % % % plot([0 0],[MIN MIN],'OK','MarkerFacecolor','k','markersize',4)

suptitle(strcat('LICK-',NOME(1:8),'-',NOME(10:12)));


cd(Results_dir)
filename = 'Licking';
print(F, '-dpdf', filename, '-r400')






