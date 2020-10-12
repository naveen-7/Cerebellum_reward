%% Code for analysing licks for reward amount


% Created by NAVEEN ON 08/31/17 at CUMC


% function CER_analysis_L_LINEAR

% % % % % % % % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % FOR REGULAR, RELEASE HERE -----------------------------------------
% % % 
% % % clc;
% % % clear 
% % % close all;
% % % 
% % % 
% % % codes_dir = fullfile('C:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW','LINEAR');
% % % data_dir  = fullfile('C:','NAVEEN_Work','Cerebellum','Data','ALL_CELLS');
% % % 
% % % cd(data_dir)
% % % disp('!!! CER_analysis_L_LINEAR has started running !!!');
% % % 
% % % 
% % % 
% % % % LOADING FILE--------------------------------------------------------
% % % disp('*************************************************')
% % % disp('*************************************************')
% % % disp('LOAD THE "M" DATA FILE')
% % % [FileName1,PathName1] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
% % % DATAfile_M = cat(2,PathName1,FileName1);
% % % disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
% % % load(DATAfile_M);
% % % cd(PathName1);
% % % 



% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % FOR BATCH, RELEASE HERE -----------------------------------------

close all;
File = DATAfile;
FileName = Tempy;

% % % % % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%
%% ONE: SIMPLE SPIKES ----------------------------------------------------


disp('*************************************************');
disp('*************************************************');
disp('!!! Analysing simple spikes !!!');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1: PSTHs --------------------

TEMP = 1./diff(Licks{1,1});
figure();
plot(TEMP)
ylim([0.9995 1])

























%%%%% probability of licking --------------------------

F = figure();


Start_T = -450;
End_T  =1050;
StartLIM_T = -400; EndLIM_T = 800;

% Start_M = -750;
% End_M  =750;
% StartLIM_M = -500; EndLIM_M = 700;


Start_M = -1050;
End_M  =1050;
StartLIM_M = -900; EndLIM_M = 900;


% seperate into reward amounts ------------------


IND_N = find(Infos(:,17)<=5);
IND_S = union(find(Infos(:,17)==50),find(Infos(:,17)==110));
IND_L = find(Infos(:,17)>200);

COLOR_N = [1 0 0];
COLOR_S = [0 1 0];
COLOR_L = [0 0 1];



START = Start_M; END = End_M;
StartLIM = StartLIM_M; EndLIM=EndLIM_M;
temp_time = START:END;


subplot(2,2,1)
hold on;

Signal = Licks(IND_N); Align = Infos(IND_N,11); Start = Start_M; End = End_M;
Cured = time2logic_n(Signal,Align,Start,End);
Percent_N = smooth(nansum(Cured)/size(Cured,1));
plot(START:END,Percent_N,'k')

Signal = Licks(IND_S); Align = Infos(IND_S,11); Start = Start_M; End = End_M;
Cured = time2logic_n(Signal,Align,Start,End);
Percent_S = smooth(nansum(Cured)/size(Cured,1));
plot(START:END,Percent_S,'r')

Signal = Licks(IND_L); Align = Infos(IND_L,11); Start = Start_M; End = End_M;
Cured = time2logic_n(Signal,Align,Start,End);
Percent_L = smooth(nansum(Cured)/size(Cured,1));
plot(START:END,Percent_L,'b')

xlim([START END])



% Percent_TOT = (Percent_N+Percent_S+Percent_L)/3;
% plot(START:END,Percent_L,'b')




suptitle(strcat('LICK-',NOME(1:8),'-',NOME(10:12)));

cd(Results_dir)
filename = 'Licking_PERCENT';
print(F, '-dpdf', filename, '-r400')






