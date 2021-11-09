

%% ANALYSE_LICK
%% Written by naveen at JLG on 12/10/18


clc;
clear all;
close all;


% Setup directories--------------------------------------------------------

codes_dir = fullfile('e:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW','LINEAR');
data_dir  = fullfile('e:','NAVEEN_Work','Cerebellum','Data','MERGED_CELLS');

cd(data_dir)

% LOADING FILE------------------------------------------------------- %% 
disp('*************************************************')
disp('*************************************************')
disp('LOAD THE CELL DATA FILE')
[FileName1,PathName1] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
DATAfile = cat(2,PathName1,FileName1);
disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
load(DATAfile);
cd(PathName1);





for i=1:size(LICK,1)
   LICK_spikes{i,1} = find(LICK{i,1}<1); 
end


F = figure()

Sigma = 20;


subplot(3,3,1)
hold on;
Start_time = -400; End_time = 800;

Signal = LICK_spikes(Infos(1:CHANGE,10)==1);
Align_time = Infos(Infos(1:CHANGE,10)==1,4);
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,[0 0 1],1,0);
xlim([-200 600])

Signal = LICK_spikes(Infos(1:CHANGE,10)==2);
Align_time = Infos(Infos(1:CHANGE,10)==2,4);
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,[1 0 0],1,0);
xlim([-200 600])


subplot(3,3,2)
hold on;
Start_time = -600; End_time = 800;

Signal = LICK_spikes(Infos(1:CHANGE,10)==1);
Align_time = Infos(Infos(1:CHANGE,10)==1,11);
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,[0 0 1],1,0);
xlim([-400 600])

Signal = LICK_spikes(Infos(1:CHANGE,10)==2);
Align_time = Infos(Infos(1:CHANGE,10)==2,11);
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,[1 0 0],1,0);
xlim([-400 600])





subplot(3,3,4)
hold on;
Start_time = -400; End_time = 800;

Signal = LICK_spikes(Infos(CHANGE+1:CHANGE+20,10)==1);
Align_time = Infos(Infos(CHANGE+1:CHANGE+20,10)==1,4);
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,[0 0 1],1,0);
xlim([-200 600])

Signal = LICK_spikes(Infos(CHANGE+1:CHANGE+20,10)==2);
Align_time = Infos(Infos(CHANGE+1:CHANGE+20,10)==2,4);
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,[1 0 0],1,0);
xlim([-200 600])


subplot(3,3,5)
hold on;
Start_time = -600; End_time = 800;

Signal = LICK_spikes(Infos(CHANGE+1:CHANGE+20,10)==1);
Align_time = Infos(Infos(CHANGE+1:CHANGE+20,10)==1,11);
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,[0 0 1],1,0);
xlim([-400 600])

Signal = LICK_spikes(Infos(CHANGE+1:CHANGE+20,10)==2);
Align_time = Infos(Infos(CHANGE+1:CHANGE+20,10)==2,11);
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,[1 0 0],1,0);
xlim([-400 600])





subplot(3,3,7)
hold on;
Start_time = -400; End_time = 800;

Signal = LICK_spikes(Infos(LEARNT-20:length(Infos),10)==1);
Align_time = Infos(Infos(LEARNT-20:length(Infos),10)==1,4);
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,[0 0 1],1,0);
xlim([-200 600])

Signal = LICK_spikes(Infos(LEARNT-20:length(Infos),10)==2);
Align_time = Infos(Infos(LEARNT-20:length(Infos),10)==2,4);
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,[1 0 0],1,0);
xlim([-200 600])


subplot(3,3,8)
hold on;
Start_time = -600; End_time = 800;

Signal = LICK_spikes(Infos(LEARNT-20:length(Infos),10)==1);
Align_time = Infos(Infos(LEARNT-20:length(Infos),10)==1,11);
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,[0 0 1],1,0);
xlim([-400 600])

Signal = LICK_spikes(Infos(LEARNT-20:length(Infos),10)==2);
Align_time = Infos(Infos(LEARNT-20:length(Infos),10)==2,11);
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,[1 0 0],1,0);
xlim([-400 600])


cd(Results_dir)
filename = 'LICK';
print(F, '-dpdf', filename, '-r400')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% AVG LICKING



F = figure()

Sigma = 20;


Start_time = -1250;
End_time = 650;
StartLIM = Start_time+150;
EndLIM = End_time-150;






% subplot(2,2,1)
% hold on;
% % Start_time = -400; End_time = 800;
% 
% Signal = LICK_spikes(Infos(:,10)==1);
% Align_time = Infos(Infos(:,10)==1,4);
% PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,[0 0 1],1,0);
% xlim([-200 600])
% 
% Signal = LICK_spikes(Infos(:,10)==2);
% Align_time = Infos(Infos(:,10)==2,4);
% PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,[1 0 0],1,0);
% xlim([-200 600])


subplot(2,2,2)
hold on;
% Start_time = -600; End_time = 800;

Signal = LICK_spikes(Infos(:,10)==1);
Align_time = Infos(Infos(:,10)==1,11);
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,[0 0 1],1,0);
xlim([StartLIM EndLIM])

Signal = LICK_spikes(Infos(:,10)==2);
Align_time = Infos(Infos(:,10)==2,11);
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,[1 0 0],1,0);
xlim([StartLIM EndLIM])




cd(Results_dir)
filename = 'LICK_AVG2';
print(F, '-dpdf', filename, '-r400')

