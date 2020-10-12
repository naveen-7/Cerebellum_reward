%% To merge the rex and spike 2 files for Lick 

% Written by naveen at cumc on 8/31/17

% CER_MERGE_LL

close all;
clear;
clc;

% Setup directories--------------------------------------------------------
codes_dir = fullfile('C:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW','LINEAR');
data_dir  = fullfile('C:','NAVEEN_Work','Cerebellum','Data','RECORDED_CELLS');

disp('!!!!!  CER_MERGE has started running  !!!!!')
cd('C:\NAVEEN_Work\Cerebellum\Data')
mkdir('MERGED_CELLS');
MERGE_dir = fullfile('C:','NAVEEN_Work','Cerebellum','Data','MERGED_CELLS');
ALLCELLS_dir = fullfile('C:','NAVEEN_Work','Cerebellum','Data','ALL_CELLS');
POP_dir = fullfile('C:','NAVEEN_Work','Cerebellum','Data','POP_CELLS');




cd(data_dir);


% LOADING FILE--------------------------------------------------------
disp('*************************************************')
disp('*************************************************')
disp('LOAD THE REX FILE')
[FileName1,PathName1] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
REXfile = cat(2,PathName1,FileName1);
disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
load(REXfile);
% Infos = load(REXfile);
cd(PathName1);


% LOADING Spike2 - Trigger FILE------------------------------------------------------- %% ch5
disp('*************************************************')
disp('*************************************************')
disp('LOAD Spk2-Trigg FILE')
[FileName2,PathName2] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
disp(strcat('!!!!!','File you entered is :',FileName2,' !!!!!'));

clear A names;
A = load(cat(2,PathName2,FileName2));
names = fieldnames(A);
Trigg = A.(names{1,1});
cd(PathName1);





% LOADING Spike2 - LICK FILE------------------------------------------------------- %% ch2
disp('*************************************************')
disp('*************************************************')
disp('LOAD Spk2-LICK FILE')
[FileName3,PathName3] = uigetfile('*.txt','File to Append');   % Open standard dialog box for retrieving files
disp(strcat('!!!!!','File you entered is :',FileName3,' !!!!!'));

LICK_DATA_IMPORT = importdata(FileName3);

% clear A names;
% A = load(cat(2,PathName3,FileName3));
% names = fieldnames(A);
% LICK_DATA_IMPORT = A.(names{1,1});
% cd(PathName1);




clear DATA_LICK
temp_LICK = LICK_DATA_IMPORT.data;



%%%% ALL LICK TIMES --------------------------------------------------
LICK_pos = temp_LICK(find(temp_LICK(:,2)==1),1);

ALL_LICKS = temp_LICK;

for i=2:2:length(LICK_pos)
    ind1 = find(temp_LICK==LICK_pos(i-1));
    ind2 = find(temp_LICK==LICK_pos(i));
    ALL_LICKS(ind1:ind2,2)=1;
end

LICK_TIMES = ALL_LICKS(find(ALL_LICKS(:,2)==0),1);
LICK.times = LICK_TIMES;
LICK.codes = ones(size(LICK_TIMES));


% % % %%%% JUST WHEN HE STARTES TO LICK TIMES --------------------------------------------------
% % % 
% % % LICK_pos = temp_LICK(find(temp_LICK(:,2)==1),1);
% % % 
% % % LICK_pos=LICK_pos(2:length(LICK_pos));
% % % LICK_pos=LICK_pos(1:2:length(LICK_pos))+1/1000;
% % % 
% % % LICK.times = LICK_pos;
% % % LICK.codes = ones(size(LICK_pos));





















disp('*************************************************')
disp('*************************************************')
disp('ALL FILES  WERE LOADED SUCCESSFULLY')




len=length(FileName1)-4;
% filenm = strcat('Data_',FileName1(1:len))
MERGE_file = strcat('C:\NAVEEN_Work\Cerebellum\Data\MERGED_CELLS\',FileName1(1:len),'\Data_',FileName1(1:len))
NOME = FileName1(1:len);

% getting the task type here -----------------
clear C;
[C,matches] = strsplit(NOME,{'_'});

if cell2mat(strfind(C(2),'L'))
    TASK_TYPE = 'L';
    disp('!!! TASK TYPE has been identified as L !!!');
else
    disp('!!! WARNING: Unidentified task type !!!');
    disp('!!! Please enter a task type to continue !!!');
    TASK_TYPE = upper(input('Here: ','s'));
end

clear C matches;

% ALL FILES HAVE BEEN ENTERED
% NOW CREATING A PROPER AND A SIMPLE STRUCTURE FOR FURTHER ANALYSES


cd(codes_dir);

[Licks,Infos] = Prelim_ADAPTEDLL_n(LICK,Trials,Trigg);



disp('!!! SAVING DATA !!!');
cd(MERGE_dir);

mkdir(FileName1(1:len));
Results_dir = strcat('C:\NAVEEN_Work\Cerebellum\Data\MERGED_CELLS\',FileName1(1:len),'\','RESULTS_LICKS_',NOME);
cd(strcat('C:\NAVEEN_Work\Cerebellum\Data\MERGED_CELLS\',FileName1(1:len)))
mkdir(strcat('RESULTS_LICKS_',NOME))


ALLCELLS_file = strcat(ALLCELLS_dir,'\DATA_',FileName1(1:len));



POP_file = strcat(POP_dir,'\POP_',FileName1(1:len));
save(POP_file,'TASK_TYPE','NOME')



save(MERGE_file,'Infos','Licks','TASK_TYPE','NOME')
save(ALLCELLS_file,'Infos','Licks','TASK_TYPE','NOME')





%%

Monkey_name= upper(NOME(1:2));
if strcmp(Monkey_name,'BR') | strcmp(Monkey_name,'SL')
    save(MERGE_file,'Monkey_name','MERGE_file','Results_dir','POP_dir','POP_file','-append')
    save(ALLCELLS_file,'Monkey_name','MERGE_file','Results_dir','POP_dir','POP_file','-append')  
elseif strcmp(Monkey_name,'PR')
    save(MERGE_file,'Monkey_name','MERGE_file','Results_dir','POP_dir','POP_file','-append')
    save(ALLCELLS_file,'Monkey_name','MERGE_file','ALLCELLS_file','Results_dir','POP_dir','POP_file','-append')
end


disp('ALL FILES WERE MERGED SUCCESSFULLY')
disp('!!! CELL READY FOR FURTHER ANALYSES !!!')


% end