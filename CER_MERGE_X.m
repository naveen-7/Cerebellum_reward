%% To merge the rex and spike 2 files for reward expectation task

% Written by naveen at cumc on 8/11/17

% CER_MERGE_X

close all;
clear;
clc;

% Setup directories--------------------------------------------------------
codes_dir = fullfile('E:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW','LINEAR');
data_dir  = fullfile('E:','NAVEEN_Work','Cerebellum','Data','RECORDED_CELLS');

disp('!!!!!  CER_MERGE has started running  !!!!!')
cd('E:\NAVEEN_Work\Cerebellum\Data')
mkdir('MERGED_CELLS');
MERGE_dir = fullfile('E:','NAVEEN_Work','Cerebellum','Data','MERGED_CELLS');
ALLCELLS_dir = fullfile('E:','NAVEEN_Work','Cerebellum','Data','ALL_CELLS');
POP_dir = fullfile('E:','NAVEEN_Work','Cerebellum','Data','POP_CELLS');




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



% LOADING Spike2 - SSpikes FILE------------------------------------------------------- %% ch2
disp('*************************************************')
disp('*************************************************')
disp('LOAD Spk2-S-Spikes FILE')
[FileName3,PathName3] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
disp(strcat('!!!!!','File you entered is :',FileName3,' !!!!!'));

clear A names;
A = load(cat(2,PathName3,FileName3));
names = fieldnames(A);
S_SPK = A.(names{1,1});
cd(PathName1);




disp('*************************************************')
disp('*************************************************')
disp('ALL FILES  WERE LOADED SUCCESSFULLY')




len=length(FileName1)-4;
% filenm = strcat('Data_',FileName1(1:len))
MERGE_file = strcat('E:\NAVEEN_Work\Cerebellum\Data\MERGED_CELLS\',FileName1(1:len),'\Data_',FileName1(1:len))
NOME = FileName1(1:len);

% getting the task type here -----------------
clear C;
[C,matches] = strsplit(NOME,{'_'});

if cell2mat(strfind(C(2),'M'))
    TASK_TYPE = 'M';
    disp('!!! TASK TYPE has been identified as M !!!');
elseif cell2mat(strfind(C(2),'X'))
    TASK_TYPE = 'X';
    disp('!!! TASK TYPE has been identified as X !!!');
else
    disp('!!! WARNING: Unidentified task type !!!');
    disp('!!! Please enter a task type to continue !!!');
    TASK_TYPE = upper(input('Here: ','s'));
end

clear C matches;

% ALL FILES HAVE BEEN ENTERED
% NOW CREATING A PROPER AND A SIMPLE STRUCTURE FOR FURTHER ANALYSES


cd(codes_dir);

[Spikes.S,Shapes.S,Infos] = Prelim_ADAPTEDX_n(S_SPK,Trials,Trigg,'Simple');



disp('!!! SAVING DATA !!!');
cd(MERGE_dir);

mkdir(FileName1(1:len));
Results_dir = strcat('E:\NAVEEN_Work\Cerebellum\Data\MERGED_CELLS\',FileName1(1:len),'\','RESULTS_CS_',NOME);
cd(strcat('E:\NAVEEN_Work\Cerebellum\Data\MERGED_CELLS\',FileName1(1:len)))
mkdir(strcat('RESULTS_CS_',NOME))


ALLCELLS_file = strcat(ALLCELLS_dir,'\DATA_',FileName1(1:len));



POP_file = strcat(POP_dir,'\POP_',FileName1(1:len));
save(POP_file,'TASK_TYPE','NOME')



save(MERGE_file,'Infos','Spikes','Shapes','TASK_TYPE','NOME')
save(ALLCELLS_file,'Infos','Spikes','Shapes','TASK_TYPE','NOME')





%%

Monkey_name= upper(NOME(1:2));
if strcmp(Monkey_name,'BR') | strcmp(Monkey_name,'SL')
    save(MERGE_file,'Monkey_name','MERGE_file','Results_dir','POP_dir','POP_file','-append')
    save(ALLCELLS_file,'Monkey_name','MERGE_file','Results_dir','POP_dir','POP_file','-append')  
elseif strcmp(Monkey_name,'PR')
    save(MERGE_file,'Monkey_name','MERGE_file','Results_dir','POP_dir','POP_file','-append')
    save(ALLCELLS_file,'Monkey_name','MERGE_file','ALLCELLS_file','Results_dir','POP_dir','POP_file','-append')
end









%% CLEAN UP CS ------------------------------------------


disp('ALL FILES WERE MERGED SUCCESSFULLY')
disp('!!! CELL READY FOR FURTHER ANALYSES !!!')


isrun = upper(input('Do you want to run PRELIM_ALL_n?  ','s'));
if strcmp(isrun,'Y')
    PRELIM_ALLX_n
end

disp('!!! Code terminated successfully !!!');
disp('!!! All data saved and ready for further processing !!!');

% end