%% To merge the rex and spike 2 files for Lick 

% Written by naveen at cumc on 8/23/17

% CER_MERGE_L

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



% LOADING Spike2 - EMPTY FILE-------------------------------------------------------- %% ch7
disp('*************************************************')
disp('*************************************************')
disp('LOAD Spk2-EMPTY FILE')
[FileName4,PathName4] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
disp(strcat('!!!!!','File you entered is :',FileName4,' !!!!!'));

clear A names;
A = load(cat(2,PathName4,FileName4));
names = fieldnames(A);
EMPTY = A.(names{1,1});
cd(PathName1);






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



LICKS = Prelim_ADAPTEDL_n(LICK_DATA_IMPORT,Trials,Trigg);
[Spikes,Shapes,Infos] = Prelim_ADAPTEDX_n(EMPTY,Trials,Trigg,'Empty');


disp('!!! SAVING DATA !!!');
cd(MERGE_dir);

mkdir(FileName1(1:len));
Results_dir = strcat('C:\NAVEEN_Work\Cerebellum\Data\MERGED_CELLS\',FileName1(1:len),'\','RESULTS_CS_',NOME);
cd(strcat('C:\NAVEEN_Work\Cerebellum\Data\MERGED_CELLS\',FileName1(1:len)))
mkdir(strcat('RESULTS_CS_',NOME))


ALLCELLS_file = strcat(ALLCELLS_dir,'\DATA_',FileName1(1:len));



POP_file = strcat(POP_dir,'\POP_',FileName1(1:len));
save(POP_file,'TASK_TYPE','NOME')



save(MERGE_file,'Infos','LICKS','TASK_TYPE','NOME')
save(ALLCELLS_file,'Infos','LICKS','TASK_TYPE','NOME')





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