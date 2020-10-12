%% Batch Code: RUN THIS FIRST


% Created by NAVEEN ON 05/29/15 at CUMC
% Modified on 04/20/16 at CUMC
% Linearized on 03/08/17 at CUMC


% function CER_analysis_LINEAR_BATCH


 
%% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % FOR REGULAR, RELEASE HERE -----------------------------------------

%clc;
%clear ;
%close all;


% codes_dir = fullfile('E:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW','LINEAR');
% data_dir  = fullfile('E:','NAVEEN_Work','Cerebellum','Data','ALL_CELLS');

% disp('!!!!!  CER_analysis_LINEAR_BATCH has started running  !!!!!')
% cd(data_dir)

% % cd('E:\NAVEEN_Work\Cerebellum\Data\MERGED_CELLS');

% % % LOADING FILE--------------------------------------------------------
% disp('*************************************************')
% disp('*************************************************')
% disp('LOAD THE CELL DATA FILE')
% [FileName1,PathName1] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
% DATAfile = cat(2,PathName1,FileName1);
% disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
% load(DATAfile);
% cd(PathName1);


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % FOR BATCH, RELEASE HERE -----------------------------------------
% % 
close all;
File = DATAfile;
FileName = Tempy;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%% TWO: Running all the codes in the batch ---------------------------

Death = 0;

try
    
% % %     Nminusn_Nplusn
    
    
 CER_analysis_03_LINEAR
% CER_analysis_TWO
%  CER_analysis_ZERO

%     CER_ISI

%     CER_dr   
%     CER_dr_W   
%     CER_dr_R
%     CER_dr_motor
 
catch
    Death = 1;
    disp('!!! Error Encountered !!!');
    disp ('!!! OPERATION TERMINATED PREMATURELY !!!');
    disp('!!! Death flag has been hoisted !!!');
end




disp(strcat( '!!!!!!  Death has been noticed as: ',num2str(Death),'  !!!!!!'));
if Death==0
    disp('!!! ALL CODES EXECUTED SUCCESSFULLY !!!');
end


















% end
