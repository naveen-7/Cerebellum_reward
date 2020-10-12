%% Batch Code: RUN THIS FIRST
%% PLUGIN code ... takes all the data automatically and runs them


% Created by NAVEEN ON 08/04/17 at CUMC



function CER_analysis_LINEAR_BATCH_BATCH

clc;
clear all;
close all;


% Setup directories--------------------------------------------------------

codes_dir = fullfile('E:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW','LINEAR');
data_dir  = fullfile('E:','NAVEEN_Work','Cerebellum','Data','ALL_CELLS');


cd(data_dir)

III=0;
LENGTH = length(dir);
ADAPTER_count=0;
PRINT_FLAG=0;
PRINT_NO=0;

for IJ=1:LENGTH
    
    if IJ==LENGTH
        PRINT_FLAG=1;
    end
    
    cd(data_dir)
    
    tempy = dir;
    TEMPY = getfield(tempy,{IJ});
    Tempy = TEMPY.name;
    
    if length(Tempy)>3
        Tempy1 = upper(Tempy);
        if strcmp(Tempy1(1:4),'DATA')
            
            III=III+1;
            clearvars -except ADAPTER_count codes_dir data_dir I III LENGTH matfile_NEW ...
                PRINT_FLAG PRINT_NO tempy Tempy TEMPY
            
            disp(strcat('!!!!!!!!!!!!!! RUNNING CELL #',num2str(III), '!!!!!!!!!!!!!!'));
            DATAfile = strcat(data_dir,'\',Tempy);
            load(DATAfile,'TASK_TYPE')
            if (strcmp(TASK_TYPE,'N') | strcmp(TASK_TYPE,'R') | strcmp(TASK_TYPE,'D') )
                
                %                  if strcmp(TASK_TYPE,'D')
                %             if strcmp(TASK_TYPE,'W')
                disp(strcat('!!!! LOADING',Tempy,'!!!!'));
                load(DATAfile)
                cd(codes_dir)
                
                try
                    %                 if strcmp(LEARNING,'Y')
                    CER_analysis_LINEAR_BATCH
                    %                 end
                end
                
                
            end
        end
    end
end

end




