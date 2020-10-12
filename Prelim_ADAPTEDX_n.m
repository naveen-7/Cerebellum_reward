%% Makes the cell plate and gives all necessary information about the code
% Created by NAVEEN ON 8/11/17 at CUMC



function [SpikeS,ShapeS,Infos] = Prelim_ADAPTEDX_n(CHAN,Trials,Trigg,String)

%
% CER_analysis_01_ADAPTEDX_n
CER_analysis_01_ADAPTED_MISMATCHX_n
% 
clearvars -except MERGE_dir Spikes Shapes C_SPK S_SPK CHtemp UNIQUE_CODES max_sspikes data_CH shapes_CH CODE_S C filenm CHAN Trials Trigg CODE PathName1 FileName1 codes_dir data_dir Results_dir CH_number FILE_name SpikeS Infos


FLAG_CODE=zeros(max_sspikes,1);


% for CODE_S = 1:max_sspikes
    disp(strcat('>>> Analysing CODE-',num2str(CODE_S),' >>>'));
    
    
    if CODE_S==0
        SpikeS = data_CH(:,1);
        ShapeS = shapes_CH(:,1);
    else
%         clear Tempy TEMPY;
        for kk=1:size(data_CH,1)
                         kk;
            
            if (size(data_CH{kk,3},1)==max_sspikes) & (size(shapes_CH{kk,2},1)==max_sspikes)
                Tempy{kk,1}= data_CH{kk,3}{CODE_S,1};
                TEMPY{kk,:}= shapes_CH{kk,2}{CODE_S,1};
                FLAG_CODE(CODE_S)=1;
            else
                Tempy{kk,1}= [];
                TEMPY{kk,:}= [];
                FLAG_CODE(CODE_S)=0;
            end
        end
        
        SpikeS = Tempy;
        ShapeS = TEMPY;
    end
    
    
% end






disp('!!! END OF CODE Prelim_ADAPTEDX_n !!!');

end