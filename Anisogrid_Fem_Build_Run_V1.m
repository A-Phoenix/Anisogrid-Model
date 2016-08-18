function[C_F]=Anisogrid_Fem_Build_Run_V1(Param_Input_Set)
clear all
 Param_Input_Set=[34,.625,1.030E+00,.005,.005,.005,.005,123,123]
 cd('C:\Users\aphoenix\Dropbox\Anisogrid_Modeling\Anisogrid_Model_Build')
format Long
Comment='Insert_comment_here';
Base_Input_File_ID='Input_Base_File_with_Mat_and_LoadingCases';
Base_Input_File_Loc='Base_File_Location';
  date2=date;
    n=ones(length(Param_Input_Set(:,1)),1)*length(Param_Input_Set(:,1));
for nn=1:length(Param_Input_Set(1,1,:))
    for i=1:length(Param_Input_Set(1:n(nn),1))
        save_File_ID=strcat(Base_Input_File_ID(1:end-4),'_','_',Comment,'_',num2str(i),'_',num2str(nn),'.dat');
         cd('Save_Location_Input')
        [cir_number(i)]=Anisogrid_Fem_Build_Function_V13(Base_Input_File_ID,Base_Input_File_Loc,Comment,Param_Input_Set(i,1,nn),Param_Input_Set(i,2,nn),Param_Input_Set(i,3,nn),Param_Input_Set(i,4,nn),Param_Input_Set(i,5,nn),Param_Input_Set(i,6,nn),Param_Input_Set(i,7,nn),Param_Input_Set(i,8,nn),Param_Input_Set(i,9,nn),i,date2);
        fclose('all') 
    end
    for i=1:length(Param_Input_Set(1:n(nn),1))
        Save_File_Location=strcat(Base_Input_File_Loc,date2);
        save_File_ID=strcat(Base_Input_File_ID(1:end-4),'_','_',Comment,'_',num2str(i),'.dat');
        PCH_File_ID=strcat(Base_Input_File_ID(1:end-4),'_','_',Comment,'_',num2str(i),'.pch');
        PCH_File_Loc_Name=strcat(Save_File_Location,'\',PCH_File_ID);
        save_File_ID=strcat(Base_Input_File_ID(1:end-4),'_','_',Comment,'_',num2str(i),'.dat');
        cd(strcat('C:\Users\aphoenix\Dropbox\Anisogrid_Modeling\Anisogrid_Model_Build\Thermal_Load\',date2))
        system(['nastran.exe Location',strcat(char(Save_File_Location),'\',char(save_File_ID))]);
    %% Read Punch File  
    startRow = 8;
    formatSpec = '%*18s%18f%18f%18f%[^\n\r]';
    fileID = fopen(char(PCH_File_Loc_Name),'r');
    PchDataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    PchInputData = [PchDataArray{1:end-1}];
    L=length(PchInputData(1:2:end,:));
    Disp_Responce(1:L,1:3,i,nn)=PchInputData(1:2:end,:);
    Disp_Responce(1:L,4:6,i,nn)=PchInputData(2:2:end,:);
    fclose('all') 
    end    
end

%      end
% save('C:\Users\aphoenix\Dropbox\Anisogrid_Modeling\Anisogrid_Model_Build\Optimzed3Param_Geom,'Cost_Function','Disp_Responce','Param_Input_Set')
%  save('C:\Users\aphoenix\Dropbox\Anisogrid_Modeling\Anisogrid_Model_Build\Max_Disp_Geom','Cost_Function','Disp_Responce','Param_Input_Set')