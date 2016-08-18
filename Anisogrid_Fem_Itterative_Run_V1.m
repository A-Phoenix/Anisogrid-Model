clear all
close all
cd('C:\Users\aphoenix\Dropbox\Anisogrid_Modeling\Anisogrid_Model_Build')
tic
format Long
Comment=''
Base_Input_File_ID='Thermal_Load_Input_Normaized_Input.dat';
Base_Input_File_Loc='C:\Users\aphoenix\Dropbox\Anisogrid_Modeling\Anisogrid_Model_Build\Thermal_Load\';
 date2=date;
% date2='29-Jul-2016'
%% Top 5 Param Sweep
% Max_Param_Input=[34,.65,3*pi+pi*90/360,.05,.05,.05,.05,123456,123456]
% Min_Param_Input=[6,.0625,pi*90/360,.001,.0025,.001,.0025,123,123]
% Nominal_Param_Input=[20,.25,pi,.01,.01,.01,.01,123,123]
Nominal_Param_Input=[20,.25,pi,.005,.005,.005,.005,123,123]
N=8
index=npermutek(linspace(1,N,N),5);
Hel_Sweep=linspace(6,34,8);
R_Sweep=[.0625,.25,0.383333333333333,0.516666666666667,.65]
Pitch_Sweep=linspace(pi*90/360,3*pi+pi*90/360,5);
Hel_T_Sweep=[.0025,.005,.02,.035,.05]
%     .linspace(.005,.05,4);
Cyl_T_Sweep=[.0025,.005,.02,.035,.05];
% Param_Input_Set=zeros(length(index),1)*Nominal_Param_Input
counter=0;

for i=1:length(index)
    if isempty(find(index(i,2:end)>=6))==1
%         index(i,2)~=8 | index(i,3)~=8 | index(i,2)~=7 | index(i,3)~=7 index(i,2)~=6 | index(i,3)~=6 | index(i,4)~=8 | index(i,5)~=8 | index(i,4)~=7 | index(i,5)~=7 |index(i,4)~=6 | index(i,5)~=6
        Param_Input_Set(i-counter,:)=Nominal_Param_Input;
        Param_Input_Set(i-counter,1:3)=[Hel_Sweep(index(i,1)),R_Sweep(index(i,2)),Pitch_Sweep(index(i,3))];
        Param_Input_Set(i-counter,5)=Hel_T_Sweep(index(i,4));
        Param_Input_Set(i-counter,7)=Hel_T_Sweep(index(i,5));
       else

        counter=counter+1;
    end
end
%%
% Param_Input_Set(3,:)=Max_Param_Input
% Param_Input_Set(2,:)=Min_Param_Input
% Param_Input_Set(1,:)=Nominal_Param_Input
   n=ones(length(Param_Input_Set(:,1)),1)*length(Param_Input_Set(:,1));
% % Anisogrid_Fem_Build_Function_V3(Base_Input_File_ID,Base_Input_File_Loc,Comment,Total_Number_of_Helicals,R,pitch,Hel_IR,Hel_OR,Cir_IR,Cir_OR,bot_Const,top_Const,i)
for nn=1:length(Param_Input_Set(1,1,:))
    for i=1:length(Param_Input_Set(1:n(nn),1))
        save_File_ID=strcat(Base_Input_File_ID(1:end-4),'_','_',Comment,'_',num2str(i),'_',num2str(nn),'.dat');
        i
        cd('C:\Users\aphoenix\Dropbox\Anisogrid_Modeling\Anisogrid_Model_Build')
        [cir_number(i)]=Anisogrid_Fem_Build_Function_V13(Base_Input_File_ID,Base_Input_File_Loc,Comment,Param_Input_Set(i,1,nn),Param_Input_Set(i,2,nn),Param_Input_Set(i,3,nn),Param_Input_Set(i,4,nn),Param_Input_Set(i,5,nn),Param_Input_Set(i,6,nn),Param_Input_Set(i,7,nn),Param_Input_Set(i,8,nn),Param_Input_Set(i,9,nn),i,date2);
        fclose('all') 
    end
    disp('All Files Created')
    toc
    for i=3101:length(Param_Input_Set(1:n(nn),1))
        Save_File_Location=strcat(Base_Input_File_Loc,date2);
        save_File_ID=strcat(Base_Input_File_ID(1:end-4),'_','_',Comment,'_',num2str(i),'.dat');
        PCH_File_ID=strcat(Base_Input_File_ID(1:end-4),'_','_',Comment,'_',num2str(i),'.pch');
        PCH_File_Loc_Name=strcat(Save_File_Location,'\',PCH_File_ID);
        save_File_ID=strcat(Base_Input_File_ID(1:end-4),'_','_',Comment,'_',num2str(i),'.dat');
        cd(strcat('C:\Users\aphoenix\Dropbox\Anisogrid_Modeling\Anisogrid_Model_Build\Thermal_Load\',date2))
        system(['C:\MSC.Software\MSC_Nastran\20160\bin\nastran.exe   ',strcat(char(Save_File_Location),'\',char(save_File_ID))]);
    %% Read Punch File  
    % filename = 'C:\Users\aphoenix\Dropbox\Anisogrid_Modeling\Anisogrid_Model_Build\Thermal_Load\24-Apr-2016\thermal_load_input_2___1.pch';
    startRow = 8;
    formatSpec = '%*18s%18f%18f%18f%[^\n\r]';
    fileID = fopen(char(PCH_File_Loc_Name),'r');
    PchDataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    PchInputData = [PchDataArray{1:end-1}];
    L=length(PchInputData(1:2:end,:));
    Disp_Responce(1:L,1:3,i,nn)=PchInputData(1:2:end,:);
    Disp_Responce(1:L,4:6,i,nn)=PchInputData(2:2:end,:);
    toc
    i
    nn
    fclose('all') 
    end
     
end
% save(strcat('Aniso_Disp_Resp_',Base_Input_File_ID(1:end-4),'_',Comment,date2),'Disp_Responce','Base_Input_File_ID','Base_Input_File_Loc','Param_Input_Set','-v7.3')
%    save(strcat('Aniso_Disp_Resp_',Base_Input_File_ID(1:end-4),'_',Comment,date2,'4745_end'),'Disp_Responce','Base_Input_File_ID','Base_Input_File_Loc','Param_Input_Set','-v7.3')
  save(strcat('Aniso_Disp_Resp_',Base_Input_File_ID(1:end-4),'_3102_end',Comment,date2),'Disp_Responce','Base_Input_File_ID','Base_Input_File_Loc','Param_Input_Set','-v7.3')
% save(strcat('Aniso_Disp_Resp_',Base_Input_File_ID(1:end-4),'_1_3101_',Comment,date2),'Disp_Responce','Base_Input_File_ID','Base_Input_File_Loc','Param_Input_Set','-v7.3')

% save(strcat('Aniso_Disp_Resp_',Base_Input_File_ID(1:end-4),'_',Comment,date2,'673_2211'),'Disp_Responce','Base_Input_File_ID','Base_Input_File_Loc','Param_Input_Set','-v7.3')
