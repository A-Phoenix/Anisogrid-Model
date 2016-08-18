 function [Cir_Number]=Anisogrid_Fem_Build_Function_V13(Base_Input_File_ID,Base_Input_File_Loc,Comment,Total_Number_of_Helicals,R,pitch,Hel_IR,Hel_Thickness,Cir_IR,Cir_Thickness,bot_Const,top_Const,i,date2)
%% debug with no funciton
%  clear all 
% 
%  Comment='Top5_Param_Sweep'
%  Base_Input_File_ID='Thermal_Load_Input_Normaized_Input.dat';
%  Base_Input_File_Loc='C:\Users\aphoenix\Dropbox\Anisogrid_Modeling\Anisogrid_Model_Build\Thermal_Load\';
%  date2=date;
% nput=[0.4001626698637          6.28319307409346                      0.05        0.0198800959110092]
% %%
% i=1
% Nominal_Param_Input(1)=34;
% R=Input(1);
% Pitch=Input(2);
% HelT=Input(3);
% CylT=Input(4);
% Nominal_Param_Input(2:9)=[R,Pitch,.01,HelT,.01,CylT,123,123];
% Param_Input_Set=Nominal_Param_Input
% % % Param_Input_Set(1,:)=Nominal_Param_Input
% Total_Number_of_Helicals=Param_Input_Set(i,1,1);
% R=Param_Input_Set(i,2,1);
% pitch=Param_Input_Set(i,3,1)
% Hel_IR=Param_Input_Set(i,4,1);
% Hel_Thickness=Param_Input_Set(i,5,1);
% Cir_IR=Param_Input_Set(i,6,1);
% Cir_Thickness=Param_Input_Set(i,7,1);
% bot_Const=Param_Input_Set(i,8,1);
% top_Const=Param_Input_Set(i,9,1);
%%
warning('off', 'MATLAB:MKDIR:DirectoryExists')
 
 
Hel_OR=Hel_IR+Hel_Thickness;
Cir_OR=Cir_IR+Cir_Thickness;
 %% non fuction inputs
%%
save_File_ID=strcat(Base_Input_File_ID(1:end-4),'_','_',Comment,'_',num2str(i),'.dat');
Base_Input_File=strcat(Base_Input_File_Loc,Base_Input_File_ID);
Save_File_Location=strcat(Base_Input_File_Loc,date2);%{['C:\Users\aphoenix\Dropbox\Anisogrid_Modeling\Anisogrid_Model_Build\Thermal_Load\',date]};
mkdir(char(Save_File_Location))
copyfile(strcat([char(Base_Input_File_Loc),char(Base_Input_File_ID)]),strcat([char(Save_File_Location),'\',char(save_File_ID)]));
cd(Save_File_Location)
format long g
%% Hard_Coded_Vars
Height=1;
Number_Of_Hel_Elements_Between_Constraints=10;
Number_Of_Cir_Elements_Between_Constraints=10;
%%
X_append=0;
Y_append=R;
Z_append=0;

Removed_Node_Counter=0;
%         pitch(k)=2*pi*2/Total_Number_of_Helicals+2*pi*pitch_input(k);
        pitch_input=(pitch-2*pi*2/Total_Number_of_Helicals)/(2*pi);
        Intersections=1+Total_Number_of_Helicals*(2/Total_Number_of_Helicals+pitch_input);
%         pitch_input=(pitch/(2*pi)-2*pi*2/Total_Number_of_Helicals)/(2*pi);
%         Intersections=1+Total_Number_of_Helicals*(2/Total_Number_of_Helicals+pitch_input);
        Number_Of_Hel_Nodes=round(2*Number_Of_Hel_Elements_Between_Constraints*Intersections-(2*Number_Of_Hel_Elements_Between_Constraints-1));

%         z(1:Number_Of_Hel_Nodes)=linspace(0,Height,Number_Of_Hel_Nodes)';

Nodes_Remainder=rem(rem(pitch/pi,1)*(10*Total_Number_of_Helicals),1);
% Nodes_Overhang=rem(rem(pitch/pi,1)*(10*Total_Number_of_Helicals),10);
Nodes_Overhang=10*rem(2*Intersections,1);

if Nodes_Remainder==0
    z=linspace(0,Height,Number_Of_Hel_Nodes)';
elseif Nodes_Overhang<1
    Intersection_Point=Height-rem(2*Intersections,1)*Height/(2*(Intersections-1));
    El_Tip_Height=(Height-Intersection_Point)/ceil(Nodes_Overhang);
    if round(El_Tip_Height*Nodes_Overhang,4)==0
        %       z=linspace(0,Height,Number_Of_Hel_Nodes)';
        Intersection_Point=Height-rem(2*Intersections,1)*Height/(2*(Intersections-1));
        El_Tip_Height=(Height-Intersection_Point)/ceil(Nodes_Overhang);
        z=linspace(0,Intersection_Point,Number_Of_Hel_Nodes)';
        z(end)=Height;
    elseif Nodes_Remainder>.5
        Intersection_Point=Height-rem(2*Intersections,1)*Height/(2*(Intersections-1));
        El_Tip_Height=(Height-Intersection_Point)/ceil(Nodes_Overhang);
        %      H_per=(Number_Of_Hel_Nodes-Nodes_Overhang)/Number_Of_Hel_Nodes;
        %     el_end_Step=Height*(1-H_per)/ceil(Nodes_Overhang);
        z=linspace(0,Intersection_Point,Number_Of_Hel_Nodes-ceil(Nodes_Overhang))';
        z(length(z)+1:length(z)+ceil(Nodes_Overhang))=linspace(Intersection_Point+El_Tip_Height,Height,ceil(Nodes_Overhang))';
%         Number_Of_Hel_Nodes=Number_Of_Hel_Nodes+1;
        %  elseif Nodes_Overhang>1
    else
        Intersection_Point=Height-rem(2*Intersections,1)*Height/(2*(Intersections-1));
        El_Tip_Height=(Height-Intersection_Point)/ceil(Nodes_Overhang);
        %      H_per=(Number_Of_Hel_Nodes-Nodes_Overhang)/Number_Of_Hel_Nodes;
        %     el_end_Step=Height*(1-H_per)/ceil(Nodes_Overhang);
        z=linspace(0,Intersection_Point,Number_Of_Hel_Nodes-ceil(Nodes_Overhang)+1)';
        z(length(z)+1:length(z)+ceil(Nodes_Overhang))=linspace(Intersection_Point+El_Tip_Height,Height,ceil(Nodes_Overhang))';
        Number_Of_Hel_Nodes=Number_Of_Hel_Nodes+1;
    end
elseif Nodes_Remainder>.5
    Intersection_Point=Height-rem(2*Intersections,1)*Height/(2*(Intersections-1));
    El_Tip_Height=(Height-Intersection_Point)/ceil(Nodes_Overhang);
    %      H_per=(Number_Of_Hel_Nodes-Nodes_Overhang)/Number_Of_Hel_Nodes;
    %     el_end_Step=Height*(1-H_per)/ceil(Nodes_Overhang);
    z=linspace(0,Intersection_Point,Number_Of_Hel_Nodes-ceil(Nodes_Overhang))';
    z(length(z)+1:length(z)+ceil(Nodes_Overhang)+1)=linspace(Intersection_Point+El_Tip_Height,Height,ceil(Nodes_Overhang)+1)';
    Number_Of_Hel_Nodes=Number_Of_Hel_Nodes+1;
elseif Nodes_Remainder<.5
    
    Intersection_Point=Height-rem(2*Intersections,1)*Height/(2*(Intersections-1));
    El_Tip_Height=(Height-Intersection_Point)/ceil(Nodes_Overhang);
    %      H_per=(Number_Of_Hel_Nodes-Nodes_Overhang)/Number_Of_Hel_Nodes;
    %     el_end_Step=Height*(1-H_per)/ceil(Nodes_Overhang);
    z=linspace(0,Intersection_Point,Number_Of_Hel_Nodes-ceil(Nodes_Overhang)+1)';
    z(length(z)+1:length(z)+ceil(Nodes_Overhang))=linspace(Intersection_Point+El_Tip_Height,Height,ceil(Nodes_Overhang))';
    Number_Of_Hel_Nodes=Number_Of_Hel_Nodes+1;
end
      
        Hel_length=Number_Of_Hel_Nodes;
        for j=1:Total_Number_of_Helicals
            if mod(j,2)==1 % even / odd check if odd=1    
                theta(1:Number_Of_Hel_Nodes,j)=pitch*z(1:Number_Of_Hel_Nodes)+ones(size(z(1:Number_Of_Hel_Nodes)))*(j-1)*(2*pi/(Total_Number_of_Helicals));
            else
                theta(1:Number_Of_Hel_Nodes,j)=-pitch*z(1:Number_Of_Hel_Nodes)+ones(size(z(1:Number_Of_Hel_Nodes)))*(j-2)*(2*pi/(Total_Number_of_Helicals));
            end
            X_append(1+(j-1)*Number_Of_Hel_Nodes:j*Number_Of_Hel_Nodes)=[R*cos(theta(1:Number_Of_Hel_Nodes,j))];
            Y_append(1+(j-1)*Number_Of_Hel_Nodes:j*Number_Of_Hel_Nodes)=[R*sin(theta(1:Number_Of_Hel_Nodes,j))];
            Z_append(1+(j-1)*Number_Of_Hel_Nodes:j*Number_Of_Hel_Nodes)=[z(1:Number_Of_Hel_Nodes)];
        end
        Hel_Nodes=[X_append(:),Y_append(:),Z_append(:)];

        Total_Number_of_Circumfrential=Intersections-1;
        n=0;% counter to append cir nodes
        for m=1:round(Total_Number_of_Circumfrential)
            %theta_Cir(:,m)=linspace(0,2*pi-2*pi/Total_Number_of_Helicals,Total_Number_of_Helicals);
            %theta_Cir(:,m)=linspace(pi/Total_Number_of_Helicals,2*pi-pi/Total_Number_of_Helicals,Total_Number_of_Helicals*Number_Of_Cir_Nodes_Between_Constraints);
            Cir_Start=pi/Total_Number_of_Helicals;
            Cir_End_2=2*pi+pi/Total_Number_of_Helicals-2*pi/Total_Number_of_Helicals/Number_Of_Cir_Elements_Between_Constraints;
            %Cir_Number=Total_Number_of_Helicals*Number_Of_Cir_Nodes_Between_Constraints+2-Number_Of_Cir_Nodes_Between_Constraints;
            Cir_Number=Total_Number_of_Helicals*Number_Of_Cir_Elements_Between_Constraints;
            theta_Cir(:,m)=linspace(Cir_Start,Cir_End_2,Cir_Number);
            x_Cir(:,m)=R*cos(theta_Cir(:,m));
            y_Cir(:,m)=R*sin(theta_Cir(:,m));
            z_Cir(:,m)=(Height/Total_Number_of_Circumfrential/2+(m-1)*Height/Total_Number_of_Circumfrential)*ones(Cir_Number,1);
            Cir_Nodes(1+n:n+Cir_Number,:)=[x_Cir(:,m),y_Cir(:,m),z_Cir(:,m)];
            n=n+Cir_Number;
        end
%         xtest(:,:)=R*cos(theta);
%         ytest(:,:)=R*sin(theta);
        %% PLot function
%         figure
%         %for l=1:length(xtest(1,:,1))  
%         title({'Pitch Input ',num2str(pitch/pi)})
%         hold on
%         plot3(X_append(:),Y_append(:),Z_append(:),'--rs','LineWidth',2)
%         plot3(xtest(:,:),ytest(:,:),z(:),'-ok')
% %     %   view(3)
%         view([1,0,0])
%         hold on
%         plot3(x_Cir(:,:),y_Cir(:,:),z_Cir(:,:),'-*r')
%         grid
%%
%     if k~=1%isempty(find(abs(Hel_Nodes(:,1,k))<=1e-12))==0
%         Hel_End=find(Hel_Nodes(:,1)==0,1)-1;
%         Cir_End=find(Cir_Nodes(:,1)==0,1)-1;
%     else
        [a,b]=size(Hel_Nodes(:,1,1));
        [aa,bb]=size(Cir_Nodes(:,1,1));
        Hel_End=a;
        Cir_End=aa;
%     end
        Element_ID=linspace(1,Hel_End+Cir_End,Hel_End+Cir_End);
    Node_ID=linspace(1,Hel_End+Cir_End,Hel_End+Cir_End)';
    Hel_Prop_ID=1;
    Cir_Prop_ID=2;
    Prop_ID=[Hel_Prop_ID*ones(Hel_End,1);Cir_Prop_ID*ones(Cir_End,1)];
    Node_Def_Cord=[0*ones(Hel_End,1);0*ones(Cir_End,1)];
    Node_Out_Cord=[0*ones(Hel_End,1);0*ones(Cir_End,1)];
    %
    Nodal_Export_1(1:Hel_End+Cir_End,:)=[Hel_Nodes(1:Hel_End,:);Cir_Nodes(1:Cir_End,:)];
    [Nodal_Export_update,Nodal_Export_Loc_Vector,Nodal_Export_ID_Vector]=unique(round(Nodal_Export_1(1:Hel_End+Cir_End,:),4),'rows','stable');%,'stable');
    Hex_ElementRemovalNodes=find(Nodal_Export_1(:,3)==Height);
    Cir_Element_Node_Update=Hel_End+linspace(Cir_End/round(Total_Number_of_Circumfrential),(round(Total_Number_of_Circumfrential))*(Cir_End/round(Total_Number_of_Circumfrential)),round(Total_Number_of_Circumfrential));
 Hex_ElementRemovalNodes=Nodal_Export_ID_Vector(Hex_ElementRemovalNodes);
New_Node_Esport=Nodal_Export_1(b,:);
    Removed_Node_Counter=0;
    for j=1:(size(Nodal_Export_ID_Vector))
         if isempty(find(Hex_ElementRemovalNodes==Nodal_Export_ID_Vector(j)))==0
             Removed_Node_Counter=Removed_Node_Counter+1;
         elseif isempty(find(j==Cir_Element_Node_Update))==0
               CBEAM_Build(j-Removed_Node_Counter)={num2str(['CBEAM,',num2str(Node_ID(j)),',',num2str(Prop_ID(j)),',',num2str(Nodal_Export_ID_Vector(j)),',',num2str(Nodal_Export_ID_Vector(j-Cir_End/round(Total_Number_of_Circumfrential)+1)),',','0.,1.,0.',',','2'])};
         else   
             CBEAM_Build(j-Removed_Node_Counter)={num2str(['CBEAM,',num2str(Node_ID(j)),',',num2str(Prop_ID(j)),',',num2str(Nodal_Export_ID_Vector(j)),',',num2str(Nodal_Export_ID_Vector(j+1)),',','0.,1.,0.',',','2'])};     
         end
    end
    for j=1:size(Nodal_Export_Loc_Vector)
        Grid_Build(j)={num2str(['GRID*,',num2str(Node_ID(j)),',',num2str(Node_Def_Cord(Nodal_Export_Loc_Vector(j))),',',num2str(Nodal_Export_1(Nodal_Export_Loc_Vector(j),1),'%8.6f'),',',num2str(Nodal_Export_1(Nodal_Export_Loc_Vector(j),2),'%8.6f'),',',num2str(Nodal_Export_1(Nodal_Export_Loc_Vector(j),3),'%8.6f'),',',num2str(Node_Out_Cord(j))])};
    end

%% Build RBE_2
Top_Nodes=find(round(Nodal_Export_1(Nodal_Export_Loc_Vector,3),5)==Height);
[length_Top_Nodes,~]=size(Top_Nodes);
Bot_Nodes=find(Nodal_Export_1(Nodal_Export_Loc_Vector,3)==0);
[length_Bot_Nodes,~]=size(Bot_Nodes);
RBE_Grid=['GRID,190000,0,0.,0.,1.,0';'GRID,190001,0,0.,0.,0.,0'];
% GRID,90000,0,0,0,1,2   
% GRID,90001,0,0,0,0,2 
% RBE2,90000,90000,123456,41,80,156
RBE_1_2_L=0;

    if length_Top_Nodes==3
        RBE_1(1,:)=strcat('RBE2,190000,190000,',num2str(top_Const),',',num2str(Top_Nodes(1)),',',num2str(Top_Nodes(2)),',',num2str(Top_Nodes(3)),',');
    elseif length_Top_Nodes==4
        RBE_1(1,:)=strcat('RBE2,190000,190000,',num2str(top_Const),',',num2str(Top_Nodes(1)),',',num2str(Top_Nodes(2)),',',num2str(Top_Nodes(3)),',',num2str(Top_Nodes(4)),',');
    elseif length_Top_Nodes>4
        RBE_1(1,:)=strcat('RBE2,190000,190000,',num2str(top_Const),',',num2str(Top_Nodes(1)),',',num2str(Top_Nodes(2)),',',num2str(Top_Nodes(3)),',',num2str(Top_Nodes(4)),',',num2str(Top_Nodes(5)),',');
    end
    if length_Top_Nodes>=13
        for i=1:(length_Top_Nodes-5)/8
                [~,RBE_1_2_L]=size(strcat(',',num2str(Top_Nodes(6+8*(i-1))),',',num2str(Top_Nodes(7+8*(i-1))),',',num2str(Top_Nodes(8+8*(i-1))),',',num2str(Top_Nodes(9+8*(i-1))),',',num2str(Top_Nodes(10+8*(i-1))),',',num2str(Top_Nodes(11+8*(i-1))),',',num2str(Top_Nodes(12+8*(i-1))),',',num2str(Top_Nodes(13+8*(i-1))),','));
                RBE_1(i+1,1:RBE_1_2_L)=strcat(',',num2str(Top_Nodes(6+8*(i-1))),',',num2str(Top_Nodes(7+8*(i-1))),',',num2str(Top_Nodes(8+8*(i-1))),',',num2str(Top_Nodes(9+8*(i-1))),',',num2str(Top_Nodes(10+8*(i-1))),',',num2str(Top_Nodes(11+8*(i-1))),',',num2str(Top_Nodes(12+8*(i-1))),',',num2str(Top_Nodes(13+8*(i-1))),',');
                L=length(RBE_1(1,:));
                RBE_1(i+1,RBE_1_2_L+1:L)=blanks(L-RBE_1_2_L);    
        end
    else
        i=0;
    end
    if rem(length_Top_Nodes-5,8)>0
        if rem(length_Top_Nodes-5,8)==1
            [~,RBE_1_2_L]=size(strcat(',',num2str(Top_Nodes(13+8*(i-1)+1)),','));
            RBE_1(i+2,1:RBE_1_2_L)=strcat(',',num2str(Top_Nodes(13+8*(i-1)+1)),',');
            L=length(RBE_1(1,:));
            RBE_1(i+2,RBE_1_2_L+1:L)=blanks(L-RBE_1_2_L);  
        elseif rem(length_Top_Nodes-5,8)==2
            [~,RBE_1_2_L]=size(strcat(',',num2str(Top_Nodes(13+8*(i-1)+1)),',',num2str(Top_Nodes(13+8*(i-1)+2)),','));
            RBE_1(i+2,1:RBE_1_2_L)=strcat(',',num2str(Top_Nodes(13+8*(i-1)+1)),',',num2str(Top_Nodes(13+8*(i-1)+2)),',');
             L=length(RBE_1(1,:));
                RBE_1(i+2,RBE_1_2_L+1:L)=blanks(L-RBE_1_2_L);  
        elseif rem(length_Top_Nodes-5,8)==3
            [~,RBE_1_2_L]=size(strcat(',',num2str(Top_Nodes(13+8*(i-1)+1)),',',num2str(Top_Nodes(13+8*(i-1)+2)),',',num2str(Top_Nodes(13+8*(i-1)+3)),','));
            RBE_1(i+2,1:RBE_1_2_L)=strcat(',',num2str(Top_Nodes(13+8*(i-1)+1)),',',num2str(Top_Nodes(13+8*(i-1)+2)),',',num2str(Top_Nodes(13+8*(i-1)+3)),',');
             L=length(RBE_1(1,:));
                RBE_1(i+2,RBE_1_2_L+1:L)=blanks(L-RBE_1_2_L);  
        elseif rem(length_Top_Nodes-5,8)==4
            [~,RBE_1_2_L]=size(strcat(',',num2str(Top_Nodes(13+8*(i-1)+1)),',',num2str(Top_Nodes(13+8*(i-1)+2)),',',num2str(Top_Nodes(13+8*(i-1)+3)),',',num2str(Top_Nodes(13+8*(i-1)+4)),','));
            RBE_1(i+2,1:RBE_1_2_L)=strcat(',',num2str(Top_Nodes(13+8*(i-1)+1)),',',num2str(Top_Nodes(13+8*(i-1)+2)),',',num2str(Top_Nodes(13+8*(i-1)+3)),',',num2str(Top_Nodes(13+8*(i-1)+4)),',');
             L=length(RBE_1(1,:));
                RBE_1(i+2,RBE_1_2_L+1:L)=blanks(L-RBE_1_2_L);  
        elseif rem(length_Top_Nodes-5,8)==5
            [~,RBE_1_2_L]=size(strcat(',',num2str(Top_Nodes(13+8*(i-1)+1)),',',num2str(Top_Nodes(13+8*(i-1)+2)),',',num2str(Top_Nodes(13+8*(i-1)+3)),',',num2str(Top_Nodes(13+8*(i-1)+4)),',',num2str(Top_Nodes(13+8*(i-1)+5)),','));
            RBE_1(i+2,1:RBE_1_2_L)=strcat(',',num2str(Top_Nodes(13+8*(i-1)+1)),',',num2str(Top_Nodes(13+8*(i-1)+2)),',',num2str(Top_Nodes(13+8*(i-1)+3)),',',num2str(Top_Nodes(13+8*(i-1)+4)),',',num2str(Top_Nodes(13+8*(i-1)+5)),',');
             L=length(RBE_1(1,:));
                RBE_1(i+2,RBE_1_2_L+1:L)=blanks(L-RBE_1_2_L);  
        elseif rem(length_Top_Nodes-5,8)==6
            [~,RBE_1_2_L]=size(strcat(',',num2str(Top_Nodes(13+8*(i-1)+1)),',',num2str(Top_Nodes(13+8*(i-1)+2)),',',num2str(Top_Nodes(13+8*(i-1)+3)),',',num2str(Top_Nodes(13+8*(i-1)+4)),',',num2str(Top_Nodes(13+8*(i-1)+5)),',',num2str(Top_Nodes(13+8*(i-1)+6)),','));
            RBE_1(i+2,1:RBE_1_2_L)=strcat(',',num2str(Top_Nodes(13+8*(i-1)+1)),',',num2str(Top_Nodes(13+8*(i-1)+2)),',',num2str(Top_Nodes(13+8*(i-1)+3)),',',num2str(Top_Nodes(13+8*(i-1)+4)),',',num2str(Top_Nodes(13+8*(i-1)+5)),',',num2str(Top_Nodes(13+8*(i-1)+6)),',');
             L=length(RBE_1(1,:));
                RBE_1(i+2,RBE_1_2_L+1:L)=blanks(L-RBE_1_2_L);  
        elseif rem(length_Top_Nodes-5,8)==7
                [~,RBE_1_2_L]=size(strcat(',',num2str(Top_Nodes(13+8*(i-1)+1)),',',num2str(Top_Nodes(13+8*(i-1)+2)),',',num2str(Top_Nodes(13+8*(i-1)+3)),',',num2str(Top_Nodes(13+8*(i-1)+4)),',',num2str(Top_Nodes(13+8*(i-1)+5)),',',num2str(Top_Nodes(13+8*(i-1)+6)),',',num2str(Top_Nodes(13+8*(i-1)+7)),','));
            RBE_1(i+2,1:RBE_1_2_L)=strcat(',',num2str(Top_Nodes(13+8*(i-1)+1)),',',num2str(Top_Nodes(13+8*(i-1)+2)),',',num2str(Top_Nodes(13+8*(i-1)+3)),',',num2str(Top_Nodes(13+8*(i-1)+4)),',',num2str(Top_Nodes(13+8*(i-1)+5)),',',num2str(Top_Nodes(13+8*(i-1)+6)),',',num2str(Top_Nodes(13+8*(i-1)+7)),',');
             L=length(RBE_1(1,:));
                RBE_1(i+2,RBE_1_2_L+1:L)=blanks(L-RBE_1_2_L);  
        else    
            disp(' you done fucked up bro')
            you did it bad bro
        end
    end
%       RBE_1(i+2,RBE_1_2_L+1:end)=blanks(length(RBE_1(1,:))-RBE_1_2_L);
     %% Bottom RBE2
RBE_2_2_L=0;

    if length_Bot_Nodes==3
        RBE_2(1,:)=strcat('RBE2,190001,190001,',num2str(bot_Const),',',num2str(Bot_Nodes(1)),',',num2str(Bot_Nodes(2)),',',num2str(Bot_Nodes(3)),',');
    elseif length_Bot_Nodes==4
        RBE_2(1,:)=strcat('RBE2,190001,190001,',num2str(bot_Const),',',num2str(Bot_Nodes(1)),',',num2str(Bot_Nodes(2)),',',num2str(Bot_Nodes(3)),',',num2str(Bot_Nodes(4)),',');
    elseif length_Bot_Nodes>4
        RBE_2(1,:)=strcat('RBE2,190001,190001,',num2str(bot_Const),',',num2str(Bot_Nodes(1)),',',num2str(Bot_Nodes(2)),',',num2str(Bot_Nodes(3)),',',num2str(Bot_Nodes(4)),',',num2str(Bot_Nodes(5)),',');
    end
    if length_Bot_Nodes>=13
        for i=1:(length_Bot_Nodes-5)/8
                [~,RBE_2_2_L]=size(strcat(',',num2str(Bot_Nodes(6+8*(i-1))),',',num2str(Bot_Nodes(7+8*(i-1))),',',num2str(Bot_Nodes(8+8*(i-1))),',',num2str(Bot_Nodes(9+8*(i-1))),',',num2str(Bot_Nodes(10+8*(i-1))),',',num2str(Bot_Nodes(11+8*(i-1))),',',num2str(Bot_Nodes(12+8*(i-1))),',',num2str(Bot_Nodes(13+8*(i-1))),','));
                RBE_2(i+1,1:RBE_2_2_L)=strcat(',',num2str(Bot_Nodes(6+8*(i-1))),',',num2str(Bot_Nodes(7+8*(i-1))),',',num2str(Bot_Nodes(8+8*(i-1))),',',num2str(Bot_Nodes(9+8*(i-1))),',',num2str(Bot_Nodes(10+8*(i-1))),',',num2str(Bot_Nodes(11+8*(i-1))),',',num2str(Bot_Nodes(12+8*(i-1))),',',num2str(Bot_Nodes(13+8*(i-1))),',');
                L=length(RBE_2(1,:));
                RBE_2(i+1,RBE_2_2_L+1:L)=blanks(L-RBE_2_2_L);    
        end
    else
        i=0;
    end
    if rem(length_Bot_Nodes-5,8)>0
        if rem(length_Bot_Nodes-5,8)==1
            [~,RBE_2_2_L]=size(strcat(',',num2str(Bot_Nodes(13+8*(i-1)+1)),','));
            RBE_2(i+2,1:RBE_2_2_L)=strcat(',',num2str(Bot_Nodes(13+8*(i-1)+1)),',');
            L=length(RBE_2(1,:));
            RBE_2(i+2,RBE_2_2_L+1:L)=blanks(L-RBE_2_2_L);  
        elseif rem(length_Bot_Nodes-5,8)==2
            [~,RBE_2_2_L]=size(strcat(',',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),','));
            RBE_2(i+2,1:RBE_2_2_L)=strcat(',',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',');
             L=length(RBE_2(1,:));
                RBE_2(i+2,RBE_2_2_L+1:L)=blanks(L-RBE_2_2_L);  
        elseif rem(length_Bot_Nodes-5,8)==3
            [~,RBE_2_2_L]=size(strcat(',',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),','));
            RBE_2(i+2,1:RBE_2_2_L)=strcat(',',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',');
             L=length(RBE_2(1,:));
                RBE_2(i+2,RBE_2_2_L+1:L)=blanks(L-RBE_2_2_L);  
        elseif rem(length_Bot_Nodes-5,8)==4
            [~,RBE_2_2_L]=size(strcat(',',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',',num2str(Bot_Nodes(13+8*(i-1)+4)),','));
            RBE_2(i+2,1:RBE_2_2_L)=strcat(',',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',',num2str(Bot_Nodes(13+8*(i-1)+4)),',');
             L=length(RBE_2(1,:));
                RBE_2(i+2,RBE_2_2_L+1:L)=blanks(L-RBE_2_2_L);  
        elseif rem(length_Bot_Nodes-5,8)==5
            [~,RBE_2_2_L]=size(strcat(',',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',',num2str(Bot_Nodes(13+8*(i-1)+4)),',',num2str(Bot_Nodes(13+8*(i-1)+5)),','));
            RBE_2(i+2,1:RBE_2_2_L)=strcat(',',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',',num2str(Bot_Nodes(13+8*(i-1)+4)),',',num2str(Bot_Nodes(13+8*(i-1)+5)),',');
             L=length(RBE_2(1,:));
                RBE_2(i+2,RBE_2_2_L+1:L)=blanks(L-RBE_2_2_L);  
        elseif rem(length_Bot_Nodes-5,8)==6
            [~,RBE_2_2_L]=size(strcat(',',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',',num2str(Bot_Nodes(13+8*(i-1)+4)),',',num2str(Bot_Nodes(13+8*(i-1)+5)),',',num2str(Bot_Nodes(13+8*(i-1)+6)),','));
            RBE_2(i+2,1:RBE_2_2_L)=strcat(',',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',',num2str(Bot_Nodes(13+8*(i-1)+4)),',',num2str(Bot_Nodes(13+8*(i-1)+5)),',',num2str(Bot_Nodes(13+8*(i-1)+6)),',');
             L=length(RBE_2(1,:));
                RBE_2(i+2,RBE_2_2_L+1:L)=blanks(L-RBE_2_2_L);  
        elseif rem(length_Bot_Nodes-5,8)==7
                [~,RBE_2_2_L]=size(strcat(',',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',',num2str(Bot_Nodes(13+8*(i-1)+4)),',',num2str(Bot_Nodes(13+8*(i-1)+5)),',',num2str(Bot_Nodes(13+8*(i-1)+6)),',',num2str(Bot_Nodes(13+8*(i-1)+7)),','));
            RBE_2(i+2,1:RBE_2_2_L)=strcat(',',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',',num2str(Bot_Nodes(13+8*(i-1)+4)),',',num2str(Bot_Nodes(13+8*(i-1)+5)),',',num2str(Bot_Nodes(13+8*(i-1)+6)),',',num2str(Bot_Nodes(13+8*(i-1)+7)),',');
             L=length(RBE_2(1,:));
                RBE_2(i+2,RBE_2_2_L+1:L)=blanks(L-RBE_2_2_L);  
        else    
            disp(' you done fucked up bro')
            you did it bad bro
        end
    end
%       RBE_2(i+2,RBE_2_2_L+1:end)=blanks(length(RBE_2(1,:))-RBE_2_2_L);
     % RBE_1_2_L=0
%     if length_Bot_Nodes==3
%         RBE_2(1,:)=strcat('RBE2,90001,90001,',num2str(bot_Const),',',num2str(Bot_Nodes(1)),',',num2str(Bot_Nodes(2)),',',num2str(Bot_Nodes(3)),',,,');
%     elseif length_Bot_Nodes==4
%         RBE_2(1,:)=strcat('RBE2,90001,90001,',num2str(bot_Const),',',num2str(Bot_Nodes(1)),',',num2str(Bot_Nodes(2)),',',num2str(Bot_Nodes(3)),',',num2str(Bot_Nodes(4)),',,');
%     elseif length_Bot_Nodes>4
%         RBE_2(1,:)=strcat('RBE2,90001,90001,',num2str(bot_Const),',',num2str(Bot_Nodes(1)),',',num2str(Bot_Nodes(2)),',',num2str(Bot_Nodes(3)),',',num2str(Bot_Nodes(4)),',',num2str(Bot_Nodes(5)),',');
%     end
%     if length_Bot_Nodes>=13
%         for i=1:(length_Bot_Nodes-5)/8
%                 [~,RBE_1_2_L]=size(strcat('+,',num2str(Bot_Nodes(6+8*(i-1))),',',num2str(Bot_Nodes(7+8*(i-1))),',',num2str(Bot_Nodes(8+8*(i-1))),',',num2str(Bot_Nodes(9+8*(i-1))),',',num2str(Bot_Nodes(10+8*(i-1))),',',num2str(Bot_Nodes(11+8*(i-1))),',',num2str(Bot_Nodes(12+8*(i-1))),',',num2str(Bot_Nodes(13+8*(i-1))),','));
%                 RBE_2(i+1,1:RBE_1_2_L)=strcat('+,',num2str(Bot_Nodes(6+8*(i-1))),',',num2str(Bot_Nodes(7+8*(i-1))),',',num2str(Bot_Nodes(8+8*(i-1))),',',num2str(Bot_Nodes(9+8*(i-1))),',',num2str(Bot_Nodes(10+8*(i-1))),',',num2str(Bot_Nodes(11+8*(i-1))),',',num2str(Bot_Nodes(12+8*(i-1))),',',num2str(Bot_Nodes(13+8*(i-1))),',');
%                 L=length(RBE_2(:,1));
%                 RBE_2(i+1,RBE_1_2_L+1:L)=blanks(L-RBE_1_2_L);
%         end
%     else
%         i=0
%     end
%     if rem(length_Bot_Nodes-5,8)>0
%         if rem(length_Bot_Nodes-5,8)==1
%             [~,RBE_1_2_L]=size(strcat('+,',num2str(Bot_Nodes(13+8*(i-1)+1)),','));
%             RBE_2(i+2,1:RBE_1_2_L)=strcat('+,',num2str(Bot_Nodes(13+8*(i-1)+1)),',');
%              L=length(RBE_2(:,1));
%                 RBE_2(i+2,RBE_1_2_L+1:L)=blanks(L-RBE_1_2_L);
%         elseif rem(length_Bot_Nodes-5,8)==2
%             [~,RBE_1_2_L]=size(strcat('+,',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),','));
%             RBE_2(i+2,1:RBE_1_2_L)=strcat('+,',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',');
%                          L=length(RBE_2(:,1));
%                 RBE_2(i+2,RBE_1_2_L+1:L)=blanks(L-RBE_1_2_L);
%         elseif rem(length_Bot_Nodes-5,8)==3
%             [~,RBE_1_2_L]=size(strcat('+,',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),','));
%             RBE_2(i+2,1:RBE_1_2_L)=strcat('+,',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',');
%                          L=length(RBE_2(:,1));
%                 RBE_2(i+2,RBE_1_2_L+1:L)=blanks(L-RBE_1_2_L);
%         elseif rem(length_Bot_Nodes-5,8)==4
%             [~,RBE_1_2_L]=size(strcat('+,',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',',num2str(Bot_Nodes(13+8*(i-1)+4)),','));
%             RBE_2(i+2,1:RBE_1_2_L)=strcat('+,',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',',num2str(Bot_Nodes(13+8*(i-1)+4)),',');
%                          L=length(RBE_2(:,1));
%                 RBE_2(i+2,RBE_1_2_L+1:L)=blanks(L-RBE_1_2_L);
%         elseif rem(length_Bot_Nodes-5,8)==5
%             [~,RBE_1_2_L]=size(strcat('+,',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',',num2str(Bot_Nodes(13+8*(i-1)+4)),',',num2str(Bot_Nodes(13+8*(i-1)+5)),','));
%             RBE_2(i+2,1:RBE_1_2_L)=strcat('+,',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',',num2str(Bot_Nodes(13+8*(i-1)+4)),',',num2str(Bot_Nodes(13+8*(i-1)+5)),',');
%                          L=length(RBE_2(:,1));
%                 RBE_2(i+2,RBE_1_2_L+1:L)=blanks(L-RBE_1_2_L);
%         elseif rem(length_Bot_Nodes-5,8)==6
%             [~,RBE_1_2_L]=size(strcat('+,',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',',num2str(Bot_Nodes(13+8*(i-1)+4)),',',num2str(Bot_Nodes(13+8*(i-1)+5)),',',num2str(Bot_Nodes(13+8*(i-1)+6)),','));
%             RBE_2(i+2,1:RBE_1_2_L)=strcat('+,',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',',num2str(Bot_Nodes(13+8*(i-1)+4)),',',num2str(Bot_Nodes(13+8*(i-1)+5)),',',num2str(Bot_Nodes(13+8*(i-1)+6)),',');
%                          L=length(RBE_2(:,1));
%                 RBE_2(i+2,RBE_1_2_L+1:L)=blanks(L-RBE_1_2_L);
%         elseif rem(length_Bot_Nodes-5,8)==7
%                 [~,RBE_1_2_L]=size(strcat('+,',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',',num2str(Bot_Nodes(13+8*(i-1)+4)),',',num2str(Bot_Nodes(13+8*(i-1)+5)),',',num2str(Bot_Nodes(13+8*(i-1)+6)),',',num2str(Bot_Nodes(13+8*(i-1)+7)),','));
%             RBE_2(i+2,1:RBE_1_2_L)=strcat('+,',num2str(Bot_Nodes(13+8*(i-1)+1)),',',num2str(Bot_Nodes(13+8*(i-1)+2)),',',num2str(Bot_Nodes(13+8*(i-1)+3)),',',num2str(Bot_Nodes(13+8*(i-1)+4)),',',num2str(Bot_Nodes(13+8*(i-1)+5)),',',num2str(Bot_Nodes(13+8*(i-1)+6)),',',num2str(Bot_Nodes(13+8*(i-1)+7)),',');
%                      L=length(RBE_2(:,1));
%                 RBE_2(i+2,RBE_1_2_L+1:L)=blanks(L-RBE_1_2_L);
%         else    
%             disp(' you done fucked up bro')
%             you did it bad bro
%         end
%     end
% %     RBE_2(end,end)=' ';
% %     for k=1:i+2
%     RBE_2(i+2,RBE_1_2_L+1:end)=blanks(length(RBE_2(1,:))-RBE_1_2_L)
%     end
%     RBE_2(2,72)=' '
%     RBE_1(2,72)=' '

%% Helical thermal Loading
Time_Function(1,:)=zeros(1,Total_Number_of_Helicals);
% var=npermutek([1 0 -1],Total_Number_of_Helicals);
% Time_Function(2:1+length(var(:,1)),:)=var;
% Time_Function(1:3,:)=[zeros(1,Total_Number_of_Helicals);ones(1,Total_Number_of_Helicals);-ones(1,Total_Number_of_Helicals)];
%% Pos and NEg inputs
 Time_Function(2:1+Total_Number_of_Helicals,:)=eye(Total_Number_of_Helicals,Total_Number_of_Helicals);
 Time_Function(Total_Number_of_Helicals+2:1+2*Total_Number_of_Helicals,:)=-eye(Total_Number_of_Helicals,Total_Number_of_Helicals);
%% Pos inputs
% Time_Function(2:1+Total_Number_of_Helicals,:)=eye(Total_Number_of_Helicals,Total_Number_of_Helicals);
 %Time_Function(Total_Number_of_Helicals+2:1+2*Total_Number_of_Helicals,:)=-eye(Total_Number_of_Helicals,Total_Number_of_Helicals);

 % Time_Function=zeros(1,Total_Number_of_Helicals);
% a=npermutek([1 0 -1],Total_Number_of_Helicals);
% Time_Function(2:length(a)+1,:)=a;


counter=linspace(0,length(Time_Function(:,1)),1+length(Time_Function(:,1)));

   
    for i=1:Total_Number_of_Helicals
        a=length(['TABLED2 ',blanks(8-length(num2str(500+i))),num2str(500+i),'      0.                                                 ']);
        Thermal_Time_Load(1,1:a,i)=['TABLED2 ',blanks(8-length(num2str(500+i))),num2str(500+i),'      0.                                                 '];
        for j=1:length(Time_Function(:,1))/4
[~,L]=size([',',num2str(counter(4*j-3)),'.,',num2str(Time_Function(4*j-3,i)),'.,',num2str(counter(4*j-2)),'.,',num2str(Time_Function(4*j-2,i)),'.,',num2str(counter(4*j-1)),'.,',num2str(Time_Function(4*j-1,i)),'.,',num2str(counter(4*j)),'.,',num2str(Time_Function(4*j,i)),'.,']);      
Thermal_Time_Load(j+1,1:L,i)=[',',num2str(counter(4*j-3)),'.,',num2str(Time_Function(4*j-3,i)),'.,',num2str(counter(4*j-2)),'.,',num2str(Time_Function(4*j-2,i)),'.,',num2str(counter(4*j-1)),'.,',num2str(Time_Function(4*j-1,i)),'.,',num2str(counter(4*j)),'.,',num2str(Time_Function(4*j,i)),'.,'];
Thermal_Time_Load(j+1,L+1:end,i)=blanks(73-L);
% TABLED2        1      0.                                                +
% +             0.      0.      1.-3.770e6      2.-3.448e6      3.-3.104e6+
% +             4.-2.740e6      5.-2.357e6      6.-1.956e6      7.-1.538e6+
% +             8.-1.104e6      9.-655493.     10.-192566.     11. 283557.+
        end
            if rem(length(Time_Function(:,1)),4)==3
                j=j+1;
                L=length([',',num2str(counter(4*j-3)),'.,',num2str(Time_Function(4*j-3,i)),'.,',num2str(counter(4*j-2)),'.,',num2str(Time_Function(4*j-2,i)),'.,',num2str(counter(4*j-1)),'.,',num2str(Time_Function(4*j-1,i)),'.,ENDT']);
                Thermal_Time_Load(j+1,1:L,i)=[',',num2str(counter(4*j-3)),'.,',num2str(Time_Function(4*j-3,i)),'.,',num2str(counter(4*j-2)),'.,',num2str(Time_Function(4*j-2,i)),'.,',num2str(counter(4*j-1)),'.,',num2str(Time_Function(4*j-1,i)),'.,ENDT'];
                Thermal_Time_Load(j+1,L+1:end,i)=blanks(73-L);
            elseif rem(length(Time_Function(:,1)),4)==2
                j=j+1;
                L=length([',',num2str(counter(4*j-3)),'.,',num2str(Time_Function(4*j-3,i)),'.,',num2str(counter(4*j-2)),'.,',num2str(Time_Function(4*j-2,i)),'.,ENDT']);
                Thermal_Time_Load(j+1,1:L,i)=[',',num2str(counter(4*j-3)),'.,',num2str(Time_Function(4*j-3,i)),'.,',num2str(counter(4*j-2)),'.,',num2str(Time_Function(4*j-2,i)),'.,ENDT'];
                Thermal_Time_Load(j+1,L+1:end,i)=blanks(73-L);
                j=j+1;
            elseif rem(length(Time_Function(:,1)),4)==1
                j=j+1;
                L=length([',',num2str(counter(4*j-3)),'.,',num2str(Time_Function(4*j-3,i)),'.,ENDT']);
                Thermal_Time_Load(j+1,1:L,i)=[',',num2str(counter(4*j-3)),'.,',num2str(Time_Function(4*j-3,i)),'.,ENDT'];     
                Thermal_Time_Load(j+1,L+1:end,i)=blanks(73-L);
            elseif rem(length(Time_Function(:,1)),4)==0
                j=j+1;
                a=length(['        ','ENDT',blanks(61)]);
                Thermal_Time_Load(j+1,1:a,i)=['        ','ENDT',blanks(61)];
            end
    end 
   Hel_Ele_Counter=(size(Hel_Nodes)-Total_Number_of_Helicals)/Total_Number_of_Helicals;
   for i=1:Total_Number_of_Helicals
       
       
%% Default Temp
TEMPD(1,:)=['TEMPD,201,0.,202,0.,203,0.,204,0.,'];
for j=1:(Total_Number_of_Helicals-4)/4
    a=length(strcat('TEMPD,',num2str(204+j*4-3),',0.,',num2str(204+j*4-2),',0.,',num2str(204+j*4-1),',0.,',num2str(204+j*4),',0.,'));
    TEMPD(j+1,1:a)=strcat('TEMPD,',num2str(204+j*4-3),',0.,',num2str(204+j*4-2),',0.,',num2str(204+j*4-1),',0.,',num2str(204+j*4),',0.,');
end
   if (Total_Number_of_Helicals-4)/4<1
       j=0;
   end
    if rem((Total_Number_of_Helicals-4),4)>0
        j=j+1;
        if rem((Total_Number_of_Helicals-4),4)==1
           a=length(strcat('TEMPD,',num2str(204+j*4-3),',0.'));
            TEMPD(j+1,1:a)=strcat('TEMPD,',num2str(204+j*4-3),',0.');
            L=length(TEMPD(1,:));
            TEMPD(j+1,a+1:L)=blanks(L-a);
        elseif rem((Total_Number_of_Helicals-4),4)==2
            a=length(strcat('TEMPD,',num2str(204+j*4-3),',0.,',num2str(204+j*4-2),',0.,'));
            TEMPD(j+1,1:a)=strcat('TEMPD,',num2str(204+j*4-3),',0.,',num2str(204+j*4-2),',0.,');
            L=length(TEMPD(1,:));
            TEMPD(j+1,a+1:L)=blanks(L-a);
        elseif rem((Total_Number_of_Helicals-4),4)==3
            a=length(strcat('TEMPD,',num2str(204+j*4-3),',0.,',num2str(204+j*4-2),',0.,',num2str(204+j*4-1),',0.,'));
            TEMPD(j+1,1:a)=strcat('TEMPD,',num2str(204+j*4-3),',0.,',num2str(204+j*4-2),',0.,',num2str(204+j*4-1),',0.,');
            L=length(TEMPD(1,:));
            TEMPD(j+1,a+1:L)=blanks(L-a);
        else    
            disp(' you done fucked up bro')
            you did it bad bro
        end
    end

%% TEMPRB       102       1      1.      1.                                +
a=length(['TEMPRB,',num2str(200+i),',',num2str((i-1)*Hel_length+1),',1.,1.,,,,,+']);
TEMPRB(1,1:a,i)=['TEMPRB,',num2str(200+i),',',num2str((i-1)*Hel_length+1),',1.,1.,,,,,+'];
TEMPRB(1,a+1:73,i)=blanks(73-a);
a=length(['+                                                                       +']);
TEMPRB(2,1:a,i)=['+                                                                       +'];
a=length(['+,',num2str((i-1)*Hel_length+2),',THRU,',num2str(i*Hel_length-1)]);
TEMPRB(3,1:a,i)=['+,',num2str((i-1)*Hel_length+2),',THRU,',num2str(i*Hel_length-1)]  ;
TEMPRB(3,a+1:73,i)=blanks(73-a);
% TLOAD1       101     102            LOAD       2
a=length(['TLOAD1,',num2str(100+i),',',num2str(200+i),',,LOAD,',num2str(500+i)]);
TLOAD1(i,1:a)=['TLOAD1,',num2str(100+i),',',num2str(200+i),',,LOAD,',num2str(500+i)];
   end
%    TSTEP          1       6     1.1       1
TSTEP=['TSTEP,1,',num2str(length(Time_Function(:,1))-1),',1.,1'];
%% DLOAD DEv

DLOAD(1,:)=['DLOAD          2      1.      1.     101      1.     102      1.     103+'] ;
for j=1:(Total_Number_of_Helicals-3)/4
    a=length([',1.,',num2str(100+4*j),',1.,',num2str(100+4*j+1),',1.,',num2str(100+4*j+2),',1.,',num2str(100+4*j+3)]);
    DLOAD(j+1,1:a)=[',1.,',num2str(100+4*j),',1.,',num2str(100+4*j+1),',1.,',num2str(100+4*j+2),',1.,',num2str(100+4*j+3)];
    DLOAD(j+1,a+1:end)=blanks(73-a);
end
if (Total_Number_of_Helicals-3)/4<1
    j=0;
end
j=j+1;
if rem((Total_Number_of_Helicals-3),4)==3
    a=length([',1.,',num2str(100+4*j),',1.,',num2str(100+4*j+1),',1.,',num2str(100+4*j+2)]);
    DLOAD(j+1,1:a)=[',1.,',num2str(100+4*j),',1.,',num2str(100+4*j+1),',1.,',num2str(100+4*j+2)];
    DLOAD(j+1,a+1:end)=blanks(73-a);
elseif rem((Total_Number_of_Helicals-3),4)==2
    a=length([',1.,',num2str(100+4*j),',1.,',num2str(100+4*j+1)]);
    DLOAD(j+1,1:a)=[',1.,',num2str(100+4*j),',1.,',num2str(100+4*j+1)];
    DLOAD(j+1,a+1:end)=blanks(73-a);
    elseif rem((Total_Number_of_Helicals-3),4)==1
    a=length([',1.,',num2str(100+4*j)]);
    DLOAD(j+1,1:a)=[',1.,',num2str(100+4*j)];
    DLOAD(j+1,a+1:end)=blanks(73-a);
else
    disp('you fucked up bro')
    something is wrong
end

%% PbeamL Helical
% PBEAML_1=['PBEAML         1       1            Tube                                +'];  
% P_length=length(['PBEAML        1       1            Tube                                +']);
% PBEAML_1(1,1:P_length)=['PBEAML        1       1            Tube                                +'];  
% 
% [~,P_length]=size(strcat('+,',num2str(Hel_OR),',',num2str(Hel_IR)));
% PBEAML_1(2,1:P_length)=strcat('+,',num2str(Hel_OR),',',num2str(Hel_IR));
% %% PbeamL circumfrential
% PBEAML_2=['PBEAML        2       1            Tube                                +'];   
% [~,P_length]=size(strcat('+,',num2str(Cir_OR),',',num2str(Cir_IR)));
% PBEAML_2(2,1:P_length)=strcat('+,',num2str(Cir_OR),',',num2str(Cir_IR));
% P_length=length(['PBEAML*       1       1            Tube                                +']);
% PBEAML_1(1,1:P_length)=['PBEAML*       1       1            Tube                                +'];  
% 
% [~,P_length]=size(strcat('+,',num2str(Hel_OR),',',num2str(Hel_IR),','));
% PBEAML_1(2,1:P_length)=strcat('+,',num2str(Hel_OR),',',num2str(Hel_IR),',');
% %% PbeamL circumfrential
% PBEAML_2=['PBEAML*       2       1            Tube                                '];   
% [~,P_length]=size(strcat(',',num2str(Cir_OR),',',num2str(Cir_IR),','));
% PBEAML_2(2,1:P_length)=strcat(',',num2str(Cir_OR),',',num2str(Cir_IR),',');
% PBEAML,1,1,MSCBML0,Tube,,,,,                                            +                                                  
% +,0.0125,0.00625,                                                       
% P_length=length(['PBEAML,1,1,MSCBML0,Tube,,,,,                                            +']);
% PBEAML_1(1,1:P_length)=['PBEAML,1,1,MSCBML0,Tube,,,,,                                            +'];  
% 
% [~,P_length]=size(strcat('+,',num2str(Hel_OR),',',num2str(Hel_IR),','));
% PBEAML_1(2,1:P_length)=strcat('+,',num2str(Hel_OR),',',num2str(Hel_IR),',');
% P_length=length(['PBEAML,1,1,MSCBML0,Tube,,,,,']);
% PBEAML_1(1,1:P_length)=['PBEAML,1,1,MSCBML0,Tube,,,,,'];  
% [~,P_length]=size(strcat(',',num2str(Hel_OR),',',num2str(Hel_IR),','));
% PBEAML_1(2,1:P_length)=strcat(',',num2str(Hel_OR),',',num2str(Hel_IR),',');

% PBEAML         2       1 MSCBML0    TUBE                                +       
% +          .0125  .00625      0.        

PBEAML_1(2,:)=['PBEAML         1       1 MSCBML0    TUBE                                +'];  
PBEAML_1(1,1:65)=['$ Femap with NX Nastran Property 1 : BEAM Property (NASTRAN Tube)'];
b=length(num2str(Hel_OR));
c=length(num2str(Hel_IR));
a=length(['+       ',blanks(8-b),num2str(Hel_OR),blanks(8-c),num2str(Hel_IR),'      .0']);
% PBEAML_1(3,1:a)=['+       ',blanks(8-b),num2str(Hel_OR),blanks(8-c),num2str(Hel_IR),'      .0']; 
b=length(['+,',num2str(Hel_OR),',',num2str(Hel_IR),',.0']);
PBEAML_1(3,1:a)=['+,',num2str(Hel_OR),',',num2str(Hel_IR),',.0',blanks(a-b)];

% [~,P_length]=size(strcat('+,',num2str(Hel_OR),',',num2str(Hel_IR),','));
% PBEAML_1(2,1:P_length)=strcat('+,',num2str(Hel_OR),',',num2str(Hel_IR),',');
% P_length=length(['PBEAML,1,1,MSCBML0,Tube,,,,,']);
% PBEAML_1(1,1:P_length)=['PBEAML,1,1,MSCBML0,Tube,,,,,'];  
% [~,P_length]=size(strcat(',',num2str(Hel_OR),',',num2str(Hel_IR),','));
% PBEAML_1(2,1:P_length)=strcat(',',num2str(Hel_OR),',',num2str(Hel_IR),',');

PBEAML_2(2,:)=['PBEAML         2       1 MSCBML0    TUBE                                +'];  
PBEAML_2(1,1:65)=['$ Femap with NX Nastran Property 2 : BEAM Property (NASTRAN Tube)'];
b=length(num2str(Cir_OR));
c=length(num2str(Cir_IR));
a=length(['+       ',blanks(8-b),num2str(Cir_OR),blanks(8-c),num2str(Cir_IR),'      .0']);
% PBEAML_2(3,1:a)=['+       ',blanks(8-b),num2str(Cir_OR),blanks(8-c),num2str(Cir_IR),'      .0'];
b=length(['+,',num2str(Cir_OR),',',num2str(Cir_IR),',.0']);
PBEAML_2(3,1:a)=['+,',num2str(Cir_OR),',',num2str(Cir_IR),',.0',blanks(a-b)]; 
% 
% [~,P_length]=size(strcat('+,',num2str(Hel_OR),',',num2str(Hel_IR),','));
% PBEAML_1(2,1:P_length)=strcat('+,',num2str(Hel_OR),',',num2str(Hel_IR),',');
% P_length=length(['PBEAML,1,1,MSCBML0,Tube,,,,,']);
% PBEAML_1(1,1:P_length)=['PBEAML,1,1,MSCBML0,Tube,,,,,'];  
% [~,P_length]=size(strcat(',',num2str(Hel_OR),',',num2str(Hel_IR),','));
% PBEAML_1(2,1:P_length)=strcat(',',num2str(Hel_OR),',',num2str(Hel_IR),',');
%% PbeamL circumfrential
% PBEAML_2=['PBEAML,2,1,MSCBML0,Tube,,,,,'];   
% [~,P_length]=size(strcat(',',num2str(Cir_OR),',',num2str(Cir_IR),','));
% PBEAML_2(2,1:P_length)=strcat(',',num2str(Cir_OR),',',num2str(Cir_IR),',');
% %% PBEAM Development
% A=pi*(Hel_OR^2-Hel_IR^2)
% I2 = pi/4*(Hel_OR^4 - Hel_IR^4), I1 = I2
% J=pi*(Hel_OR^4-Hel_IR^4)/2
% PBEAM_1=['PBEAM*,1,1,',num2str(A),',',num2str(I1)]
% a=length(strcat('*,',num2str(I2),',0.,',num2str(J),',0.'))
% PBEAM_1(2,1:a)=['*,',num2str(I2),',0.,',num2str(J),',0.']
% a=length(strcat('*,',num2str(Hel_OR),',0.,0.,',num2str(Hel_OR)))
% PBEAM_1(3,1:a)=['*,',num2str(Hel_OR),',0.,0.,',num2str(Hel_OR)]
% a=length(strcat('*,0.,0.,',num2str(-Hel_OR),',0.'))
% PBEAM_1(4,1:a)=['*,0.,0.,',num2str(-Hel_OR),',0.']
% a=length(strcat('*,0.,',num2str(-Hel_OR),',0.,0.'))
% PBEAM_1(5,1:a)=['*,0.,',num2str(-Hel_OR)]
% 
% 
% PBEAM  *               1               1  3.68155335E-04  1.79762942E-08
% *         1.79762925E-08  0.00000000E+00  3.59525885E-08  0.00000000E+00
% *         1.25000002E-02  0.00000000E+00  0.00000000E+00  1.25000002E-02
% *        -1.25000002E-02  0.00000000E+00  0.00000000E+00 -1.25000002E-02
% *         5.88560820E-01  5.88558912E-01  0.00000000E+00  0.00000000E+00
% *         0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
% *         0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
% *         0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
% PBEAM_1
%      PBEAML_1(1,73:74)='++'
%      PBEAML_2(1,73:75)='+++'
%% Dlmwrite 
dlmwrite(save_File_ID,'              ','newline','pc','delimiter','','-append') 
dlmwrite(save_File_ID,'$ Updated Section','newline','pc','delimiter','','-append') 
% dlmwrite(save_File_ID,PBEAML_1,'newline','pc','delimiter','','-append') 
% dlmwrite(save_File_ID,PBEAML_1,'newline','pc','delimiter','','-append') 
% dlmwrite(save_File_ID,PBEAML_2,'newline','pc','delimiter','','-append') 
dlmwrite(save_File_ID,PBEAML_1(1,1:65),'newline','pc','delimiter','','-append') 
dlmwrite(save_File_ID,PBEAML_1(2,1:end),'newline','pc','delimiter','','-append') 
dlmwrite(save_File_ID,PBEAML_1(3,1:32),'newline','pc','delimiter','','-append') 
dlmwrite(save_File_ID,PBEAML_2(1,1:65),'newline','pc','delimiter','','-append') 
dlmwrite(save_File_ID,PBEAML_2(2,1:end),'newline','pc','delimiter','','-append') 
dlmwrite(save_File_ID,PBEAML_2(3,1:32),'newline','pc','delimiter','','-append') 
for i=1:length(TEMPRB(1,1,:))
    dlmwrite(save_File_ID,TEMPRB(:,:,i),'newline','pc','delimiter','','-append') 
    dlmwrite(save_File_ID,Thermal_Time_Load(:,:,i),'newline','pc','delimiter','','-append') 
end
dlmwrite(save_File_ID,TEMPD,'newline','pc','delimiter','','-append')
dlmwrite(save_File_ID,TLOAD1,'newline','pc','delimiter','','-append')
dlmwrite(save_File_ID,TSTEP,'newline','pc','delimiter','','-append') 
dlmwrite(save_File_ID,DLOAD,'newline','pc','delimiter','','-append')
dlmwrite(save_File_ID,RBE_1,'newline','pc','delimiter','','-append') 
dlmwrite(save_File_ID,RBE_2,'newline','pc','delimiter','','-append')
dlmwrite(save_File_ID,RBE_Grid,'newline','pc','delimiter','','-append')
dlmwrite(save_File_ID,char(CBEAM_Build(:)),'newline','pc','delimiter','','-append') 
dlmwrite(save_File_ID,char(Grid_Build(:)),'newline','pc','delimiter','','-append') 
dlmwrite(save_File_ID,'ENDDATA','newline','pc','delimiter','','-append') 