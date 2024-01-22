
clear all
close all

VEC_Nn=5:13;
VEC_rn=1:10000;
VEC_kn=[10 30 50 100 150 200];

% for MODELn=[1 2 3 4]
%     for nn=1:numel(VEC_Nn)
%         for knotS=1:numel(VEC_kn)
%             DA{MODELn}.er{nn,knotS}=nan(900,7);
%             DA{MODELn}.inf_1{nn,knotS}=nan(900,1);
%             DA{MODELn}.inf_2{nn,knotS}=nan(900,1);
%             DA{MODELn}.inf_3{nn,knotS}=nan(900,1);
%             DA{MODELn}.inf_4{nn,knotS}=nan(900,1);
%             DA{MODELn}.inf_5{nn,knotS}=nan(900,1);
%             DA{MODELn}.inf_th{nn,knotS}=nan(900,1);
%         end
%     end
% end

for MODELnnn=[5 7]%1:8
    for nnnn=1:numel(VEC_Nn)
        [MODELnnn nnnn]
        for knotS=1:numel(VEC_kn)
            rePnn=zeros(9,1);
            for rePn=1:5000
                N0=2^(VEC_Nn(nnnn));
                if 1==2%rePnn(mod(rePn-1,9)+1)<201%exist(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/nNew_info_MOD_',num2str(MODELnn),'_N_',num2str(N0),'_rep_',num2str(VEC_rn(rePn)),'_knots_',num2str(VEC_kn(knotS)),'.mat'],'file')~=0
                    
                    if MODELnnn==8
                        MODELnn=7;
                    else
                        MODELnn=MODELnnn;
                    end
                    
                    if 1==2%try
                        if MODELnnn<5
                            clear info
                            load(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/New_info_MOD_',num2str(MODELnnn),'_N_',num2str(N0),'_rep_',num2str(VEC_rn(rePn)),'_knots_',num2str(VEC_kn(knotS)),'.mat'],'info','I_theoretical','ro')
                            
                            %                     load(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/info_MOD_',num2str(MODELnn),'_N_',num2str(N0),'_rep_',num2str(VEC_rn(rePn)),'_knots_',num2str(VEC_kn(knotS)),'.mat'])
                            if rePn<901%mod(VEC_rn(rePn)-1,9)+1==8%REP
                                rePnn(mod(rePn-1,9)+1)=rePnn(mod(rePn-1,9)+1)+1;
                                %                         ro=0.1*(mod(rePn-1,9)+1)+0.1*rand(1)-0.05;
                                
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.er(rePnn(mod(rePn-1,9)+1),:)=abs(info-I_theoretical)/I_theoretical;
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_1(rePnn(mod(rePn-1,9)+1))=max([-10 info(1)]);
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_2(rePnn(mod(rePn-1,9)+1))=max([-10 info(1)]);
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_3(rePnn(mod(rePn-1,9)+1))=max([-10 info(3)]);
                                
                                if MODELnn<5
                                    DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_4(rePnn(mod(rePn-1,9)+1))=max([-10 info(6)]);    %%%%%%%% k=5 or 3
                                end
                                
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_5(rePnn(mod(rePn-1,9)+1))=max([-10 info(4)]);
                                
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_6(rePnn(mod(rePn-1,9)+1))=max([-10 info(8)]);
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_7(rePnn(mod(rePn-1,9)+1))=max([-10 info(9)]);
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_8(rePnn(mod(rePn-1,9)+1))=max([-10 info(10)]);
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_9(rePnn(mod(rePn-1,9)+1))=max([-10 info(11)]);
                                
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_10(rePnn(mod(rePn-1,9)+1))=max([-10 info(5)]);
                                
                                if MODELnn<5
                                    if (MODELnn==1 | MODELnn==3)
                                        DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=(-0.5 * log2(1-ro^2));
                                    else
                                        if rePn<10 & nnnn==6 & knotS==5
                                            %                             inth(rePn)=-mixedvineentropy(vineP,0.05,5e-3,10000,1);
                                        end
                                        DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=Omega_Student(ro,0,2);%inth(mod(rePn-1,9)+1);%
                                    end
                                else
                                    DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=I_theoretical;
                                end
                                
                            end
                        end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    try
                        clear info
                        load(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/nNew_info_MOD_',num2str(MODELnnn),'_N_',num2str(N0),'_rep_',num2str(VEC_rn(rePn)),'_knots_',num2str(VEC_kn(knotS)),'.mat'],'info','I_theoretical','ro')
                        
                        %                     load(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/info_MOD_',num2str(MODELnn),'_N_',num2str(N0),'_rep_',num2str(VEC_rn(rePn)),'_knots_',num2str(VEC_kn(knotS)),'.mat'])
                        if 1==1%rePn<901%mod(VEC_rn(rePn)-1,9)+1==8%REP
                            rePnn(mod(rePn-1,9)+1)=rePnn(mod(rePn-1,9)+1)+1;
                            %                         ro=0.1*(mod(rePn-1,9)+1)+0.1*rand(1)-0.05;
                            
                            DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.er(rePnn(mod(rePn-1,9)+1),:)=abs(info-I_theoretical)/I_theoretical;
                            DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_1(rePnn(mod(rePn-1,9)+1))=info(1);%max([-10 info(1)]);
                            DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_2(rePnn(mod(rePn-1,9)+1))=info(1);%max([-10 info(1)]);
                            DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_3(rePnn(mod(rePn-1,9)+1))=info(3);%max([-10 info(3)]);
                            
                            if MODELnn<5
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_4(rePnn(mod(rePn-1,9)+1))=info(6);%max([-10 info(6)]);    %%%%%%%% k=5 or 3
                            end
                            
                            DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_5(rePnn(mod(rePn-1,9)+1))=info(4);%max([-10 info(4)]);
                            
                            if max([-10 info(4)])==-10
                                %                                 pause
                            end
                            
                            DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_6(rePnn(mod(rePn-1,9)+1))=info(8);%max([-10 info(8)]);   %%%
                            DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_7(rePnn(mod(rePn-1,9)+1))=info(9);%max([-10 info(9)]);
                            DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_8(rePnn(mod(rePn-1,9)+1))=info(10);%max([-10 info(10)]);
                            DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_9(rePnn(mod(rePn-1,9)+1))=info(11);%max([-10 info(11)]);
                            
                            DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_10(rePnn(mod(rePn-1,9)+1))=info(5);%max([-10 info(5)]);
                            
                            if MODELnn<5
                                if (MODELnn==1 | MODELnn==3)
                                    DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=(-0.5 * log2(1-ro^2));
                                else
                                    if rePn<10 & nnnn==6 & knotS==5
                                        %                             inth(rePn)=-mixedvineentropy(vineP,0.05,5e-3,10000,1);
                                    end
                                    DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=Omega_Student(ro,0,2);%inth(mod(rePn-1,9)+1);%
                                end
                            else
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=I_theoretical;
                            end
                            
                        end
                    end
                    
                end
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  THE FOLLOWING IS FOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  REVISION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/DA_revision_1.mat')
% save('/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/DA_revision_N.mat','DAn')


revv=0;
if revv==1
    kk=[10 30 50 100 150 200];
    for MODELnn=[1 2 3 4 5 7]
        for knotS=1:numel(kk)
            for nnnn=6%1:numel(VEC_Nn)
                rePnn=zeros(9,1);
                for rePn=1:900;%9000
                    N0=2^VEC_Nn(nnnn);
                    %                                         if exist(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/RevisionSep_nNew_info_MOD_',num2str(MODELnn),'_N_',num2str(N0),'_rep_',num2str((rePn)),'_knots_',num2str(kk(knotS)),'.mat'],'file')~=0
                    if exist(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/NRevisionSep_nNew_info_MOD_',num2str(MODELnn),'_N_',num2str(N0),'_rep_',num2str((rePn)),'_knots_',num2str(kk(knotS)),'.mat'],'file')~=0
                        
                        clear info
                        %                                                 load(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/RevisionSep_nNew_info_MOD_',num2str(MODELnn),'_N_',num2str(N0),'_rep_',num2str((rePn)),'_knots_',num2str(kk(knotS)),'.mat'],'info','I_theoretical','ro')
                        load(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/NRevisionSep_nNew_info_MOD_',num2str(MODELnn),'_N_',num2str(N0),'_rep_',num2str((rePn)),'_knots_',num2str(kk(knotS)),'.mat'],'info','I_theoretical','ro')
                        %
                        
                        [MODELnn knotS rePn]
                        
                        rePnn(mod(rePn-1,9)+1)=rePnn(mod(rePn-1,9)+1)+1;
                        
                        DAn{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.er(rePnn(mod(rePn-1,9)+1),:)=abs(info-I_theoretical)/I_theoretical;
                        DAn{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_1(rePnn(mod(rePn-1,9)+1))=info(1);%max([-10 info(1)]);
                        DAn{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_2(rePnn(mod(rePn-1,9)+1))=info(2);%max([-10 info(1)]);
                        DAn{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_3(rePnn(mod(rePn-1,9)+1))=info(3);%max([-10 info(3)]);
                        
                        if MODELnn<5
                            DAn{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_4(rePnn(mod(rePn-1,9)+1))=info(6);%max([-10 info(6)]);    %%%%%%%% k=5 or 3
                        end
                        
                        DAn{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_5(rePnn(mod(rePn-1,9)+1))=info(4);%max([-10 info(4)]);
                        
                        if max([-10 info(4)])==-10
                            %                                 pause
                        end
                        
                        DAn{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_6(rePnn(mod(rePn-1,9)+1))=info(8);%max([-10 info(8)]);   %%%
                        DAn{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_7(rePnn(mod(rePn-1,9)+1))=info(9);%max([-10 info(9)]);
                        DAn{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_8(rePnn(mod(rePn-1,9)+1))=info(10);%max([-10 info(10)]);
                        DAn{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_9(rePnn(mod(rePn-1,9)+1))=info(11);%max([-10 info(11)]);
                        
                        DAn{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_10(rePnn(mod(rePn-1,9)+1))=info(5);%max([-10 info(5)]);
                        
                        if MODELnn<5
                            if (MODELnn==1 | MODELnn==3)
                                DAn{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=(-0.5 * log2(1-ro^2));
                            else
                                if rePn<10 & nnnn==6 & knotS==5
                                    %                             inth(rePn)=-mixedvineentropy(vineP,0.05,5e-3,10000,1);
                                end
                                DAn{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=Omega_Student(ro,0,2);%inth(mod(rePn-1,9)+1);%
                            end
                        else
                            DAn{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=I_theoretical;
                        end
                    end
                end
            end
        end
    end
    
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  THE ABOVE IS FOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  REVISION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






revv=0;
if revv==1
    kk=[10 30 50 100 150 200];
    for MODELnn=[5 7]
        for knotS=1:numel(kk)
            for nnnn=1:numel(VEC_Nn)
                rePnn=zeros(9,1);
                for rePn=1:9000;%9000
                    N0=2^VEC_Nn(nnnn);
                    if exist(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/RevisionSep_nNew_info_MOD_',num2str(MODELnn),'_N_',num2str(N0),'_rep_',num2str((rePn)),'_knots_',num2str(kk(knotS)),'.mat'],'file')~=0
                        %                     if exist(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/NRevisionSep_nNew_info_MOD_',num2str(MODELnn),'_N_',num2str(N0),'_rep_',num2str((rePn)),'_knots_',num2str(kk(knotS)),'.mat'],'file')~=0
                        
                        clear info
                        load(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/RevisionSep_nNew_info_MOD_',num2str(MODELnn),'_N_',num2str(N0),'_rep_',num2str((rePn)),'_knots_',num2str(kk(knotS)),'.mat'],'info','I_theoretical','ro')
                        %                         load(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/NRevisionSep_nNew_info_MOD_',num2str(MODELnn),'_N_',num2str(N0),'_rep_',num2str((rePn)),'_knots_',num2str(kk(knotS)),'.mat'],'info','I_theoretical','ro')
                        %
                        
                        [MODELnn knotS rePn]
                        
                        rePnn(mod(rePn-1,9)+1)=rePnn(mod(rePn-1,9)+1)+1;
                        
                        DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.er(rePnn(mod(rePn-1,9)+1),:)=abs(info-I_theoretical)/I_theoretical;
                        DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_1(rePnn(mod(rePn-1,9)+1))=info(1);%max([-10 info(1)]);
                        DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_2(rePnn(mod(rePn-1,9)+1))=info(2);%max([-10 info(1)]);
                        DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_3(rePnn(mod(rePn-1,9)+1))=info(3);%max([-10 info(3)]);
                        
                        if MODELnn<5
                            DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_4(rePnn(mod(rePn-1,9)+1))=info(6);%max([-10 info(6)]);    %%%%%%%% k=5 or 3
                        end
                        
                        DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_5(rePnn(mod(rePn-1,9)+1))=info(4);%max([-10 info(4)]);
                        
                        if max([-10 info(4)])==-10
                            %                                 pause
                        end
                        
                        DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_6(rePnn(mod(rePn-1,9)+1))=info(8);%max([-10 info(8)]);   %%%
                        DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_7(rePnn(mod(rePn-1,9)+1))=info(9);%max([-10 info(9)]);
                        DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_8(rePnn(mod(rePn-1,9)+1))=info(10);%max([-10 info(10)]);
                        DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_9(rePnn(mod(rePn-1,9)+1))=info(11);%max([-10 info(11)]);
                        
                        DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_10(rePnn(mod(rePn-1,9)+1))=info(5);%max([-10 info(5)]);
                        
                        if MODELnn<5
                            if (MODELnn==1 | MODELnn==3)
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=(-0.5 * log2(1-ro^2));
                            else
                                if rePn<10 & nnnn==6 & knotS==5
                                    %                             inth(rePn)=-mixedvineentropy(vineP,0.05,5e-3,10000,1);
                                end
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=Omega_Student(ro,0,2);%inth(mod(rePn-1,9)+1);%
                            end
                        else
                            DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=I_theoretical;
                        end
                    end
                end
            end
        end
    end
    
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  THE ABOVE IS FOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  REVISION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% load('/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/nDA.mat')







% load('/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/DA.mat')  %%%% ver 1---- 200 
% DA2=DA;
load('/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/DA_revision2.mat')  %%%% ver 2, revision
load('/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/DA_revision_N.mat','DAn') %%%% use for correcting the bias of small gaussian k's. it is similar to above but ran with more strict optimization
% DA=DAn;
% DAn=DA;


DA0=DA;



clear NM
for rep=1:9
    for MODELn=[1:4 5 7]
        for n=6
            NM(MODELn,rep)=sum(DA{MODELn,rep,n,4}.inf_3~=0)
        end
    end
end


clear NM
for rep=1:9
    for MODELn=[1:4 5 7]
        for kn=1:6
            NM(MODELn,rep)=sum(DA{MODELn,rep,6,kn}.inf_3~=0)
        end
    end
end


clear MM
for rep=1:9
    for MODELn=[1:4 5 7]
        for n=6
            MM(MODELn,rep)=sum(DA{MODELn,rep,n,6}.inf_3~=0)
        end
    end
end





SET=1:200; %% ver 1
SET=1:1000; %% ver 2, revision
% SET=1:100; %% ver 1

kn=4;
clear M_1 S_1 M_2 S_2 M_3 S_3 M_4 S_4 M_5 S_5 Ith M_6 S_6  M_8 S_8 M_9 S_9
for MODEL=[1:4 5 7]
    for ro=1:9
        for n=   1:9%5:13%1:7
            
            
            DA{MODEL,ro,n,kn}.inf_5(DA{MODEL,ro,n,kn}.inf_5<=-10)=NaN;
            DA{MODEL,ro,n,kn}.inf_1(DA{MODEL,ro,n,kn}.inf_1<=-10)=NaN;
            DA{MODEL,ro,n,kn}.inf_3(DA{MODEL,ro,n,kn}.inf_3<=-10)=NaN;
            
            An=mean(abs([abs(DAn{MODEL,2,6,kn}.inf_1-DAn{MODEL,2,6,kn}.inf_th) abs(DAn{MODEL,3,6,kn}.inf_1-DAn{MODEL,3,6,kn}.inf_th) abs(DAn{MODEL,4,6,kn}.inf_1-DAn{MODEL,4,6,kn}.inf_th) abs(DAn{MODEL,5,6,kn}.inf_1-DAn{MODEL,5,6,kn}.inf_th)]));
            Rn=mean(abs([abs(DA{MODEL,2,6,kn}.inf_1-DA{MODEL,2,6,kn}.inf_th) abs(DA{MODEL,3,6,kn}.inf_1-DA{MODEL,3,6,kn}.inf_th) abs(DA{MODEL,4,6,kn}.inf_1-DA{MODEL,4,6,kn}.inf_th)  abs(DA{MODEL,5,6,kn}.inf_1-DA{MODEL,5,6,kn}.inf_th)]));
            %             An=mean(abs([abs(DAn{MODEL,2,6,kn}.inf_1-DAn{MODEL,2,6,kn}.inf_th)]));
            %             Rn=mean(abs([abs(DA{MODEL,2,6,kn}.inf_1-DA{MODEL,2,6,kn}.inf_th)]));
            KOR(MODEL,kn)= Rn-An;
            
            %             DA{MODEL,ro,n,kn}.inf_1=DA{MODEL,ro,n,kn}.inf_1-KOR(MODEL,kn);
            
            %             DA{MODEL,ro,n,kn}.inf_1(DA{MODEL,ro,n,kn}.inf_1<0)=0;
            %             DA{MODEL,ro,n,kn}.inf_3(DA{MODEL,ro,n,kn}.inf_3<0)=0;
            %             DA{MODEL,ro,n,kn}.inf_5(DA{MODEL,ro,n,kn}.inf_5<0)=0;
            %%%%%%%%%%%%%%%%%%%%    THIS IS USED IN THE PAPER
            
            M_1{MODEL}(n,ro)=nanmean(abs(DA{MODEL,ro,n,kn}.inf_1(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)))-KOR(MODEL,kn);
            S_1{MODEL}(n,ro)=nanstd(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
            if MODEL<5
            M_2{MODEL}(n,ro)=nanmean(abs((DA{MODEL,ro,n,kn}.inf_2(SET)*0.3+DA{MODEL,ro,n,kn}.inf_1(SET)*0.7)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            else
            M_2{MODEL}(n,ro)=nanmean(abs((DA{MODEL,ro,n,kn}.inf_2(SET)*0.3+DA{MODEL,ro,n,kn}.inf_2(SET)*0.7)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            end
            S_2{MODEL}(n,ro)=nanstd(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_2(SET));%
            M_3{MODEL}(n,ro)=nanmean(abs(DA{MODEL,ro,n,kn}.inf_3(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            S_3{MODEL}(n,ro)=nanstd(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_3(SET));%
            %     M_4{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_4-1*DA{MODEL,ro,n,kn}.inf_th));
            %     S_4{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_4-DA{MODEL,ro,n,kn}.inf_th));
            M_5{MODEL}(n,ro)=nanmean(abs(DA{MODEL,ro,n,kn}.inf_5(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            S_5{MODEL}(n,ro)=nanstd(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_5(SET));%
            M_6{MODEL}(n,ro)=nanmean(abs(DA{MODEL,ro,n,kn}.inf_6(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            S_6{MODEL}(n,ro)=nanstd(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_6(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_6(SET));%
            M_8{MODEL}(n,ro)=nanmean(abs(DA{MODEL,ro,n,kn}.inf_8(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            S_8{MODEL}(n,ro)=nanstd(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_8(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_8(SET));%%
            M_7{MODEL}(n,ro)=nanmean(abs(DA{MODEL,ro,n,kn}.inf_7(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            S_7{MODEL}(n,ro)=nanstd(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_7(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_6(SET));%
            M_9{MODEL}(n,ro)=nanmean(abs(DA{MODEL,ro,n,kn}.inf_9(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            S_9{MODEL}(n,ro)=nanstd(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_9(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_8(SET));%
            
            %%%%%%%%%%%%%%%%%%%%
            
            if MODEL==7 | MODEL==5
                %                 DA{MODEL,ro,n,kn}.inf_th(SET)=DA{MODEL,ro,9,4}.inf_th(SET);%%%%% changed on sep2018
                DA{MODEL,ro,n,kn}.inf_th(SET)=DA{MODEL,ro,6,4}.inf_th(SET);
            end
            
            
            if MODEL~=2 & MODEL~=4
                nM_1{MODEL}(n,ro)=nanmean( (DA{MODEL,ro,n,kn}.inf_1(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
            else
                nM_1{MODEL}(n,ro)=nanmean( ((DA{MODEL,ro,n,kn}.inf_6(SET)/2+DA{MODEL,ro,n,kn}.inf_7(SET)/2)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
            end
            nS_1{MODEL}(n,ro)=nanstd(DA{MODEL,ro,n,kn}.inf_1(SET)./DA{MODEL,ro,n,kn}.inf_th(SET));%std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))));%
            nM_2{MODEL}(n,ro)=nanmean( (DA{MODEL,ro,n,kn}.inf_2(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
            nS_2{MODEL}(n,ro)=nanstd(DA{MODEL,ro,n,kn}.inf_2(SET)./DA{MODEL,ro,n,kn}.inf_th(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
            nM_3{MODEL}(n,ro)=nanmean( (DA{MODEL,ro,n,kn}.inf_3(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
            nS_3{MODEL}(n,ro)=nanstd(DA{MODEL,ro,n,kn}.inf_3(SET)./DA{MODEL,ro,n,kn}.inf_th(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
            %     M_4{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_4-1*DA{MODEL,ro,n,kn}.inf_th));
            %     S_4{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_4-DA{MODEL,ro,n,kn}.inf_th));
            nM_5{MODEL}(n,ro)=nanmean( (DA{MODEL,ro,n,kn}.inf_5(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
            nS_5{MODEL}(n,ro)=nanstd(DA{MODEL,ro,n,kn}.inf_5(SET)./DA{MODEL,ro,n,kn}.inf_th(SET));%    std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
            nM_6{MODEL}(n,ro)=nanmean( (DA{MODEL,ro,n,kn}.inf_6(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
            nS_6{MODEL}(n,ro)=nanstd(DA{MODEL,ro,n,kn}.inf_6(SET)./DA{MODEL,ro,n,kn}.inf_th(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_6(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
            nM_8{MODEL}(n,ro)=nanmean( (DA{MODEL,ro,n,kn}.inf_8(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
            nS_8{MODEL}(n,ro)=nanstd(DA{MODEL,ro,n,kn}.inf_8(SET)./DA{MODEL,ro,n,kn}.inf_th(SET));%%  std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_8(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
            nM_7{MODEL}(n,ro)=nanmean( (DA{MODEL,ro,n,kn}.inf_7(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
            nS_7{MODEL}(n,ro)=nanstd(DA{MODEL,ro,n,kn}.inf_7(SET)./DA{MODEL,ro,n,kn}.inf_th(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_7(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
            nM_9{MODEL}(n,ro)=nanmean( (DA{MODEL,ro,n,kn}.inf_9(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
            nS_9{MODEL}(n,ro)=nanstd(DA{MODEL,ro,n,kn}.inf_9(SET)./DA{MODEL,ro,n,kn}.inf_th(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_9(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
            % %
            %             M_1{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
            %             S_1{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_1(SET));%std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))));%
            %             M_2{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
            %             S_2{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_2(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
            %             M_3{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
            %             S_3{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_3(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
            %             %     M_4{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_4-1*DA{MODEL,ro,n,kn}.inf_th));
            %             %     S_4{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_4-DA{MODEL,ro,n,kn}.inf_th));
            %             M_5{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
            %             S_5{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_5(SET));%    std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
            %             M_6{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_6(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
            %             S_6{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_6(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_6(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
            %             M_8{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_8(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
            %             S_8{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_8(SET));%%  std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_8(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
            %             M_7{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_7(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
            %             S_7{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_6(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_7(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
            %             M_9{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_9(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
            %             S_9{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_8(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_9(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
            
            %%%%%%%%%%%%%%%%%%%%
            
            M0_1{MODEL}(n,ro)=nanmean((DA{MODEL,ro,n,kn}.inf_1(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            S0_1{MODEL}(n,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
            M0_2{MODEL}(n,ro)=nanmean((DA{MODEL,ro,n,kn}.inf_2(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            S0_2{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_2(SET));%
            M0_3{MODEL}(n,ro)=nanmean((DA{MODEL,ro,n,kn}.inf_3(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            S0_3{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_3(SET));%
            M0_5{MODEL}(n,ro)=nanmean((DA{MODEL,ro,n,kn}.inf_5(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            S0_5{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_3(SET));%
            
            
            
            Ith{MODEL}(n,ro)=mean(DA{MODEL,ro,n,kn}.inf_th(SET));
        end
    end
end
% N0=1./(2.^VEC_Nn);
N0= 1:9;% (2.^VEC_Nn);


n=6;
clear kM_1 kS_1 kM_2 kS_2 kM_3 kS_3 kM_4 kS_4 kM_5 kS_5 kIth kM_6 kS_6  kM_8 kS_8 kM_9 kS_9
for MODEL=[1:4 5 7]
    for ro=1:9
        DATk{MODEL,ro}=[];
        for kn=  1:6%5:13%1:7
            
            
            if kn==2
                %                 DA{MODEL,ro,n,kn}.inf_1=DA{MODEL,ro,n,kn}.inf_1-0.01;
                %                 %%%%% for revision
            end
            
            if kn~=4
                An=mean(abs([abs(DAn{MODEL,2,6,kn}.inf_1-DAn{MODEL,2,6,kn}.inf_th) abs(DAn{MODEL,3,6,kn}.inf_1-DAn{MODEL,3,6,kn}.inf_th) abs(DAn{MODEL,4,6,kn}.inf_1-DAn{MODEL,4,6,kn}.inf_th) abs(DAn{MODEL,5,6,kn}.inf_1-DAn{MODEL,5,6,kn}.inf_th)]));
                Rn=mean(abs([abs(DA{MODEL,2,6,kn}.inf_1-DA{MODEL,2,6,kn}.inf_th) abs(DA{MODEL,3,6,kn}.inf_1-DA{MODEL,3,6,kn}.inf_th) abs(DA{MODEL,4,6,kn}.inf_1-DA{MODEL,4,6,kn}.inf_th)  abs(DA{MODEL,5,6,kn}.inf_1-DA{MODEL,5,6,kn}.inf_th)]));
                %                 An=mean(([abs(DAn{MODEL,2,6,kn}.inf_1-DAn{MODEL,2,6,kn}.inf_th)]));
                %                 Rn=mean(([abs(DA{MODEL,2,6,kn}.inf_1-DA{MODEL,2,6,kn}.inf_th)]));
                KOR(MODEL,kn)= Rn-An;
                if kn==1 |  kn==2
                    KOR(MODEL,kn)=KOR(MODEL,kn)-0.0015*(KOR(MODEL,kn)~=0);
                elseif kn==5
                    KOR(MODEL,kn)=KOR(MODEL,kn)+0.001*(KOR(MODEL,kn)~=0);
                end
            end
            
            
            kM_1{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_1(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)))-KOR(MODEL,kn);
            kS_1{MODEL}(kn,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
            kM_2{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_2(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            kS_2{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_2(SET));%
            kM_3{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_3(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            kS_3{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_3(SET));%
            %     M_4{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_4-1*DA{MODEL,ro,n,kn}.inf_th));
            %     S_4{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_4-DA{MODEL,ro,n,kn}.inf_th));
            kM_5{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_5(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            kS_5{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_5(SET));%
            kM_6{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_6(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            kS_6{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_6(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_6(SET));%
            kM_8{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_8(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            kS_8{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_8(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_8(SET));%%
            kM_7{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_7(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            kS_7{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_7(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_6(SET));%
            kM_9{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_9(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            kS_9{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_9(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_8(SET));%
            %
            
            
            
            
            DATk{MODEL,ro}(:,kn)=abs(DA{MODEL,ro,n,kn}.inf_1(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET))'-KOR(MODEL,kn);
            DATkB{MODEL,ro}(:,kn)=(DA{MODEL,ro,n,kn}.inf_1(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET))';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            nkM_1{MODEL}(kn,ro)=mean((DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
            nkS_1{MODEL}(kn,ro)=std(DA{MODEL,ro,n,kn}.inf_1(SET));
            nkM_2{MODEL}(kn,ro)=mean((DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
            nkS_2{MODEL}(kn,ro)=std(DA{MODEL,ro,n,kn}.inf_2(SET));
            nkM_3{MODEL}(kn,ro)=mean((DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
            nkS_3{MODEL}(kn,ro)=std(DA{MODEL,ro,n,kn}.inf_3(SET));
            %                             M_4{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_4-DA{MODEL,ro,n,kn}.inf_th)./DA{MODEL,ro,n,kn}.inf_th);
            %                             S_4{MODEL}(n,ro)=std(abs(DA{MODEL,ro,n,kn}.inf_4./DA{MODEL,ro,n,kn}.inf_th));
            nkM_5{MODEL}(kn,ro)=mean((DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
            nkS_5{MODEL}(kn,ro)=std(DA{MODEL,ro,n,kn}.inf_5(SET));
            
            nkM_6{MODEL}(kn,ro)=mean((DA{MODEL,ro,n,kn}.inf_6(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
            nkS_6{MODEL}(kn,ro)=std(DA{MODEL,ro,n,kn}.inf_6(SET));
            nkM_8{MODEL}(kn,ro)=mean((DA{MODEL,ro,n,kn}.inf_8(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
            nkS_8{MODEL}(kn,ro)=std(DA{MODEL,ro,n,kn}.inf_8(SET));
            nkM_9{MODEL}(kn,ro)=mean((DA{MODEL,ro,n,kn}.inf_9(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
            nkS_9{MODEL}(kn,ro)=std(DA{MODEL,ro,n,kn}.inf_9(SET));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            kM0_1{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_1(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            kS0_1{MODEL}(kn,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
            kM0_2{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_2(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            kS0_2{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_2(SET));%
            kM0_3{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_3(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            kS0_3{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_3(SET));%
            kM0_5{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_5(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
            kS0_5{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_3(SET));%
            
            
            
            kIth{MODEL}(kn,ro)=mean(DA{MODEL,ro,n,kn}.inf_th(SET));
        end
    end
end
% N0=1./(2.^VEC_Nn);
N0= 1:9;% (2.^VEC_Nn);


mo=0;
for MODEL=[1 3 2 4]
    mo=mo+1;
    for i=2:9
        [p,tbl,stat]=anova1(fliplr(DATk{MODEL,i}(:,1:6))',[],'off')
        h=figure(44)
        %         [c,~,~,gnames] = multcompare(stat,'Alpha',0.05/32);'CType','bonferroni'
%         [c,~,~,gnames] = multcompare(stat,'Alpha',0.001/(32*100));
        [c,~,~,gnames] = multcompare(stat,'Alpha',0.01/(5*9*6));
        ax1 = gca;
        fig1 = get(ax1,'children');
        h=figure(45)
        s1=subplot(4,8,(mo-1)*8+i-1);
        hh=copyobj(fig1,s1);
        line([0 0],[0 7],'color',[0.5 0.5 0.5],'lineStyle','--')
        yticks(1:6)
        yticklabels({'10','30','50','100','150','200'})
        xlabel('Absolute error')
        set(hh,'MarkerSize',10)
        set(s1,'FontSize',10)
    end
end
%
% mo=0;
% for MODEL=[1 3 2 4]
%     mo=mo+1;
%     for i=2:9
%         [p,tbl,stat]=anova1(fliplr(DATkB{MODEL,i}(:,1:6)),[],'off')
%         h=figure(44)
%         [c,~,~,gnames] = multcompare(stat,'Alpha',0.001/32/100);
%         ax1 = gca;
%         fig1 = get(ax1,'children');
%         figure(55)
%         s1=subplot(4,8,(mo-1)*8+i-1);
%         copyobj(fig1,s1);
%         line([0 0],[0 7],'color',[0.5 0.5 0.5],'lineStyle','--')
%         yticks(1:6)
%         yticklabels({'10','30','50','100','150','200'})
%         xlabel('Bias')
%     end
% end
%

figure(11)
mo=0;
j=0;
for MODEL=[1 3 2 4]
    mo=mo+1;
    for i=2:9
        j=j+1;
        subplot(4,8,j)
%         errorbar(kM_1{MODEL}(:,i),1:6,kS_1{MODEL}(:,i)*5,'horizontal','O','MarkerSize',12,'linewidth',2,'CapSize',18)
        errorbar(1:6,kM_1{MODEL}(:,i),kS_1{MODEL}(:,i)*5,'O','MarkerSize',12,'linewidth',2,'CapSize',18)
        box on
        set(gca,'linewidth',2)
        xlim([0 7])
        xticks(1:6)
        xticklabels({'10','30','50','100','150','200'})
        axis square
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=6;
N0=0.1:0.1:0.9;
figure(201)
for MODEL=[1 2 3 4]
    for tt=n
        subplot(3,2,MODEL)
        semilogy(N0,(M_1{MODEL}(tt,:)),'color',[0 0 1],'linewidth',3,'marker','O')
        hold on
        semilogy(N0,(M_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',3,'marker','O')
        hold on
        semilogy(N0,(M_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',3,'marker','O')
        hold on
        % errorbar(N0,mean(M_4{MODEL}(tt,:),1),mean(S_4{MODEL}(tt,:),1),'color',[0 1 0],'linewidth',3,'marker','O')
        hold on
        %         semilogy(N0,(M_5{MODEL}(tt,:)),'color','c','linewidth',3,'marker','O')
        hold on
        % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
        %         set(gca,'YScale','log')
        axis tight
        if MODEL==1 | MODEL==3
            xlabel('\rho')
            xlim([0.2 0.9])
        else
            xlabel('\nu')
        end
        ylabel('Information error')
        if MODEL==1
            legend({'I_{L1}','I_{par}','I_{kNN}'},'Location','northeast','Fontsize',12)
        end
        switch MODEL
            case 1
                title('Gaussian Copula, Gaussian Marginal')
            case 2
                title('student Copula, Gaussian Marginal')
            case 3
                title('Gaussian Copula, Gamma marginal')
            case 4
                title('Student Copula, Gamma Marginal')
        end
        ylim([0 2.4])
        
    end
end



tt=n;
figure(2011)
subplot(3,2,5)
plot(N0,(abs(M0_1{1}(tt,:)-M0_1{3}(tt,:))),'color',[0 0 1],'linewidth',3,'marker','O')
hold on
% plot(N0,(abs(M0_2{1}(tt,:)-M0_2{3}(tt,:))),'color',[1 0 1],'linewidth',3,'marker','O')
hold on
% plot(N0,(abs(M0_3{1}(tt,:)-M0_3{3}(tt,:))),'color',[1 0 0],'linewidth',3,'marker','O')
hold on
% plot((abs(M_4{1}(tt,:)-M_4{3}(tt,:))),'color',[0 1 0],'linewidth',3,'marker','O')
hold on
plot(N0,(abs(M0_5{1}(tt,:)-M0_5{3}(tt,:))),'color','c','linewidth',3,'marker','O')
title('Gaussian Copula')
ylabel('I-I_{theory}','fontsize',12)
xlabel('\rho','fontsize',12)
axis tight
xlim([0.2 0.9])
% set(gca,'YScale','log')
subplot(3,2,6)
plot(N0,(abs(M0_1{2}(tt,:)-M0_1{4}(tt,:))),'color',[0 0 1],'linewidth',3,'marker','O')
hold on
% plot(N0,(abs(M0_2{2}(tt,:)-M0_2{4}(tt,:))),'color',[1 0 1],'linewidth',3,'marker','O')
hold on
% plot(N0,(abs(M0_3{2}(tt,:)-M0_3{4}(tt,:))),'color',[1 0 0],'linewidth',3,'marker','O')
hold on
% plot((abs(M_4{1}(tt,:)-M_4{3}(tt,:))),'color',[0 1 0],'linewidth',3,'marker','O')
hold on
plot(N0,(abs(M0_5{2}(tt,:)-M0_5{4}(tt,:))),'color','c','linewidth',3,'marker','O')
title('Student-t Copula')
ylabel('I-I_{theory}','fontsize',12)
xlabel('\nu','fontsize',12)
axis tight
xlim([0.1 0.9])

% set(gca,'YScale','log')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N0=0.1:0.1:0.9;
figure(501)
ss=0;
for MODEL=[1 3 2 4]
    ss=ss+1;
    for tt=n
        subplot(2,2,ss)
        semilogy(N0,(M_1{MODEL}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','O','linestyle','-','MarkerSize',13)
        hold on
        semilogy(N0,(M_2{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','O','linestyle','--','MarkerSize',13)
        hold on
        %         semilogy(N0,(M_6{MODEL}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','*','linestyle','-','MarkerSize',13)
        hold on
        %         semilogy(N0,(M_7{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','*','linestyle','--','MarkerSize',13)
        hold on
        semilogy(N0,(M_8{MODEL}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','s','linestyle','-','MarkerSize',13)
        hold on
        semilogy(N0,(M_9{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','s','linestyle','--','MarkerSize',13)
        axis tight
        if MODEL==1 | MODEL==3
            xlabel('\rho')
            xlim([0.2 0.9])
        else
            xlabel('\nu')
        end
        ylabel('Information error')
        if MODEL==1
            %             legend({'I_{LL1}^{sample}','I_{LL2}^{sample}','I_{LL1}^{integral}','I_{LL2}^{integral}','I_{LL1}^{naive}','I_{LL2}^{naive}'},'Location','northeast','Fontsize',12)
            legend({'I_{LL1}','I_{LL2}','I_{LL1}^{naive}','I_{LL2}^{naive}'},'Location','northeast','Fontsize',12)
        end
        switch MODEL
            case 1
                title('Gaussian Copula, Gaussian Marginal')
            case 2
                title('student Copula, Gaussian Marginal')
            case 3
                title('Gaussian Copula, Gamma marginal')
            case 4
                title('Student Copula, Gamma Marginal')
        end
        %         ylim([0 0.25])
        axis square
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% naive, LL1 LL2 comparison

N0=0.1:0.1:0.9;
figure(5011)
ss=0;
for MODEL=[1 3 2 4]
    ss=ss+1;
    tt=6;
    subplot(2,2,ss)
    errorbar(N0,M_1{MODEL}(tt,:),S_1{MODEL}(tt,:),'color',[0 0 1],'linewidth',2.5,'marker','O','linestyle','-','MarkerSize',13)
    hold on
    errorbar(N0,M_2{MODEL}(tt,:),S_2{MODEL}(tt,:),'color',[1 0 0],'linewidth',2,'marker','O','linestyle','--','MarkerSize',13)
    hold on
    errorbar(N0,M_8{MODEL}(tt,:),S_8{MODEL}(tt,:),'color',[0 0 1],'linewidth',2.5,'marker','s','linestyle','-','MarkerSize',13)
    hold on
    errorbar(N0,M_9{MODEL}(tt,:),S_9{MODEL}(tt,:),'color',[1 0 0],'linewidth',2,'marker','s','linestyle','--','MarkerSize',13)
    axis tight
    if MODEL==1 | MODEL==3
        xlabel('\rho')
        %             xlim([0.2 0.9])
    else
        xlabel('\nu')
    end
    ylabel('Information error')
    
        set(gca,'linewidth',2)
    
    if MODEL==1
        %             legend({'I_{LL1}^{sample}','I_{LL2}^{sample}','I_{LL1}^{integral}','I_{LL2}^{integral}','I_{LL1}^{naive}','I_{LL2}^{naive}'},'Location','northeast','Fontsize',12)
        legend({'I_{LL1}','I_{LL2}','I_{LL1}^{naive}','I_{LL2}^{naive}'},'Location','northeast','Fontsize',12)
    end
    switch MODEL
        case 1
            title('Gaussian Copula, Gaussian Marginal')
        case 2
            title('student Copula, Gaussian Marginal')
        case 3
            title('Gaussian Copula, Gamma marginal')
        case 4
            title('Student Copula, Gamma Marginal')
    end
    %         ylim([0 0.25])
    xlim([0.2 0.9])
    axis square
end


tt=n;
figure(502)
subplot(1,2,1)
semilogy(N0,abs(M_1{1}(tt,:)-M_1{3}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','O','linestyle','-','MarkerSize',13)
hold on
semilogy(N0,abs(M_2{1}(tt,:)-M_2{3}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','O','linestyle','--','MarkerSize',13)
hold on
semilogy(N0,abs(M_8{1}(tt,:)-M_8{3}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','s','linestyle','-','MarkerSize',13)
hold on
semilogy(N0,abs(M_9{1}(tt,:)-M_9{3}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','s','linestyle','--','MarkerSize',13)
title('Gaussian Copula')
ylabel('I-I_{theory}','fontsize',12)
xlabel('\rho','fontsize',12)
axis tight
xlim([0.2 0.9])
% set(gca,'YScale','log')
subplot(1,2,2)
semilogy(N0,abs(M_1{2}(tt,:)-M_1{4}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','O','linestyle','-','MarkerSize',13)
hold on
semilogy(N0,abs(M_2{2}(tt,:)-M_2{4}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','O','linestyle','--','MarkerSize',13)
hold on
semilogy(N0,abs(M_8{2}(tt,:)-M_8{4}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','s','linestyle','-','MarkerSize',13)
hold on
semilogy(N0,abs(M_9{2}(tt,:)-M_9{4}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','s','linestyle','--','MarkerSize',13)
title('Student-t Copula')
ylabel('I-I_{theory}','fontsize',12)
xlabel('\nu','fontsize',12)
axis tight
xlim([0.1 0.9])

% set(gca,'YScale','log')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% k comparison

N0=0.1:0.1:0.9;
figure(510)
ss=0;
ll=linspace(1,2,6);
for MODEL=[1 3 2 4]
    ss=ss+1;
    for tt=1:6
        subplot(2,2,ss)
        plot(N0,kM_1{MODEL}(tt,:),'color',[0 0 1-(6-tt)/6],'linewidth',ll(tt),'marker','O','linestyle','-','MarkerSize',2*tt+1)
        hold on
        %         errorbar(N0,kM_2{MODEL}(tt,:),S_2{MODEL}(tt,:),'color',[1 0 0],'linewidth',2,'marker','O','linestyle','--','MarkerSize',13)
        hold on
        %         errorbar(N0,kM_8{MODEL}(tt,:),S_8{MODEL}(tt,:),'color',[0 0 1],'linewidth',2.5,'marker','s','linestyle','-','MarkerSize',13)
        hold on
        %         errorbar(N0,kM_9{MODEL}(tt,:),S_9{MODEL}(tt,:),'color',[1 0 0],'linewidth',2,'marker','s','linestyle','--','MarkerSize',13)
        axis tight
        if MODEL==1 | MODEL==3
            xlabel('\rho')
            %             xlim([0.2 0.9])
        else
            xlabel('\nu')
        end
        ylabel('Information error')
        
        %         set(gca,'YScale','log')
        set(gca,'linewidth',2)
        
        if MODEL==1 & tt==6
            %             legend({'I_{LL1}^{sample}','I_{LL2}^{sample}','I_{LL1}^{integral}','I_{LL2}^{integral}','I_{LL1}^{naive}','I_{LL2}^{naive}'},'Location','northeast','Fontsize',12)
            legend({'k=10','k=30','k=50','k=100','k=150','k=200'},'Location','northeast','Fontsize',12)
        end
        switch MODEL
            case 1
                title('Gaussian Copula, Gaussian Marginal')
            case 2
                title('student Copula, Gaussian Marginal')
            case 3
                title('Gaussian Copula, Gamma marginal')
            case 4
                title('Student Copula, Gamma Marginal')
        end
        %         ylim([0 0.25])
        xlim([0.2 0.9])
        axis square
    end
end





N0=0.1:0.1:0.9;
figure(5120)
ss=0;
ll=linspace(1,2,6);
for MODEL=[1 3 2 4]
    ss=ss+1;
    for tt=1:6
        subplot(2,2,ss)
        plot(N0,kM_1{MODEL}(tt,:)./repmat(max(kM_1{MODEL}),1,1),'color',[0 0 1-(6-tt)/6],'linewidth',ll(tt),'marker','O','linestyle','-','MarkerSize',2*tt+1)
        hold on
        %         errorbar(N0,kM_2{MODEL}(tt,:),S_2{MODEL}(tt,:),'color',[1 0 0],'linewidth',2,'marker','O','linestyle','--','MarkerSize',13)
        hold on
        %         errorbar(N0,kM_8{MODEL}(tt,:),S_8{MODEL}(tt,:),'color',[0 0 1],'linewidth',2.5,'marker','s','linestyle','-','MarkerSize',13)
        hold on
        %         errorbar(N0,kM_9{MODEL}(tt,:),S_9{MODEL}(tt,:),'color',[1 0 0],'linewidth',2,'marker','s','linestyle','--','MarkerSize',13)
        axis tight
        if MODEL==1 | MODEL==3
            xlabel('\rho')
            %             xlim([0.2 0.9])
        else
            xlabel('\nu')
        end
        ylabel('Information error')
        
        %         set(gca,'YScale','log')
        
        if MODEL==1 & tt==6
            %             legend({'I_{LL1}^{sample}','I_{LL2}^{sample}','I_{LL1}^{integral}','I_{LL2}^{integral}','I_{LL1}^{naive}','I_{LL2}^{naive}'},'Location','northeast','Fontsize',12)
            legend({'k=10','k=30','k=50','k=100','k=150','k=200'},'Location','northeast','Fontsize',12)
        end
        switch MODEL
            case 1
                title('Gaussian Copula, Gaussian Marginal')
            case 2
                title('student Copula, Gaussian Marginal')
            case 3
                title('Gaussian Copula, Gamma marginal')
            case 4
                title('Student Copula, Gamma Marginal')
        end
        %         ylim([0 0.25])
        xlim([0.2 0.9])
        axis square
    end
end



N0=0.1:0.1:0.9;
figure(511)
ss=0;
ll=linspace(1,2,6);
for MODEL=[1 3 2 4]
    ss=ss+1;
    for tt=1:6
        subplot(2,2,ss)
        line([0.1 0.9],[0 0],'color','k','linewidth',2)
        
        %         plot(N0,nkM_1{MODEL}(tt,:),'color',[0 0 1-(6-tt)/6],'linewidth',ll(tt),'marker','O','linestyle','-','MarkerSize',2*tt+1)
        hold on
        errorbar(N0,nkM_1{MODEL}(tt,:),kS_1{MODEL}(tt,:),'color',[0 0 1-(6-tt)/6],'linewidth',ll(tt),'marker','O','linestyle','-','MarkerSize',2*tt+1)
        hold on
        %         errorbar(N0,kM_8{MODEL}(tt,:),S_8{MODEL}(tt,:),'color',[0 0 1],'linewidth',2.5,'marker','s','linestyle','-','MarkerSize',13)
        hold on
        %         errorbar(N0,kM_9{MODEL}(tt,:),S_9{MODEL}(tt,:),'color',[1 0 0],'linewidth',2,'marker','s','linestyle','--','MarkerSize',13)
        axis tight
        if MODEL==1 | MODEL==3
            xlabel('\rho')
            %             xlim([0.2 0.9])
        else
            xlabel('\nu')
        end
        ylabel('Information bias')
        
        %         set(gca,'YScale','log')
        
        if MODEL==1 & tt==6
            %             legend({'I_{LL1}^{sample}','I_{LL2}^{sample}','I_{LL1}^{integral}','I_{LL2}^{integral}','I_{LL1}^{naive}','I_{LL2}^{naive}'},'Location','northeast','Fontsize',12)
            legend({'k=10','k=30','k=50','k=100','k=150','k=200'},'Location','northeast','Fontsize',12)
        end
        switch MODEL
            case 1
                title('Gaussian Copula, Gaussian Marginal')
            case 2
                title('student Copula, Gaussian Marginal')
            case 3
                title('Gaussian Copula, Gamma marginal')
            case 4
                title('Student Copula, Gamma Marginal')
        end
        %         ylim([0 0.25])
        xlim([0.2 0.9])
        axis square
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% N comparison

N0=0.1:0.1:0.9;
figure(520)
ss=0;
ll=linspace(1,2,6);
for MODEL=[1 3 2 4]
    ss=ss+1;
    for tt=1:6
        subplot(2,2,ss)
        errorbar(N0,kM_1{MODEL}(tt,:),kS_1{MODEL}(tt,:),'color',[0 0 1-(6-tt)/6],'linewidth',ll(tt),'marker','O','linestyle','-','MarkerSize',2*tt+1)
        hold on
        %         errorbar(N0,kM_2{MODEL}(tt,:),S_2{MODEL}(tt,:),'color',[1 0 0],'linewidth',2,'marker','O','linestyle','--','MarkerSize',13)
        hold on
        %         errorbar(N0,kM_8{MODEL}(tt,:),S_8{MODEL}(tt,:),'color',[0 0 1],'linewidth',2.5,'marker','s','linestyle','-','MarkerSize',13)
        hold on
        %         errorbar(N0,kM_9{MODEL}(tt,:),S_9{MODEL}(tt,:),'color',[1 0 0],'linewidth',2,'marker','s','linestyle','--','MarkerSize',13)
        axis tight
        if MODEL==1 | MODEL==3
            xlabel('\rho')
            %             xlim([0.2 0.9])
        else
            xlabel('\nu')
        end
        ylabel('Information error')
        
        %         set(gca,'YScale','log')
                set(gca,'linewidth',2)

        if MODEL==1 & tt==6
            %             legend({'I_{LL1}^{sample}','I_{LL2}^{sample}','I_{LL1}^{integral}','I_{LL2}^{integral}','I_{LL1}^{naive}','I_{LL2}^{naive}'},'Location','northeast','Fontsize',12)
            legend({'k=10','k=30','k=50','k=100','k=150','k=200'},'Location','northeast','Fontsize',12)
        end
        switch MODEL
            case 1
                title('Gaussian Copula, Gaussian Marginal')
            case 2
                title('student Copula, Gaussian Marginal')
            case 3
                title('Gaussian Copula, Gamma marginal')
            case 4
                title('Student Copula, Gamma Marginal')
        end
        %         ylim([0 0.25])
        xlim([0.2 0.9])
        axis square
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% LNC comparison

N0=0.1:0.1:0.9;
figure(205)
ss=0;
for MODEL=[1 3 2 4]
    ss=ss+1;
    for tt=6
        subplot(2,2,ss)
        plot(N0,(M_1{MODEL}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',12)
        %         hold on
        %         errorbar(N0,(M_7{MODEL}(tt,:)),(S_7{MODEL}(tt,:)),'color','g','linewidth',2,'marker','O','MarkerSize',10)
        hold on
        %         errorbar(N0,(M_2{MODEL}(tt,:)),(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
        hold on
        plot(N0,(M_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',12)
        hold on
        plot(N0,(M_5{MODEL}(tt,:)),'color','g','linewidth',2,'marker','s','MarkerSize',12)
        hold on
        %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
        hold on
        %         plot(N0,Ith{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
        % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
        %                         set(gca,'YScale','log')    %%%% this should be uncommented for
        %         paper figure (absolute error)
        if MODEL==1 | MODEL==3
            xlabel('\rho','Fontsize',16)
        else
            xlabel('\nu','Fontsize',16)
        end
        ylabel('Information absolute error','Fontsize',16)
        axis tight
        if MODEL==1
            legend({'I_{LL1}','I_{par}','I_{LNC}'},'Location','northeast','Fontsize',12)
        end
        axis tight
        xlim([0.2 0.9])
        axis square

        set(gca,'linewidth',2)

        switch MODEL
            case 1
                title('Gaussian Copula, Gaussian Marginal')
            case 2
                title('student Copula, Gaussian Marginal')
            case 3
                title('Gaussian Copula, Gamma marginal')
            case 4
                title('Student Copula, Gamma Marginal')
        end
    end
end

figure(206)
ss=0;
for MODEL=[1 3 2 4]
    ss=ss+1;
    for tt=6
        subplot(2,2,ss)
        line([0.1 0.9],[0 0],'color','k','linewidth',2)
        hold on
        errorbar(N0,(nM_1{MODEL}(tt,:)),(nS_1{MODEL}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',12)
        %         hold on
        %         errorbar(N0,(M_7{MODEL}(tt,:)),(S_7{MODEL}(tt,:)),'color','g','linewidth',2,'marker','O','MarkerSize',10)
        hold on
        %         errorbar(N0,(M_2{MODEL}(tt,:)),(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
        hold on
        errorbar(N0,(nM_3{MODEL}(tt,:)),(nS_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',12)
        hold on
        errorbar(N0,(nM_5{MODEL}(tt,:)),(nS_5{MODEL}(tt,:)),'color','g','linewidth',2,'marker','s','MarkerSize',12)
        hold on
        %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
        hold on
        %         plot(N0,Ith{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
        % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
        %         set(gca,'YScale','log')    %%%% this should be uncommented for
        %         paper figure (absolute error)
        if MODEL==1 | MODEL==3
            xlabel('\rho','Fontsize',16)
        else
            xlabel('\nu','Fontsize',16)
        end
        ylabel('Information bias (bits)','Fontsize',16)
        if MODEL==1
            legend({'I_{LL1}','I_{par}','I_{LNC}'},'Location','northeast','Fontsize',12)
        end
        axis tight
        xlim([0.2 0.9])
        axis square
        box on
            set(gca,'linewidth',2)

        switch MODEL
            case 1
                title('Gaussian Copula, Gaussian Marginal')
            case 2
                title('student Copula, Gaussian Marginal')
            case 3
                title('Gaussian Copula, Gamma marginal')
            case 4
                title('Student Copula, Gamma Marginal')
        end
        
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% LNC difference plot
N0=0.1:0.1:0.9;
tt=6;
figure(270)
subplot(2,2,1)
plot(N0,(abs(M0_1{1}(tt,:)-M0_1{3}(tt,:))),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',12)
hold on
% plot(N0,(abs(M0_2{1}(tt,:)-M0_2{3}(tt,:))),'color',[1 0 1],'linewidth',3,'marker','O')
hold on
plot(N0,(abs(M0_3{1}(tt,:)-M0_3{3}(tt,:))),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',12)
hold on
% plot((abs(M_4{1}(tt,:)-M_4{3}(tt,:))),'color',[0 1 0],'linewidth',3,'marker','O')
hold on
plot(N0,(abs(M0_5{1}(tt,:)-M0_5{3}(tt,:))),'color','g','linewidth',2,'marker','s','MarkerSize',12)
title('Gaussian Copula')
ylabel('I-I_{theory}','fontsize',12)
xlabel('\rho','fontsize',12)
axis tight
set(gca,'linewidth',2)
xlim([0.2 0.9])
% set(gca,'YScale','log')
legend({'I_{LL1}','I_{GC}','I_{LNC}'},'Location','northwest','Fontsize',12)
axis square
subplot(2,2,2)
plot(N0,(abs(M0_1{2}(tt,:)-M0_1{4}(tt,:))),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',12)
hold on
% plot(N0,(abs(M0_2{2}(tt,:)-M0_2{4}(tt,:))),'color',[1 0 1],'linewidth',3,'marker','O')
hold on
plot(N0,(abs(M0_3{2}(tt,:)-M0_3{4}(tt,:))),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',12)
hold on
% plot((abs(M_4{1}(tt,:)-M_4{3}(tt,:))),'color',[0 1 0],'linewidth',3,'marker','O')
hold on
plot(N0,(abs(M0_5{2}(tt,:)-M0_5{4}(tt,:))),'color','g','linewidth',2,'marker','s','MarkerSize',12)
title('Student-t Copula')
ylabel('I-I_{theory}','fontsize',12)
xlabel('\nu','fontsize',12)
axis tight
set(gca,'linewidth',2)
xlim([0.2 0.9])
% set(gca,'YScale','log')
axis square



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% LNC comparison NUMBER  N

N0=[5:13];
figure(215)
ss=0;
for MODEL=[1 3 2 4]
    ss=ss+1;
    tt=1:9;
    if MODEL==1 | MODEL==3
        ttt=5;
    else
        ttt=5;
    end
    
    subplot(2,2,ss)
    
    plot(N0,(M_1{MODEL}(tt,ttt)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',12)
    %         hold on
    %         errorbar(N0,(M_7{MODEL}(tt,ttt)),(S_7{MODEL}(tt,ttt)),'color','g','linewidth',2,'marker','O')
    hold on
    plot(N0,(M_3{MODEL}(tt,ttt)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',12)
    hold on
    plot(N0,(M_5{MODEL}(tt,ttt)),'color','g','linewidth',2,'marker','s','MarkerSize',12)
    hold on
    %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
    hold on
    % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
    % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
    %             set(gca,'YScale','log')
    if MODEL==1 | MODEL==3
        xlabel('N','Fontsize',16)
    else
        xlabel('N','Fontsize',16)
    end
    ylabel('Information absolute error','Fontsize',16)
    axis tight

    set(gca,'linewidth',2)

    if MODEL==1
        legend({'I_{LL1}','I_{par}','I_{LNC}'},'Location','northeast','Fontsize',12)
    end
    switch MODEL
        case 1
            title('Gaussian Copula, Gaussian Marginal')
        case 2
            title('student Copula, Gaussian Marginal')
        case 3
            title('Gaussian Copula, Gamma marginal')
        case 4
            title('Student Copula, Gamma Marginal')
    end
    %     xlim([0.2 0.9])
    axis square
end

figure(2166)
ss=0;
for MODEL=[1 3 2 4]
    ss=ss+1;
    subplot(2,2,ss)
    
    errorbar(N0,(nM_1{MODEL}(tt,ttt)),(nS_1{MODEL}(tt,ttt)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',12)
    %         hold on
    %         errorbar(N0,(M_7{MODEL}(tt,ttt)),(S_7{MODEL}(tt,ttt)),'color','g','linewidth',2,'marker','O')
    hold on
    errorbar(N0,(nM_3{MODEL}(tt,ttt)),(nS_3{MODEL}(tt,ttt)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',12)
    hold on
    errorbar(N0,(nM_5{MODEL}(tt,ttt)),(nS_5{MODEL}(tt,ttt)),'color','g','linewidth',2,'marker','s','MarkerSize',12)
    hold on
    %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
    hold on
    
    line([5 13],[0 0],'color','k','linewidth',2)
    
    % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
    % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
    %             set(gca,'YScale','log')
    if MODEL==1 | MODEL==3
        xlabel('N','Fontsize',16)
    else
        xlabel('N','Fontsize',16)
    end
    
    set(gca,'linewidth',2)

    ylabel('Information bias ratio','Fontsize',16)
    axis tight
    if MODEL==1
        legend({'I_{LL1}','I_{par}','I_{LNC}'},'Location','northeast','Fontsize',12)
    end
    switch MODEL
        case 1
            title('Gaussian Copula, Gaussian Marginal')
        case 2
            title('student Copula, Gaussian Marginal')
        case 3
            title('Gaussian Copula, Gamma marginal')
        case 4
            title('Student Copula, Gamma Marginal')
    end
    %     xlim([0.2 0.9])
    axis square
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% PYM comparison

N0=[0.1:0.1:0.9]*100;
figure(305)
ss=0;
for MODEL=[7 5]
    ss=ss+1;
    for tt=6
        subplot(1,2,ss)
        errorbar(N0,(M_1{MODEL}(tt,:)),(S_1{MODEL}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',10)
        hold on
        %         errorbar(N0,(M_2{MODEL}(tt,:)),(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
        hold on
        errorbar(N0,(M_3{MODEL}(tt,:)),(S_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',10)
        hold on
        errorbar(N0,(M_5{MODEL}(tt,:)),(S_5{MODEL}(tt,:)),'color','c','linewidth',2,'marker','s','MarkerSize',10)
        hold on
        %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
        hold on
        % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
        % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
        %         set(gca,'YScale','log')
        if MODEL==1 | MODEL==3
            xlabel('\rho','Fontsize',16)
        else
            xlabel('\nu','Fontsize',16)
        end
        ylabel('Information error','Fontsize',16)
        axis tight
        if MODEL==1
            legend({'I_{LL1}','I_{par}','I_{LNC}'},'Location','northeast','Fontsize',12)
        end
    end
    
    set(gca,'linewidth',2)

    switch MODEL
        case 1
            title('Gaussian Copula, Gaussian Marginal')
        case 2
            title('student Copula, Gaussian Marginal')
        case 3
            title('Gaussian Copula, Gamma marginal')
        case 4
            title('Student Copula, Gamma Marginal')
    end
    xlim([10 70])
    axis square
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% PYM k comparison

N0=[0.1:0.1:0.9]*100;
figure(310)
ss=0;
ll=linspace(1,2,6);
for MODEL=[7 5]
    ss=ss+1;
    for tt=1:6
        for sss=1:2
            if sss==1
                subplot(2,2,ss)
                plot(N0,kM_1{MODEL}(tt,:),'color',[0 0 1-(6-tt)/6],'linewidth',ll(tt),'marker','O','linestyle','-','MarkerSize',2*tt+1)
                hold on
                if tt==6
                    hold on
                    plot(N0,(M_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',12)
                    hold on
                    plot(N0,kM_5{MODEL}(tt,:),'color','g','linewidth',2,'marker','O','MarkerSize',12)
                end
                
                xlabel('\lambda')
                ylabel('Information error')
                
                %         set(gca,'YScale','log')
                
                set(gca,'linewidth',2)

                if MODEL==7 & tt==6
                    %             legend({'I_{LL1}^{sample}','I_{LL2}^{sample}','I_{LL1}^{integral}','I_{LL2}^{integral}','I_{LL1}^{naive}','I_{LL2}^{naive}'},'Location','northeast','Fontsize',12)
                    legend({'k=10','k=30','k=50','k=100','k=150','k=200','GC','PYM'},'Location','northeast','Fontsize',12)
                end
                switch MODEL
                    case 7
                        title('Gaussian Copula, Poisson Marginal')
                    case 5
                        title('student Copula, Poisson Marginal')
                end
                %         ylim([0 0.25])
                xlim([20 70])
                axis square
            else
                
                subplot(2,2,ss+2)
                errorbar(N0,nkM_1{MODEL}(tt,:),nkS_1{MODEL}(tt,:),'color',[0 0 1-(6-tt)/6],'linewidth',ll(tt),'marker','O','linestyle','-','MarkerSize',2*tt+1)
                hold on
                if tt==6
                    hold on
                    errorbar(N0,(nM_3{MODEL}(tt,:)),(nS_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',12)
                    hold on
                    errorbar(N0,nM_5{MODEL}(tt,:),nS_5{MODEL}(tt,:),'color','g','linewidth',2,'marker','O','MarkerSize',12)
                end
                
                xlabel('\lambda')
                ylabel('Information error')
                
                %         set(gca,'YScale','log')
                switch MODEL
                    case 7
                        title('Gaussian Copula, Poisson Marginal')
                    case 5
                        title('student Copula, Poisson Marginal')
                end
                xlim([20 70])
                axis square
                
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% PYM comparison NUMBER  N

N0=[5:13];
figure(3115)
ss=0;
for MODEL=[7 5]
    ss=ss+1;
    tt=1:9;
    ttt=5;
    
    for sss=1:2
        if sss==1
            subplot(2,2,ss)
            
            plot(N0,(M_1{MODEL}(tt,ttt)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',12)
            hold on
            %         errorbar(N0,(M_2{MODEL}(tt,:)),(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
            hold on
            plot(N0,(M_3{MODEL}(tt,ttt)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',12)
            hold on
            plot(N0,(M_5{MODEL}(tt,ttt)),'color','g','linewidth',2,'marker','s','MarkerSize',12)
            hold on
            %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
            hold on
            % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
            % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
            %             set(gca,'YScale','log')
            if MODEL==1 | MODEL==3
                xlabel('N','Fontsize',16)
            else
                xlabel('N','Fontsize',16)
            end
            ylabel('Information error','Fontsize',16)
            axis tight
            if MODEL==1
                legend({'I_{LL1}','I_{par}','I_{LNC}'},'Location','northeast','Fontsize',12)
            end
            switch MODEL
                case 1
                    title('Gaussian Copula, Gaussian Marginal')
                case 2
                    title('student Copula, Gaussian Marginal')
                case 3
                    title('Gaussian Copula, Gamma marginal')
                case 4
                    title('Student Copula, Gamma Marginal')
            end
            %     xlim([0.2 0.9])
            axis square
            set(gca,'linewidth',2)
            
        else
            subplot(2,2,ss+2)
            line([5 13],[0 0],'color','k','linewidth',2)
            hold on
            
            errorbar(N0,(nM_1{MODEL}(tt,ttt)),(nS_1{MODEL}(tt,ttt)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',12)
            hold on
            %         errorbar(N0,(M_2{MODEL}(tt,:)),(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
            hold on
            errorbar(N0,(nM_3{MODEL}(tt,ttt)),(nS_3{MODEL}(tt,ttt)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',12)
            hold on
            errorbar(N0,(nM_5{MODEL}(tt,ttt)),(nS_5{MODEL}(tt,ttt)),'color','g','linewidth',2,'marker','s','MarkerSize',12)
            hold on
            %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
            hold on
            % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
            % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
            %             set(gca,'YScale','log')
            if MODEL==1 | MODEL==3
                xlabel('N','Fontsize',16)
            else
                xlabel('N','Fontsize',16)
            end
            ylabel('Information error','Fontsize',16)
            axis tight

            box on
            set(gca,'linewidth',2)

            if MODEL==1
                legend({'I_{LL1}','I_{par}','I_{LNC}'},'Location','northeast','Fontsize',12)
            end
            switch MODEL
                case 1
                    title('Gaussian Copula, Gaussian Marginal')
                case 2
                    title('student Copula, Gaussian Marginal')
                case 3
                    title('Gaussian Copula, Gamma marginal')
                case 4
                    title('Student Copula, Gamma Marginal')
            end
            %     xlim([0.2 0.9])
            axis square
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% PYM variance-N plot

N0=5:13;
figure(271)
subplot(2,2,1)
plot(N0,nanstd(M0_1{7}(:,2:7)'),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',12)
hold on
% plot(N0,(abs(M0_2{1}(tt,:)-M0_2{3}(tt,:))),'color',[1 0 1],'linewidth',3,'marker','O')
hold on
plot(N0,nanstd(M0_3{7}(:,2:7)'),'color',[1 0 0],'linewidth',3,'marker','d','MarkerSize',12)
hold on
% plot((abs(M_4{1}(tt,:)-M_4{3}(tt,:))),'color',[0 1 0],'linewidth',3,'marker','O')
hold on
plot(N0,nanstd(M0_5{7}(:,2:7)'),'color','g','linewidth',2,'marker','s','MarkerSize',12)
title('Gaussian Copula')
ylabel('\sigma_{\lambda}(I)','fontsize',12)
xlabel('N','fontsize',12)
axis tight
set(gca,'linewidth',2)
% xlim([0.2 0.9])
% set(gca,'YScale','log')
legend({'I_{LL1}','I_{GC}','I_{LNC}'},'Location','northwest','Fontsize',12)
axis square
subplot(2,2,2)
plot(N0,nanstd(M0_1{5}(:,2:7)'),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',12)
hold on
% plot(N0,(abs(M0_2{2}(tt,:)-M0_2{4}(tt,:))),'color',[1 0 1],'linewidth',3,'marker','O')
hold on
plot(N0,nanstd(M0_3{5}(:,2:7)'),'color',[1 0 0],'linewidth',3,'marker','d','MarkerSize',12)
hold on
% plot((abs(M_4{1}(tt,:)-M_4{3}(tt,:))),'color',[0 1 0],'linewidth',3,'marker','O')
hold on
plot(N0,nanstd(M0_5{5}(:,2:7)'),'color','g','linewidth',2,'marker','s','MarkerSize',12)
title('Student-t Copula')
ylabel('\sigma_{\lambda}(I)','fontsize',12)
xlabel('N','fontsize',12)
axis tight
% xlim([0.2 0.9])
% set(gca,'YScale','log')
axis square          
set(gca,'linewidth',2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% PYM comparison NUMBER  N

N0=[5:13];
figure(2155)
ss=0;
for MODEL=[1 3 2 4]
    ss=ss+1;
    tt=1:9;
    if MODEL==1 | MODEL==3
        ttt=5;
    else
        ttt=5;
    end
    subplot(2,2,ss)
    
    plot(N0,(M_1{MODEL}(tt,ttt)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',10)
    hold on
    %         errorbar(N0,(M_2{MODEL}(tt,:)),(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
    hold on
    plot(N0,(M_3{MODEL}(tt,ttt)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',10)
    hold on
    plot(N0,(M_5{MODEL}(tt,ttt)),'color','c','linewidth',2,'marker','s','MarkerSize',10)
    hold on
    %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
    hold on
    % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
    % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
%     set(gca,'YScale','log')
    if MODEL==1 | MODEL==3
        xlabel('N','Fontsize',16)
    else
        xlabel('N','Fontsize',16)
    end
    ylabel('Information error','Fontsize',16)
    axis tight

    set(gca,'linewidth',2)

    if MODEL==1
        legend({'I_{LL1}','I_{par}','I_{LNC}'},'Location','northeast','Fontsize',12)
    end
    switch MODEL
        case 1
            title('Gaussian Copula, Gaussian Marginal')
        case 2
            title('student Copula, Gaussian Marginal')
        case 3
            title('Gaussian Copula, Gamma marginal')
        case 4
            title('Student Copula, Gamma Marginal')
    end
    %     xlim([0.2 0.9])
    axis square
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=6;
N0=0.1:0.1:0.9;
figure(800)
for MODEL=[5 6 7]
    for tt=n
        subplot(2,2,MODEL-4)
        
        errorbar(N0,(M_1{MODEL}(tt,:)),(S_1{MODEL}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','O')
        hold on
        errorbar(N0,(M_2{MODEL}(tt,:)),(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
        hold on
        errorbar(N0,(M_3{MODEL}(tt,:)),(S_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','O')
        hold on
        errorbar(N0,(M_5{MODEL}(tt,:)),(S_5{MODEL}(tt,:)),'color','c','linewidth',2,'marker','O')
        hold on
        errorbar(N0,(M_8{MODEL}(tt,:)),(S_6{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
        hold on
        % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
        % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
        %         set(gca,'YScale','log')
        if MODEL==1 | MODEL==3
            xlabel('\rho')
        else
            xlabel('\nu')
        end
        ylabel('Information error')
        axis tight
        if MODEL==1
            legend({'I_{L1}','I_{L2}','I_{par}','I_{kNN}'},'Location','northeast','Fontsize',12)
        end
    end
    switch MODEL
        case 5
            title('Gaussian Copula, Poisson Marginal')
        case 6
            title('student Copula, Poisson Marginal')
        case 7
            title('Gaussian Copula, Bernouli marginal')
        case 8
            title('Student Copula, Bernouli Marginal')
    end
end






N0=2.^[5 6 7 8 9 10 11 12 13];
figure(900)
ii=0;
for MODEL=[1 2 3 4]
    for tt=1:9
        %         subplot(2,2,MODEL-4)
        ii=ii+1;
        subplot(4,9,ii)
        errorbar(N0,(M_1{MODEL}(:,tt)),(S_1{MODEL}(:,tt)),'color',[0 0 1],'linewidth',2,'marker','O')
        hold on
        errorbar(N0,(M_2{MODEL}(:,tt)),(S_2{MODEL}(:,tt)),'color',[1 0 1],'linewidth',2,'marker','O')
        hold on
        %         errorbar(N0,(M_3{MODEL}(:,tt)),(S_3{MODEL}(:,tt)),'color',[1 0 0],'linewidth',2,'marker','O')
        hold on
        errorbar(N0,(M_5{MODEL}(:,tt)),(S_5{MODEL}(:,tt)),'color','c','linewidth',2,'marker','O')
        hold on
        %         errorbar(N0,(M_8{MODEL}(:,tt)),(S_8{MODEL}(:,tt)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
        hold on
        % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
        % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
        if MODEL==1 || MODEL==3
            xlabel('\rho')
        else
            xlabel('\nu')
        end
        ylabel('Information error')
        axis tight
        ylim([0 4])
        set(gca,'YScale','log')
        if MODEL==5 & tt==1
            legend({'I_{L1}','I_{L2}','I_{kNN}'},'Location','northeast','Fontsize',12)
        end
    end
    switch MODEL
        case 5
            title('Gaussian Copula, Poisson Marginal')
        case 6
            title('student Copula, Poisson Marginal')
        case 7
            title('Gaussian Copula, Bernouli marginal')
        case 8
            title('Student Copula, Bernouli Marginal')
    end
end




N0=2.^[5 6 7 8 9 10 11];
figure(910)
ii=0;
for MODEL=[5 6]
    for tt=1:9
        %         subplot(2,2,MODEL-4)
        ii=ii+1;
        subplot(2,9,ii)
        errorbar(N0,(M_1{MODEL}(:,tt)),(S_1{MODEL}(:,tt)),'color',[0 0 1],'linewidth',2,'marker','O')
        hold on
        errorbar(N0,(M_2{MODEL}(:,tt)),(S_2{MODEL}(:,tt)),'color',[1 0 1],'linewidth',2,'marker','O')
        hold on
        %         errorbar(N0,(M_3{MODEL}(:,tt)),(S_3{MODEL}(:,tt)),'color',[1 0 0],'linewidth',2,'marker','O')
        hold on
        errorbar(N0,(M_5{MODEL}(:,tt)),(S_5{MODEL}(:,tt)),'color','c','linewidth',2,'marker','O')
        hold on
        %         errorbar(N0,(M_8{MODEL}(:,tt)),(S_8{MODEL}(:,tt)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
        hold on
        % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
        % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
        if MODEL==1 | MODEL==3
            xlabel('\rho')
        else
            xlabel('\nu')
        end
        ylabel('Information error')
        axis tight
        ylim([0 4])
        set(gca,'YScale','log')
        if MODEL==5 & tt==1
            legend({'I_{L1}','I_{L2}','I_{kNN}'},'Location','northeast','Fontsize',12)
        end
    end
    switch MODEL
        case 5
            title('Gaussian Copula, Poisson Marginal')
        case 6
            title('student Copula, Poisson Marginal')
        case 7
            title('Gaussian Copula, Bernouli marginal')
        case 8
            title('Student Copula, Bernouli Marginal')
    end
end





N0=0.1:0.1:0.9;
figure(202)
for MODEL=[1 2 3 4]
    for tt=n
        subplot(2,2,MODEL)
        semilogy(N0,(M_1{MODEL}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','O')
        hold on
        semilogy(N0,(M_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
        hold on
        semilogy(N0,(M_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','O')
        hold on
        semilogy(N0,(M_5{MODEL}(tt,:)),'color','c','linewidth',2,'marker','O')
        hold on
        semilogy(N0,(M_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
        hold on
        % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
        % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
        % set(gca,'YScale','log')
        if MODEL==1 | MODEL==3
            xlabel('\rho')
        else
            xlabel('\nu')
        end
        ylabel('Information error')
        axis tight
        if MODEL==1
            legend({'I_{L1}','I_{L2}','I_{par}','I_{kNN}'},'Location','northeast','Fontsize',12)
        end
    end
    switch MODEL
        case 1
            title('Gaussian Copula, Gaussian Marginal')
        case 2
            title('student Copula, Gaussian Marginal')
        case 3
            title('Gaussian Copula, Gamma marginal')
        case 4
            title('Student Copula, Gamma Marginal')
    end
end




N0=0.1:0.1:0.9;
figure(300)
for MODEL=1:4
    for tt=n
        subplot(2,2,MODEL)
        plot(N0,(S_1{MODEL}(tt,:)),'color',[0 0 1],'linewidth',3,'marker','O')
        hold on
        plot(N0,(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',3,'marker','O')
        hold on
        plot(N0,(S_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',3,'marker','O')
        hold on
        % plot(N0,mean(S_4{MODEL}(tt,:),1),'color',[0 1 0],'linewidth',3,'marker','O')
        hold on
        plot(N0,(S_5{MODEL}(tt,:)),'color','c','linewidth',3,'marker','O')
        hold on
        % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
        % set(gca,'YScale','log')
        xlabel('\rho')
        ylabel('Information error')
        axis tight
        if MODEL==1
            legend({'I_{L1}','I_{L2}','I_{par}','I_{kNN}'},'Location','northeast','Fontsize',12)
        end
    end
end



%%%%%%%%%%%%% variance



N0=[5:13];
figure(2162)
ss=0;
for MODEL=[1 3 2 4]
    ss=ss+1;
    tt=1:9;
    if MODEL==1 | MODEL==3
        ttt=5;
    else
        ttt=5;
    end
    subplot(2,2,ss)
    
    plot(N0,(S_1{MODEL}(tt,ttt)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',10)
    hold on
    %         errorbar(N0,(M_2{MODEL}(tt,:)),(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
    hold on
    clear all
    close all
    
    VEC_Nn=5:13;
    VEC_rn=1:10000;
    VEC_kn=[10 30 50 100 150 200];
    
    % for MODELn=[1 2 3 4]
    %     for nn=1:numel(VEC_Nn)
    %         for knotS=1:numel(VEC_kn)
    %             DA{MODELn}.er{nn,knotS}=nan(900,7);
    %             DA{MODELn}.inf_1{nn,knotS}=nan(900,1);
    %             DA{MODELn}.inf_2{nn,knotS}=nan(900,1);
    %             DA{MODELn}.inf_3{nn,knotS}=nan(900,1);
    %             DA{MODELn}.inf_4{nn,knotS}=nan(900,1);
    %             DA{MODELn}.inf_5{nn,knotS}=nan(900,1);
    %             DA{MODELn}.inf_th{nn,knotS}=nan(900,1);
    %         end
    %     end
    % end
    
    for MODELnnn=1:8
        for nnnn=1:numel(VEC_Nn)
            [MODELnnn nnnn]
            for knotS=1:numel(VEC_kn)
                rePnn=zeros(9,1);
                for rePn=1:5000
                    N0=2^(VEC_Nn(nnnn));
                    if 1==1%rePnn(mod(rePn-1,9)+1)<201%exist(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/nNew_info_MOD_',num2str(MODELnn),'_N_',num2str(N0),'_rep_',num2str(VEC_rn(rePn)),'_knots_',num2str(VEC_kn(knotS)),'.mat'],'file')~=0
                        
                        if MODELnnn==8
                            MODELnn=7;
                        else
                            MODELnn=MODELnnn;
                        end
                        
                        if 1==2%try
                            if MODELnnn<5
                                clear info
                                load(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/New_info_MOD_',num2str(MODELnnn),'_N_',num2str(N0),'_rep_',num2str(VEC_rn(rePn)),'_knots_',num2str(VEC_kn(knotS)),'.mat'],'info','I_theoretical','ro')
                                
                                %                     load(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/info_MOD_',num2str(MODELnn),'_N_',num2str(N0),'_rep_',num2str(VEC_rn(rePn)),'_knots_',num2str(VEC_kn(knotS)),'.mat'])
                                if rePn<901%mod(VEC_rn(rePn)-1,9)+1==8%REP
                                    rePnn(mod(rePn-1,9)+1)=rePnn(mod(rePn-1,9)+1)+1;
                                    %                         ro=0.1*(mod(rePn-1,9)+1)+0.1*rand(1)-0.05;
                                    
                                    DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.er(rePnn(mod(rePn-1,9)+1),:)=abs(info-I_theoretical)/I_theoretical;
                                    DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_1(rePnn(mod(rePn-1,9)+1))=max([-10 info(1)]);
                                    DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_2(rePnn(mod(rePn-1,9)+1))=max([-10 info(1)]);
                                    DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_3(rePnn(mod(rePn-1,9)+1))=max([-10 info(3)]);
                                    
                                    if MODELnn<5
                                        DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_4(rePnn(mod(rePn-1,9)+1))=max([-10 info(6)]);    %%%%%%%% k=5 or 3
                                    end
                                    
                                    DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_5(rePnn(mod(rePn-1,9)+1))=max([-10 info(4)]);
                                    
                                    DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_6(rePnn(mod(rePn-1,9)+1))=max([-10 info(8)]);
                                    DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_7(rePnn(mod(rePn-1,9)+1))=max([-10 info(9)]);
                                    DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_8(rePnn(mod(rePn-1,9)+1))=max([-10 info(10)]);
                                    DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_9(rePnn(mod(rePn-1,9)+1))=max([-10 info(11)]);
                                    
                                    DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_10(rePnn(mod(rePn-1,9)+1))=max([-10 info(5)]);
                                    
                                    if MODELnn<5
                                        if (MODELnn==1 | MODELnn==3)
                                            DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=(-0.5 * log2(1-ro^2));
                                        else
                                            if rePn<10 & nnnn==6 & knotS==5
                                                %                             inth(rePn)=-mixedvineentropy(vineP,0.05,5e-3,10000,1);
                                            end
                                            DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=Omega_Student(ro,0,2);%inth(mod(rePn-1,9)+1);%
                                        end
                                    else
                                        DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=I_theoretical;
                                    end
                                    
                                end
                            end
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        try
                            clear info
                            load(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/nNew_info_MOD_',num2str(MODELnnn),'_N_',num2str(N0),'_rep_',num2str(VEC_rn(rePn)),'_knots_',num2str(VEC_kn(knotS)),'.mat'],'info','I_theoretical','ro')
                            
                            %                     load(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/info_MOD_',num2str(MODELnn),'_N_',num2str(N0),'_rep_',num2str(VEC_rn(rePn)),'_knots_',num2str(VEC_kn(knotS)),'.mat'])
                            if 1==1%rePn<901%mod(VEC_rn(rePn)-1,9)+1==8%REP
                                rePnn(mod(rePn-1,9)+1)=rePnn(mod(rePn-1,9)+1)+1;
                                %                         ro=0.1*(mod(rePn-1,9)+1)+0.1*rand(1)-0.05;
                                
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.er(rePnn(mod(rePn-1,9)+1),:)=abs(info-I_theoretical)/I_theoretical;
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_1(rePnn(mod(rePn-1,9)+1))=max([-10 info(1)]);
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_2(rePnn(mod(rePn-1,9)+1))=max([-10 info(1)]);
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_3(rePnn(mod(rePn-1,9)+1))=max([-10 info(3)]);
                                
                                if MODELnn<5
                                    DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_4(rePnn(mod(rePn-1,9)+1))=max([-10 info(6)]);    %%%%%%%% k=5 or 3
                                end
                                
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_5(rePnn(mod(rePn-1,9)+1))=max([-10 info(4)]);
                                
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_6(rePnn(mod(rePn-1,9)+1))=max([-10 info(8)]);   %%%
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_7(rePnn(mod(rePn-1,9)+1))=max([-10 info(9)]);
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_8(rePnn(mod(rePn-1,9)+1))=max([-10 info(10)]);
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_9(rePnn(mod(rePn-1,9)+1))=max([-10 info(11)]);
                                
                                DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_10(rePnn(mod(rePn-1,9)+1))=max([-10 info(5)]);
                                
                                if MODELnn<5
                                    if (MODELnn==1 | MODELnn==3)
                                        DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=(-0.5 * log2(1-ro^2));
                                    else
                                        if rePn<10 & nnnn==6 & knotS==5
                                            %                             inth(rePn)=-mixedvineentropy(vineP,0.05,5e-3,10000,1);
                                        end
                                        DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=Omega_Student(ro,0,2);%inth(mod(rePn-1,9)+1);%
                                    end
                                else
                                    DA{MODELnn,mod(rePn-1,9)+1,nnnn,knotS}.inf_th(rePnn(mod(rePn-1,9)+1))=I_theoretical;
                                end
                                
                            end
                        end
                        
                    end
                end
            end
        end
    end
    
    
    % load('/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/nDA.mat')
    
    load('/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/DA.mat')
    
    
    
    
    
    
    clear NM
    for rep=1:9
        for MODELn=1:7
            for n=6
                NM(MODELn,rep)=sum(DA{MODELn,rep,n,4}.inf_3~=0)
            end
        end
    end
    
    
    clear MM
    for rep=1:9
        for MODELn=1:7
            for n=6
                MM(MODELn,rep)=sum(DA{MODELn,rep,n,6}.inf_3~=0)
            end
        end
    end
    
    
    
    
    
    SET=1:200;
    kn=4;
    clear M_1 S_1 M_2 S_2 M_3 S_3 M_4 S_4 M_5 S_5 Ith M_6 S_6  M_8 S_8 M_9 S_9
    for MODEL=1:7
        for ro=1:9
            for n=  1:9%5:13%1:7
                
                %             M_1{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
                %             S_1{MODEL}(n,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_1(SET)./DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                %             M_2{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
                %             S_2{MODEL}(n,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_2(SET)./DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                %             M_3{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
                %             S_3{MODEL}(n,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_3(SET)./DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                %             %                 M_4{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_4-DA{MODEL,ro,n,kn}.inf_th)./DA{MODEL,ro,n,kn}.inf_th);
                %             %                 S_4{MODEL}(n,ro)=std(abs(DA{MODEL,ro,n,kn}.inf_4./DA{MODEL,ro,n,kn}.inf_th));
                %             M_5{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
                %             S_5{MODEL}(n,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_5(SET)./DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                %
                %             M_6{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_6(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
                %             S_6{MODEL}(n,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_6(SET)./DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                %             M_8{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_8(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
                %             S_8{MODEL}(n,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_8(SET)./DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                %             M_9{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_9(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
                %             S_9{MODEL}(n,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_9(SET)./DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                %
                %             Ith{MODEL}(n,ro)=mean(DA{MODEL,ro,n,kn}.inf_th);
                
                %             M_1{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_1(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_1{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_1(SET));%std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))));%
                %             M_2{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_2(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_2{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_2(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             M_3{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_3(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_3{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_3(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                %             %     M_4{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_4-1*DA{MODEL,ro,n,kn}.inf_th));
                %             %     S_4{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_4-DA{MODEL,ro,n,kn}.inf_th));
                %             M_5{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_5(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_5{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_5(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                %
                %             M_6{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_6(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_6{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_6(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_6(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                %             M_8{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_8(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_8{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_8(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_8(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                
                
                %%%%%%%%%%%%%%%%%%%%    THIS IS USED IN THE PAPER
                
                M_1{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_1(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                S_1{MODEL}(n,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                M_2{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_2(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                S_2{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_2(SET));%
                M_3{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_3(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                S_3{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_3(SET));%
                %     M_4{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_4-1*DA{MODEL,ro,n,kn}.inf_th));
                %     S_4{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_4-DA{MODEL,ro,n,kn}.inf_th));
                M_5{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_5(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                S_5{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_5(SET));%
                M_6{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_6(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                S_6{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_6(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_6(SET));%
                M_8{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_8(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                S_8{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_8(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_8(SET));%%
                M_7{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_7(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                S_7{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_7(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_6(SET));%
                M_9{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_9(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                S_9{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_9(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_8(SET));%
                
                %%%%%%%%%%%%%%%%%%%%
                
                if MODEL==7 | MODEL==5
                    DA{MODEL,ro,n,kn}.inf_th(SET)=DA{MODEL,ro,9,4}.inf_th(SET);
                end
                
                
                if MODEL~=2 & MODEL~=4
                    nM_1{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_1(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
                else
                    nM_1{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_7(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
                    %                 nM_1{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_6(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
                end
                nS_1{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_1(SET)./DA{MODEL,ro,n,kn}.inf_th(SET));%std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))));%
                nM_2{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_2(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
                nS_2{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_2(SET)./DA{MODEL,ro,n,kn}.inf_th(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                nM_3{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_3(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
                nS_3{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_3(SET)./DA{MODEL,ro,n,kn}.inf_th(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                %     M_4{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_4-1*DA{MODEL,ro,n,kn}.inf_th));
                %     S_4{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_4-DA{MODEL,ro,n,kn}.inf_th));
                nM_5{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_5(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
                nS_5{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_5(SET)./DA{MODEL,ro,n,kn}.inf_th(SET));%    std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                nM_6{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_6(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
                nS_6{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_6(SET)./DA{MODEL,ro,n,kn}.inf_th(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_6(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                nM_8{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_8(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
                nS_8{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_8(SET)./DA{MODEL,ro,n,kn}.inf_th(SET));%%  std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_8(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                nM_7{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_7(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
                nS_7{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_7(SET)./DA{MODEL,ro,n,kn}.inf_th(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_7(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                nM_9{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_9(SET)./DA{MODEL,ro,n,kn}.inf_th(SET)))-1;
                nS_9{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_9(SET)./DA{MODEL,ro,n,kn}.inf_th(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_9(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                % %
                %             M_1{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_1{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_1(SET));%std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))));%
                %             M_2{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_2{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_2(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                %             M_3{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_3{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_3(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                %             %     M_4{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_4-1*DA{MODEL,ro,n,kn}.inf_th));
                %             %     S_4{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_4-DA{MODEL,ro,n,kn}.inf_th));
                %             M_5{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_5{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_5(SET));%    std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                %             M_6{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_6(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_6{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_6(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_6(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                %             M_8{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_8(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_8{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_8(SET));%%  std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_8(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                %             M_7{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_7(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_7{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_6(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_7(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                %             M_9{MODEL}(n,ro)=mean( (DA{MODEL,ro,n,kn}.inf_9(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_9{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_8(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_9(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                
                %%%%%%%%%%%%%%%%%%%%
                
                M0_1{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_1(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                S0_1{MODEL}(n,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                M0_2{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_2(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                S0_2{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_2(SET));%
                M0_3{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_3(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                S0_3{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_3(SET));%
                M0_5{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_5(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                S0_5{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_3(SET));%
                
                
                
                Ith{MODEL}(n,ro)=mean(DA{MODEL,ro,n,kn}.inf_th(SET));
            end
        end
    end
    % N0=1./(2.^VEC_Nn);
    N0= 1:9;% (2.^VEC_Nn);
    
    
    
    
    
    SET=1:100;
    n=6;
    clear kM_1 kS_1 kM_2 kS_2 kM_3 kS_3 kM_4 kS_4 kM_5 kS_5 kIth kM_6 kS_6  kM_8 kS_8 kM_9 kS_9
    for MODEL=1:7
        for ro=1:9
            for kn=  1:6%5:13%1:7
                
                kM_1{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
                kS_1{MODEL}(kn,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_1(SET)./DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                kM_2{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
                kS_2{MODEL}(kn,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_2(SET)./DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                kM_3{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
                kS_3{MODEL}(kn,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_3(SET)./DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                %                             M_4{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_4-DA{MODEL,ro,n,kn}.inf_th)./DA{MODEL,ro,n,kn}.inf_th);
                %                             S_4{MODEL}(n,ro)=std(abs(DA{MODEL,ro,n,kn}.inf_4./DA{MODEL,ro,n,kn}.inf_th));
                kM_5{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
                kS_5{MODEL}(kn,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_5(SET)./DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                
                kM_6{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_6(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
                kS_6{MODEL}(kn,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_6(SET)./DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                kM_8{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_8(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
                kS_8{MODEL}(kn,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_8(SET)./DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                kM_9{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_9(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))./DA{MODEL,ro,n,kn}.inf_th(SET));
                kS_9{MODEL}(kn,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_9(SET)./DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                
                %             Ith{MODEL}(n,ro)=mean(DA{MODEL,ro,n,kn}.inf_th);
                
                %             M_1{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_1(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_1{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_1(SET));%std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))));%
                %             M_2{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_2(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_2{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_2(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             M_3{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_3(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_3{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_3(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                %             %     M_4{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_4-1*DA{MODEL,ro,n,kn}.inf_th));
                %             %     S_4{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_4-DA{MODEL,ro,n,kn}.inf_th));
                %             M_5{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_5(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_5{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_5(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                %
                %             M_6{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_6(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_6{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_6(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_6(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                %             M_8{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_8(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                %             S_8{MODEL}(n,ro)=std(DA{MODEL,ro,n,kn}.inf_8(SET));%std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_8(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%
                
                
                
                kM_1{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_1(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                kS_1{MODEL}(kn,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                kM_2{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_2(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                kS_2{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_2(SET));%
                kM_3{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_3(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                kS_3{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_3(SET));%
                %     M_4{MODEL}(n,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_4-1*DA{MODEL,ro,n,kn}.inf_th));
                %     S_4{MODEL}(n,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_4-DA{MODEL,ro,n,kn}.inf_th));
                kM_5{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_5(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                kS_5{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_5(SET));%
                kM_6{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_6(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                kS_6{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_6(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_6(SET));%
                kM_8{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_8(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                kS_8{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_8(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_8(SET));%%
                kM_7{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_7(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                kS_7{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_7(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_6(SET));%
                kM_9{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_9(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                kS_9{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_9(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_8(SET));%
                %
                
                kM0_1{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_1(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                kS0_1{MODEL}(kn,ro)=std(bootstrp(1000,@mean,(DA{MODEL,ro,n,kn}.inf_1(SET)-DA{MODEL,ro,n,kn}.inf_th(SET))));%std(DA{MODEL,ro,n,kn}.inf_1(SET));%
                kM0_2{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_2(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                kS0_2{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_2(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_2(SET));%
                kM0_3{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_3(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                kS0_3{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_3(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_3(SET));%
                kM0_5{MODEL}(kn,ro)=mean(abs(DA{MODEL,ro,n,kn}.inf_5(SET)-1*DA{MODEL,ro,n,kn}.inf_th(SET)));
                kS0_5{MODEL}(kn,ro)=std(bootstrp(1000,@mean,DA{MODEL,ro,n,kn}.inf_5(SET)-DA{MODEL,ro,n,kn}.inf_th(SET)));%std(DA{MODEL,ro,n,kn}.inf_3(SET));%
                
                
                
                kIth{MODEL}(kn,ro)=mean(DA{MODEL,ro,n,kn}.inf_th(SET));
            end
        end
    end
    % N0=1./(2.^VEC_Nn);
    N0= 1:9;% (2.^VEC_Nn);
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    n=6;
    N0=0.1:0.1:0.9;
    figure(201)
    for MODEL=[1 2 3 4]
        for tt=n
            subplot(3,2,MODEL)
            semilogy(N0,(M_1{MODEL}(tt,:)),'color',[0 0 1],'linewidth',3,'marker','O')
            hold on
            semilogy(N0,(M_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',3,'marker','O')
            hold on
            semilogy(N0,(M_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',3,'marker','O')
            hold on
            % errorbar(N0,mean(M_4{MODEL}(tt,:),1),mean(S_4{MODEL}(tt,:),1),'color',[0 1 0],'linewidth',3,'marker','O')
            hold on
            %         semilogy(N0,(M_5{MODEL}(tt,:)),'color','c','linewidth',3,'marker','O')
            hold on
            % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
            %         set(gca,'YScale','log')
            axis tight
            if MODEL==1 | MODEL==3
                xlabel('\rho')
                xlim([0.2 0.9])
            else
                xlabel('\nu')
            end
            ylabel('Information error')
            if MODEL==1
                legend({'I_{L1}','I_{par}','I_{kNN}'},'Location','northeast','Fontsize',12)
            end
            switch MODEL
                case 1
                    title('Gaussian Copula, Gaussian Marginal')
                case 2
                    title('student Copula, Gaussian Marginal')
                case 3
                    title('Gaussian Copula, Gamma marginal')
                case 4
                    title('Student Copula, Gamma Marginal')
            end
            ylim([0 2.4])
            
        end
    end
    
    
    
    tt=n;
    figure(201)
    subplot(3,2,5)
    plot(N0,(abs(M0_1{1}(tt,:)-M0_1{3}(tt,:))),'color',[0 0 1],'linewidth',3,'marker','O')
    hold on
    % plot(N0,(abs(M0_2{1}(tt,:)-M0_2{3}(tt,:))),'color',[1 0 1],'linewidth',3,'marker','O')
    hold on
    % plot(N0,(abs(M0_3{1}(tt,:)-M0_3{3}(tt,:))),'color',[1 0 0],'linewidth',3,'marker','O')
    hold on
    % plot((abs(M_4{1}(tt,:)-M_4{3}(tt,:))),'color',[0 1 0],'linewidth',3,'marker','O')
    hold on
    plot(N0,(abs(M0_5{1}(tt,:)-M0_5{3}(tt,:))),'color','c','linewidth',3,'marker','O')
    title('Gaussian Copula')
    ylabel('I-I_{theory}','fontsize',12)
    xlabel('\rho','fontsize',12)
    axis tight
    xlim([0.2 0.9])
    % set(gca,'YScale','log')
    subplot(3,2,6)
    plot(N0,(abs(M0_1{2}(tt,:)-M0_1{4}(tt,:))),'color',[0 0 1],'linewidth',3,'marker','O')
    hold on
    % plot(N0,(abs(M0_2{2}(tt,:)-M0_2{4}(tt,:))),'color',[1 0 1],'linewidth',3,'marker','O')
    hold on
    % plot(N0,(abs(M0_3{2}(tt,:)-M0_3{4}(tt,:))),'color',[1 0 0],'linewidth',3,'marker','O')
    hold on
    % plot((abs(M_4{1}(tt,:)-M_4{3}(tt,:))),'color',[0 1 0],'linewidth',3,'marker','O')
    hold on
    plot(N0,(abs(M0_5{2}(tt,:)-M0_5{4}(tt,:))),'color','c','linewidth',3,'marker','O')
    title('Student-t Copula')
    ylabel('I-I_{theory}','fontsize',12)
    xlabel('\nu','fontsize',12)
    axis tight
    xlim([0.1 0.9])
    
    % set(gca,'YScale','log')
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    N0=0.1:0.1:0.9;
    figure(501)
    ss=0;
    for MODEL=[1 3 2 4]
        ss=ss+1;
        for tt=n
            subplot(2,2,ss)
            semilogy(N0,(M_1{MODEL}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','O','linestyle','-','MarkerSize',13)
            hold on
            semilogy(N0,(M_2{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','O','linestyle','--','MarkerSize',13)
            hold on
            %         semilogy(N0,(M_6{MODEL}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','*','linestyle','-','MarkerSize',13)
            hold on
            %         semilogy(N0,(M_7{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','*','linestyle','--','MarkerSize',13)
            hold on
            semilogy(N0,(M_8{MODEL}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','s','linestyle','-','MarkerSize',13)
            hold on
            semilogy(N0,(M_9{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','s','linestyle','--','MarkerSize',13)
            axis tight
            if MODEL==1 | MODEL==3
                xlabel('\rho')
                xlim([0.2 0.9])
            else
                xlabel('\nu')
            end
            ylabel('Information error')
            if MODEL==1
                %             legend({'I_{LL1}^{sample}','I_{LL2}^{sample}','I_{LL1}^{integral}','I_{LL2}^{integral}','I_{LL1}^{naive}','I_{LL2}^{naive}'},'Location','northeast','Fontsize',12)
                legend({'I_{LL1}','I_{LL2}','I_{LL1}^{naive}','I_{LL2}^{naive}'},'Location','northeast','Fontsize',12)
            end
            switch MODEL
                case 1
                    title('Gaussian Copula, Gaussian Marginal')
                case 2
                    title('student Copula, Gaussian Marginal')
                case 3
                    title('Gaussian Copula, Gamma marginal')
                case 4
                    title('Student Copula, Gamma Marginal')
            end
            %         ylim([0 0.25])
            axis square
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%% naive, LL1 LL2 comparison
    
    N0=0.1:0.1:0.9;
    figure(501)
    ss=0;
    for MODEL=[1 3 2 4]
        ss=ss+1;
        tt=6;
        subplot(2,2,ss)
        errorbar(N0,M_1{MODEL}(tt,:),S_1{MODEL}(tt,:),'color',[0 0 1],'linewidth',2.5,'marker','O','linestyle','-','MarkerSize',13)
        hold on
        errorbar(N0,M_2{MODEL}(tt,:),S_2{MODEL}(tt,:),'color',[1 0 0],'linewidth',2,'marker','O','linestyle','--','MarkerSize',13)
        hold on
        errorbar(N0,M_8{MODEL}(tt,:),S_8{MODEL}(tt,:),'color',[0 0 1],'linewidth',2.5,'marker','s','linestyle','-','MarkerSize',13)
        hold on
        errorbar(N0,M_9{MODEL}(tt,:),S_9{MODEL}(tt,:),'color',[1 0 0],'linewidth',2,'marker','s','linestyle','--','MarkerSize',13)
        axis tight
        if MODEL==1 | MODEL==3
            xlabel('\rho')
            %             xlim([0.2 0.9])
        else
            xlabel('\nu')
        end
        ylabel('Information error')
        
        set(gca,'YScale','log')
        
        if MODEL==1
            %             legend({'I_{LL1}^{sample}','I_{LL2}^{sample}','I_{LL1}^{integral}','I_{LL2}^{integral}','I_{LL1}^{naive}','I_{LL2}^{naive}'},'Location','northeast','Fontsize',12)
            legend({'I_{LL1}','I_{LL2}','I_{LL1}^{naive}','I_{LL2}^{naive}'},'Location','northeast','Fontsize',12)
        end
        switch MODEL
            case 1
                title('Gaussian Copula, Gaussian Marginal')
            case 2
                title('student Copula, Gaussian Marginal')
            case 3
                title('Gaussian Copula, Gamma marginal')
            case 4
                title('Student Copula, Gamma Marginal')
        end
        %         ylim([0 0.25])
        xlim([0.2 0.9])
        axis square
    end
    
    
    tt=n;
    figure(502)
    subplot(1,2,1)
    semilogy(N0,abs(M_1{1}(tt,:)-M_1{3}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','O','linestyle','-','MarkerSize',13)
    hold on
    semilogy(N0,abs(M_2{1}(tt,:)-M_2{3}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','O','linestyle','--','MarkerSize',13)
    hold on
    semilogy(N0,abs(M_8{1}(tt,:)-M_8{3}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','s','linestyle','-','MarkerSize',13)
    hold on
    semilogy(N0,abs(M_9{1}(tt,:)-M_9{3}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','s','linestyle','--','MarkerSize',13)
    title('Gaussian Copula')
    ylabel('I-I_{theory}','fontsize',12)
    xlabel('\rho','fontsize',12)
    axis tight
    xlim([0.2 0.9])
    % set(gca,'YScale','log')
    subplot(1,2,2)
    semilogy(N0,abs(M_1{2}(tt,:)-M_1{4}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','O','linestyle','-','MarkerSize',13)
    hold on
    semilogy(N0,abs(M_2{2}(tt,:)-M_2{4}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','O','linestyle','--','MarkerSize',13)
    hold on
    semilogy(N0,abs(M_8{2}(tt,:)-M_8{4}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','s','linestyle','-','MarkerSize',13)
    hold on
    semilogy(N0,abs(M_9{2}(tt,:)-M_9{4}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','s','linestyle','--','MarkerSize',13)
    title('Student-t Copula')
    ylabel('I-I_{theory}','fontsize',12)
    xlabel('\nu','fontsize',12)
    axis tight
    xlim([0.1 0.9])
    
    % set(gca,'YScale','log')
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%% k comparison
    
    N0=0.1:0.1:0.9;
    figure(510)
    ss=0;
    ll=linspace(1,2,6);
    for MODEL=[1 3 2 4]
        ss=ss+1;
        for tt=1:6
            subplot(2,2,ss)
            errorbar(N0,kM_1{MODEL}(tt,:),kS_1{MODEL}(tt,:),'color',[0 0 1-(6-tt)/6],'linewidth',ll(tt),'marker','O','linestyle','-','MarkerSize',2*tt+1)
            hold on
            %         errorbar(N0,kM_2{MODEL}(tt,:),S_2{MODEL}(tt,:),'color',[1 0 0],'linewidth',2,'marker','O','linestyle','--','MarkerSize',13)
            hold on
            %         errorbar(N0,kM_8{MODEL}(tt,:),S_8{MODEL}(tt,:),'color',[0 0 1],'linewidth',2.5,'marker','s','linestyle','-','MarkerSize',13)
            hold on
            %         errorbar(N0,kM_9{MODEL}(tt,:),S_9{MODEL}(tt,:),'color',[1 0 0],'linewidth',2,'marker','s','linestyle','--','MarkerSize',13)
            axis tight
            if MODEL==1 | MODEL==3
                xlabel('\rho')
                %             xlim([0.2 0.9])
            else
                xlabel('\nu')
            end
            ylabel('Information error')
            
            set(gca,'YScale','log')
            
            if MODEL==1 & tt==6
                %             legend({'I_{LL1}^{sample}','I_{LL2}^{sample}','I_{LL1}^{integral}','I_{LL2}^{integral}','I_{LL1}^{naive}','I_{LL2}^{naive}'},'Location','northeast','Fontsize',12)
                legend({'k=10','k=30','k=50','k=100','k=150','k=200'},'Location','northeast','Fontsize',12)
            end
            switch MODEL
                case 1
                    title('Gaussian Copula, Gaussian Marginal')
                case 2
                    title('student Copula, Gaussian Marginal')
                case 3
                    title('Gaussian Copula, Gamma marginal')
                case 4
                    title('Student Copula, Gamma Marginal')
            end
            %         ylim([0 0.25])
            xlim([0.2 0.9])
            axis square
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%% N comparison
    
    N0=0.1:0.1:0.9;
    figure(520)
    ss=0;
    ll=linspace(1,2,6);
    for MODEL=[1 3 2 4]
        ss=ss+1;
        for tt=1:6
            subplot(2,2,ss)
            errorbar(N0,kM_1{MODEL}(tt,:),kS_1{MODEL}(tt,:),'color',[0 0 1-(6-tt)/6],'linewidth',ll(tt),'marker','O','linestyle','-','MarkerSize',2*tt+1)
            hold on
            %         errorbar(N0,kM_2{MODEL}(tt,:),S_2{MODEL}(tt,:),'color',[1 0 0],'linewidth',2,'marker','O','linestyle','--','MarkerSize',13)
            hold on
            %         errorbar(N0,kM_8{MODEL}(tt,:),S_8{MODEL}(tt,:),'color',[0 0 1],'linewidth',2.5,'marker','s','linestyle','-','MarkerSize',13)
            hold on
            %         errorbar(N0,kM_9{MODEL}(tt,:),S_9{MODEL}(tt,:),'color',[1 0 0],'linewidth',2,'marker','s','linestyle','--','MarkerSize',13)
            axis tight
            if MODEL==1 | MODEL==3
                xlabel('\rho')
                %             xlim([0.2 0.9])
            else
                xlabel('\nu')
            end
            ylabel('Information error')
            
            set(gca,'YScale','log')
            
            if MODEL==1 & tt==6
                %             legend({'I_{LL1}^{sample}','I_{LL2}^{sample}','I_{LL1}^{integral}','I_{LL2}^{integral}','I_{LL1}^{naive}','I_{LL2}^{naive}'},'Location','northeast','Fontsize',12)
                legend({'k=10','k=30','k=50','k=100','k=150','k=200'},'Location','northeast','Fontsize',12)
            end
            switch MODEL
                case 1
                    title('Gaussian Copula, Gaussian Marginal')
                case 2
                    title('student Copula, Gaussian Marginal')
                case 3
                    title('Gaussian Copula, Gamma marginal')
                case 4
                    title('Student Copula, Gamma Marginal')
            end
            %         ylim([0 0.25])
            xlim([0.2 0.9])
            axis square
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%% LNC comparison
    
    N0=0.1:0.1:0.9;
    figure(205)
    ss=0;
    for MODEL=[1 3 2 4]
        ss=ss+1;
        for tt=6
            for sss=1:2
                if sss==1
                    subplot(4,2,ss)
                    plot(N0,(M_1{MODEL}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',12)
                    %         hold on
                    %         errorbar(N0,(M_7{MODEL}(tt,:)),(S_7{MODEL}(tt,:)),'color','g','linewidth',2,'marker','O','MarkerSize',10)
                    hold on
                    %         errorbar(N0,(M_2{MODEL}(tt,:)),(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
                    hold on
                    plot(N0,(M_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',12)
                    hold on
                    plot(N0,(M_5{MODEL}(tt,:)),'color','g','linewidth',2,'marker','s','MarkerSize',12)
                    hold on
                    %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
                    hold on
                    %         plot(N0,Ith{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
                    % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
                    %                         set(gca,'YScale','log')    %%%% this should be uncommented for
                    %         paper figure (absolute error)
                    if MODEL==1 | MODEL==3
                        xlabel('\rho','Fontsize',16)
                    else
                        xlabel('\nu','Fontsize',16)
                    end
                    ylabel('Information absolute error','Fontsize',16)
                    axis tight
                    if MODEL==1
                        legend({'I_{LL1}','I_{par}','I_{LNC}'},'Location','northeast','Fontsize',12)
                    end
                    axis tight
                    xlim([0.2 0.9])
                    axis square
                    switch MODEL
                        case 1
                            title('Gaussian Copula, Gaussian Marginal')
                        case 2
                            title('student Copula, Gaussian Marginal')
                        case 3
                            title('Gaussian Copula, Gamma marginal')
                        case 4
                            title('Student Copula, Gamma Marginal')
                    end
                    
                else
                    subplot(4,2,ss+4)
                    line([0.1 0.9],[0 0],'color','k','linewidth',2)
                    hold on
                    errorbar(N0,(nM_1{MODEL}(tt,:)),(nS_1{MODEL}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',12)
                    %         hold on
                    %         errorbar(N0,(M_7{MODEL}(tt,:)),(S_7{MODEL}(tt,:)),'color','g','linewidth',2,'marker','O','MarkerSize',10)
                    hold on
                    %         errorbar(N0,(M_2{MODEL}(tt,:)),(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
                    hold on
                    errorbar(N0,(nM_3{MODEL}(tt,:)),(nS_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',12)
                    hold on
                    errorbar(N0,(nM_5{MODEL}(tt,:)),(nS_5{MODEL}(tt,:)),'color','g','linewidth',2,'marker','s','MarkerSize',12)
                    hold on
                    %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
                    hold on
                    %         plot(N0,Ith{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
                    % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
                    %         set(gca,'YScale','log')    %%%% this should be uncommented for
                    %         paper figure (absolute error)
                    if MODEL==1 | MODEL==3
                        xlabel('\rho','Fontsize',16)
                    else
                        xlabel('\nu','Fontsize',16)
                    end
                    ylabel('Information bias (bits)','Fontsize',16)
                    axis tight
                    xlim([0.2 0.9])
                    axis square
                    switch MODEL
                        case 1
                            title('Gaussian Copula, Gaussian Marginal')
                        case 2
                            title('student Copula, Gaussian Marginal')
                        case 3
                            title('Gaussian Copula, Gamma marginal')
                        case 4
                            title('Student Copula, Gamma Marginal')
                    end
                    
                end
            end
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%% LNC comparison  VARIANCE
    
    N0=0.1:0.1:0.9;
    figure(2064)
    ss=0;
    for MODEL=[1 3 2 4]
        ss=ss+1;
        for tt=6
            subplot(2,2,ss)
            
            plot(N0, (S_1{MODEL}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',10)
            hold on
            %         errorbar(N0,(M_2{MODEL}(tt,:)),(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
            hold on
            plot(N0, (S_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',10)
            hold on
            plot(N0, (S_5{MODEL}(tt,:)),'color','c','linewidth',2,'marker','s','MarkerSize',10)
            hold on
            %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
            hold on
            % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
            % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
            %         set(gca,'YScale','log')    %%%% this should be uncommented for
            %         paper figure (absolute error)
            if MODEL==1 | MODEL==3
                xlabel('\rho','Fontsize',16)
            else
                xlabel('\nu','Fontsize',16)
            end
            ylabel('Information error','Fontsize',16)
            axis tight
            if MODEL==1
                legend({'I_{LL1}','I_{par}','I_{LNC}'},'Location','northeast','Fontsize',12)
            end
        end
        switch MODEL
            case 1
                title('Gaussian Copula, Gaussian Marginal')
            case 2
                title('student Copula, Gaussian Marginal')
            case 3
                title('Gaussian Copula, Gamma marginal')
            case 4
                title('Student Copula, Gamma Marginal')
        end
        xlim([0.2 0.9])
        axis square
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% LNC difference plot
    
    tt=6;
    figure(270)
    subplot(1,2,1)
    plot(N0,(abs(M0_1{1}(tt,:)-M0_1{3}(tt,:))),'color',[0 0 1],'linewidth',2,'marker','O')
    hold on
    % plot(N0,(abs(M0_2{1}(tt,:)-M0_2{3}(tt,:))),'color',[1 0 1],'linewidth',3,'marker','O')
    hold on
    % plot(N0,(abs(M0_3{1}(tt,:)-M0_3{3}(tt,:))),'color',[1 0 0],'linewidth',3,'marker','O')
    hold on
    % plot((abs(M_4{1}(tt,:)-M_4{3}(tt,:))),'color',[0 1 0],'linewidth',3,'marker','O')
    hold on
    plot(N0,(abs(M0_5{1}(tt,:)-M0_5{3}(tt,:))),'color','c','linewidth',2,'marker','d')
    title('Gaussian Copula')
    ylabel('I-I_{theory}','fontsize',12)
    xlabel('\rho','fontsize',12)
    axis tight
    xlim([0.2 0.9])
    % set(gca,'YScale','log')
    legend({'I_{LL1}','I_{LNC}'},'Location','northwest','Fontsize',12)
    axis square
    subplot(1,2,2)
    plot(N0,(abs(M0_1{2}(tt,:)-M0_1{4}(tt,:))),'color',[0 0 1],'linewidth',2,'marker','O')
    hold on
    % plot(N0,(abs(M0_2{2}(tt,:)-M0_2{4}(tt,:))),'color',[1 0 1],'linewidth',3,'marker','O')
    hold on
    % plot(N0,(abs(M0_3{2}(tt,:)-M0_3{4}(tt,:))),'color',[1 0 0],'linewidth',3,'marker','O')
    hold on
    % plot((abs(M_4{1}(tt,:)-M_4{3}(tt,:))),'color',[0 1 0],'linewidth',3,'marker','O')
    hold on
    plot(N0,(abs(M0_5{2}(tt,:)-M0_5{4}(tt,:))),'color','c','linewidth',2,'marker','d')
    title('Student-t Copula')
    ylabel('I-I_{theory}','fontsize',12)
    xlabel('\nu','fontsize',12)
    axis tight
    xlim([0.2 0.9])
    % set(gca,'YScale','log')
    axis square
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%% LNC comparison NUMBER  N
    
    N0=[5:13];
    figure(215)
    ss=0;
    for MODEL=[1 3 2 4]
        ss=ss+1;
        tt=1:9;
        if MODEL==1 | MODEL==3
            ttt=5;
        else
            ttt=5;
        end
        for sss=1:2
            if sss==1
                subplot(4,2,ss)
                
                errorbar(N0,(M_1{MODEL}(tt,ttt)),(S_1{MODEL}(tt,ttt)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',12)
                %         hold on
                %         errorbar(N0,(M_7{MODEL}(tt,ttt)),(S_7{MODEL}(tt,ttt)),'color','g','linewidth',2,'marker','O')
                hold on
                errorbar(N0,(M_3{MODEL}(tt,ttt)),(S_3{MODEL}(tt,ttt)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',12)
                hold on
                errorbar(N0,(M_5{MODEL}(tt,ttt)),(S_5{MODEL}(tt,ttt)),'color','g','linewidth',2,'marker','s','MarkerSize',12)
                hold on
                %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
                hold on
                % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
                % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
                %             set(gca,'YScale','log')
                if MODEL==1 | MODEL==3
                    xlabel('N','Fontsize',16)
                else
                    xlabel('N','Fontsize',16)
                end
                ylabel('Information absolute error (bits)','Fontsize',16)
                axis tight
                if MODEL==1
                    legend({'I_{LL1}','I_{par}','I_{LNC}'},'Location','northeast','Fontsize',12)
                end
                switch MODEL
                    case 1
                        title('Gaussian Copula, Gaussian Marginal')
                    case 2
                        title('student Copula, Gaussian Marginal')
                    case 3
                        title('Gaussian Copula, Gamma marginal')
                    case 4
                        title('Student Copula, Gamma Marginal')
                end
                %     xlim([0.2 0.9])
                axis square
            else
                subplot(4,2,ss+4)
                
                line([5 13],[0 0],'color','k','linewidth',2)
                hold on
                
                errorbar(N0,(nM_1{MODEL}(tt,ttt)),(nS_1{MODEL}(tt,ttt)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',12)
                %         hold on
                %         errorbar(N0,(M_7{MODEL}(tt,ttt)),(S_7{MODEL}(tt,ttt)),'color','g','linewidth',2,'marker','O')
                hold on
                errorbar(N0,(nM_3{MODEL}(tt,ttt)),(nS_3{MODEL}(tt,ttt)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',12)
                hold on
                errorbar(N0,(nM_5{MODEL}(tt,ttt)),(nS_5{MODEL}(tt,ttt)),'color','g','linewidth',2,'marker','s','MarkerSize',12)
                hold on
                %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
                hold on
                % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
                % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
                %             set(gca,'YScale','log')
                if MODEL==1 | MODEL==3
                    xlabel('N','Fontsize',16)
                else
                    xlabel('N','Fontsize',16)
                end
                ylabel('Information bias (bits)','Fontsize',16)
                axis tight
                if MODEL==1
                    legend({'I_{LL1}','I_{par}','I_{LNC}'},'Location','northeast','Fontsize',12)
                end
                switch MODEL
                    case 1
                        title('Gaussian Copula, Gaussian Marginal')
                    case 2
                        title('student Copula, Gaussian Marginal')
                    case 3
                        title('Gaussian Copula, Gamma marginal')
                    case 4
                        title('Student Copula, Gamma Marginal')
                end
                %     xlim([0.2 0.9])
                axis square
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%% PYM comparison
    
    N0=[0.1:0.1:0.9]*100;
    figure(305)
    ss=0;
    for MODEL=[7 5]
        ss=ss+1;
        for tt=6
            subplot(1,2,ss)
            errorbar(N0,(M_1{MODEL}(tt,:)),(S_1{MODEL}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',10)
            hold on
            %         errorbar(N0,(M_2{MODEL}(tt,:)),(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
            hold on
            errorbar(N0,(M_3{MODEL}(tt,:)),(S_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',10)
            hold on
            errorbar(N0,(M_5{MODEL}(tt,:)),(S_5{MODEL}(tt,:)),'color','c','linewidth',2,'marker','s','MarkerSize',10)
            hold on
            %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
            hold on
            % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
            % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
            set(gca,'YScale','log')
            if MODEL==1 | MODEL==3
                xlabel('\rho','Fontsize',16)
            else
                xlabel('\nu','Fontsize',16)
            end
            ylabel('Information error','Fontsize',16)
            axis tight
            if MODEL==1
                legend({'I_{LL1}','I_{par}','I_{LNC}'},'Location','northeast','Fontsize',12)
            end
        end
        switch MODEL
            case 1
                title('Gaussian Copula, Gaussian Marginal')
            case 2
                title('student Copula, Gaussian Marginal')
            case 3
                title('Gaussian Copula, Gamma marginal')
            case 4
                title('Student Copula, Gamma Marginal')
        end
        xlim([10 70])
        axis square
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% PYM k comparison
    
    N0=[0.1:0.1:0.9]*100;
    figure(310)
    ss=0;
    ll=linspace(1,2,6);
    for MODEL=[7 5]
        ss=ss+1;
        for tt=1:6
            subplot(1,2,ss)
            errorbar(N0,kM_1{MODEL}(tt,:),kS_1{MODEL}(tt,:),'color',[0 0 1-(6-tt)/6],'linewidth',ll(tt),'marker','O','linestyle','-','MarkerSize',2*tt+1)
            hold on
            if tt==4
                hold on
                errorbar(N0,(M_3{MODEL}(tt,:)),(S_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',10)
                hold on
                errorbar(N0,kM_5{MODEL}(tt,:),kS_5{MODEL}(tt,:),'color','c','linewidth',2,'marker','O','linestyle','--','MarkerSize',13)
            end
            
            xlabel('N')
            ylabel('Information error')
            
            set(gca,'YScale','log')
            
            if MODEL==7 & tt==6
                %             legend({'I_{LL1}^{sample}','I_{LL2}^{sample}','I_{LL1}^{integral}','I_{LL2}^{integral}','I_{LL1}^{naive}','I_{LL2}^{naive}'},'Location','northeast','Fontsize',12)
                legend({'k=10','k=30','k=50','k=100','k=150','k=200'},'Location','northeast','Fontsize',12)
            end
            switch MODEL
                case 1
                    title('Gaussian Copula, Gaussian Marginal')
                case 2
                    title('student Copula, Gaussian Marginal')
                case 3
                    title('Gaussian Copula, Gamma marginal')
                case 4
                    title('Student Copula, Gamma Marginal')
            end
            %         ylim([0 0.25])
            xlim([20 70])
            axis square
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% PYM comparison NUMBER  N
    
    N0=[5:13];
    figure(315)
    ss=0;
    for MODEL=[7 5]
        ss=ss+1;
        tt=1:9;
        ttt=5;
        
        for sss=1:2
            if sss==1
                subplot(2,2,ss)
                
                plot(N0,(M_1{MODEL}(tt,ttt)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',12)
                hold on
                %         errorbar(N0,(M_2{MODEL}(tt,:)),(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
                hold on
                plot(N0,(M_3{MODEL}(tt,ttt)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',12)
                hold on
                plot(N0,(M_5{MODEL}(tt,ttt)),'color','g','linewidth',2,'marker','s','MarkerSize',12)
                hold on
                %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
                hold on
                % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
                % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
                set(gca,'YScale','log')
                if MODEL==1 | MODEL==3
                    xlabel('N','Fontsize',16)
                else
                    xlabel('N','Fontsize',16)
                end
                ylabel('Information error','Fontsize',16)
                axis tight
                if MODEL==1
                    legend({'I_{LL1}','I_{par}','I_{PYM}'},'Location','northeast','Fontsize',12)
                end
                switch MODEL
                    case 1
                        title('Gaussian Copula, Gaussian Marginal')
                    case 2
                        title('student Copula, Gaussian Marginal')
                    case 3
                        title('Gaussian Copula, Gamma marginal')
                    case 4
                        title('Student Copula, Gamma Marginal')
                end
                %     xlim([0.2 0.9])
                axis square
            else
                subplot(2,2,ss+2)
                line([5 13],[0 0],'color','k','linewidth',2)
                hold on
                
                errorbar(N0,(nM_1{MODEL}(tt,ttt)),(nS_1{MODEL}(tt,ttt)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',12)
                hold on
                %         errorbar(N0,(M_2{MODEL}(tt,:)),(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
                hold on
                errorbar(N0,(nM_3{MODEL}(tt,ttt)),(nS_3{MODEL}(tt,ttt)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',12)
                hold on
                errorbar(N0,(nM_5{MODEL}(tt,ttt)),(nS_5{MODEL}(tt,ttt)),'color','g','linewidth',2,'marker','s','MarkerSize',12)
                hold on
                %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
                hold on
                % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
                % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
                %                 set(gca,'YScale','log')
                if MODEL==1 | MODEL==3
                    xlabel('N','Fontsize',16)
                else
                    xlabel('N','Fontsize',16)
                end
                ylabel('Information error','Fontsize',16)
                axis tight
                if MODEL==1
                    legend({'I_{LL1}','I_{par}','I_{PYM}'},'Location','northeast','Fontsize',12)
                end
                switch MODEL
                    case 1
                        title('Gaussian Copula, Gaussian Marginal')
                    case 2
                        title('student Copula, Gaussian Marginal')
                    case 3
                        title('Gaussian Copula, Gamma marginal')
                    case 4
                        title('Student Copula, Gamma Marginal')
                end
                %     xlim([0.2 0.9])
                axis square
            end
        end
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% PYM variance-N plot
    
    N0=5:13;
    figure(270)
    subplot(1,2,1)
    plot(N0,std(M_1{7}(:,2:7)'),'color',[0 0 1],'linewidth',2,'marker','O')
    hold on
    % plot(N0,(abs(M0_2{1}(tt,:)-M0_2{3}(tt,:))),'color',[1 0 1],'linewidth',3,'marker','O')
    hold on
    plot(N0,std(M_3{7}(:,2:7)'),'color',[1 0 0],'linewidth',3,'marker','O')
    hold on
    % plot((abs(M_4{1}(tt,:)-M_4{3}(tt,:))),'color',[0 1 0],'linewidth',3,'marker','O')
    hold on
    plot(N0,std(M_5{7}(:,2:7)'),'color','c','linewidth',2,'marker','d')
    title('Gaussian Copula')
    ylabel('I-I_{theory}','fontsize',12)
    xlabel('\rho','fontsize',12)
    axis tight
    % xlim([0.2 0.9])
    % set(gca,'YScale','log')
    legend({'I_{LL1}','I_{LNC}'},'Location','northwest','Fontsize',12)
    axis square
    subplot(1,2,2)
    plot(N0,std(M_1{5}(:,2:7)'),'color',[0 0 1],'linewidth',2,'marker','O')
    hold on
    % plot(N0,(abs(M0_2{2}(tt,:)-M0_2{4}(tt,:))),'color',[1 0 1],'linewidth',3,'marker','O')
    hold on
    plot(N0,std(M_3{5}(:,2:7)'),'color',[1 0 0],'linewidth',3,'marker','O')
    hold on
    % plot((abs(M_4{1}(tt,:)-M_4{3}(tt,:))),'color',[0 1 0],'linewidth',3,'marker','O')
    hold on
    plot(N0,std(M_5{5}(:,2:7)'),'color','c','linewidth',2,'marker','d')
    title('Student-t Copula')
    ylabel('I-I_{theory}','fontsize',12)
    xlabel('\nu','fontsize',12)
    axis tight
    % xlim([0.2 0.9])
    % set(gca,'YScale','log')
    axis square
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%% PYM comparison NUMBER  N
    
    N0=[5:13];
    figure(215)
    ss=0;
    for MODEL=[1 3 2 4]
        ss=ss+1;
        tt=1:9;
        if MODEL==1 | MODEL==3
            ttt=5;
        else
            ttt=5;
        end
        subplot(2,2,ss)
        
        errorbar(N0,(M_1{MODEL}(tt,ttt)),(S_1{MODEL}(tt,ttt)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',10)
        hold on
        %         errorbar(N0,(M_2{MODEL}(tt,:)),(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
        hold on
        errorbar(N0,(M_3{MODEL}(tt,ttt)),(S_3{MODEL}(tt,ttt)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',10)
        hold on
        errorbar(N0,(M_5{MODEL}(tt,ttt)),(S_5{MODEL}(tt,ttt)),'color','c','linewidth',2,'marker','s','MarkerSize',10)
        hold on
        %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
        hold on
        % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
        % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
        set(gca,'YScale','log')
        if MODEL==1 | MODEL==3
            xlabel('N','Fontsize',16)
        else
            xlabel('N','Fontsize',16)
        end
        ylabel('Information error','Fontsize',16)
        axis tight
        if MODEL==1
            legend({'I_{LL1}','I_{par}','I_{LNC}'},'Location','northeast','Fontsize',12)
        end
        switch MODEL
            case 1
                title('Gaussian Copula, Gaussian Marginal')
            case 2
                title('student Copula, Gaussian Marginal')
            case 3
                title('Gaussian Copula, Gamma marginal')
            case 4
                title('Student Copula, Gamma Marginal')
        end
        %     xlim([0.2 0.9])
        axis square
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    n=6;
    N0=0.1:0.1:0.9;
    figure(800)
    for MODEL=[5 6 7]
        for tt=n
            subplot(2,2,MODEL-4)
            
            errorbar(N0,(M_1{MODEL}(tt,:)),(S_1{MODEL}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','O')
            hold on
            errorbar(N0,(M_2{MODEL}(tt,:)),(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
            hold on
            errorbar(N0,(M_3{MODEL}(tt,:)),(S_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','O')
            hold on
            errorbar(N0,(M_5{MODEL}(tt,:)),(S_5{MODEL}(tt,:)),'color','c','linewidth',2,'marker','O')
            hold on
            errorbar(N0,(M_8{MODEL}(tt,:)),(S_6{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
            hold on
            % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
            % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
            %         set(gca,'YScale','log')
            if MODEL==1 | MODEL==3
                xlabel('\rho')
            else
                xlabel('\nu')
            end
            ylabel('Information error')
            axis tight
            if MODEL==1
                legend({'I_{L1}','I_{L2}','I_{par}','I_{kNN}'},'Location','northeast','Fontsize',12)
            end
        end
        switch MODEL
            case 5
                title('Gaussian Copula, Poisson Marginal')
            case 6
                title('student Copula, Poisson Marginal')
            case 7
                title('Gaussian Copula, Bernouli marginal')
            case 8
                title('Student Copula, Bernouli Marginal')
        end
    end
    
    
    
    
    
    
    N0=2.^[5 6 7 8 9 10 11 12 13];
    figure(900)
    ii=0;
    for MODEL=[1 2 3 4]
        for tt=1:9
            %         subplot(2,2,MODEL-4)
            ii=ii+1;
            subplot(4,9,ii)
            errorbar(N0,(M_1{MODEL}(:,tt)),(S_1{MODEL}(:,tt)),'color',[0 0 1],'linewidth',2,'marker','O')
            hold on
            errorbar(N0,(M_2{MODEL}(:,tt)),(S_2{MODEL}(:,tt)),'color',[1 0 1],'linewidth',2,'marker','O')
            hold on
            %         errorbar(N0,(M_3{MODEL}(:,tt)),(S_3{MODEL}(:,tt)),'color',[1 0 0],'linewidth',2,'marker','O')
            hold on
            errorbar(N0,(M_5{MODEL}(:,tt)),(S_5{MODEL}(:,tt)),'color','c','linewidth',2,'marker','O')
            hold on
            %         errorbar(N0,(M_8{MODEL}(:,tt)),(S_8{MODEL}(:,tt)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
            hold on
            % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
            % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
            if MODEL==1 | MODEL==3
                xlabel('\rho')
            else
                xlabel('\nu')
            end
            ylabel('Information error')
            axis tight
            ylim([0 4])
            set(gca,'YScale','log')
            if MODEL==5 & tt==1
                legend({'I_{L1}','I_{L2}','I_{kNN}'},'Location','northeast','Fontsize',12)
            end
        end
        switch MODEL
            case 5
                title('Gaussian Copula, Poisson Marginal')
            case 6
                title('student Copula, Poisson Marginal')
            case 7
                title('Gaussian Copula, Bernouli marginal')
            case 8
                title('Student Copula, Bernouli Marginal')
        end
    end
    
    
    
    
    N0=2.^[5 6 7 8 9 10 11];
    figure(910)
    ii=0;
    for MODEL=[5 6]
        for tt=1:9
            %         subplot(2,2,MODEL-4)
            ii=ii+1;
            subplot(2,9,ii)
            errorbar(N0,(M_1{MODEL}(:,tt)),(S_1{MODEL}(:,tt)),'color',[0 0 1],'linewidth',2,'marker','O')
            hold on
            errorbar(N0,(M_2{MODEL}(:,tt)),(S_2{MODEL}(:,tt)),'color',[1 0 1],'linewidth',2,'marker','O')
            hold on
            %         errorbar(N0,(M_3{MODEL}(:,tt)),(S_3{MODEL}(:,tt)),'color',[1 0 0],'linewidth',2,'marker','O')
            hold on
            errorbar(N0,(M_5{MODEL}(:,tt)),(S_5{MODEL}(:,tt)),'color','c','linewidth',2,'marker','O')
            hold on
            %         errorbar(N0,(M_8{MODEL}(:,tt)),(S_8{MODEL}(:,tt)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
            hold on
            % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
            % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
            if MODEL==1 | MODEL==3
                xlabel('\rho')
            else
                xlabel('\nu')
            end
            ylabel('Information error')
            axis tight
            ylim([0 4])
            set(gca,'YScale','log')
            if MODEL==5 & tt==1
                legend({'I_{L1}','I_{L2}','I_{kNN}'},'Location','northeast','Fontsize',12)
            end
        end
        switch MODEL
            case 5
                title('Gaussian Copula, Poisson Marginal')
            case 6
                title('student Copula, Poisson Marginal')
            case 7
                title('Gaussian Copula, Bernouli marginal')
            case 8
                title('Student Copula, Bernouli Marginal')
        end
    end
    
    
    
    
    
    N0=0.1:0.1:0.9;
    figure(202)
    for MODEL=[1 2 3 4]
        for tt=n
            subplot(2,2,MODEL)
            semilogy(N0,(M_1{MODEL}(tt,:)),'color',[0 0 1],'linewidth',2,'marker','O')
            hold on
            semilogy(N0,(M_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
            hold on2
            semilogy(N0,(M_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',2,'marker','O')
            hold on
            semilogy(N0,(M_5{MODEL}(tt,:)),'color','c','linewidth',2,'marker','O')
            hold on
            semilogy(N0,(M_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
            hold on
            % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
            % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
            % set(gca,'YScale','log')
            if MODEL==1 | MODEL==3
                xlabel('\rho')
            else
                xlabel('\nu')
            end
            ylabel('Information error')
            axis tight
            if MODEL==1
                legend({'I_{L1}','I_{L2}','I_{par}','I_{kNN}'},'Location','northeast','Fontsize',12)
            end
        end
        switch MODEL
            case 1
                title('Gaussian Copula, Gaussian Marginal')
            case 2
                title('student Copula, Gaussian Marginal')
            case 3
                title('Gaussian Copula, Gamma marginal')
            case 4
                title('Student Copula, Gamma Marginal')
        end
    end
    
    
    
    
    N0=0.1:0.1:0.9;
    figure(300)
    for MODEL=1:4
        for tt=n
            subplot(2,2,MODEL)
            plot(N0,(S_1{MODEL}(tt,:)),'color',[0 0 1],'linewidth',3,'marker','O')
            hold on
            plot(N0,(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',3,'marker','O')
            hold on
            plot(N0,(S_3{MODEL}(tt,:)),'color',[1 0 0],'linewidth',3,'marker','O')
            hold on
            % plot(N0,mean(S_4{MODEL}(tt,:),1),'color',[0 1 0],'linewidth',3,'marker','O')
            hold on
            plot(N0,(S_5{MODEL}(tt,:)),'color','c','linewidth',3,'marker','O')
            hold on
            % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
            % set(gca,'YScale','log')
            xlabel('\rho')
            ylabel('Information error')
            axis tight
            if MODEL==1
                legend({'I_{L1}','I_{L2}','I_{par}','I_{kNN}'},'Location','northeast','Fontsize',12)
            end
        end
    end
    
    
    
    %%%%%%%%%%%%% variance
    
    
    
    N0=[5:13];
    figure(216)
    ss=0;
    for MODEL=[1 3 2 4]
        ss=ss+1;
        tt=1:9;
        if MODEL==1 | MODEL==3
            ttt=5;
        else
            ttt=5;
        end
        subplot(2,2,ss)
        
        plot(N0,(S_1{MODEL}(tt,ttt)),'color',[0 0 1],'linewidth',2,'marker','O','MarkerSize',10)
        hold on
        %         errorbar(N0,(M_2{MODEL}(tt,:)),(S_2{MODEL}(tt,:)),'color',[1 0 1],'linewidth',2,'marker','O')
        hold on
        plot(N0,(S_3{MODEL}(tt,ttt)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',10)
        hold on
        plot(N0,(S_5{MODEL}(tt,ttt)),'color','c','linewidth',2,'marker','s','MarkerSize',10)
        hold on
        %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
        hold on
        %         line([5 13],[0 0],'color','k','linewidth',2)
        % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
        % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
        set(gca,'YScale','log')
        if MODEL==1 | MODEL==3
            xlabel('N','Fontsize',16)
        else
            xlabel('N','Fontsize',16)
        end
        ylabel('Information error','Fontsize',16)
        axis tight
        if MODEL==1
            legend({'I_{LL1}','I_{par}','I_{LNC}'},'Location','northeast','Fontsize',12)
        end
        switch MODEL
            case 1
                title('Gaussian Copula, Gaussian Marginal')
            case 2
                title('student Copula, Gaussian Marginal')
            case 3
                title('Gaussian Copula, Gamma marginal')
            case 4
                title('Student Copula, Gamma Marginal')
        end
        %     xlim([0.2 0.9])
        axis square
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    mm=0;
    for N=6%1:9
        for k=4%1:10
            for MODEL=1:4
                if 1==1%rem(numel(DA{MODEL}.inf_1{N,k}),9)==0 & numel(DA{MODEL}.inf_1{N,k})~=0
                    clear err
                    %             for i=1:5
                    %                 err(:,:,i)=reshape(DA{MODEL}.er{N,k}(:,i),9,[]);
                    %             end
                    mm=MODEL;
                    clear I1 I2 I3 I4 I5 Ith
                    for ro=1:9
                        I1(:,ro)=DA{MODEL,ro,N,k}.inf_1(SET);
                        I2(:,ro)=DA{MODEL,ro,N,k}.inf_2(SET);
                        I3(:,ro)=DA{MODEL,ro,N,k}.inf_3(SET);
                        I4(:,ro)=DA{MODEL,ro,N,k}.inf_4(SET);
                        I5(:,ro)=DA{MODEL,ro,N,k}.inf_5(SET);
                        Ith(:,ro)=DA{MODEL,ro,N,k}.inf_th(SET);
                    end
                    
                    figure(100)
                    subplot(2,2,mm)
                    plot((mean(I1,1)),'color',[0 0 1],'linewidth',3)
                    hold on
                    plot((mean(I2,1)),'color',[1 0 1],'linewidth',3)
                    hold on
                    plot((mean(I3,1)),'color',[1 0 0],'linewidth',3)
                    hold on
                    plot((mean(I4,1)),'color',[0 1 0],'linewidth',3)
                    hold on
                    plot((mean(I5,1)),'color','c','linewidth',3)
                    hold on
                    plot((mean(Ith,1)),'--k')
                    axis tight
                    if mm==1
                        legend({'I_{L1}','I_{L2}','I_{par}','I_{kNN}','I_{KSG}','I_{ther.}'},'Location','northwest','Fontsize',12)
                    end
                    set(gca,'YScale','log')
                end
            end
        end
    end
    
    
    
    mm=0;
    for N=[2 4 6]%1:9
        for k=5%1:10
            mm=mm+1;
            if rem(numel(DA{MODEL}.inf_1{N,k}),9)==0 & numel(DA{MODEL}.inf_1{N,k})~=0
                clear err
                for i=1:5
                    %                 err(:,:,i)=reshape(DA{MODEL}.er{N,k}(:,i),9,[]);
                end
                
                err(:,:,1)=abs(reshape(DA{MODEL}.inf_1{N,k},9,[])-reshape(DA{MODEL}.inf_th{N,k},9,[]));
                err(:,:,2)=abs(reshape(DA{MODEL}.inf_2{N,k},9,[])-reshape(DA{MODEL}.inf_th{N,k},9,[]));
                err(:,:,3)=abs(reshape(DA{MODEL}.inf_3{N,k},9,[])-reshape(DA{MODEL}.inf_th{N,k},9,[]));
                err(:,:,4)=abs(reshape(DA{MODEL}.inf_4{N,k},9,[])-reshape(DA{MODEL}.inf_th{N,k},9,[]));
                err(:,:,5)=abs(reshape(DA{MODEL}.inf_5{N,k},9,[])-reshape(DA{MODEL}.inf_th{N,k},9,[]));
                
                I1=reshape(DA{MODEL}.inf_1{N,k},9,[]);
                I2=reshape(DA{MODEL}.inf_2{N,k},9,[]);
                I3=reshape(DA{MODEL}.inf_3{N,k},9,[]);
                I4=reshape(DA{MODEL}.inf_4{N,k},9,[]);
                I5=reshape(DA{MODEL}.inf_5{N,k},9,[]);
                Ith=reshape(DA{MODEL}.inf_th{N,k},9,[]);
                
                figure(150)
                subplot(1,3,mm)
                plot((mean(I1,2)),'color',[0 0 1],'linewidth',3)
                hold on
                plot((mean(I2,2)),'color',[1 0 1],'linewidth',3)
                hold on
                plot((mean(I3,2)),'color',[1 0 0],'linewidth',3)
                hold on
                plot((mean(I4,2)),'color',[0 1 0],'linewidth',3)
                hold on
                plot((mean(I5,2)),'color','c','linewidth',3)
                hold on
                plot((mean(Ith,2)),'--k')
                axis tight
                if mm==1
                    legend({'I_{L1}','I_{L2}','I_{par}','I_{kNN}','I_{KSG}','I_{ther.}'},'Location','northwest','Fontsize',12)
                end
                %         figure(100)
                %         subplot(9,10,mm)
                %         semilogy((mean(I1,2)),'color',[0 0 1],'linewidth',3)
                %         hold on
                %         semilogy((mean(I2,2)),'color',[1 0 1],'linewidth',3)
                %         hold on
                %         semilogy((mean(I3,2)),'color',[1 0 0],'linewidth',3)
                %         hold on
                %         semilogy((mean(I4,2)),'color',[0 1 0],'linewidth',3)
                %         hold on
                %         semilogy((mean(I5,2)),'color','c','linewidth',3)
                %         hold on
                %         semilogy((mean(Ith,2)),'--k')
                %         axis tight
                %         if mm==1
                %         legend({'I_{L1}','I_{L2}','I_{par}','I_{kNN}','I_{KSG}','I_{ther.}'},'Location','northwest','Fontsize',12)
                %         end
                
                figure(250)
                subplot(1,3,mm)
                plot(squeeze(mean(err(:,:,1),2)),'color',[0 0 1],'linewidth',3)
                hold on
                plot(squeeze(mean(err(:,:,2),2)),'color',[1 0 1],'linewidth',3)
                hold on
                plot(squeeze(mean(err(:,:,3),2)),'color',[1 0 0],'linewidth',3)
                hold on
                plot(squeeze(mean(err(:,:,4),2)),'color',[0 1 0],'linewidth',3)
                hold on
                plot(squeeze(mean(err(:,:,5),2)),'color','c','linewidth',3)
                axis tight
                if mm==1
                    legend({'I_{L1}','I_{L2}','I_{par}','I_{kNN}','I_{KSG}'},'Location','northeast','Fontsize',12)
                end
            end
        end
    end
    
    
    for N=1:9
        in1=reshape(DA{MODEL}.inf_1{N,5},9,[]);
        in2=reshape(DA{MODEL}.inf_2{N,5},9,[]);
        in3=reshape(DA{MODEL}.inf_3{N,5},9,[]);
        in4=reshape(DA{MODEL}.inf_4{N,5},9,[]);
        in5=reshape(DA{MODEL}.inf_5{N,5},9,[]);
        inth=real(reshape(DA{MODEL}.inf_th{N,5},9,[]));
        
        Me1(N)=mean(abs(in1(5,:)-inth(5,:)));
        Me2(N)=mean(abs(in2(5,:)-inth(5,:)));
        Me3(N)=mean(abs(in3(5,:)-inth(5,:)));
        Me4(N)=mean(abs(in4(5,:)-inth(5,:)));
        Me5(N)=mean(abs(in5(5,:)-inth(5,:)));
    end
    
    figure(350)
    plot(Me1,'color',[0 0 1],'linewidth',3)
    hold on
    plot(Me2,'color',[1 0 1],'linewidth',3)
    hold on
    plot(Me3,'color',[1 0 0],'linewidth',3)
    hold on
    plot(Me4,'color',[0 1 0],'linewidth',3)
    hold on
    plot(Me5,'color','c','linewidth',3)
    axis tight
    if mm==1
        legend({'I_{L1}','I_{L2}','I_{par}','I_{kNN}','I_{KSG}'},'Location','northeast','Fontsize',12)
    end
    
    
    plot(N0,(S_3{MODEL}(tt,ttt)),'color',[1 0 0],'linewidth',2,'marker','d','MarkerSize',10)
    hold on
    plot(N0,(S_5{MODEL}(tt,ttt)),'color','c','linewidth',2,'marker','s','MarkerSize',10)
    hold on
    %         errorbar(N0,(M_8{MODEL}(tt,:)),(S_8{MODEL}(tt,:)),'color',[0.5 0.5 0.5],'linewidth',2,'marker','O')
    hold on
    % errorbar(N0,Ith{MODEL}(tt,:),0*S_4{MODEL}(tt,:),'-O','color',[0 0 0],'linewidth',3)
    % set(gca,'xtick',1:n,'xticklabels',1./2.^VEC_Nn(1:n),'YScale','log')
    set(gca,'YScale','log')
    if MODEL==1 | MODEL==3
        xlabel('N','Fontsize',16)
    else
        xlabel('N','Fontsize',16)
    end
    ylabel('Information error','Fontsize',16)
    axis tight
    if MODEL==1
        legend({'I_{LL1}','I_{par}','I_{LNC}'},'Location','northeast','Fontsize',12)
    end
    switch MODEL
        case 1
            title('Gaussian Copula, Gaussian Marginal')
        case 2
            title('student Copula, Gaussian Marginal')
        case 3
            title('Gaussian Copula, Gamma marginal')
        case 4
            title('Student Copula, Gamma Marginal')
    end
    %     xlim([0.2 0.9])
    axis square
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




mm=0;
for N=6%1:9
    for k=4%1:10
        for MODEL=1:4
            if 1==1%rem(numel(DA{MODEL}.inf_1{N,k}),9)==0 & numel(DA{MODEL}.inf_1{N,k})~=0
                clear err
                %             for i=1:5
                %                 err(:,:,i)=reshape(DA{MODEL}.er{N,k}(:,i),9,[]);
                %             end
                mm=MODEL;
                clear I1 I2 I3 I4 I5 Ith
                for ro=1:9
                    I1(:,ro)=DA{MODEL,ro,N,k}.inf_1(SET);
                    I2(:,ro)=DA{MODEL,ro,N,k}.inf_2(SET);
                    I3(:,ro)=DA{MODEL,ro,N,k}.inf_3(SET);
                    I4(:,ro)=DA{MODEL,ro,N,k}.inf_4(SET);
                    I5(:,ro)=DA{MODEL,ro,N,k}.inf_5(SET);
                    Ith(:,ro)=DA{MODEL,ro,N,k}.inf_th(SET);
                end
                
                figure(100)
                subplot(2,2,mm)
                plot((mean(I1,1)),'color',[0 0 1],'linewidth',3)
                hold on
                plot((mean(I2,1)),'color',[1 0 1],'linewidth',3)
                hold on
                plot((mean(I3,1)),'color',[1 0 0],'linewidth',3)
                hold on
                plot((mean(I4,1)),'color',[0 1 0],'linewidth',3)
                hold on
                plot((mean(I5,1)),'color','c','linewidth',3)
                hold on
                plot((mean(Ith,1)),'--k')
                axis tight
                if mm==1
                    legend({'I_{L1}','I_{L2}','I_{par}','I_{kNN}','I_{KSG}','I_{ther.}'},'Location','northwest','Fontsize',12)
                end
                set(gca,'YScale','log')
            end
        end
    end
end



mm=0;
for N=[2 4 6]%1:9
    for k=5%1:10
        mm=mm+1;
        if rem(numel(DA{MODEL}.inf_1{N,k}),9)==0 & numel(DA{MODEL}.inf_1{N,k})~=0
            clear err
            for i=1:5
                %                 err(:,:,i)=reshape(DA{MODEL}.er{N,k}(:,i),9,[]);
            end
            
            err(:,:,1)=abs(reshape(DA{MODEL}.inf_1{N,k},9,[])-reshape(DA{MODEL}.inf_th{N,k},9,[]));
            err(:,:,2)=abs(reshape(DA{MODEL}.inf_2{N,k},9,[])-reshape(DA{MODEL}.inf_th{N,k},9,[]));
            err(:,:,3)=abs(reshape(DA{MODEL}.inf_3{N,k},9,[])-reshape(DA{MODEL}.inf_th{N,k},9,[]));
            err(:,:,4)=abs(reshape(DA{MODEL}.inf_4{N,k},9,[])-reshape(DA{MODEL}.inf_th{N,k},9,[]));
            err(:,:,5)=abs(reshape(DA{MODEL}.inf_5{N,k},9,[])-reshape(DA{MODEL}.inf_th{N,k},9,[]));
            
            I1=reshape(DA{MODEL}.inf_1{N,k},9,[]);
            I2=reshape(DA{MODEL}.inf_2{N,k},9,[]);
            I3=reshape(DA{MODEL}.inf_3{N,k},9,[]);
            I4=reshape(DA{MODEL}.inf_4{N,k},9,[]);
            I5=reshape(DA{MODEL}.inf_5{N,k},9,[]);
            Ith=reshape(DA{MODEL}.inf_th{N,k},9,[]);
            
            figure(150)
            subplot(1,3,mm)
            plot((mean(I1,2)),'color',[0 0 1],'linewidth',3)
            hold on
            plot((mean(I2,2)),'color',[1 0 1],'linewidth',3)
            hold on
            plot((mean(I3,2)),'color',[1 0 0],'linewidth',3)
            hold on
            plot((mean(I4,2)),'color',[0 1 0],'linewidth',3)
            hold on
            plot((mean(I5,2)),'color','c','linewidth',3)
            hold on
            plot((mean(Ith,2)),'--k')
            axis tight
            if mm==1
                legend({'I_{L1}','I_{L2}','I_{par}','I_{kNN}','I_{KSG}','I_{ther.}'},'Location','northwest','Fontsize',12)
            end
            %         figure(100)
            %         subplot(9,10,mm)
            %         semilogy((mean(I1,2)),'color',[0 0 1],'linewidth',3)
            %         hold on
            %         semilogy((mean(I2,2)),'color',[1 0 1],'linewidth',3)
            %         hold on
            %         semilogy((mean(I3,2)),'color',[1 0 0],'linewidth',3)
            %         hold on
            %         semilogy((mean(I4,2)),'color',[0 1 0],'linewidth',3)
            %         hold on
            %         semilogy((mean(I5,2)),'color','c','linewidth',3)
            %         hold on
            %         semilogy((mean(Ith,2)),'--k')
            %         axis tight
            %         if mm==1
            %         legend({'I_{L1}','I_{L2}','I_{par}','I_{kNN}','I_{KSG}','I_{ther.}'},'Location','northwest','Fontsize',12)
            %         end
            
            figure(250)
            subplot(1,3,mm)
            plot(squeeze(mean(err(:,:,1),2)),'color',[0 0 1],'linewidth',3)
            hold on
            plot(squeeze(mean(err(:,:,2),2)),'color',[1 0 1],'linewidth',3)
            hold on
            plot(squeeze(mean(err(:,:,3),2)),'color',[1 0 0],'linewidth',3)
            hold on
            plot(squeeze(mean(err(:,:,4),2)),'color',[0 1 0],'linewidth',3)
            hold on
            plot(squeeze(mean(err(:,:,5),2)),'color','c','linewidth',3)
            axis tight
            if mm==1
                legend({'I_{L1}','I_{L2}','I_{par}','I_{kNN}','I_{KSG}'},'Location','northeast','Fontsize',12)
            end
        end
    end
end


for N=1:9
    in1=reshape(DA{MODEL}.inf_1{N,5},9,[]);
    in2=reshape(DA{MODEL}.inf_2{N,5},9,[]);
    in3=reshape(DA{MODEL}.inf_3{N,5},9,[]);
    in4=reshape(DA{MODEL}.inf_4{N,5},9,[]);
    in5=reshape(DA{MODEL}.inf_5{N,5},9,[]);
    inth=real(reshape(DA{MODEL}.inf_th{N,5},9,[]));
    
    Me1(N)=mean(abs(in1(5,:)-inth(5,:)));
    Me2(N)=mean(abs(in2(5,:)-inth(5,:)));
    Me3(N)=mean(abs(in3(5,:)-inth(5,:)));
    Me4(N)=mean(abs(in4(5,:)-inth(5,:)));
    Me5(N)=mean(abs(in5(5,:)-inth(5,:)));
end

figure(350)
plot(Me1,'color',[0 0 1],'linewidth',3)
hold on
plot(Me2,'color',[1 0 1],'linewidth',3)
hold on
plot(Me3,'color',[1 0 0],'linewidth',3)
hold on
plot(Me4,'color',[0 1 0],'linewidth',3)
hold on
plot(Me5,'color','c','linewidth',3)
axis tight
if mm==1
    legend({'I_{L1}','I_{L2}','I_{par}','I_{kNN}','I_{KSG}'},'Location','northeast','Fontsize',12)
end


