
%%%%%% Fitting kernel vine copula and computing the joint density funciton


function [p,pdata,Copula,PDF,p_copula] = Fit_vCopula(vine,x,TL,METH_fit,fit,lfit,sort_m,vineest,knots)

COP_BR=0;
if fit==-2
    COP_BR=1;
    fit=-1;
end



COP_BR=1;
if fit==-3
    COP_BR=0;
    fit=-1;
end

%  parallel.defaultClusterProfile('local')
% ClusterInfo.setQueueName('mpi')
% ClusterInfo.setWallTime('24:00:00')
% ClusterInfo.setMemUsage('500M')

% parpool('SpmdEnabled',false)


% c=parcluster;
% c.AdditionalProperties.WallTime= '12:00:00';
% c.AdditionalProperties.QueueName='mpi';
% c.AdditionalProperties.AddisionalSubmitArgs= '--mem-per-cpu=2G'
% c.saveprofile;

if nargin < 9
    knots = lfit{1,2}.fit.knots;%sqrt(size(lfit{1,2}.fit.F.pdf.Points,1));%size(lfit{1,2}.C_grid,1);
    % for m=1:size(lfit,1)
    % knots(m)=size(lfit{1,m}.C_grid,2);
    % end
end

d = length(vine.margins);

if nargin < 1
    error('Fit_vCopula: Usage [p,logp] = mixedvinepdf(vine,x)');
end
if ~isstruct(vine)
    error('Fit_vCopula: Argument "vine" must be a struct');
end
if ~ismatrix(x)
    error('Fit_vCopula: Argument "x" must be a matrix');
end
if size(x,2) ~= d
    error('Fit_vCopula: Second dimension of u must match number of margins');
end
for i=1:d
    PDF{i}=[];
end
cases = size(x,1);
newnode = true(d);
Fp = zeros(cases,d,d);
% logf(s,i,j) represents the s-th sample of log(f(j|i))
logf = zeros(cases,d);
logfdata=zeros(size(vine.theta{1,1},1),d);
tr=[];

no=0;
if fit==3
    fit=0;
    no=1;
end


[pts,GRID_u]= mk_grid(knots,'');

for i = 1:d  %%%parfor
    
    if fit~=0
        parG{i}.fit=fit;
        parS{i}.fit=fit;
        parT{i}.fit=fit;
        parP{i}.fit=fit;
        
        parG{i}.max=max(GRID_u(:));
        parG{i}.min=min(GRID_u(:));
        parT{i}.max=max(GRID_u(:));
        parT{i}.min=min(GRID_u(:));
    else
        parG{i}.fit=0;
        parG{i}.p=lfit{i,1}.MarginG.p;
        parG{i}.s=lfit{i,1}.MarginG.s;
        parS{i}.fit=0;
        parS{i}.p=lfit{i,1}.MarginS.p;
        parS{i}.s=lfit{i,1}.MarginS.s;
        parT{i}.fit=0;
        parT{i}.p=lfit{i,1}.MarginG.p;
        parT{i}.s=lfit{i,1}.MarginG.s;
        parP{i}.fit=0;
        parP{i}.p=lfit{i,1}.MarginS.p;
        parP{i}.s=lfit{i,1}.MarginS.s;
        
        parG{i}.max=max(GRID_u(:));
        parG{i}.min=min(GRID_u(:));
        parT{i}.max=max(GRID_u(:));
        parT{i}.min=min(GRID_u(:));
    end
    
    if isfield(vine.margins{i},'MAX')
        parS{i}.MAX=vine.margins{i}.MAX;
        parP{i}.MAX=vine.margins{i}.MAX;
    else
        parS{i}.MAX=[];
        parP{i}.MAX=[];
    end
    
end


for i = 1:d  %%%parfor
    %     disp(['pdf ',num2str(i)])
    
    %     G0{i} = margincdf(vine.margins{i},x(:,i));
    [G0{i},Mar_G{i}]=kernelcdf(vine.margins{i}.ker,x(:,i),parG{i});
    %     S0{i} = marginpdf(vine.margins{i},x(:,i));
    [S0{i},Mar_S{i}]=kernelpdf(vine.margins{i}.ker,vine.margins{i}.iscont,x(:,i),parS{i});
    %     T0{i} = margincdf(vine.margins{i},vine.margins{i}.ker);
    [T0{i},Mar_T{i}]=kernelcdf(vine.margins{i}.ker,vine.margins{i}.ker,parT{i});
    %     P0{i} = marginpdf(vine.margins{i},vine.margins{i}.ker);
    [P0{i},Mar_P{i}]=kernelpdf(vine.margins{i}.ker,vine.margins{i}.iscont,vine.margins{i}.ker,parP{i});
    if fit==0 | fit==-1
        PDF{i}.margins=S0{i};
    end
end

for i = 1:d
    if fit~=-4
    Fp(:,1,i)=G0{i};
    else
    Fp(:,1,i)=x(:,i);
    end
    logf(:,1)=logf(:,1)+log(S0{i});
    logfdata(:,1)=logfdata(:,1)+log(P0{i});
    vine.theta{1,i}=T0{i};
end
if fit==-4
    fit=-1;
end

clear G0 T0 S0

if COP_BR
    switch lower(vine.type)
        case 'c-vine'
            
            DD=min([TL d-1]);
            
            for tr=1:DD
                
                if ~strcmp(vine.families{tr,2},'ind')
                    
                    clear co logftr G1p G2p G3p G4p CopC Copf Copul CCC
                    
                    if fit==1
                        
                        for str=2:d- (tr-1)    %%%parfor
                            tic
                            [G1p{str},G2p{str},CopC{str},Copf{str},G3p{str},G4p{str},Copul{str}]=KernelCop(knots,[vine.theta{tr,1} vine.theta{tr,str}],[vine.theta{tr,1}(1) vine.theta{tr,str}(1)],vine.METH{tr,str},METH_fit,0,1);
                            
                            if vine.condition==1 && tr>1
                                [bw3]=KernelCop1D(knots,[vine.theta{tr-1,1}],[vine.theta{tr-1,1}(1)],METH,Copul{str},1);
                                Copul{str}.fit.bw=[Copul{str}.fit.bw bw3*5];
                                [G1p{str},G2p{str},CopC{str},Copf{str},G3p{str},G4p{str},Copul{str}]=KernelCop3D(knots,[vine.theta{tr,1} vine.theta{tr,str} vine.theta{tr-1,1}],[vine.theta{tr,1}(1) vine.theta{tr,str}(1) vine.theta{tr-1,1}(1)],METH,Copul{str},-1,Copul{str}.fit.bw);
                            end
                            logftr{str}=log(G2p{str});
                            logfda{str}=log(G4p{str});
                            disp(['Fit copula -> ',num2str(tr),',',num2str(str),' , Time= ',num2str(toc)])
                            
                        end
                        
                    else
                        
                        for str=2:d- (tr-1)   %%%parfor
                            Lfit=lfit{tr,str}.fit;
                            if vine.condition~=1 || tr==1
                                [G1p{str},G2p{str},CopC{str},Copf{str},G3p{str},G4p{str},LF{str}]=KernelCop(knots,[vine.theta{tr,1} vine.theta{tr,str}],[Fp(:,tr,1) Fp(:,tr,str)],vine.METH{tr,str},METH_fit,Lfit,fit);
                                
                                
                                % u1=unique(GRID_u);
                                % figure
                                % surf(u1,u1,Copf{str})
                                % hold on
                                % plot3(vine.theta{tr,1}(618),vine.theta{tr,str}(618),G4p{str}(618),'O')
                                % hold on
                                % plot3(vine.theta{tr,str}(618),vine.theta{tr,1}(618),G4p{str}(618),'O')
                                % % colorbar
                                %                         figure
                                %                         gg=reshape(G2p{str},200,[]);
                                %                         plot(gg(:,618),'-O')
                                % title([num2str(tr),' , ',num2str(str),',',num2str(G4p{str}(618))])
                                
                            else
                                [G1p{str},G2p{str},CopC{str},Copf{str},G3p{str},G4p{str},LF{str}]=KernelCop3D(knots,[vine.theta{tr,1} vine.theta{tr,str} vine.theta{tr-1,1}],[Fp(:,tr,1) Fp(:,tr,str) Fp(:,tr-1,str)],vine.METH{tr,str},Lfit,fit,lfit{tr,str}.fit.bw);
                            end
                            logftr{str}=log(G2p{str});
                            logfda{str}=log(G4p{str});
                            if no~=1
                                %                         disp(['Eval copula -> ',num2str(tr),',',num2str(str)])
                            end
                        end
                    end
                    
                    if fit==0
                        if vine.condition==1
                            error  %% this case is not needed now. i might develope later
                        end
                    end
                    
                    for str=2:d- (tr-1)
                        if fit~=0
                            par1{str}.fit=1;
                            par2{str}.fit=1;
                            
                            par1{str}.max=max(GRID_u(:));
                            par1{str}.min=min(GRID_u(:));
                            par2{str}.max=max(GRID_u(:));
                            par2{str}.min=min(GRID_u(:));
                        else
                            par1{str}.fit=0;
                            par1{str}.p=lfit{tr,str}.Margin1.p;
                            par1{str}.s=lfit{tr,str}.Margin1.s;
                            par2{str}.fit=0;
                            par2{str}.p=lfit{tr,str}.Margin2.p;
                            par2{str}.s=lfit{tr,str}.Margin2.s;
                            
                            par1{str}.max=max(GRID_u(:));
                            par1{str}.min=min(GRID_u(:));
                            par2{str}.max=max(GRID_u(:));
                            par2{str}.min=min(GRID_u(:));
                        end
                    end
                    
                    
                    
                    clear Y1 Y2 Y3 Y4
                    for str=2:d- (tr-1)
                        Y1(:,str)=NaN*zeros(size(G1p{str}{1}));
                        Y2(:,str)=NaN*zeros(size(G3p{str}{1}));
                        Y3(:,str)=NaN*zeros(size(G4p{str}));
                        
                        if vine.condition==1 && tr>1
                            %                     Gr=GRID_Bands(vine.theta{tr,1},knots,1);
                            %                     [~,gg1]=histc(G1p{str},Gr.u);
                            %                     [~,gg3]=histc(G3p{str},Gr.u);
                            %                     gg1=gg1+1;gg3=gg3+1;
                            %                     for jc=1:size(Gr.u,1)
                            %                         if sum(gg1==jc)~=0 & sum(gg3==jc)~=0
                            %                             [Y1(gg1==jc,str),Mar_1{str}] = kernelcdf(G3p{str}(gg3==jc),G1p{str}(gg1==jc),par1{str});
                            %                             [Y2(gg3==jc,str),Mar_2{str}] = kernelcdf(G3p{str}(gg3==jc),G3p{str}(gg3==jc),par2{str});
                            %                             Y3(gg3==jc,str) = G4p{str}(gg3==jc);
                            %                         else
                            %                             Y1(:,str) = G1p{str};
                            %                             Y2(:,str) = G3p{str};
                            %                             Y3(:,str) = G4p{str};
                            %                         end
                            %                     end
                            
                            
                            [Y1(:,str),Mar_1{str}] = kernelcdf(G3p{str},G1p{str},par1{str});
                            [Y2(:,str),Mar_2{str}] = kernelcdf(G3p{str},G3p{str},par2{str});
                            Y3(:,str) = G4p{str};
                            
                            %                     Y1(:,str) = G1p{str};
                            %                     Y2(:,str) = G3p{str};
                            %                     Y3(:,str) = G4p{str};
                        else
                            [Y1(:,str),Mar_1{str}] = kernelcdf(G3p{str}{1},G1p{str}{1},par1{str});
                            [Y2(:,str),Mar_2{str}] = kernelcdf(G3p{str}{1},G3p{str}{1},par2{str});
                            Y3(:,str) = G4p{str};
                            %                                         Y1(:,str) = G1p{str};
                            %                                         Y2(:,str) = G3p{str};
                            %                                         Y3(:,str) = G4p{str};
                        end
                    end
                    
                    
                    if strcmp(sort_m,'sort') && fit==1
                        if tr~=DD
                            XX=Y2;
                            clear C
                            for i=1:d- tr
                                for j=1:d- tr
                                    if i~=j
                                        C(i,j)=abs(corr(XX(:,i+1),XX(:,j+1),'type','kendall'));
                                    else
                                        C(i,j)=NaN;
                                    end
                                end
                            end
                            
                            clear ord
                            [~,ord(1)]=max(nanmean(C));
                            for i=2:d- tr-1
                                ss=setdiff(1:d- tr,ord);
                                [m n]=max(C(ord(i-1),ss));
                                ord(i)=ss(n);
                            end
                            ord(d- tr)=setdiff(1:d- tr,ord);
                            ord=ord+1;
                            Y2(:,2:numel(ord)+1)=Y2(:,ord);
                        else
                            ord=2;
                        end
                    end
                    
                    if strcmp(sort_m,'rand')
                        ord=1:size(Y2,2)-1;
                    end
                    
                    
                    if strcmp(sort_m,'sort') & fit==0
                        Y2(:,2:d- (tr-1))=Y2(:,lfit{tr,2}.ord);
                    end
                    
                    for str=2:d- (tr-1)
                        
                        if fit==1
                            Copula{tr,str}=Copul{str};
                            Copula{tr,str}.fit.F='not saved here';
                            Copula{tr,str}.fit.knots=knots;
                            Copula{tr,str}.C_grid='not saved here';%CopC{str};
                            Copula{tr,str}.f_grid='not saved here';%Copf{str};
                            Copula{tr,str}.ord=ord;
                            %                     Copula{tr,str}.Margin1=Mar_1{str};
                            %                     Copula{tr,str}.Margin2=Mar_2{str};
                        elseif fit==0
                            Copula=NaN;
                        elseif fit==-1
                            Copula{tr,str}=LF{str};
                            Copula{tr,str}.fit=lfit{tr,str}.fit;
                            Copula{tr,str}.fit.F='not saved here';%LF{str}.fit.F;
                            Copula{tr,str}.C_grid=CopC{str};
                            Copula{tr,str}.f_grid=Copf{str};
                            Copula{tr,str}.ord=ord;
                            Copula{tr,str}.Margin1=[];%Mar_1{str};
                            Copula{tr,str}.Margin2=[];%Mar_2{str};
                            Copula{tr,str}.knots=size(Copula{tr,str}.f_grid,1);
                            
                            if COP_BR==1
                                Copula{tr,str}.co_br=G2p{str};
                            end
                            
                        end
                        
                        Fp(:,tr+1,str-1) = Y1(:,str);
                        vine.theta{tr+1,str-1} = Y2(:,str);
                        vine.c{tr+1,str-1} = Y3(:,str);
                        logf(:,tr+1) = logf(:,tr+1) + logftr{str};
                        logfdata(:,tr+1) = logfdata(:,tr+1) + logfda{str};
                        
                    end
                    if tr==DD
                        %             disp(['------------------ tree = ' , num2str(tr)])
                    end
                    
                else
                    
                    %%%%%%% puy the independent copula here
                    
                    
                end
            end
            
            logp = logf(:,1);
            logp_copula = 0;
            for i = 2:d
                if ~strcmp(vine.families{i-1,2},'ind')
                    logp = logp + logf(:,i);
                    logp_copula = logp_copula + logf(:,i);
                end
            end
            
            logdata = logfdata(:,1);
            for i = 2:d
                if ~strcmp(vine.families{i-1,2},'ind')
                    logdata = logdata + logfdata(:,i);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        case 'd-vine'
            
            DD=min([TL d-1]);
            
            for tr=1:DD
                if ~strcmp(vine.families{tr,2},'ind')
                    
                    clear co logftr G1p G2p G3p G4p CopC Copf Copul CCC
                    
                    if fit==1
                        
                        for str=2:d- (tr-1)    %%%parfor
                            tic
                            
%                             [G1p{str},G2p{str},CopC{str},Copf{str},G3p{str},G4p{str},Copul{str}]=KernelCop(knots,[vine.theta{tr,str-1} vine.theta{tr,str}],[vine.theta{tr,str-1}(1) vine.theta{tr,str}(1)],vine.METH{tr,str},METH_fit,0,1);
                        if tr==1
                            vine.theta2{tr,str}=vine.theta{tr,str};
                        end
                              [G1p{str},G2p{str},CopC{str},Copf{str},G3p{str},G4p{str},Copul{str}]=KernelCop(knots,[vine.theta{tr,str-1} vine.theta2{tr,str}],[vine.theta{tr,str-1}(1) vine.theta2{tr,str}(1)],vine.METH{tr,str},METH_fit,0,1);
                           
                            if vine.condition==1 && tr>1
                                [bw3]=KernelCop1D(knots,[vine.theta{tr-1,str-1}],[vine.theta{tr-1,str-1}(1)],METH,Copul{str},1);
                                Copul{str}.fit.bw=[Copul{str}.fit.bw bw3*5];
                                [G1p{str},G2p{str},CopC{str},Copf{str},G3p{str},G4p{str},Copul{str}]=KernelCop3D(knots,[vine.theta{tr,str-1} vine.theta{tr,str} vine.theta{tr-1,str-1}],[vine.theta{tr,str-1}(1) vine.theta{tr,str}(1) vine.theta{tr-1,str-1}(1)],METH,Copul{str},-1,Copul{str}.fit.bw);
                            end
                            logftr{str}=log(G2p{str});
                            logfda{str}=log(G4p{str});
                            disp(['Fit copula -> ',num2str(tr),',',num2str(str),' , Time= ',num2str(toc)])
                        end
                        
                    else
                        
                        for str=2:d- (tr-1)   %%%parfor
                            
                            if tr==1
                            vine.theta2{tr,str}=vine.theta{tr,str};    
                            Fp2(:,tr,str)=Fp(:,tr,str);
                            end
                        end
                        
                        for str=2:d- (tr-1)   %%%parfor
                                                    
                            Lfit=lfit{tr,str}.fit;
                            if vine.condition~=1 || tr==1
%                                 [G1p{str},G2p{str},CopC{str},Copf{str},G3p{str},G4p{str},LF{str}]=KernelCop(knots,[vine.theta{tr,str-1} vine.theta{tr,str}],[Fp(:,tr,str-1) Fp(:,tr,str)],vine.METH{tr,str},METH_fit,Lfit,fit);
                                [G1p{str},G2p{str},CopC{str},Copf{str},G3p{str},G4p{str},LF{str}]=KernelCop(knots,[vine.theta{tr,str-1} vine.theta2{tr,str}],[Fp(:,tr,str-1) Fp2(:,tr,str)],vine.METH{tr,str},METH_fit,Lfit,fit);
                                
                                
                                % u1=unique(GRID_u);
                                % figure
                                % surf(u1,u1,Copf{str})
                                % hold on
                                % plot3(vine.theta{tr,1}(618),vine.theta{tr,str}(618),G4p{str}(618),'O')
                                % hold on
                                % plot3(vine.theta{tr,str}(618),vine.theta{tr,1}(618),G4p{str}(618),'O')
                                % % colorbar
                                %                         figure
                                %                         gg=reshape(G2p{str},200,[]);
                                %                         plot(gg(:,618),'-O')
                                % title([num2str(tr),' , ',num2str(str),',',num2str(G4p{str}(618))])
                                
                            else
                                [G1p{str},G2p{str},CopC{str},Copf{str},G3p{str},G4p{str},LF{str}]=KernelCop3D(knots,[vine.theta{tr,str-1} vine.theta{tr,str} vine.theta{tr-1,str-1}],[Fp(:,tr,str-1) Fp(:,tr,str) Fp(:,tr-1,str)],vine.METH{tr,str},Lfit,fit,lfit{tr,str}.fit.bw);
                            end
                            logftr{str}=log(G2p{str});
                            logfda{str}=log(G4p{str});
                            if no~=1
                                %                         disp(['Eval copula -> ',num2str(tr),',',num2str(str)])
                            end
                        end
                    end
                    
                    if fit==0
                        if vine.condition==1
                            error  %% this case is not needed now. i might develope later
                        end
                    end
                    
                    for str=2:d- (tr-1)
                        if fit~=0
                            par1{str}.fit=1;
                            par2{str}.fit=1;
                            
                            par1{str}.max=max(GRID_u(:));
                            par1{str}.min=min(GRID_u(:));
                            par2{str}.max=max(GRID_u(:));
                            par2{str}.min=min(GRID_u(:));
                        else
                            par1{str}.fit=0;
                            par1{str}.p=lfit{tr,str}.Margin1.p;
                            par1{str}.s=lfit{tr,str}.Margin1.s;
                            par2{str}.fit=0;
                            par2{str}.p=lfit{tr,str}.Margin2.p;
                            par2{str}.s=lfit{tr,str}.Margin2.s;
                            
                            par1{str}.max=max(GRID_u(:));
                            par1{str}.min=min(GRID_u(:));
                            par2{str}.max=max(GRID_u(:));
                            par2{str}.min=min(GRID_u(:));
                        end
                    end
                    
                    clear Y1 Y2 Y3 Y4
                    for str=2:d- (tr-1)
                        Y1(:,str)=NaN*zeros(size(G1p{str}{1}));
                        Y2(:,str)=NaN*zeros(size(G3p{str}{1}));
                        Y3(:,str)=NaN*zeros(size(G4p{str}));
                        Z1(:,str)=NaN*zeros(size(G1p{str}{2}));
                        Z2(:,str)=NaN*zeros(size(G3p{str}{2}));
                        
                        if vine.condition==1 && tr>1
                            %                     Gr=GRID_Bands(vine.theta{tr,1},knots,1);
                            %                     [~,gg1]=histc(G1p{str},Gr.u);
                            %                     [~,gg3]=histc(G3p{str},Gr.u);
                            %                     gg1=gg1+1;gg3=gg3+1;
                            %                     for jc=1:size(Gr.u,1)
                            %                         if sum(gg1==jc)~=0 & sum(gg3==jc)~=0
                            %                             [Y1(gg1==jc,str),Mar_1{str}] = kernelcdf(G3p{str}(gg3==jc),G1p{str}(gg1==jc),par1{str});
                            %                             [Y2(gg3==jc,str),Mar_2{str}] = kernelcdf(G3p{str}(gg3==jc),G3p{str}(gg3==jc),par2{str});
                            %                             Y3(gg3==jc,str) = G4p{str}(gg3==jc);
                            %                         else
                            %                             Y1(:,str) = G1p{str};
                            %                             Y2(:,str) = G3p{str};
                            %                             Y3(:,str) = G4p{str};
                            %                         end
                            %                     end
                            
                            [Y1(:,str),Mar_1{str}] = kernelcdf(G3p{str}{1},G1p{str}{1},par1{str});
                            [Y2(:,str),Mar_2{str}] = kernelcdf(G3p{str}{1},G3p{str}{1},par2{str});
                            [Z1(:,str),zMar_1{str}] = kernelcdf(G3p{str}{2},G1p{str}{2},par1{str});
                            [Z2(:,str),zMar_2{str}] = kernelcdf(G3p{str}{2},G3p{str}{2},par2{str});
                            Y3(:,str) = G4p{str};
                            
                            %                     Y1(:,str) = G1p{str};
                            %                     Y2(:,str) = G3p{str};
                            %                     Y3(:,str) = G4p{str};
                        else
                            [Y1(:,str),Mar_1{str}] = kernelcdf(G3p{str}{1},G1p{str}{1},par1{str});
                            [Y2(:,str),Mar_2{str}] = kernelcdf(G3p{str}{1},G3p{str}{1},par2{str});
                            [Z1(:,str),zMar_1{str}] = kernelcdf(G3p{str}{2},G1p{str}{2},par1{str});
                            [Z2(:,str),zMar_2{str}] = kernelcdf(G3p{str}{2},G3p{str}{2},par2{str});
                            Y3(:,str) = G4p{str};
                            %                                         Y1(:,str) = G1p{str};
                            %                                         Y2(:,str) = G3p{str};
                            %                                         Y3(:,str) = G4p{str};
                        end
                    end
                    
%                     if strcmp(sort_m,'sort') && fit==1
%                         if tr~=DD
%                             XX=Y2;
%                             clear C
%                             for i=1:d- tr
%                                 for j=1:d- tr
%                                     if i~=j
%                                         C(i,j)=abs(corr(XX(:,i+1),XX(:,j+1),'type','kendall'));
%                                     else
%                                         C(i,j)=NaN;
%                                     end
%                                 end
%                             end
%                             
%                             clear ord
%                             [~,ord(1)]=max(nanmean(C));
%                             for i=2:d- tr-1
%                                 ss=setdiff(1:d- tr,ord);
%                                 [m n]=max(C(ord(i-1),ss));
%                                 ord(i)=ss(n);
%                             end
%                             ord(d- tr)=setdiff(1:d- tr,ord);
%                             ord=ord+1;
%                             Y2(:,2:numel(ord)+1)=Y2(:,ord);
%                         else
%                             ord=2;
%                         end
%                     end
                    
                    if strcmp(sort_m,'rand')
                        ord=1:size(Y2,2)-1;
                    end
                    
                    
                    if strcmp(sort_m,'sort') & fit==0
                        Y2(:,2:d- (tr-1))=Y2(:,lfit{tr,2}.ord);
                    end
                    
                    for str=2:d- (tr-1)
                        
                        if fit==1
                            Copula{tr,str}=Copul{str};
                            Copula{tr,str}.fit.F='not saved here';
                            Copula{tr,str}.fit.knots=knots;
                            Copula{tr,str}.C_grid='not saved here';%CopC{str};
                            Copula{tr,str}.f_grid='not saved here';%Copf{str};
                            Copula{tr,str}.ord=ord;
                            %                     Copula{tr,str}.Margin1=Mar_1{str};
                            %                     Copula{tr,str}.Margin2=Mar_2{str};
                        elseif fit==0
                            Copula=NaN;
                        elseif fit==-1
                            Copula{tr,str}=LF{str};
                            Copula{tr,str}.fit=lfit{tr,str}.fit;
                            Copula{tr,str}.fit.F='not saved here';%LF{str}.fit.F;
                            Copula{tr,str}.C_grid=CopC{str};
                            Copula{tr,str}.f_grid=Copf{str};
                            Copula{tr,str}.ord=ord;
                            Copula{tr,str}.Margin1=[];%Mar_1{str};
                            Copula{tr,str}.Margin2=[];%Mar_2{str};
                            Copula{tr,str}.knots=size(Copula{tr,str}.f_grid,1);
                            
                            if COP_BR==1
                                Copula{tr,str}.co_br=G2p{str};
                            end
                            
                        end
                        
                        Fp(:,tr+1,str-1) = Y1(:,str);
                        vine.theta{tr+1,str-1} = Y2(:,str);
                        Fp2(:,tr+1,str-1) = Z1(:,str);
                        vine.theta2{tr+1,str-1} = Z2(:,str);

                        vine.c{tr+1,str-1} = Y3(:,str);
                        logf(:,tr+1) = logf(:,tr+1) + logftr{str};
                        logfdata(:,tr+1) = logfdata(:,tr+1) + logfda{str};
                        
                    end
                    if tr==DD
                        %             disp(['------------------ tree = ' , num2str(tr)])
                    end
                end
            end
            
            logp = logf(:,1);
            logp_copula = 0;
            for i = 2:d
                if ~strcmp(vine.families{i-1,2},'ind')
                    logp = logp + logf(:,i);
                    logp_copula = logp_copula + logf(:,i);
                end
            end
            
            logdata = logfdata(:,1);
            for i = 2:d
                if ~strcmp(vine.families{i-1,2},'ind')
                    logdata = logdata + logfdata(:,i);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        otherwise
            error(['Unknown vine type "' vine.type '"']);
    end
    
    
    if fit==-1
        for i=1:d
            Copula{i,1}.MarginG=Mar_G{i};
            Copula{i,1}.MarginS=Mar_S{i};
        end
    end
    
    
    if fit==1
        [nor]=1;%densitynorm(vine,Copula,5e-2,50000,vineest);
        Copula{1,1}.norm=nor;
    elseif fit==0
        if no==1
            nor=1;
        else
            nor=lfit{1,1}.norm;%*prod(abs([min(vine.range')-max(vine.range')]));
        end
    elseif fit==-1
        Copula{1,1}.norm=1;%lfit{1,1}.norm;
        nor=Copula{1,1}.norm;
    end
    
    
    
    
    logp = reshape(logp,[cases 1]);
    % Correct numerical inaccuracies
    logp = real(logp);
    % Massless intervals can result in NaNs; set to 0
    logp(isnan(logp)) = -inf;
    p = exp(round(logp,6));%/nor;
    
    logp_copula = reshape(logp_copula,[],1);
    % Correct numerical inaccuracies
    logp_copula = real(logp_copula);
    % Massless intervals can result in NaNs; set to 0
    logp_copula(isnan(logp_copula)) = -inf;
    p_copula = exp(round(logp_copula,6));%/nor;
    
    
    if ~isempty(tr)
%         logdata = reshape(logdata,[size(logp,1) 1]);
        logdata = reshape(logdata,[],1);
        % Correct numerical inaccuracies
        logdata = real(logdata);
        % Massless intervals can result in NaNs; set to 0
        logdata(isnan(logdata)) = -inf;
        pdata = exp(round(logdata,6));%/nor;
    else
        pdata=[];
    end
    
else
    
    p=NaN;
    pdata=NaN;
    Copula=NaN;
    p_copula=NaN;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% we can define the copulas as a list C_{1:i}= {C{tr=1,str=1},......ith copula}
%%%%%% then the copula density can be defined as the max(Prod C_{1:i})

% delete(gcp('nocreate'))

end





