% #' Evaluate the density of the transformation log likelihood estimator


function [pd_data,ccdf_data,ccdf2_data,pd_grid,ccdf_grid,ccdf2_grid,pd_points,ccdf_points,ccdf2_points,F,norma] = func_tll(lfit,Grid,points,data,fit,short,NORM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(Grid.X,2)==1
    kn=round((size(Grid.u,1)));
end
if size(Grid.X,2)==2
    kn=round(sqrt(size(Grid.u,1)));
end
if size(Grid.X,2)==3
    kn=round((size(Grid.u,1)^(1/3)));
end

if size(Grid.X,2)==1
    x1=unique(Grid.u(:,1));
    xd1=[diff(x1)];xd1=[xd1;xd1(1)];
end
if size(Grid.X,2)==2
    x1=unique(Grid.u(:,1));
    xd1=[diff(x1)];xd1=[xd1;xd1(1)];
    x2=unique(Grid.u(:,2));
    xd2=[diff(x2)];xd2=[xd2;xd2(1)];
end
if size(Grid.X,2)==3
    x1=unique(Grid.u(:,1));
    xd1=[diff(x1)];xd1=[xd1;xd1(1)];
    x2=unique(Grid.u(:,2));
    xd2=[diff(x2)];xd2=[xd2;xd2(1)];
    x3=unique(Grid.u(:,3));
    xd3=[diff(x3)];xd3=[xd3;xd3(1)];
end

if numel(x1)~=numel(xd1)
    error
end

if isfield(lfit,'Kergrid')
    
    if size(Grid.X,2)==1
        t1=(lfit.Kergrid);
        F1=griddedInterpolant(Grid.S,double(t1),'linear','none');
    end
    
    if size(Grid.X,2)==2
        t1=reshape(lfit.Kergrid,kn,kn);%./NORM;
        F1=scatteredInterpolant(Grid.S,double(t1(:)),'linear','none');
    end
    
    if size(Grid.X,2)==3
        error('correct this')
        t1=reshape(lfit.Kergrid,kn,kn,kn);
        I2=sum(xd2.*squeeze(sum(xd1.*permute(t1,[1 2 3]),1)),1);
        I2=reshape(I2,1,1,numel(I2));
        tcon=t1./I2;
        gg=Grid.S(1:kn*kn,1:2);
%         gg=Grid.u(1:kn*kn,1:2);
        for j=1:kn
            pp=tcon(:,:,j);
            F1{j}=scatteredInterpolant(gg,double(pp(:)),'linear','none');
        end
    end
    
else
    LF = loclik_fit(lfit.bw,data,Grid);
    if size(Grid.X,2)==1
        t1=reshape(LF.Kergrid,kn);
    end
    if size(Grid.X,2)==2
        t1=reshape(LF.Kergrid,kn,kn);
    end
    if size(Grid.X,2)==3
        t1=reshape(LF.Kergrid,kn,kn,kn);
    end
    F1=scatteredInterpolant(Grid.S,double(t1(:)),'linear','none');
end


if size(Grid.X,2)==1
    if short==0
        pd_grid=t1(:);
        tnorm=pd_grid;
        
        for n=1
            I2=sum(xd1.*tnorm,1);
            K=I2;
            tnorm=tnorm./K;
        end
        pd_grid=tnorm;
    end
end


%%%%% Normalisation
%%%%% f(S,T)/phi(S)phi(T)--> c(u,v)  

% P1=normpdf(s1);
% P2=normpdf(s2);
% PHI=P1*P2';
% t1=t1./PHI;

if short==1
    NUM=50;
else
    NUM=1000;
end

% t1=t1+1e-50;

% if max(t1(:))<1e-3
%     t1=NORM;%ones(kn,kn);
% end
    
if size(Grid.X,2)==2
    if  short~=1
        pd_grid=t1(:);
        t1=reshape(pd_grid,kn,kn);
        tnorm=t1./NORM;
        for n=1:NUM
%             if n==1
%             tnorm(tnorm<1e-300)=1e-300;
%             end

            I1=sum(xd2'.*tnorm,2);
            I2=sum(xd1.*tnorm,1);
            
%             I1=trapz(x2,tnorm,2);
%             I2=trapz(x1,tnorm,1);
            
            K=I1.*I2;
            tnorm=tnorm./K;
            
        end
            II=sum(xd1.*sum(xd2'.*tnorm,2));
%             II=trapz(x2,trapz(x1,tnorm,1));

            tnorm=tnorm/II;           
        pd_grid=tnorm;
    end
end


if short==2 & fit==0
        t1=pd_grid.*NORM;
        F1=scatteredInterpolant(Grid.S,double(t1(:)),'linear','none');
        short=1;
end

if size(Grid.X,2)==3
    error('check this before USE')
    for j=1:kn
        tc=squeeze(tcon(:,:,j));
        pd_grid=tc(:);
        t1=reshape(pd_grid,kn,kn);
        tnorm=t1;
        for n=1:1000
            I2=sum(xd1.*tnorm,1);
            I1=sum(xd2.*tnorm',1);
            K=I1'*I2;
            tnorm=tnorm./K;
            tnorm(K==0)=0;
            II=sum(sum(tnorm.*xd1).*xd2');
            tnorm=tnorm/II;
        end
        pd_grid_cond(:,:,j)=tnorm;
    end
    pd_grid=pd_grid_cond;
end

norma=1;


if short==0 | fit==-1
    
    if fit~=0
        if size(Grid.X,2)==1
%             F.pdf=griddedInterpolant(Grid.u,double(pd_grid(:)),'linear','none');
%             F.pdf=griddedInterpolant(Grid.u,double(pd_grid(:)),'linear','nearest');
            F.pdf=griddedInterpolant(Grid.S,double(pd_grid(:)),'linear','none');
        end
        if size(Grid.X,2)==2
%             F.pdf=scatteredInterpolant(Grid.u,double(pd_grid(:)),'linear','none');
%             F.pdf=scatteredInterpolant(Grid.u,double(pd_grid(:)),'linear','nearest');
            F.pdf=scatteredInterpolant(Grid.S,double(pd_grid(:)),'nearest','none');
        end
        if size(Grid.X,2)==3
            error('correct this')
%             gg=Grid.u(1:kn*kn,1:2);
            gg=Grid.X(1:kn*kn,1:2);
            for j=1:kn
                pp=pd_grid(:,:,j);
                F{j}.pdf=scatteredInterpolant(gg,double(pp(:)),'linear','nearest');
            end
        end
        
    else
        
        if size(Grid.X,2)==2 || size(Grid.X,2)==1
            F.pdf=lfit.F.pdf;
        end
        if size(Grid.X,2)==3
            for j=1:kn
                F{j}.pdf=lfit.F{j}.pdf;
            end
        end
        
    end
    
    
    if size(Grid.X,2)==2 || size(Grid.X,2)==1
        pd_points=F.pdf(points);
%         pd_data=F.pdf(data.u);
        pd_data=F.pdf(data.S);
        pd_data(pd_data<0)=0;
        pd_points(pd_points<0)=0;
        pd_points(isnan(pd_points))=0;
        pd_data(isnan(pd_data))=0;
    end
    
    if size(Grid.X,2)==3
        error('correct this')
        [~,pn3]=histc(points(:,3),x3);
        [~,dn3]=histc(data.u(:,3),x3);
%         [~,dn3]=histc(data.X(:,3),x3);
        pd_points=nan(size(points,1),1);
        pd_data=nan(size(data.u,1),1);
        for nn=1:max(pn3)
            if sum(pn3==nn)~=0
                pd_points(find(pn3==nn))=F{nn}.pdf(points(find(pn3==nn),1:2));
            end
            if sum(dn3==nn)~=0
                pd_data(find(dn3==nn))=F{nn}.pdf(data.u(find(dn3==nn),1:2));
%                 pd_data(find(dn3==nn))=F{nn}.pdf(data.X(find(dn3==nn),1:2));
            end
        end
        pd_data(pd_data<0)=eps;
        pd_points(pd_points<0)=eps;
        pd_points(isnan(pd_points))=eps;
        pd_data(isnan(pd_data))=eps;
    end
    
else
    
    if size(Grid.X,2)==2 | size(Grid.X,2)==1
        pd_grid=t1;
        pd_points=NaN;
        pd_data=F1(data.S);
%         pd_data=F1(data.u);
        pd_data(pd_data<0)=eps;
        pd_data(isnan(pd_data))=eps;
    end
    
    if size(Grid.X,2)==3
        error('correct this too')
        pd_grid=tcon;
        pd_points=NaN;
        [~,dn3]=histc(data.u(:,3),x3);
        for nn=1:max(dn3)
            if sum(dn3==nn)~=0
                pd_data(find(dn3==nn))=F1{nn}(data.S(find(dn3==nn),1:2));
%                 pd_data(find(dn3==nn))=F1{nn}(data.u(find(dn3==nn),1:2));
            end
        end
        
        %         pd_data=F1(data.X);
        
        pd_data(pd_data<0)=eps;
        pd_data(isnan(pd_data))=eps;
    end
    
end




if short==0
    t1=pd_grid;    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if size(Grid.X,2)==1
        h=t1*NaN;        
    end

    if size(Grid.X,2)==2
        h=t1*NaN;
        for i2=1:size(t1,2)
            for i1=1:size(t1,1)
                if i2==1
                    h(i1,1)=t1(i1,1)*xd2(1)/sum(t1(i1,1:end)'.*xd2(1:end));                    
                else
                    if sum(t1(i1,1:end)'.*xd2(1:end))~=0
                        h(i1,i2)=sum(t1(i1,1:i2)'.*xd2(1:i2))/sum(t1(i1,1:end)'.*xd2(1:end));
                    else
                        h(i1,i2)=0;
                    end
                end
            end
        end        
    end
    

    if size(Grid.X,2)==3
        error('check this before use')
        hh=nan(kn,kn,kn);
        h=nan(kn,kn);        
        for j=1:kn            
            t1=squeeze(pd_grid(:,:,j));            
            for i2=1:size(t1,2)
                for i1=1:size(t1,1)
                    if i2==1
                        h(i1,1)=t1(i1,1)*xd2(1);
                    else
                        if sum(t1(i1,1:end)'.*xd2(1:end))~=0
                            h(i1,i2)=sum(t1(i1,1:i2)'.*xd2(1:i2))/sum(t1(i1,1:end)'.*xd2(1:end));
                        else
                            h(i1,i2)=0;
                        end
                    end
                end
            end   
            hh(:,:,j)=h;            
        end
        h=hh;        
    end
    
    [~,GRID_u]= mk_grid(kn,'');

        h(h<=min(Grid.u(:)))=min(GRID_u(:));%1/(size(data.u,1)+1);
        h(h>=max(Grid.u(:)))=max(GRID_u(:));%1-1/(size(data.u,1)+1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% h2
    if size(Grid.X,2)==1
        h2=t1*NaN;        
    end

    if size(Grid.X,2)==2
        h2=t1*NaN;
        for i1=1:size(t1,1)
            for i2=1:size(t1,2)
                if i1==1
                    h2(1,i2)=t1(1,i2)*xd1(1)/sum(t1(1:end,i2).*xd1(1:end));
                else
                    if sum(t1(1:end,i2).*xd1(1:end))~=0
                        h2(i1,i2)=sum(t1(1:i1,i2).*xd1(1:i1))/sum(t1(1:end,i2).*xd1(1:end));
                    else
                        h2(i1,i2)=0;
                    end
                end
            end
        end
    end
    

%     if size(Grid.X,2)==3
%         error('check this before use')
%         hh2=nan(kn,kn,kn);
%         h=nan(kn,kn);        
%         for j=1:kn            
%             t1=squeeze(pd_grid(:,:,j));            
%             for i2=1:size(t1,2)
%                 for i1=1:size(t1,1)
%                     if i2==1
%                         h(i1,1)=t1(i1,1)*xd2(1);
%                     else
%                         if sum(t1(i1,1:end)'.*xd2(1:end))~=0
%                             h(i1,i2)=sum(t1(i1,1:i2)'.*xd2(1:i2))/sum(t1(i1,1:end)'.*xd2(1:end));
%                         else
%                             h(i1,i2)=0;
%                         end
%                     end
%                 end
%             end   
%             hh(:,:,j)=h;            
%         end
%         h=hh;        
%     end
    

        h2(h2<=min(Grid.u(:)))=min(GRID_u(:));%1/(size(data.u,1)+1);
        h2(h2>=max(Grid.u(:)))=max(GRID_u(:));%1-1/(size(data.u,1)+1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if fit~=0
        if size(Grid.X,2)==2
%             F.ccdf=scatteredInterpolant(Grid.u,double(h(:)),'linear','nearest');
            F.ccdf=scatteredInterpolant(Grid.S,double(h(:)),'linear','none');
            F.ccdf2=scatteredInterpolant(Grid.S,double(h2(:)),'linear','none');
%             F.ccdf=scatteredInterpolant(Grid.u,double(h(:)),'linear','nearest');
        end
        if size(Grid.X,2)==3
            error('correct this')
            gg=Grid.X(1:kn*kn,1:2);
            for j=1:kn
                pp=hh(:,:,j);
%                 F{j}.ccdf=scatteredInterpolant(gg,double(pp(:)),'linear','nearest');
                F{j}.ccdf=scatteredInterpolant(gg,double(pp(:)),'linear','nearest');
            end
        end
        
    else
        
        if size(Grid.X,2)==2
            F.ccdf=lfit.F.ccdf;
            F.ccdf2=lfit.F.ccdf2;
        end
        
        if size(Grid.X,2)==3
            for j=1:kn
                F{j}.ccdf=lfit.F{j}.ccdf;
            end
        end
        
    end
    
    
    if size(Grid.X,2)==1
        ccdf_points=NaN;
        ccdf_data=NaN;
        ccdf_grid=NaN;        
    end
    
    if size(Grid.X,2)==2
        ccdf_points=F.ccdf(points);
%         ccdf_data=F.ccdf(data.u);
        ccdf_data=F.ccdf(data.S);
        ccdf_grid=h;
        ccdf_points(isnan(ccdf_points))=0;
        ccdf_data(isnan(ccdf_data))=0;
        ccdf_data(ccdf_data<=min(GRID_u(:)))=min(GRID_u(:));%1/(size(data.u,1)+1);
        ccdf_data(ccdf_data>=max(GRID_u(:)))=max(GRID_u(:));%1-1/(size(data.u,1)+1);
        ccdf_points(ccdf_points<=min(GRID_u(:)))=min(GRID_u(:));%1/(size(data.u,1)+1);
        ccdf_points(ccdf_points>=max(GRID_u(:)))=max(GRID_u(:));%1-1/(size(data.u,1)+1);
        
        ccdf2_points=F.ccdf2(points);
%         ccdf_data=F.ccdf(data.u);
        ccdf2_data=F.ccdf2(data.S);
        ccdf2_grid=h2;
        ccdf2_points(isnan(ccdf2_points))=0;
        ccdf2_data(isnan(ccdf2_data))=0;
        ccdf2_data(ccdf2_data<=min(GRID_u(:)))=min(GRID_u(:));%1/(size(data.u,1)+1);
        ccdf2_data(ccdf2_data>=max(GRID_u(:)))=max(GRID_u(:));%1-1/(size(data.u,1)+1);
        ccdf2_points(ccdf2_points<=min(GRID_u(:)))=min(GRID_u(:));%1/(size(data.u,1)+1);
        ccdf2_points(ccdf2_points>=max(GRID_u(:)))=max(GRID_u(:));%1-1/(size(data.u,1)+1);
        
    end
    
    if size(Grid.X,2)==3
        error('correct this')
        [~,pn3]=histc(points(:,3),x3);
        [~,dn3]=histc(data.u(:,3),x3);
%         [~,dn3]=histc(data.X(:,3),x3);
        ccdf_points=nan(size(points,1),1);
        ccdf_data=nan(size(data.u,1),1);
        for nn=1:max(pn3)
            if sum(pn3==nn)~=0
                ccdf_points(find(pn3==nn))=F{nn}.ccdf(points(find(pn3==nn),1:2));
            end
            if sum(dn3==nn)~=0
                ccdf_data(find(dn3==nn))=F{nn}.ccdf(data.u(find(dn3==nn),1:2));
%                 ccdf_data(find(dn3==nn))=F{nn}.ccdf(data.X(find(dn3==nn),1:2));
            end
        end
        ccdf_points(isnan(ccdf_points))=0;
        ccdf_data(isnan(ccdf_data))=0;
        ccdf_grid=h;
        ccdf_data(ccdf_data<=0)=1/(size(data.u,1)+1);
        ccdf_data(ccdf_data>=1)=1-1/(size(data.u,1)+1);
        ccdf_points(ccdf_points<=0)=1/(size(data.u,1)+1);
        ccdf_points(ccdf_points>=1)=1-1/(size(data.u,1)+1);
        
    end
    
    
else
    ccdf_grid=NaN;
    ccdf_data=NaN;
    ccdf_points=NaN;
    F=NaN;
    ccdf2_grid=NaN;
    ccdf2_data=NaN;
    ccdf2_points=NaN;
 end





