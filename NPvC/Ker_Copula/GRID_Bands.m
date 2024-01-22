
function Grid=GRID_Bands(dat,knots,method)


data.u=dat;
data.S=norminv(data.u,0,1);    %%%%norminv
% qrs=linsolve(bw,zdata');

if size(dat,2)~=1
    [COEFF,data.X(:,1:2),~,~,~,mu] = pca(data.S(:,1:2));
    data.X(:,1:2)=data.S(:,1:2) * COEFF - repmat(mu,size(data.S(:,1:2),1),1) * COEFF;
else
    data.X=data.S;
end

if size(dat,2)==3
    data.X(:,3)=data.S(:,3);
end

if size(dat,2)==1
    [pts,GRID_u]= mk_grid(knots,method);
    GRID_S=norminv(GRID_u);
    GRID_X=GRID_S;
end

if size(dat,2)==2
    [pts,GRID_u]= mk_grid(knots,method);
    GRID_S=norminv(GRID_u);
%     GRID_X(:,1:2)=GRID_S*(COEFF);
    GRID_X(:,1:2)=GRID_S * COEFF - repmat(mu,size(GRID_S,1),1) * COEFF;
    
end

if size(dat,2)==3
    [pts,GRID_u]= mk_grid(knots,method);
    [pts3,G]= mk_grid(knots);
    if size(pts,2)==1
        [ex1 ex2 ex3] = meshgrid(pts,pts,pts3(:,1));
    end
    if size(pts,2)==2
        [ex1 ex2 ex3] = meshgrid(pts(:,1),pts(:,2),pts3(:,1));
    end
    [ex1 ex2 ex3] = meshgrid(pts(:,1),pts(:,2),pts3(:,1));
    
    GRID_u=[ex1(:) ex2(:) ex3(:)];
    GRID_S=norminv(GRID_u);
%     GRID_X(:,1:2)=GRID_S(:,1:2)*(COEFF);
    GRID_X(:,1:2)=GRID_S(:,1:2) * COEFF - repmat(mu,size(GRID_S(:,1:2),1),1) * COEFF;
    
    GRID_X(:,3)=GRID_S(:,3);
end

Grid.S=double(GRID_S);
Grid.X=double(GRID_X);
Grid.u=double(GRID_u);

