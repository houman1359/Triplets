function Ker_grid = DenseNaive(B,data,grid)


% for i=1:size(grid,2)
% [a]=DenseNaive_1D(grid(:,i),B,data);
% Ker_grid{1}(i)=a(1);
% Ker_grid{2}(i)=a(2);
% Ker_grid{3}(i)=a(3);
% Ker_grid{4}(i)=a(4);
% Ker_grid{5}(i)=a(5);
% end


if 1==1%size(data,1)==2

M=10;%%% this was 10 I changed on May 2019 to 5

%%%% for the mex file:
S=1:size(grid,2);
Sm=size(grid,2)/M;

if rem(size(S,2),M)~=0
error('choose knots a product of 5!  or change DneseNaive code')    
end

if size(data,1)==1
Ker_gridS=zeros(3,size(grid,2));
end
if size(data,1)==2
Ker_gridS=zeros(5,size(grid,2));
end
if size(data,1)==3
Ker_gridS=zeros(7,size(grid,2));
end


for i=1:M
    SS=(i-1)*Sm+1:i*Sm;
if i==1
    if size(data,1)==1
    Ker_gridS = DenseNaiveMat1D(B,data,grid(:,SS));    
    end
    
    if size(data,1)==2
        if ispc
%             Ker_gridS = DenseNaiveMatwin(B,data,grid(:,SS));
            Ker_gridS =NPC_DenseNaiveMatwin(B,data,grid(:,SS));%%%%alessandro
        else
%             Ker_gridS = DenseNaiveMat(B,data,grid(:,SS));     
            Ker_gridS =NPC_DenseNaiveMatwin(B,data,grid(:,SS));%%%%alessandro     %%%%%% CHANGED on May2020
        end
    end
    
    if size(data,1)==3
    Ker_gridS = DenseNaiveMat3D(B,data,grid(:,SS));    
    end
else
    if size(data,1)==1
    Ker_gridS = cat(2,Ker_gridS,DenseNaiveMat1D(B,data,grid(:,SS)));
    end
    if size(data,1)==2
    
        if ispc
%             Ker_gridS = cat(2,Ker_gridS,DenseNaiveMatwin(B,data,grid(:,SS)));
            Ker_gridS = cat(2,Ker_gridS,NPC_DenseNaiveMatwin(B,data,grid(:,SS)));
        else
%             Ker_gridS = cat(2,Ker_gridS,DenseNaiveMat(B,data,grid(:,SS)));
            Ker_gridS = cat(2,Ker_gridS,NPC_DenseNaiveMatwin(B,data,grid(:,SS)));
        end
        
    end
    if size(data,1)==3
    Ker_gridS = cat(2,Ker_gridS,DenseNaiveMat3D(B,data,grid(:,SS)));
    end
end
end

% Ker_grid = Ker_gridS;%DenseNaiveMat(B,data,grid);
% Ker_grid = DenseNaiveMat(B,data,grid);
    
if size(data,1)==1
Ker_grid = mat2cell(Ker_gridS,ones(1,3),size(grid,2))';
end
if size(data,1)==2
Ker_grid = mat2cell(Ker_gridS,ones(1,5),size(grid,2))';
end
if size(data,1)==3
Ker_grid = mat2cell(Ker_gridS,ones(1,7),size(grid,2))';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

