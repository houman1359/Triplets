function u=kerncopccdfinv(copula,v,method)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(copula,'knots')
knots = copula.knots;%sqrt(size(copula.fit.lfit.grid,1));
else
knots = copula.fit.knots;%sqrt(size(copula.fit.lfit.grid,1));
end
[gr,GRID_u]= mk_grid(knots,method);
copula.Grid.u=double(GRID_u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cc=reshape(copula.C_grid,knots,knots);
% for i=1:size(v,1)
% F=interp1(copula.Grid.u(:,2),cc(:,i),,'linear','none');
% end
% u=F(v);



S=copula.Grid.u(:,1);
SS=reshape(S,sqrt(size(S,1)),sqrt(size(S,1)));
X=SS(:,1);
S=copula.Grid.u(:,2);
SS=reshape(S,sqrt(size(S,1)),sqrt(size(S,1)));
Y=SS(1,:);

G=copula.C_grid;


%%%% here is the old version

ww=[];
for i=1:size(v,1)

% vv=find(X>=v(i,1)); 
% [~,m1]=min(abs(X(vv)-v(i,1))); 
% g=G(vv(m1),:);

[~,m1]=min(abs(X-v(i,1))); 
g=G(m1,:);

gg=find(g>=v(i,2));
ggg=g(gg);%+1e-10*rand(size(g(gg)));

% [~,m2]=min(abs(v(i,2)-g));
[~,m22]=min(abs(v(i,2)-ggg)); %%%%% changed may 2018
m2=gg(m22);

u(i)=Y(m2);


% u(i)=interp1(g,Y,v(i,2),'linear','extrap');
end


return





ww=[];
for i=1:size(v,1)

% vv=find(X>=v(i,1)); 
% [~,m1]=min(abs(X(vv)-v(i,1))); 
% g=G(vv(m1),:);

[~,m1]=min(abs(X-v(i,1))); 
g=G(m1,:);

mag=max(X(:))+1e-10;
mig=min(X(:))-1e-10;
g=(mag-mig)*(g-min(g))./(max(g)-min(g))+mig;

gg=find(g>=v(i,2));
ggg=g(gg);%+1e-10*rand(size(g(gg)));

if numel(gg)~=numel(g) & gg(1)~=1 %%%% august 2018 added
    ggg=g(gg(1)-1:gg(end));%+1e-10*rand(size(g(gg)));
    gg=[gg(1)-1 gg];
end

[~,m22]=min(abs(v(i,2)-ggg)); %%%%% changed may 2018
m2=gg(m22);

u(i)=Y(m2);

[x n]=unique(ggg); %%august 2018
u(i)=interp1(x,Y(gg(n)),v(i,2),'linear');

if u(i)>max(X(:))
    u(i)=max(X(:))-1e-10;
end
if u(i)<min(X(:))
    u(i)=min(X(:));
end


end

return


