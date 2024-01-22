function u=kerncopccdf(copula,v,method)

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

u=[];
for i=1:size(v,1)

[~,m1]=min(abs(X-v(i,1))); 
[~,m2]=min(abs(Y-v(i,2))); 
g=G(m1,m2);


u(i)=g;


end


