
% #' Custom function for local likelihood fitting

function lfit=loclik_fit(bw,data,Grid)


[Ker_grid_Naive]=DenseNaive(bw,(data.X)',(Grid.X)');
[Ker_grid]=Kern_LL(Ker_grid_Naive,bw);

% lfit_grid_mesh=GRID.X;
lfit.grid=Grid.X;
lfit.Kergrid=Ker_grid;
lfit.KergridNaive=Ker_grid_Naive;








