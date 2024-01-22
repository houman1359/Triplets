function Ker_grid=NPC_DenseNaiveMatwin(B,data,grid)
%This function is related to equation (26) of the paper. It represents an
% % intermediate passage to compute the local likelihood function (LLF).
% 
% %%% INPUT:
% % - B    --> It represents the bandwidth of the copula
% % - data --> It represents the data on which the LLF is being computed
% % - grid --> It represents the point of the grid on which the LLF is being
% %            computed, i.e. the (p,q) space.
%
% %%% OUTPUT:
% % - Ker_grid ---> It is a matrix (5 x length of the data) which row
% %                 correspons respectively to f_naive, f1, f2, f3 and f4.
%
g_n = size(grid,2); %n of columns grid
d_n = size(data,2); %n of columns data
Ker_grid = zeros(5, g_n);

if size(B,1) == 1 || size(B,2) == 1
    b1 = B(1);
    b2 = B(2);
else
    for i = 1:d_n
        b1 = B(1,i);
        b2 = B(2,i);
    end
end
for k = 1:g_n
    c1 = grid(1,k) - data(1,:);
    c2 = grid(2,k) - data(2,:);
    a = exp(-c1.^2 / (2*b1^2)) .* exp(-c2.^2 / (2*b2^2)) / (2*pi*b1*b2*d_n);
    
    Ker_grid(1,k) = sum(a);
    Ker_grid(2,k) = sum(a.*c1);
    Ker_grid(3,k) = sum(a.*c2);
    Ker_grid(4,k) = sum(a.*c1.^2);
    Ker_grid(5,k) = sum(a.*c2.^2);
end
end