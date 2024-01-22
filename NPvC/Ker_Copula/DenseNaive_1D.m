function [Ker_grid]=DenseNaive_1D(grid,B,data)

%%%% grid is only one point

C1=-data(1,:)+grid(1,1);
C2=-data(2,:)+grid(2,1);

if numel(B)~=2
A=exp((-C1.^2)./(2*B(1,:).^2)).*exp((-C2.^2)./(2*B(2,:).^2))./(2*pi*B(1,:).*B(2,:)*size(data,2));
else
A=exp((-C1.^2)/(2*B(1).^2)).*exp((-C2.^2)/(2*B(2).^2))./(2*pi*B(1).*B(2)*size(data,2));
end

Ker_grid(1)=sum( A );
Ker_grid(2)=sum( A .*C1 );
Ker_grid(3)=sum( A .*C2 );
Ker_grid(4)=sum( A .*(C1.^2) );
Ker_grid(5)=sum( A .*(C2.^2) );


