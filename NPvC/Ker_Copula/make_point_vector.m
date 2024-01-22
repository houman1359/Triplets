function points=make_point_vector(dim,y_vector,vvt)


% points=zeros(numel(y_vector)*(size(vvt,1)),size(vvt,2));
% 
% ne=0;
% for j=1:size(vvt,1)
%     for i=1:numel(y_vector)
%         ne=ne+1;
%         points(ne,dim)=y_vector(i);
%         points(ne,setdiff(1:size(vvt,2),dim))=vvt(j,setdiff(1:size(vvt,2),dim));
%     end
% end


yy=repmat(y_vector',size(vvt,1),1);
vv=repelem(vvt,numel(y_vector),1);

points=zeros(numel(y_vector)*(size(vvt,1)),size(vvt,2));

points(:,dim)=yy;
points(:,setdiff(1:size(vvt,2),dim))=vv(:,setdiff(1:size(vvt,2),dim));