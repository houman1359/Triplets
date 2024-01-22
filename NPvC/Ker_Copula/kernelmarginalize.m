function [vine,copula]=kernelmarginalize(dims,vine_s,copul_s)

GO=0;
if dims(1)~=1 & max(dims)~=numel(dims)
    if numel(dims)~=2 
        error
    else
        
        if dims(1)~=1
            error
        end
        copula{1,2}=copul_s{1,dims(2)};
        copula{1,1}=copul_s{1,1};
        copula{2,1}=copul_s{dims(2),1};
        
        vine.type='c-vine';
       
        vine.condition=vine_s.condition;

        vine.theta=cell(2,2);
        for i=1:numel(dims)
            for j=1:numel(dims)
                vine.families{i,j}=vine_s.families{dims(i),dims(j)};
                vine.METH{i,j}=vine_s.METH{dims(i),dims(j)};            
            end
        end
        for i=1:numel(dims)
            vine.theta{i,i}=vine_s.theta{dims(i),dims(i)};
            vine.METH{i,i}=vine_s.METH{dims(i),dims(i)};
            vine.margins{i}=vine_s.margins{dims(i)};
            vine.range(i,:)=vine_s.range(dims(i),:);
        end
        
        GO=1;
        
    end
end

% 2D marginalization over dims
if GO==0
    clear vine
    
    vine.type='c-vine';
    vine.condition=vine_s.condition;
    
    vine.theta=cell(2,2);
    for i=1:numel(dims)
        for j=1:numel(dims)
            vine.families{i,j}=vine_s.families{dims(i),dims(j)};
            vine.METH{i,j}=vine_s.METH{dims(i),dims(j)};            
        end
    end
    for i=1:numel(dims)
        vine.theta{i,i}=vine_s.theta{dims(i),dims(i)};
        vine.METH{i,i}=vine_s.METH{dims(i),dims(i)};
        vine.margins{i}=vine_s.margins{dims(i)};
        vine.range(i,:)=vine_s.range(dims(i),:);
    end
    
    
    % if dims(1)==1
    % copula{1,2}=copul_s{1,dims(2)};
    % copula{1,1}=copul_s{dims(1),1};
    % copula{2,1}=copul_s{dims(2),1};
    % elseif dims(2)==1
    % copula{1,2}=copul_s{1,dims(1)};
    % copula{1,1}=copul_s{dims(2),1};
    % copula{2,1}=copul_s{dims(1),1};
    % else
    % error('kernel marginalize error: one of the dimensions should be 1')
    % end
    
    D=max(dims);
    
    for i=1:numel(dims)
        copula{i,1}=copul_s{dims(i),1};
    end
    
    for i=1:D
        for j=2:D- (i-1)
            copula{i,j}=copul_s{i,j};
            copula{i,j}.fit.knots=copul_s{i,j}.fit.knots;
        end
    end
end
