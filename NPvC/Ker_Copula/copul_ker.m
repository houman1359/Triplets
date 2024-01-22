

function [c1 c2 c3 c4 co] = copul_ker(family,u,theta,dim,bw,METH)

% Argument checks
if nargin < 3 & ~strcmp(family,'kercop')
    error('copulaccdf: Usage c = copulaccdf(family,u,theta,dim)');
end
if ~ischar(family)
    error('copulaccdf: Argument "family" must be a string');
end
if ~ismatrix(u)
    error('copulaccdf: Argument "u" must be a matrix');
end
if size(u,2) ~= 2
    error('copulaccdf: Second dimension of u must be 2');
end
if ~ismatrix(theta)
    error('copulaccdf: Argument "theta" must be a matrix');
end
if nargin < 4
    dim = 2;
elseif ~isscalar(dim) || (dim ~= 1 && dim ~= 2)
    error('copulaccdf: dim must be 1 or 2');
end

family = lower(family);
if dim == 1
    if strcmp(family,'claytonrot090')
        family = 'claytonrot270';
    elseif strcmp(family,'claytonrot270')
        family = 'claytonrot090';
    end
%     c = copulaccdf(family,[u(:,1) u(:,2)],[theta(:,1) theta(:,2)],1,bw);
%     return;
end

% u(u < 0) = 0;
% u(u > 1) = 1;

switch family
    case 'ind'
        c = u(:,1);
    case 'gaussian'
        if ~isscalar(theta) || theta < -1 || theta > 1
            error('copulaccdf: For Gaussian, theta must be in [-1, 1]');
        end
        x = norminv(u);
        c = normcdf((x(:,1) - theta .* x(:,2)) ./ sqrt(1-theta.^2));
    case 'student'
        if numel(theta) ~= 2
            error('copulaccdf: For Student, theta must be a vector with two elements');
        end
        if theta(1) < -1 || theta(1) > 1
            error('copulaccdf: For Student, theta(1) must be in [-1, 1]');
        end
        if theta(2) <= 0
            error('copulaccdf: For Student, theta(2) must be positive');
        end
        x = tinv(u,theta(2));
        c = tcdf(sqrt((theta(2)+1) ./ (theta(2)+x(:,2).^2)) .* (x(:,1) - theta(1) .* x(:,2)) ./ (sqrt(1-theta(1).^2)),theta(2)+1);
    case 'clayton'
        if ~isscalar(theta) || theta < 0
            error('copulaccdf: For Clayton, theta must be in [0, inf)');
        end
        if theta == 0
            c = u(:, 1);
        else
            c = max(u(:,2).^(-1-theta) .* (u(:,1).^(-theta) + u(:,2).^(-theta) - 1).^(-1-1./theta),0);
        end
    case 'claytonrot090'
        if ~isscalar(theta) || theta < 0
            error('copulaccdf: For ClaytonRot090, theta must be in [0, inf)');
        end
        if theta == 0
            c = u(:,1);
        else
            c = max((1-u(:,2)).^(-1-theta) .* (u(:,1).^(-theta) + (1-u(:,2)).^(-theta) - 1).^(-1-1./theta),0);
        end
    case 'claytonrot180'
        if ~isscalar(theta) || theta < 0
            error('copulaccdf: For ClaytonRot180, theta must be in [0, inf)');
        end
        if theta == 0
            c = u(:,1);
        else
            c = 1 - max((1-u(:,2)).^(-1-theta) .* ((1-u(:,1)).^(-theta) + (1-u(:,2)).^(-theta) - 1).^(-1-1./theta),0);
        end
    case 'claytonrot270'
        if ~isscalar(theta) || theta < 0
            error('copulaccdf: For ClaytonRot270, theta must be in [0, inf)');
        end
        if theta == 0
            c = u(:,1);
        else
            c = 1 - max(u(:,2).^(-1-theta) .* ((1-u(:,1)).^(-theta) + u(:,2).^(-theta) - 1).^(-1-1./theta),0);
        end
    case 'kercop'
            [c1 c2 c3 c4 co] = kernel_all(u,theta,dim,bw,METH);
    otherwise
        error(['copulaccdf: Unknown family "' family '"']);
end

end
