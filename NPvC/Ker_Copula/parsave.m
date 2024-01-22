



function parsave(varargin)

savefile = varargin{1}; % first input argument

% for i = 2:nargin
%     savevar.(inputname(i)) = varargin{i}; % other input arguments
% end

% save(savefile,'-struct','savevar')


for i = 2:nargin
savevar = varargin{i}; % other input arguments
save(savefile,'savevar')
end

