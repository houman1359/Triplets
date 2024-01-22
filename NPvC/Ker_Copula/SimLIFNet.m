function [spk NetParams V] = SimLIFNet(W,varargin)
% Easily numerically integrate a network of LIF neurons
%
% [spk NetParams V] = SimLIFNet(W);
% [spk NetParams V] = SimLIFNet(W,'parameter',value,...);
%
% Easily simulate a network of leaky integrate and fire (LIF) neurons.
%
% Required input: W
%   W is a weight matrix specifying the magnitude of the effect neuron i
% has on neuron j, where i are rows and j are columns of W, which must
% therefore be a square matrix. Weights can be positive (excitatory),
% negative (inhibitory), or zero (no connection). Essentially, each row
% represents a neuron, and the column values for that row represent the
% connection strength the row neuron has on the column neuron.
%
% Alternate syntax:
% [...] = SimLIFNet( NetParams );
%
% Where NetParams is the structure output by a previous call to SimLIFNet
% (or manually constructed by the user). This allows a simulation to be
% rerun easily, or to compactly store simulation parameters for later use.
%
%
% === Optional Parameter Inputs ===
% simTime : The amount of simulated time over which to integrate the model,
%   in units of membrane time constants. {100}
% tstep : Simulated time step-size. The algorithm for numerical integration
%   is the fixed step size Euler method. {0.01}
% initialConditions : An Nx1 vector specifying the initial membrane voltage
%   of each neuron, where N=length(W) and where the nth element of the
%   array corrisponds to the value of the nth cell.
% refractoryTime : An Nx1 vector specifying the refractory period of each
%   neuron, the amount of time the membrane potential is fixed at rest
%   following a spike, where N=length(W). {0}
% offsetCurrents : An Nx1 vector specifying a constant bias current to be
%   injected continuously into each cell. {0}
% synapticDensity : An NxN matrix specifying the density of the synaptic
%   connection from the ith neuron to the jth neuron, where i is the row of
%   the input, j is the column of the input, and N=length(W). {4}
%   Alternitively, an Nx1 array is an acceptable input, where it is assumed
%   that the synaptic connection from the ith neuron to all other neurons
%   is the ith element of the input. This input is the 'a' parameter for
%   the Isynapse equation below. The higher the 'a' parameter the faster
%   the full charge from the pre-synaptic neuron is delivered to the post-
%   synaptic neuron. The total delivered charge is then scaled by the
%   corrisponding parameter in the W matrix.
% forcingFunctions : An Xx2 cell array with the first column being an array
%   of function handles specifying an arbitrary time-varying forcing
%   function of the neuron where the second column contains a scalar
%   indicating the neuron to which the forcing is applied. X must be
%   less than or equal to length(W), and indicies in the second column
%   cannot be repeated. Neurons without a specified forcing function have
%   zero forcing. The function must take the current time as the input
%   and return a scalar output.
% noiseAmplitude : An Nx1 array specifying the amplitude of a zero-mean
%   gaussian noise source of current added to the corrisponding neuron,
%   where N=length(W). If not specified, noise amplitude is set to zero.
% displayProgress : Display a progress bar of the simulation. [0 {1}]
% plotResults : Display a plot of simulation results. [0 {1}] The display
%   includes a plot of all LIF neuron membrane voltages (if all 3 output
%   variables are requested), a raster indicating the spike times of each
%   neuron, and two graphical representations of the neural network. The
%   first is a diagram of connected nodes with each circle representing a
%   neuron and each line representing the connection between neurons
%   (colored by connection strength). Black nodes indicate no external
%   forcing and red nodes indicates external focring is applied. The 'x'
%   markes the neuron 1, and neurons continue in acending order counter
%   clockwise. The second representation is a heatmap of the W input, the
%   color scale of which also corrisponds to the connections on the node
%   graph.
% 
%
% === Outputs ===
% spk : A cell array of spike times, each row corrisponding to a neuron and
%   each cell containing the vector of times at which each spike occured.
% NetParams : A structure containing information and parameters about the
%   simulation.
% V : An NxK matrix where each row is the time-series membrane voltage of
%   neuron and each column is a successive time point. If the simulataion
%   is large, there may not be sufficient memory to store all voltages at
%   all times. If this output is omitted, the membrane voltages will not be
%   stored.
%
%
% === Discussion ===
% The equations are presented in dimensionless form (as outlined in Lewis &
% Rinzel 2003, equation 5), where the membrane potential is scaled by the
% ratio of the voltage minus basline to the threshold minus baseline, and
% time is rescaled by the ratio of the membrane time constant to the
% conductance. Therefore, the derivative of the membrane potential is:
%
% dVi/dt = -Vi + I_applied + sum_k[ Wki*Isynapse ]
%
%                     2   -a*(t - Kj) 
% Isynapse = sum_j[  a * e          *(t - Kj) ] * Wji
%
%   where Kj is time of the vector of spike times from neuron j, t is the
% current time, a is the synaptic density parameter, and Wji is the
% coupling strength from the jth neuron to the ith. The expression
% describes the dynamics of current input to neuron i based on the
% time-history of spikes from neuron j. I_applied are all external forcing
% currents, baises and noise. All parameters have been nondimensionalized
% and scaled such that the spiking threshold is 1 and resting potential is
% 0.
% 
% Lewis TJ & Rinzel J (2003) Dynamics of spiking neurons connected by both
% inhibitory and electrical coupling. J Comput Neurosci 14(3):283-309.
%
% The method of numerical integration is a fixed step size Euler method.
%
% The code run-time increases with the number of time-steps, the number of
% neurons and, importnatly, the number of non-zero synaptic connections
% between neurons. If simulations are slow, consider setting small weights
% to zero.
%
%
%
% === Examples ===
% 1) Reproduce figure 1 from Lewis & Rinzel 2003:
% Two neurons coupled with inhibition and synaptic strength of -0.2
% >> W = [0 -0.2; -0.2 0];     
% Asynchronous behaviour at low driving currents (1.1) and synaptic
% densities of 3, and initial membrane potentials of 0.4 and 0
% >> [spkFull NetParams V] = SimLIFNet(W,'simTime',35,'tstep',1e-2,... 
%           'offsetCurrents',1.1*[1; 1],'synapticDensity',[3; 3],...
%           'initialConditions',[0.4; 0]);
% Now observe synchronous behavior with stronger forcing, 1.6
%  >> [spkFull NetParams V] = SimLIFNet(W,'simTime',35,'tstep',1e-2,... 
%           'offsetCurrents',1.6*[1; 1],'synapticDensity',[3; 3],...
%           'initialConditions',[0.4; 0]);
% 
% 2) Investigate random network architectures
% A random sparse network with noise
% >> W = 1*rand(8)-0.5;
% >> W(randperm(numel(W))>round(0.25*numel(W))) = 0;
% >> [spk NetParams V] = SimLIFNet(W,'simTime',100,'tstep',1e-2,...
%           'offsetCurrents',0.8*ones(length(W),1),...
%           'noiseAmplitude',0.2*ones(length(W),1));
% A deterministic medium-sized highly inhibitory network initialized
% randomly
% >> W = log(abs(randn(12)));
% >> [spk NetParams V] = SimLIFNet(W,'simTime',35,'tstep',1e-2,...
%           'offsetCurrents',1.1*ones(length(W),1));
%
% 3) An example using a neuron with a refractory period:
% All neurons will drive neuron 3 at high rate, but neuron 3 will have a
% long refractory period
% >> W = [0 0 0.5; 0 0 0.5; 0 0 0.2];
% >> [spk NetParams V] = SimLIFNet(W,'simTime',50,'tstep',1e-2,...
%           'offsetCurrents',[1.5 1.5 0.5]','refractoryTime',[0 0 5]',...
%           'initialConditions',[0.2 0.4 0]');
%
% 4) Examine effects of synaptic density:
% Connect a driver neuron (1) to others via varied synaptic densities
% >> W = [0 ones(1,3)*0.4; zeros(3,4)];
% >> [spk NetParams V] = SimLIFNet(W,'simTime',35,'tstep',1e-2,...
%           'offsetCurrents',[1.6 0.6 0.6 0.6]',...
%           'synapticDensity',[0 1 4 7; zeros(3,4)]);
%
% 5) Apply forcing functions:
% Neuron 1 excites 3, neuron 2 inhibits 3.
% >> W = [0 0 0.5; 0 0 -0.5; 0 0 0];
% Neuron 1 has constant forcing starting at t=10, neuron 2 has sinusoidal
% forcing throughout
% >> Ffcns = {@(t) 1.5*heaviside(t-10), 1; @sin, 2};
% >> [spk NetParams V] = SimLIFNet(W,'simTime',35,'tstep',1e-2,...
%           'offsetCurrents',[0.8 0.8 0.8]','forcingFunctions',Ffcns);
%
% 6) Forcing with noise:
% Apply a constant current to a single neuron from t=10 to 50 in the
% presence of noise (heaviside is slow for larger networks)
% >> [spk NetParams V] = SimLIFNet(0,'forcingFunctions', ... 
%           {@(t) 1.6*(heaviside(t-10)-heaviside(t-50)), 1}, ... 
%           'noiseAmplitude',0.3);
%
% 7) Building logical operators with LIF networks (exploiting adjustable
% synaptic densities):
% Build a network with modular logical elements. The last (10th) neuron
% will serve as the output. Neurons 1 and 2 are inputs to an OR operator;
% neurons 3 and 4 are inputs to an AND operator; neuron 6 is the input to a
% NOT opperator while input 5 is an interneuron providing chronic
% stimulation to the output neuron, neurons 8 and 9 are inputs to an XOR
% operator with slow synaptic connections while neuron 7 is a fast
% inhibitory interneuron.
% >> W = zeros(10);
% >> W([99 98 97 96 95 94 93 92 91 69 68])=[3 3 -15 -1.5 2 1 1 2 2 1 1];
% >> A = 4*ones(10); A([69 68 99 98 97])=[8 8 1 1 8];
% >> F = {@(u) 1.4*(ge(u,5)-ge(u,15)+ge(u,35)-ge(u,45)),1; ...
%         @(u) 1.4*(ge(u,20)-ge(u,30)+ge(u,35)-ge(u,45)),2; ...
%         @(u) 1.4*(ge(u,55)-ge(u,65)+ge(u,70)-ge(u,80)),3; ...
%         @(u) 1.4*(ge(u,70)-ge(u,80)+ge(u,85)-ge(u,95)),4; ...
%         @(u) 1.4*(ge(u,100)-ge(u,140)),5; ...
%         @(u) 1.4*(ge(u,120)-ge(u,130)),6; ...
%         @(u) 1.4*(ge(u,160)-ge(u,170)+ge(u,190)-ge(u,210)),8; ...
%         @(u) 1.4*(ge(u,175)-ge(u,185)+ge(u,190)-ge(u,210)),9};
% >>  [spk NetParams V] = SimLIFNet(W,'simTime',220,'tstep',1e-2, ... 
%       'forcingFunctions',F,'refractoryTime',ones(length(W),1)/2, ...
%       'synapticDensity',A);
% >>  set(NetParams.handles.sub1,'ylim',[-0.4 1.9])
%
% 8) Errors due to time-step sizes using multi-part forcing
% Build a network with complex forcing
% >> W = [0 -0.2; -0.1 0];
% >> F = { @(u) (u>3 & u<5)*0.8 + (u>=5 & u<10.5)*sin(2*pi/2*u) + ...
%               (u>=10.5)*exp(-(u-10.5)/3),1};
% Run the simulation using a range of time steps.
% >> dt = [logspace(-4,-1,4) 0.5];
% >> V = cell(size(dt));
% >> pr = [1 zeros(1,length(dt)-1)];
% >> for k=1:length(dt)
%       [spk NetParams V{k}] = SimLIFNet(W,'simTime',25,'tstep',dt(k), ... 
%           'forcingFunctions',F,'offsetCurrents',[1; 1.2],'plotResults',pr(k));
%    end
% Plot the difference between the the smallest time step simulation and
% each of the other larger time step simulations (the error) for the second
% neuron in the network.
% >> figure; hold on
% >> cmp = hsv(length(dt));
% >> for k=1:length(dt)
%       plot( linspace(0,25,length(V{k})), V{1}(2,1:dt(k)/dt(1):end)-V{k}(2,:),'color',cmp(k,:),'linewidth',2 )
%    end
% >> ylabel('Difference from Smallest Timestep'); xlabel('time')
% >> legend(num2str(dt'),'location','SW')      
%
%
% %%% ZC Danziger March 2015 %%%
%

% Bugs:
% - full V vs temporary V slight results discrepency (numerical error?)


%% Initialize

if isstruct(W)
    % This code only for special case when user calls the function with a
    % structure of parameters
    
    paramList = fieldnames(W);                      % list of parameters set by user
    valueList = struct2cell(W);                     % list of corresponding field values
    % make adjustments to parameters to conform to input syntax
    Wloc = strcmp(paramList,'W');
    Wnew = valueList{Wloc};                         % find W, assign it to variable
    % remove unusable inputs
    badLocs = cellfun(@(u) any(strcmp(u,{'W','elapsedTime','simDate','handles'})),paramList);
    paramList(badLocs)=[];                             % remove W from other inputs
    valueList(badLocs)=[];
    % locate forcing functions for syntax updates
    FFloc = strcmp(paramList,'forcingFunctions');   
    if any(FFloc)
        valueList{FFloc} = [valueList{FFloc} num2cell((1:length(valueList{FFloc}))')];
        % remove null functions
        removeNullsIX = cellfun(@(x) isequal(func2str(x),'@(u)0'),valueList{FFloc}(:,1));
        if all(removeNullsIX)
            paramList(FFloc) = [];
            valueList(FFloc) = [];
        else
            valueList{FFloc}(removeNullsIX,:) = [];
        end
    end
    % combine lists for input
    combinedList = [paramList, valueList]';
    
    % re-call the function with the appropriate argument syntax
    if nargout == 3
        [spk, NetParams, V] = SimLIFNet(Wnew,combinedList{:});
    else
        [spk, NetParams] = SimLIFNet(Wnew,combinedList{:});
    end
    return;
    
else
    
    % this is a normal call, initialize variables
    e1 = clock;
    N = length(W);
end



%% Parse inputs
P = inputParser;
P.addRequired('W', @(u) isnumeric(u) && isreal(u) && diff(size(W))==0);
P.addParamValue('simTime', 100, @(u)isnumeric(u) && numel(u)==1);
P.addParamValue('tstep', 1e-2, @(u)isnumeric(u) && numel(u)==1);
P.addParamValue('initialConditions', zeros(N,1), @(u)isnumeric(u) && size(u,1)==N);
P.addParamValue('refractoryTime', zeros(N,1), @(u)isnumeric(u) && size(u,1)==N);
P.addParamValue('offsetCurrents', zeros(N,1), @(u)isnumeric(u) && size(u,1)==N);
P.addParamValue('displayProgress', true, @(u)u==1 || u==0);
P.addParamValue('plotResults', true, @(u)u==1 || u==0);
P.addParamValue('forcingFunctions', {@(u) 0, 1}, ...
    @(u)iscell(u) && size(u,2)==2 && ...
    all(cellfun(@(v) strcmp(class(v),'function_handle'),u(:,1))) && ...
    max([u{:,2}])<=N && ~(length(unique([u{:,2}]))<length([u{:,2}])) );
P.addParamValue('synapticDensity', ones(N,1)*4, @(u) isnumeric(u) && ( all(size(u)==[N 1]) || all(size(u)==[N N]) ));
P.addParamValue('noiseAmplitude', zeros(N,1), @(u) isnumeric(u) && all(size(u)==[N 1]));
% do not allow changes to these parameters to avoid violation of the
% nondimensionalization constructions
% P.addParamValue('thresholds', ones(N,1), @(u)isnumeric(u) && size(u,1)==N);   
% P.addParamValue('restVoltage', zeros(N,1), @(u)isnumeric(u) && size(u,1)==N);

P.parse(W,varargin{:});
NetParams = P.Results;



%% Set up simulation
% initialize state variables
ix = repmat(1:N,[N 1]);             % network index (row:from col:to)
spk = num2cell(-1e10*ones(1,N));    % cell array of spikes

% network parameters
vt = 1;
vr = 0;
tr = max([NetParams.refractoryTime repmat(NetParams.tstep,[N 1])],[],2);
Ia = NetParams.offsetCurrents;
V0 = NetParams.initialConditions;
nAmp = NetParams.noiseAmplitude;

% build up forcing functions
forceFcns = cell(N,1);
forceFcns([NetParams.forcingFunctions{:,2}]) = NetParams.forcingFunctions(:,1);
forceFcns(setxor([NetParams.forcingFunctions{:,2}],1:N)) = {@(u) 0};


% build up synaptic density array
if size(NetParams.synapticDensity,2)==1
    % each cell has constant 'a' when delivering pulses to other neurons
    a  = NetParams.synapticDensity;
    a = a(ix);
else
    % all 'a' between each cell are specified and possibly different
    a  = NetParams.synapticDensity;
end

% build up sparse version synaptic connection weights and densities to
% avoid calculating synaptic transmission of zero-weighted connections
[Wrow Wcol] = find(W);
Wsparse = [Wrow Wcol];      % a from->to array of only non-zero W elements
Aix = find(W);              % as with Wsparse but with linear indexing
if isempty(Wsparse)         % special case when there are no connections
    Wsparse = [1 1];
    Aix = 1;
end
% synaptic connection function
G = @(u,v,tG) sum(v.^2.*(tG-u).*exp(-v.*(tG-u)));


    
% simulation parameters
t0 = 0;                     % start time
tf = NetParams.simTime;     % stop time
dt = NetParams.tstep;       % time incrememnt
if nargout==3
    % case where user requests all membrane data
    try
        t = t0:dt:tf;           % vector of simulation times
        K = length(t);          % number of total iterations in simulation
        V = nan(N,K);           % membrane voltages
        V(:,1) = V0;            % apply initial conditions
    catch err
        disp(err.message)
        error('Unable to Store Full State Variables: please re-try without 3rd output argument')
    end
else
    t = [nan t0];           % t(1) is a placeholder only, to make indexing convenient
    K = tf/dt+1;
    V = [nan(N,1) V0];
end



%% Run simulation
vis = NetParams.displayProgress;
if vis, hw = waitbar(0,'Integrating Network...'); e2=clock; end
for k=2:K
    % index to data
    if nargout==3
        id=k;                       % case where membrane voltage stored
    else
        id=2;                       % case where only spikes are stored
        t(id) = (k-1)*dt;           % increment time
        V(:,1) = V(:,2);            % overwrite old state with new state
    end
    
    % Synapse contribution to membrane derivative:
    % calculate all non-zero synaptic transmissions
    % -- actually, we just want speed so put this in a for loop --
    spikeEffects = zeros(size(Wsparse,1),1);
    for sl=1:size(Wsparse,1)
        spikeEffects(sl) = G(spk{Wsparse(sl,1)},a(Aix(sl)),t(id)) * W(Aix(sl));
    end
    % sum all transmissions incoming to each neuron
    synapses = accumarray(Wsparse(:,2),spikeEffects,[],@sum,0);
    
    
    % Euler integration
    V(:,id) = V(:,id-1) + dt*( -V(:,id-1) + ...                     % difference of prior membrane voltage from rest
                                Ia + ...                            % bias current
                                dt^(-0.5)*nAmp.*randn(N,1) + ...    % addition of simple white noise (noise scaling b/c we are inside dt)
                                synapses + ...                      % synaptic contribution
                                cellfun(@(u) u(t(id)),forceFcns) ); % forcing functions
    
    % refraction
    isRef = t(id) - cellfun(@max,spk)' < tr;
    if any(isRef)
        V(isRef,id) = vr;
    end
    
    % spiking
    isSpk = V(:,id)>=vt;
    if any(isSpk)
        % append spike times to spk
        spk(isSpk) = cellfun( @(u,v) [u v],spk(isSpk),num2cell(isSpk(isSpk)'*t(id)),'unif',false );
        % reset membrane voltage
        V(isSpk,id) = vr;
    end

    if mod(k,round(K/50))==0 && vis
        e3=clock;
        estTimeRem = min(etime(e3,e2)/(K/50)*(K-k),86400);
        e2=e3;
        waitbar(k/K,hw,...
            ['Integrating. Est. time remaining: ' datestr(estTimeRem/86400,'HH:MM:SS')]);
    end
end
if vis, close(hw); end

% remove place-holder spikes
spk = cellfun(@(u) u(2:end),spk,'unif',false)';


% Extra simulation information
NetParams.forcingFunctions = forceFcns;
NetParams.elapsedTime = etime(clock,e1);
NetParams.simDate = datestr(now);




%% Plotting
if NetParams.plotResults;
    figure('name',mfilename,'color','w')
    
    if nargout==3
        spkHeight = 2;
        Vplot = V;
        for k=1:N
            Vplot( k, round(spk{k}./dt+1) ) = spkHeight;
        end
        NetParams.handles.sub1 = subplot(5,2,[1 2]);
        hmem = plot(t,Vplot);
        hold on
        % plot(zeros(2,N),repmat(vt',[2 1]),'>','linewidth',2)
        line([0 tf],vt*[1 1],'color','k','linestyle','--')
        set(gca,'box','off','ylim',[min(min(Vplot(:)),vr) spkHeight*0.95],...
            'xlim',[0 tf])
        ylabel({'Membrane'; 'Potential'})
        NetParams.handles.membrane = hmem;
    end
    
    NetParams.handles.sub2 = subplot(5,2,3:6); hold on
    hrast = zeros(N,1);
    for ip=1:N
        if ~isempty(spk{ip})
            hrast(ip) = plot(spk{ip},ones(1,length(spk{ip}))*ip,'.k');
        end
    end
    ylim([0.5 N+0.5]); xlim([0 tf])
    ylabel('Raster')
    NetParams.handles.raster = hrast;
    
    % network nodes
    hnet = subplot(5,2,[7 9]); hold on
    th = 2*pi/N:2*pi/N:2*pi;        % theta for node locations
    eps = 0.15;                     % node-size parameter (radii)
    verts = [sin(th); cos(th)];
    for k=1:size(verts,2)
        if strcmp(func2str(NetParams.forcingFunctions{k}),'@(u)0')
            nodeCol = 'k';
        else
            nodeCol = 'r';
        end
        hverts(k) = rectangle('position',[verts(1,k)-eps verts(2,k)-eps 2*eps 2*eps], ...
            'curvature',[1 1],'linewidth',2,'edgecolor',nodeCol,'facecolor','w');
    end
    set(hnet,'xlim',[-1.5 1.5],'ylim',[-1.5 1.5],'visible','off')
    axis equal; 
    % place weighted connections
    colI = 12;   % color gradations
    cmp = [interp1([1 colI],[0.1 0.6 1; 1 1 1],1:colI); interp1([1 colI],[1 1 1; 1 0.7 0.1],2:colI)];
    if any(W(:))
        cmpVec = linspace(-max(abs(W(:))),max(abs(W(:))),2*colI-1);
    else
        cmpVec = [inf(1,colI-1), 0, inf(1,colI-1)];   % no connections case
    end
    hcons = zeros(size(W));
    for r=1:N
    for c=1:N
        [minVal, wcol] = min(abs(W(r,c)-cmpVec));
        if r~=c && wcol~=colI
            shiftVec = [verts(1,ix(r,c))-verts(1,ix(c,r)); ...
                        verts(2,ix(c,r))-verts(2,ix(r,c))];
            shiftVec = eps*shiftVec/norm(shiftVec);
            hcons(r,c) = plot([verts(1,ix(r,c)) verts(1,ix(c,r))] + shiftVec(2),...
                              [verts(2,ix(r,c)) verts(2,ix(c,r))] + shiftVec(1),...
                              'color',cmp(wcol,:),'linewidth',2);
            
        elseif wcol~=colI
            % self connection plotting
            hcons(r,c) = rectangle('position',[verts(1,ix(r,c))*1.2-eps verts(2,ix(r,c))*1.2-eps 2*eps 2*eps], ...
                'curvature',[1 1],'linewidth',2,'edgecolor',cmp(wcol,:));
        end
    end
    end
    uistack(hverts,'top')
    plot(verts(1,1),verts(2,1),'kx','markersize',8,'linewidth',2);
    NetParams.handles.nodes = hverts;
    NetParams.handles.cons = hcons;
    NetParams.handles.sub3 = hnet;
    
    NetParams.handles.sub4 = subplot(5,2,[8 10]);
    imagesc(W)
    colormap(cmp)
    colorbar
    caxis(max(abs(W(:)))*[-1 1]);
    ylabel('Network Connectivity')
    
    
end
