P=[];
for i = 1:19
    P=[P;readmatrix(sprintf('Pl_%d.txt', i))];    % Extract chunk
end

[Rb,Modules,Costs,Outer,Inner,Att,ID_P,ppg] = Modularization_Facade(P,'Regularize',0,'Plot_1',1);

% Generating Random Location for Source
Rb.k = 3; Rb.m = 5; % Choice of # of sources and # of Robots
Rb.X_k=[5,-3;37,-6;62,-4]; % Choice of source locations

% Delivery and Assembly Variables 
ms=1; Ms=5; % standardization of weights
Rb.tau_go = pdist2(Rb.X_k,Rb.X_n,"squaredeuclidean"); ng=min(min(Rb.tau_go)); Ng=max(max(Rb.tau_go));
Rb.tau_go = (Rb.tau_go-ng)/(Ng-ng)*(Ms-ms)+ms;
Rb.tau_back = pdist2(Rb.X_n,Rb.X_k); nb=min(min(Rb.tau_back)); Nb=max(max(Rb.tau_back));
Rb.tau_back = (Rb.tau_back-nb)/(Nb-nb)*(sqrt(Ms)-sqrt(ms))+sqrt(ms);
ns=min(Rb.s); Ns=max(Rb.s);
Rb.s = (Rb.s-nb)/(Ns-ns)*(Ms*1.5-ms)+ms;
Ss=min(Rb.tau_go'+Rb.s+Rb.tau_back,[],2);
Rb.Q=full(double(adjacency(graph(Rb.E(:,1),Rb.E(:,2)))));

%% Hard Connectivity
% M-ATSP no Connectivity and no Constraints
data.E=Rb.Q; data.Tgo=Rb.tau_go; data.Tback=Rb.tau_back; data.S=Rb.s; 
data.destConn='SCF';
data.n=size(Rb.s,1); data.k=Rb.k; data.m=Rb.m;
opts.verbose   = 1;
opts.useGurobi = true; % set true if you have Gurobi installed

resN = solve_matsp_minimal(data, opts);
%% Connectivity Constraint with Source Change within Sub-routes
% Minimizing Maximum Cost
resM = solve_minmax_connected_sources(Rb.tau_go, Rb.tau_back, Rb.s, Rb.Q, Rb.m,'TimeLimit',200);

% Minimizing Sum of Absolute Deviation
resMB = solve_balanced_connected_sources(Rb.tau_go, Rb.tau_back, Rb.s, Rb.Q, Rb.m,'TimeLimit',200);

%% Connectivity Constraint with Single Source within Sub-routes

% Minimizing Maximum Cost
resMS = solve_minmax_connected_fixedsources(Rb.tau_go, Rb.tau_back, Rb.s, Rb.Q, Rb.m,'TimeLimit',200);

% Minimizing Sum of Absolute Deviation
resMSB = solve_balanced_connected_fixedsources(Rb.tau_go, Rb.tau_back, Rb.s, Rb.Q, Rb.m,'TimeLimit',200);

%% Visualization
optsVisual = struct('SXY', Rb.X_k, 'DXY', Rb.X_n, ...
    'ShowNodeIDs', false, ...
    'ShowStartEndText', true, ...
    'ArrowHeadAt', 'middle', ...       % put the arrowhead in the middle
    'ArrowHeadLenFrac', 0.01, ...      % same head size for all arrows
    'ArrowHeadWidFrac', 0.01, ...
    'ArrowWidth', 1.5, ...
    'FigTitle', 'Minimizing Makespan – Flexible Source – No Connectivity', ...
    'MarkerSizeD', 50, 'Colors', rand(Rb.m,3),'Layout', 'auto');
vizN = parse_solution_viz_new(resN, Rb.tau_go, Rb.tau_back, Rb.s, Rb.m, Rb.E, optsVisual); optsVisual.FigTitle='Minimizing Makespan – Flexible Source';
vizM = parse_solution_viz_new(resM, Rb.tau_go, Rb.tau_back, Rb.s, Rb.m, Rb.E, optsVisual); optsVisual.FigTitle='Minimizing Sum of Absolute Deviation – Flexible Source';
vizMB = parse_solution_viz_new(resMB, Rb.tau_go, Rb.tau_back, Rb.s, Rb.m, Rb.E, optsVisual); optsVisual.FigTitle='Minimizing Makespan – Fixed Source';
vizMS = parse_solution_viz_new(resMS, Rb.tau_go, Rb.tau_back, Rb.s, Rb.m, Rb.E, optsVisual); optsVisual.FigTitle='Minimizing Sum of Absolute Deviation – Fixed Source';
vizMSB = parse_solution_viz_new(resMSB, Rb.tau_go, Rb.tau_back, Rb.s, Rb.m, Rb.E, optsVisual);

%% Result Table
Rb.C=Rb.tau_go'+Rb.s+Rb.tau_back;
Rt=[resN.routeCost;resM.routeCost;resMB.routeCost;resMS.routeCost;resMSB.routeCost];
T=ResTable(Rt,Rb)