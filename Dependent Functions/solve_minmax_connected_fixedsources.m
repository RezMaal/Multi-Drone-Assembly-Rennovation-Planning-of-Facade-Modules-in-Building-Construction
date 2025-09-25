function [res,sol] = solve_minmax_connected_fixedsources(TGo, TBack, S, Adj, m, varargin)
% Minimize max per-salesman workload with connected districts and explicit sources
% Balanced SCF connectivity with N_r and McCormick linearization of t_{ir}=N_r*g_{ir}
%
% Inputs:
%   TGo   : k x n, cost source s -> city j
%   TBack : n x k or k x n, cost city j -> source s
%   S     : n x 1 (or 1 x n), service times
%   Adj   : n x n adjacency (symmetricized)
%   m     : # salesmen
% Options (Name,Value):
%   'RequireNonEmpty' (true) : force >=1 city per salesman
%   'OutputFlag'      (0)
%   'TimeLimit'       ([] seconds)
%
% Output 'sol' with fields: status, objval, x (n x m), y (n x k), root (n x m),
%   v (n x m x k if small), workloads (1 x m), components{m}, params, model.

% --------- Parse options ---------
p = inputParser;
p.addParameter('RequireNonEmpty', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('OutputFlag', 1, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('TimeLimit', 240, @(x)isempty(x)||(isnumeric(x)&&isscalar(x)));
p.parse(varargin{:});
RequireNonEmpty = p.Results.RequireNonEmpty;
OutputFlag = p.Results.OutputFlag;
TimeLimit  = p.Results.TimeLimit;

% --------- Build costs ---------
[k, n] = size(TGo);
S = S(:)';                                   % 1 x n
if size(TBack,1)==n && size(TBack,2)==k
    TBack_kxn = TBack';
elseif size(TBack,1)==k && size(TBack,2)==n
    TBack_kxn = TBack;
else
    error('TBack must be n x k or k x n to match TGo (k x n).');
end
C = TGo + TBack_kxn + repmat(S, k, 1);       % k x n : trip cost per (s,j)
Tc = sum(min(C, [], 1)) / m;       % global target / valid lower bound for T
cmin = min(C, [], 1);                % 1 x n per-city min cost (used for LB rows)


% --------- Graph arcs ---------
if ~isequal(size(Adj), [n n]), error('Adj must be n x n'); end
Adj = Adj~=0; Adj = Adj | Adj'; Adj(1:n+1:end)=0;
[iu, ju] = find(triu(Adj,1));
A_tails = [iu; ju];                           % directed arcs both ways
A_heads = [ju; iu];
a = numel(A_tails);                           % # directed arcs = 2|E|
Mflow = n - 1;                                % capacity big-M

% --------- Variable layout ---------
nm = n*m; nk = n*k; nmk = n*m*k; ma = m*a;
mk = m*k;   % NEW: salesman-source vars (yRS)

% Order: [x(nm), y(nk), v(nmk), g(nm), f(ma), N(m), t(nm), T(1)]
off.x = 0;
off.y = off.x + nm;
off.v = off.y + nk;
off.g = off.v + nmk;
off.f = off.g + nm;
off.N = off.f + ma;
off.t = off.N + m;
off.T = off.t + nm;
off.yRS = off.T + 1;               % NEW block for yRS (m x k)
numVars = off.yRS + mk;

idxX = @(j,r) off.x + (r-1)*n + j;
idxY = @(j,s) off.y + (s-1)*n + j;
idxV = @(j,r,s) off.v + ((r-1)*n + (j-1))*k + s;
idxG = @(j,r) off.g + (r-1)*n + j;
idxF = @(r,arc) off.f + (r-1)*a + arc;
idxN = @(r)     off.N + r;
idxt = @(j,r)   off.t + (r-1)*n + j;
idxT = @()      off.T + 1;

idxYRS = @(r,s) off.yRS + (s-1)*m + r;   % NEW: salesman-source choice
% --------- Objective & var types ---------
model.modelsense = 'min';
obj = zeros(numVars,1); obj(idxT()) = 1;      % minimize T

vtype = repmat('B', numVars, 1);
% NEW block yRS is binary by default
vtype(off.f+1 : off.f+ma) = 'C';              % flows
vtype(off.N+1 : off.N+m)  = 'C';              % N_r
vtype(off.t+1 : off.t+nm) = 'C';              % t_{ir}
vtype(idxT()) = 'C';                           % T

lb = zeros(numVars,1); ub = inf(numVars,1);
ub(1:off.f) = 1;                               ub(off.yRS+1 : off.yRS+mk) = 1;    % NEW: cap yRS ≤ 1
% binaries ≤1
ub(off.f+1:off.f+ma) = Mflow;                  % f ≤ M
ub(off.N+1:off.N+m)  = n;                      % 0 ≤ N_r ≤ n
% t has no preset UB beyond constraints

% --------- Constraint buffers ---------
I=[]; J=[]; Vals=[]; sense=[]; rhs=[]; row=0;
    function addRow(cols, vals, sg, bb)
        I=[I; repmat(row+1,numel(cols),1)]; %#ok<AGROW>
        J=[J; cols(:)]; Vals=[Vals; vals(:)]; %#ok<AGROW>
        sense=[sense; sg]; rhs=[rhs; bb]; row=row+1; %#ok<AGROW>
    end

% 1) City -> exactly one salesman
for j=1:n
    addRow(arrayfun(@(r) idxX(j,r),1:m), ones(1,m), '=', 1);
end

% 2) City -> exactly one source
for j=1:n
    addRow(arrayfun(@(s) idxY(j,s),1:k), ones(1,k), '=', 1);
end

% 3) v linking: sum_s v(j,r,s) = x(j,r)
for j=1:n, for r=1:m
    cols=[arrayfun(@(s) idxV(j,r,s),1:k), idxX(j,r)];
    vals=[ones(1,k), -1];
    addRow(cols, vals, '=', 0);
end, end
% 3b) NEW: Single-source per salesman
% (i) Exactly one source per salesman: sum_s yRS(r,s) = 1
for r = 1:m
    addRow(arrayfun(@(s) idxYRS(r,s), 1:k), ones(1,k), '=', 1);
end

% (ii) Gate v by chosen source: v(j,r,s) <= yRS(r,s)
for r = 1:m
  for s = 1:k
    for j = 1:n
      addRow([idxV(j,r,s), idxYRS(r,s)], [1, -1], '<', 0);
    end
  end
end


% 4) v linking: sum_r v(j,r,s) = y(j,s)
for j=1:n, for s=1:k
    cols=[arrayfun(@(r) idxV(j,r,s),1:m), idxY(j,s)];
    vals=[ones(1,m), -1];
    addRow(cols, vals, '=', 0);
end, end

% 5) Makespan per salesman: sum_{j,s} C(s,j)*v(j,r,s) <= T
for r=1:m
    cols=[]; vals=[];
    for s=1:k, for j=1:n
        cols(end+1)=idxV(j,r,s); %#ok<AGROW>
        vals(end+1)=C(s,j);      %#ok<AGROW>
    end, end
    cols=[cols, idxT()]; vals=[vals, -1];
    addRow(cols, vals, '<', 0);
end


% 5b) Valid lower bounds to tighten makespan T
% Global:  T >= Tc  (i.e., -T <= -Tc)
addRow(idxT(), -1, '<', -Tc);

% Per-salesman: sum_j cmin_j * x(j,r) <= T
for r = 1:m
    cols = arrayfun(@(j) idxX(j,r), 1:n);
    vals = cmin(:)';
    cols = [cols, idxT()]; vals = [vals, -1];
    addRow(cols, vals, '<', 0);
end

% 6) Roots: one per salesman; root must be assigned
for r=1:m
    addRow(arrayfun(@(j) idxG(j,r),1:n), ones(1,n), '=', 1);
end
for j=1:n, for r=1:m
    addRow([idxG(j,r), idxX(j,r)], [1,-1], '<', 0); % g <= x
end, end

% 7) N_r definition: N_r - sum_j x(j,r) = 0
for r=1:m
    cols=[idxN(r), arrayfun(@(j) idxX(j,r),1:n)];
    vals=[1, -ones(1,n)];
    addRow(cols, vals, '=', 0);
end
% (Optional) Require non-empty: N_r >= 1
if RequireNonEmpty
  for r=1:m
      addRow(idxN(r), 1, '>', 1);
  end
end

% 8) McCormick linearization of t_{ir} = N_r * g_{ir}, with 0<=N_r<=n
for r=1:m, for j=1:n
    % t_ir <= N_r
    addRow([idxt(j,r), idxN(r)], [1, -1], '<', 0);
    % t_ir <= n * g_ir
    addRow([idxt(j,r), idxG(j,r)], [1, -n], '<', 0);
    % t_ir >= N_r - n*(1 - g_ir)  ->  -t_ir + N_r - n + n*g_ir <= 0
    addRow([idxt(j,r), idxN(r), idxG(j,r)], [-1, 1, n], '<', n);
    % t_ir >= 0 is already ensured by lb=0
end, end

% 9) Balanced SCF flow conservation: in - out = x - t
inArcs  = cell(n,1); outArcs = cell(n,1);
for arc=1:a
    outArcs{A_tails(arc)}(end+1)=arc; %#ok<AGROW>
    inArcs{A_heads(arc)}(end+1)=arc;  %#ok<AGROW>
end
for r=1:m, for i=1:n
    cols=[]; vals=[];
    if ~isempty(inArcs{i})
        fa=arrayfun(@(arc) idxF(r,arc), inArcs{i});
        cols=[cols, fa]; vals=[vals,  ones(1,numel(fa))];
    end
    if ~isempty(outArcs{i})
        fa=arrayfun(@(arc) idxF(r,arc), outArcs{i});
        cols=[cols, fa]; vals=[vals, -ones(1,numel(fa))];
    end
    cols=[cols, idxX(i,r), idxt(i,r)];
    vals=[vals, -1, 1];                 % in - out - x + t = 0
    addRow(cols, vals, '=', 0);
end, end

% 10) Flow capacity: f(r,i->j) <= M*x(i,r) and <= M*x(j,r)
for r=1:m 
    for arc=1:a
    tail=A_tails(arc); head=A_heads(arc);
    addRow([idxF(r,arc), idxX(tail,r)], [1, -Mflow], '<', 0);
    addRow([idxF(r,arc), idxX(head,r)], [1, -Mflow], '<', 0);
    end
end

% --------- Finalize & solve ---------
model.A     = sparse(I,J,Vals,row,numVars);
model.rhs   = double(rhs);
model.sense = char(sense(:));
model.obj   = obj;
model.vtype = vtype;
model.lb    = lb;
model.ub    = ub;

params.OutputFlag = OutputFlag;
if ~isempty(TimeLimit), params.TimeLimit = TimeLimit; end
    params.MIPFocus = 1;
    params.Heuristics = 0.6;
    params.RINS = 20;
    params.NoRelHeurTime = min(100,params.TimeLimit/2);
    params.PumpPasses = 10; params.ZeroObjNodes = 1000;
    params.OutFlag=1; 
    params.Cuts=3; params.ConcurrentMIP = 4;
result = gurobi(model, params);

% --------- Extract solution ---------
sol = struct('status', result.status, 'params', params, 'model', model);
if isfield(result,'x')
    x = result.x;
    X = zeros(n,m); Y = zeros(n,k); G = zeros(n,m);
    for r=1:m, X(:,r) = x(off.x + (r-1)*n + (1:n)); end
    for s=1:k, Y(:,s) = x(off.y + (s-1)*n + (1:n)); end
    for r=1:m, G(:,r) = x(off.g + (r-1)*n + (1:n)); end
    sol.x = round(X); sol.y = round(Y); sol.root = round(G);
    sol.yRS = reshape(round(x(off.yRS + (1:mk))), m, k);  % NEW: salesman-source choices
sol.objval = result.objval;

    % Workloads
    workloads = zeros(1,m);
    useV = (nmk <= 1e6);
    if useV
        V = zeros(n,m,k);
        for r=1:m
          base_r = off.v + (r-1)*n*k;
          for j=1:n
            base_j = base_r + (j-1)*k;
            V(j,r,:) = x(base_j + (1:k));
          end
        end
        sol.v = round(V);
        for r=1:m
            Vr = squeeze(V(:,r,:));              % n x k
            workloads(r) = sum(sum((C').*Vr));
        end
    else
        [~, srcIdx] = max(sol.y,[],2);           % n x 1
        for r=1:m
            cities = find(sol.x(:,r)>0.5);
            for jj = cities'
                workloads(r) = workloads(r) + C(srcIdx(jj), jj);
            end
        end
    end
    sol.workloads = workloads;

    % Components (reporting)
    Gadj = graph(Adj);
    sol.components = cell(1,m);
    for r=1:m
        nodes = find(sol.x(:,r)>0.5);
        if isempty(nodes)
            sol.components{r} = {};
        else
            H = subgraph(Gadj, nodes);
            bins = conncomp(H);
            comps = arrayfun(@(b) nodes(bins==b), 1:max(bins), 'uni', 0);
            sol.components{r} = comps;
        end
    end
    res = format_routes(sol, TGo, TBack, S, Adj);
end
end

function res = format_routes(sol, TGo, TBack, S, Adj)
% FORMAT_ROUTES  Translate solver output 'sol' into route-level results
% with a source "timeline": startSource, bridgeSrc (length q-1),
% endSource, so total sources = q+1 where q = #cities in ordersD.
%
% res fields:
%   res.routeCost{r}   : 1 x q  cost of each sub-route (source->city->source)
%   res.startSource(r) : scalar source ID before the 1st city
%   res.bridgeSrc{r}   : 1 x (q-1) source IDs before trips 2..q (timeline)
%   res.endSource(r)   : scalar source ID after the last city
%   res.ordersD{r}     : 1 x q  city IDs in BFS order
%
% Note: In the current model each city uses a single source for both legs,
% so the last bridge source equals endSource (duplicate by design).

[k, n] = size(TGo);
S = S(:)';

% Normalize TBack to k x n
if size(TBack,1)==n && size(TBack,2)==k
    TBack_kxn = TBack';
elseif size(TBack,1)==k && size(TBack,2)==n
    TBack_kxn = TBack;
else
    error('TBack must be n x k or k x n to match TGo (k x n).');
end

if ~isfield(sol,'x') || ~isfield(sol,'y')
    error('sol must contain fields x (n x m) and y (n x k).');
end
[nx, m] = size(sol.x);
[ny, ky] = size(sol.y);
if nx~=n || ny~=n || ky~=k
    error('Dimension mismatch: x is %dx%d, y is %dx%d, but TGo is %dx%d.', nx, m, ny, ky, k, n);
end

Adj = Adj~=0; Adj = Adj | Adj'; Adj(1:n+1:end)=0;

% Pick the source per city from y
[~, srcOfCity] = max(sol.y, [], 2);   % n x 1
srcOfCity = srcOfCity(:)';

% Graph for BFS ordering
G = graph(Adj);

% Outputs
res.routeCost   = nan(1,m);
res.startSource = nan(1,m);
res.bridgeSrc   = cell(1,m);
res.endSource   = nan(1,m);
res.ordersD     = cell(1,m);

for r = 1:m
    cities = find(sol.x(:,r) > 0.5)';   % cities for salesman r
    if isempty(cities)
        res.routeCost(r) = NaN;
        res.bridgeSrc{r} = [];
        res.ordersD{r}   = [];
        res.startSource(r) = NaN;
        res.endSource(r)   = NaN;
        continue;
    end

    % Start city: use given root if present; otherwise first city
    if isfield(sol,'root') && any(sol.root(:,r) > 0.5)
        startCity = find(sol.root(:,r) > 0.5, 1, 'first');
        if ~ismember(startCity, cities), startCity = cities(1); end
    else
        startCity = cities(1);
    end

    % BFS order on the induced subgraph
    H = subgraph(G, cities);
    startLocal = find(cities == startCity, 1, 'first');
    if isempty(startLocal), startLocal = 1; end
    bfsLocal = bfsearch(H, startLocal);
    order = cities(bfsLocal);      % city IDs in timeline order
    q = numel(order);

    % Source used for each city's trip (start==end for that city)
    srcSeq = arrayfun(@(j) srcOfCity(j), order);

    % Costs per sub-route (source->city->source)
    costs = arrayfun(@(t) TGo(srcSeq(t), order(t)) + S(order(t)) + TBack_kxn(srcSeq(t), order(t)), 1:q);

    % --- Timeline outputs ---
    res.ordersD{r}     = order;            % q cities
    res.routeCost(r)   = sum(costs);            % q costs
    res.startSource(r) = srcSeq(1);        % before city 1
    if q >= 2
        % sources before trips 2..q  → length q-1
        res.bridgeSrc{r} = srcSeq(2:end);
    else
        res.bridgeSrc{r} = [];
    end
    % After the last city (same as its start source in this model)
    res.endSource(r)   = srcSeq(end);
end
end