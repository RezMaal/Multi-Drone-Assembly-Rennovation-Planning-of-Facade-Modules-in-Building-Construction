function [res,sol] = solve_minmax_connected_sources(TGo, TBack, S, Adj, m, varargin)
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

% Order: [x(nm), y(nk), v(nmk), g(nm), f(ma), N(m), t(nm), T(1)]
off.x = 0;
off.y = off.x + nm;
off.v = off.y + nk;
off.g = off.v + nmk;
off.f = off.g + nm;
off.N = off.f + ma;
off.t = off.N + m;
off.T = off.t + nm;
numVars = off.T + 1;

idxX = @(j,r) off.x + (r-1)*n + j;
idxY = @(j,s) off.y + (s-1)*n + j;
idxV = @(j,r,s) off.v + ((r-1)*n + (j-1))*k + s;
idxG = @(j,r) off.g + (r-1)*n + j;
idxF = @(r,arc) off.f + (r-1)*a + arc;
idxN = @(r)     off.N + r;
idxt = @(j,r)   off.t + (r-1)*n + j;
idxT = @()      off.T + 1;

% --------- Objective & var types ---------
model.modelsense = 'min';
obj = zeros(numVars,1); obj(idxT()) = 1;      % minimize T

vtype = repmat('B', numVars, 1);
vtype(off.f+1 : off.f+ma) = 'C';              % flows
vtype(off.N+1 : off.N+m)  = 'C';              % N_r
vtype(off.t+1 : off.t+nm) = 'C';              % t_{ir}
vtype(idxT()) = 'C';                           % T

lb = zeros(numVars,1); ub = inf(numVars,1);
ub(1:off.f) = 1;                               % binaries ≤1
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
    for i=1:m
        ff=res.ordersD{1,i};
        Tg=TGo(:,ff); Tb=TBack(ff,:); ss=S(ff);
        OutA = ATSP_Im(Tg+ss, Tb);
        res.routeCost(1,i)=OutA.cost;
        res.startSource(1,i)=OutA.startSource;
        res.bridgeSrc{1,i}=OutA.bridgeSrc;
        res.endSource(1,i)=OutA.endSource;
        res.ordersD{1,i}=ff(OutA.orderD);
    end
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


function out = ATSP_Im(G, R, params)
% ATSP_fast  High-performance ATSP-path solver with heuristic warm start,
%            early stopping, and optional precomputed reduced costs.
%
%   out = ATSP_fast(G, R, params)
%
% INPUTS
%   G : k-by-n matrix, cost from source i to destination j  (e.g., TGo)
%   R : n-by-k matrix, cost from destination j to source i  (e.g., TBack)
%   params : (optional) struct with fields (all optional)
%       % --- Solver steering (overrides respected) ------------------------
%       .TimeLimit        (default 2)       % seconds
%       .MIPGap           (default 0.02)    % relative gap
%       .Heuristics       (default 0.8)
%       .MIPFocus         (default 1)
%       .Cuts             (default 2)
%       .Presolve         (default 2)
%       .OutputFlag       (default 0)
%       % --- Heuristic control -------------------------------------------
%       .HeuristicStart   (default true)    % build warm-start tour
%       .HeuristicOnly    (default false)   % skip MILP; return heuristic
%       .Max2OptIters     (default 50)      % 2-opt moves cap
%       .CutoffSlack      (default 1e-6)    % Cutoff = heur_cost*(1+slack)
%       .WarmStartOrder   (default [])      % existing tour to refine & use
%       % --- Precompute reduced costs (skip per-call build) ---------------
%       .alpha            (1-by-n)          % start -> j reduced cost
%       .beta             (1-by-n)          % j -> end reduced cost
%       .C                (n-by-n)          % j -> h reduced cost (diag Inf)
%
% OUTPUT struct 'out' (same shape as your original ATSP)
%   out.orderD      : 1-by-n visiting order of destinations (indices 1..n)
%   out.startSource : scalar index i0 of starting source
%   out.bridgeSrc   : 1-by-(n-1), sources used between consecutive cities
%   out.endSource   : scalar index of ending source
%   out.cost        : objective evaluated on the expanded S–D–S sequence
%   out.result      : raw Gurobi result struct (missing if HeuristicOnly)
%
% Notes
%   - Keeps the same scoring model as your ATSP: reduced costs use G and R.
%   - All speed improvements are annotated with '%% IMPROVED:' comments.
%
% -------------------------------------------------------------------------
% Validate & orient
[k, n] = size(G);
if size(R,1) ~= n || size(R,2) ~= k
    error('Size mismatch: G must be k-by-n and R must be n-by-k.');
end
G = double(G); R = double(R);
if nargin < 3 || isempty(params), params = struct; end

% Defaults (only set if user hasn't provided)
def = @(f,v) (isfield(params,f) && ~isempty(params.(f))) ;
if ~def('OutputFlag'),     params.OutputFlag   = 1;     end
if ~def('HeuristicStart'), params.HeuristicStart = true; end
if ~def('HeuristicOnly'),  params.HeuristicOnly  = false;end
if ~def('Max2OptIters'),   params.Max2OptIters   = 50;   end
if ~def('CutoffSlack'),    params.CutoffSlack    = 1e-6; end
% Aggressive but safe solver defaults for speed (overridden if provided)
if ~def('TimeLimit'),      params.TimeLimit   = 60;    end
if ~def('MIPGap'),         params.MIPGap      = 0.0001; end
if ~def('Heuristics'),     params.Heuristics  = 0.5;  end
if ~def('MIPFocus'),       params.MIPFocus    = 1;    end
if ~def('Cuts'),           params.Cuts        = 2;    end
if ~def('Presolve'),       params.Presolve    = 2;    end
if ~def('RINS'),           params.RINS = 20;    end
if ~def('NoRelHeurTime'),  params.NoRelHeurTime = params.TimeLimit/5;    end
    

% Early exits
if n == 0
    out = struct('orderD',[], 'startSource',[], 'bridgeSrc',[], 'endSource',[], 'cost',0);
    if ~params.HeuristicOnly, out.result = struct('status','OPTIMAL','objval',0); end
    return
elseif n == 1
    [~, i0] = min(G(:,1));
    [~, in] = min(R(1,:));
    out.orderD = 1;
    out.startSource = i0;
    out.bridgeSrc = [];
    out.endSource = in;
    out.cost = G(i0,1) + R(1,in);
    if ~params.HeuristicOnly, out.result = struct('status','OPTIMAL','objval',out.cost); end
    return
end

% -------------------------------------------------------------------------
% Reduced costs: alpha(j), beta(j), C(j,h)
% %% IMPROVED: optionally accept precomputed alpha/beta/C (skip heavy step)
havePre = isfield(params,'alpha') && isfield(params,'beta') && isfield(params,'C') ...
          && isequal(size(params.alpha),[1 n]) && isequal(size(params.beta),[1 n]) ...
          && isequal(size(params.C),[n n]);
if havePre
    alpha = params.alpha; beta = params.beta; C = params.C;
else
    % %% IMPROVED: vectorized, cache-friendly precomputation
    alpha = min(G, [], 1);          % 1 x n
    beta  = min(R, [], 2)';         % 1 x n
    C = inf(n,n);
    for i = 1:k                     % accumulate elementwise min over sources
        cand = R(:,i) + G(i,:);
        C = min(C, cand);
    end
    C(1:n+1:end) = inf;             % forbid self-arcs
end

% Guard: infeasible if both in- and out- are blocked for some node
if any(~isfinite(alpha) & all(~isfinite(C),1))
    error('ATSP_fast: infeasible start into some destination (alpha & row(C) Inf).');
end
if any(~isfinite(beta) & all(~isfinite(C),2)')
    error('ATSP_fast: infeasible exit from some destination (beta & col(C) Inf).');
end

% -------------------------------------------------------------------------
% Heuristic warm start (greedy best-insertion + 2-opt)
% %% IMPROVED: fast tour to seed model.start and cutoff
start_order = [];
if params.HeuristicStart
    if isfield(params,'WarmStartOrder') && ~isempty(params.WarmStartOrder)
        start_order = params.WarmStartOrder(:)';
        start_order = start_order(ismember(start_order,1:n)); % sanitize
        if numel(start_order) ~= n
            % Repair by adding missing nodes with best insertion
            start_order = best_insertion_complete(start_order, alpha, C, beta);
        end
    else
        % Build from scratch
        % Start at best alpha; insert remaining by best increase
        [~, j0] = min(alpha);
        start_order = best_insertion_complete(j0, alpha, C, beta);
    end
    % 2-opt improvement
    start_order = two_opt_improve(start_order, C, alpha, beta, params.Max2OptIters);
end

% If heuristic-only requested, just expand to sources and return
if params.HeuristicOnly
    [startSource, bridgeSrc, endSource, heur_cost] = expand_to_sources(start_order, G, R);
    out.orderD = start_order;
    out.startSource = startSource;
    out.bridgeSrc = bridgeSrc;
    out.endSource = endSource;
    out.cost = heur_cost;
    return
end

% -------------------------------------------------------------------------
% Build MILP (MTZ path ATSP)
% x0j  : n start arcs; xjh : n*(n-1) internal arcs; xjt : n end arcs; u : n
m1 = n; m2 = n*(n-1); m3 = n; mX = m1+m2+m3; Nvar = mX+n;

idx.x0  = @(j) j;
map_jh = @(j,h) ((j-1)*(n-1) + (h - (h>j)));
idx.xjh = @(j,h) m1 + map_jh(j,h);
idx.xjt = @(j)   m1 + m2 + j;
idx.u   = @(j)   mX + j;

obj = zeros(Nvar,1);
obj(1:m1) = alpha(:);
ptr = m1;
for j = 1:n
    for h = 1:n
        if h==j, continue; end
        ptr = ptr + 1;
        obj(ptr) = C(j,h);
    end
end
obj(m1+m2+1 : m1+m2+m3) = beta(:);

vtype = repmat('B', Nvar, 1);
lb = zeros(Nvar,1); ub = ones(Nvar,1);
for j = 1:n, vtype(idx.u(j)) = 'C'; lb(idx.u(j)) = 1; ub(idx.u(j)) = n; end

% Disallow infinite arcs
mask_alpha = ~isfinite(alpha); ub(arrayfun(idx.x0, find(mask_alpha))) = 0;
mask_beta = ~isfinite(beta);   ub(m1+m2 + find(mask_beta)) = 0;
for j = 1:n
    bad = ~isfinite(C(j,:)); bad(j) = false;
    if any(bad)
        cols = arrayfun(@(h) idx.xjh(j,h), find(bad));
        ub(cols) = 0;
    end
end

% Constraints
nRows = 2 + 2*n + n*(n-1);
estNnz = 2*n + 2*n*n + 3*n*(n-1);
I = zeros(estNnz,1); J = zeros(estNnz,1); V = zeros(estNnz,1);
rhs = zeros(nRows,1); sense = repmat('=', nRows,1);
ri = 0; nz = 0;
addRow = @(cols,vals,sg,b) add_row(cols,vals,sg,b);
function add_row(cols, vals, sg, b)
    nonzeros = numel(cols);
    if nz + nonzeros > numel(I)
        grow = max(nonzeros, ceil(0.2*numel(I)));
        I(end+grow) = 0; J(end+grow) = 0; V(end+grow) = 0; %#ok<AGROW>
    end
    I(nz+(1:nonzeros)) = ri + 1;
    J(nz+(1:nonzeros)) = cols(:);
    V(nz+(1:nonzeros)) = vals(:);
    nz = nz + nonzeros;
    rhs(ri+1) = b;
    sense(ri+1) = sg;
    ri = ri + 1;
end

% Start once
addRow(1:m1, ones(1,m1), '=', 1);
% End once
addRow(m1+m2+1 : m1+m2+m3, ones(1,m3), '=', 1);
% In-degree
for j = 1:n
    cols = idx.x0(j); vals = 1;
    for h = 1:n
        if h==j, continue; end
        cols(end+1,1) = idx.xjh(h,j); vals(end+1,1) = 1; %#ok<AGROW>
    end
    addRow(cols, vals, '=', 1);
end
% Out-degree
for j = 1:n
    cols = idx.xjt(j); vals = 1;
    for h = 1:n
        if h==j, continue; end
        cols(end+1,1) = idx.xjh(j,h); vals(end+1,1) = 1; %#ok<AGROW>
    end
    addRow(cols, vals, '=', 1);
end
% MTZ
for j = 1:n
    for h = 1:n
        if j==h, continue; end
        addRow([idx.u(j), idx.u(h), idx.xjh(j,h)], [1, -1, n], '<', n-1);
    end
end
I = I(1:nz); J = J(1:nz); V = V(1:nz);
model.A     = sparse(I,J,V,ri,Nvar);
model.rhs   = rhs; model.sense = sense; model.obj = obj;
model.vtype = vtype; model.lb = lb; model.ub = ub;
model.modelsense = 'min';

% ---- Heuristic warm start as MIP start + Cutoff ------------------------
if params.HeuristicStart
    if isempty(start_order)
        [~, j0] = min(alpha);
        start_order = best_insertion_complete(j0, alpha, C, beta);
        start_order = two_opt_improve(start_order, C, alpha, beta, params.Max2OptIters);
    end
    [x0s, xjhs, xjts, heur_cost] = encode_tour(start_order, idx, n);
    model.start = nan(Nvar,1);
    model.start(1:m1) = x0s(:);
    % fill xjh in block order
    ptr = m1;
    for j = 1:n
        for h = 1:n
            if h==j, continue; end
            ptr = ptr + 1;
            model.start(ptr,1) = xjhs(j,h);
        end
    end
    model.start(m1+m2+1 : m1+m2+m3) = xjts(:);
    % %% IMPROVED: set Cutoff slightly above heuristic cost
    if ~isfield(params,'Cutoff')
        params.Cutoff = (1+params.CutoffSlack) * heur_cost;
    end
end

% ---- Solve --------------------------------------------------------------
result = gurobi(model, params);

% ---- Decode solution ----------------------------------------------------
x = result.x;
% start city
[~, j0] = max(x(1:m1));
order = zeros(1,n); order(1) = j0;
% adjacency
XJH = reshape(x(m1+1:m1+m2), n-1, n)';  % j-row with n-1 entries
for t = 1:n-1
    j = order(t);
    row = XJH(j,:);
    pos = 0; nxt = -1;
    for h = 1:n
        if h==j, continue; end
        pos = pos + 1;
        if row(pos) > 0.5, nxt = h; break; end
    end
    order(t+1) = nxt;
end

% Expand to sources
[startSource, bridgeSrc, endSource, exact_cost] = expand_to_sources(order, G, R);

% Output
out.orderD      = order;
out.startSource = startSource;
out.bridgeSrc   = bridgeSrc;
out.endSource   = endSource;
out.cost        = exact_cost;
out.result      = result;

% ========================= Helper functions ==============================
function order = best_insertion_complete(seed, alpha, C, beta)
% Build a complete order from initial list 'seed' using best insertion.
    n = numel(alpha);
    if isempty(seed)
        [~, seed] = min(alpha);
    end
    order = seed(:)';
    left = setdiff(1:n, order);
    % path cost function over reduced costs
    function inc = ins_cost(path, pos, v)
        % insert v between path(pos) -> path(pos+1), with start/end handling
        if isempty(path)
            inc = alpha(v) + beta(v); return
        end
        if pos==0
            % start -> v -> path(1)
            inc = alpha(v) + C(v, path(1)) - alpha(path(1));
        elseif pos==numel(path)
            % path(end) -> v -> end
            inc = C(path(end), v) + beta(v) - beta(path(end));
        else
            % path(pos) -> v -> path(pos+1) instead of path(pos)->path(pos+1)
            inc = C(path(pos), v) + C(v, path(pos+1)) - C(path(pos), path(pos+1));
        end
    end
    while ~isempty(left)
        bestInc = inf; bestV = left(1); bestPos = 0;
        % evaluate all insertions (n^2 worst-case; small n_r is fine)
        for v = left
            L = numel(order);
            % try positions 0..L (0 means before first)
            inc0 = ins_cost(order, 0, v);
            incL = ins_cost(order, L, v);
            if inc0 < bestInc, bestInc = inc0; bestV = v; bestPos = 0; end
            if incL < bestInc, bestInc = incL; bestV = v; bestPos = L; end
            for pos = 1:L-1
                inc = ins_cost(order, pos, v);
                if inc < bestInc, bestInc = inc; bestV = v; bestPos = pos; end
            end
        end
        % apply
        if bestPos==0
            order = [bestV, order];
        elseif bestPos==numel(order)
            order = [order, bestV];
        else
            order = [order(1:bestPos), bestV, order(bestPos+1:end)];
        end
        left = setdiff(left, bestV);
    end
end

function order = two_opt_improve(order, C, alpha, beta, maxIters)
% Simple 2-opt on the ATSP-path with start/end wrappers.
    if maxIters<=0, return; end
    n = numel(order);
    % helper for path delta with endpoints
    function d = delta(i, k)
        % reverse segment (i..k), compute delta in reduced cost
        a = order; j = a; %#ok<NASGU>
        d = 0;
        if i==1
            % start -> order(1) becomes start -> order(k)
            d = d + (alpha(order(k)) - alpha(order(1)));
        else
            d = d + (C(order(i-1), order(k)) - C(order(i-1), order(i)));
        end
        if k==n
            % order(n) -> end becomes order(i) -> end
            d = d + (beta(order(i)) - beta(order(n)));
        else
            d = d + (C(order(i), order(k+1)) - C(order(k), order(k+1)));
        end
        % inside segment current edges replaced by reversed ones
        for t = i:k-1
            d = d + (C(order(t+1), order(t)) - C(order(t), order(t+1)));
        end
    end
    it = 0; improved = true;
    while improved && it < maxIters
        improved = false; it = it + 1;
        for i = 1:n-1
            for k = i+1:n
                d = delta(i,k);
                if d < -1e-12
                    order(i:k) = order(k:-1:i);
                    improved = true;
                end
            end
        end
    end
end

function [x0s, xjhs, xjts, total] = encode_tour(order, idx, n)
% Encode a path 'order' into x0, xjh, xjt binaries, and compute reduced cost.
    x0s = zeros(n,1); xjhs = zeros(n,n); xjts = zeros(n,1);
    x0s(order(1)) = 1;
    for t = 1:numel(order)-1
        j = order(t); %#ok<NASGU>
    end
    for t = 1:numel(order)-1
        j = order(t); h = order(t+1);
        xjhs(j,h) = 1;
    end
    xjts(order(end)) = 1;
    % Reduced cost is alpha(first) + sum C(j->h) + beta(last)
    total = 0;
    total = total + alpha(order(1));
    for t = 1:numel(order)-1
        total = total + C(order(t), order(t+1));
    end
    total = total + beta(order(end));
end

function [i0, bridge, in_src, cost] = expand_to_sources(order, G, R)
% Expand reduced-cost path into concrete sources for start/bridges/end.
    [~, i0] = min(G(:, order(1)));
    nloc = numel(order);
    bridge = zeros(1, nloc-1);
    for t = 1:nloc-1
        j = order(t); h = order(t+1);
        [~, it] = min(R(j,:) + G(:,h)');
        bridge(t) = it;
    end
    [~, in_src] = min(R(order(end),:));
    % Compute exact cost
    cost = G(i0, order(1));
    for t = 1:nloc-1
        j = order(t); h = order(t+1); it = bridge(t);
        cost = cost + R(j, it) + G(it, h);
    end
    cost = cost + R(order(end), in_src);
end

end