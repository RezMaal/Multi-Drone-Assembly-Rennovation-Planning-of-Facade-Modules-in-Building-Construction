function [Rb,Modules,Costs,Outer,Inner,Att,ID_P,pg]=Modularization_Facade(P,varargin)
%% Input Parameters
% % Boundary Detection Inputs (below are default parameters):
% dd=[4,4,4,4,6,8,8]/100; % Pixel Size
% rr=[5,5,5,5,10,15,15]/100; % Feature Search Radius
% PlaneIter=50; % Robust Plane Fitting Iteration
% 
% % Regularization Parameters
% Regularize=1; % Regularization True
% Tm_Re=60; % if Regularize True, Optimization Time Limit
% Out_Re=1; % if Regularize True, Optimization Output Flag
% Gap_Re=10^-3; % if Regularize True, Optimization Gap Tolerance
% 
% % Plot 1: Detected Boundaries Parameter Inputs (default parameters below):
% Plot_1=1; % Plot 1 True
% sz_1=20; % Markersize
% off_1=0.5; % Figure side Offset
% 
% % Modularization Inputs (Default Parameters)
% wm=2; % Minimum Width of Module
% wM=2; % Maximum Width of Module
% hm=1.5; % Minimum Height of Module
% hM=1.5; % Maximum Height of Module
% dx=0.2; % Modularization Precision in X
% dy=0.2; % Accuracy of Modularization in Y
% 
% % Modularization Optimization Parameters:
% Tm_Mod=240; % Modularization Optimization Time Limit
% Out_Mod=1; % Modularization Optimization Output Flag
% Gap_Mod=10^-10; % Modularization Optimization Gap Tolerance
% 
% % Plot 2: Modularized Shape Parameter Inputs (default parameters below):
% Plot_2=1; % Plot 2 True
% sz_2=20; % Boundary Markersize
% off_2=0.5; % Figure Side Offset
% L_sz=2; % Module Line Weight
% P_sz=1; % Module Point Weight

% --- Defaults Values -------------------------
defaults = struct( ...
    'dd',       [4,4,4,4,6,8,8]/100, ...        % Pixel Size
    'rr',       [5,5,5,5,10,15,15]/100, ...     % Feature Search Radius
    'PlaneIter',50, ...                          % Robust Plane Fitting Iteration
    ...
    'Regularize',1, ...                          % Regularization True
    'Tm_Re',    60, ...                          % Optimization Time Limit (Reg)
    'Out_Re',   0, ...                           % Optimization Output Flag (Reg)
    'Gap_Re',   1e-3, ...                        % Optimization Gap Tolerance (Reg)
    ...
    'Plot_1',   1, ...                           % Plot 1 True
    'sz_1',     5, ...                          % Markersize
    'off_1',    0.5, ...                         % Figure side Offset
    ...
    'wm',       0.9, ...                        % Min Module Width
    'wM',       2.1, ...                        % Max Module Width
    'hm',       0.9, ...                         % Min Module Height
    'hM',       2.1, ...                        % Max Module Height
    'dx',       0.2, ...                         % Precision X
    'dy',       0.2, ...                         % Accuracy Y
    ...
    'Tm_Mod',   240, ...                         % Mod. Optimization Time Limit
    'Out_Mod',  0, ...                           % Mod. Optimization Output Flag
    'Gap_Mod',  1e-10, ...                       % Mod. Optimization Gap Tolerance
    ...
    'Plot_2',   1, ...                           % Plot 2 True
    'sz_2',     5, ...                          % Boundary Markersize
    'off_2',    0.5, ...                         % Figure Side Offset
    'L_sz',     2, ...                           % Module Line Weight
    'P_sz',     1, ...
    ...
    'Plot_3',   1, ...                           % Plot 3 True
    'sz_3',     35);                          % Min Size of Graph Node);                              % Module Point Weight

% --- Parse/validate overrides --------------------------------
opts = parseParams(defaults, varargin{:});

% (Optional) normalize logical-ish fields
opts.Plot_1    = logical(opts.Plot_1);
opts.Plot_2    = logical(opts.Plot_2);
opts.Regularize= logical(opts.Regularize);

struct2vars(opts);

% --- Regularization solve ---
params_re = struct();
params_re.TimeLimit  = double(Tm_Re(1));          % not Tm_Reg
params_re.OutputFlag = double(Out_Re(1) ~= 0);
params_re.MIPGap     = double(Gap_Re(1));

% --- Modularization solve ---
params_mod = struct();
params_mod.TimeLimit  = double(Tm_Mod(1));
params_mod.OutputFlag = double(Out_Mod(1) ~= 0);  % <-- key fix
params_mod.MIPGap     = double(Gap_Mod(1));


%% Plane and Boundary Detection
mp=max(P(:,5)); ID_P=cell(mp+1,1); Att=zeros(mp+1,4); RMSE=zeros(mp+1,1); Cent=zeros(mp+1,3);

Inner=cell(mp+1,1); Outer=Inner; pg=cell(mp+1,1); pin=pg;
startTime = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
    
% Display start message
fprintf('[%s] Start of Robust Plane Fitting...\n', startTime);

for i=1:mp+1
    i1= find(P(:,5)==i-1); [Att(i,:),RMSE(i,1),ID1]=RobustConcentrationPlaneFit(P(i1,1:3),PlaneIter);
    ID_P{i,1}=i1(ID1); Cent(i,:)=mean(P(i1(ID1),1:3));
end
Rr = bestRotation_flattenZ_svd(Att(:,1:3), max(RMSE)./RMSE);
Att(:,1:3)=Att(:,1:3)*Rr; Att(:,4)=-dot(Att(:,1:3),Cent*Rr,2);

%% Surface Ordering
ppx=[];
for i=1:mp+1
    Rt=Rodrigues(Att(i,1:3),[1,0,0]);
    pp=P(ID_P{i,1},1:3)*Rr*Rt;
    px1=(min(pp(:,2)):rr(i):max(pp(:,2)))';
    px=[repmat(mean(pp(:,1)),length(px1),1),px1,repmat(mean(pp(:,3)),length(px1),1)]*Rt';
    ppx=[ppx;px];
end
iix = mstPreorderOrder(ppx(:,1:2));
ic=knnsearch(ppx(iix,1:2),Cent(:,1:2),'k',1);
[~,ib]=sort(ic);
cc=Cent(ib(2:end),:)-Cent(ib(1:end-1),:); cc=[cc(1,:);cc];
aa=cross([Att(ib,1:2),zeros(mp+1,1)],[cc(:,1:2),zeros(mp+1,1)]);
Att(aa(:,3)<0,:)=-Att(aa(:,3)<0,:);

startTime = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');   
fprintf('[%s] Start of Plane Facade Boundary Detection...\n', startTime);
for i=1:mp+1
    Rt=Rodrigues(Att(i,1:3),[1,0,0]);
    pp=P(ID_P{i,1},1:3)*Rr*Rt;
    [rectx,recty] = minboundrect(pp(:,2),pp(:,3)); lr=max(max(rectx)-min(rectx),max(recty)-min(recty));
    num=ceil(lr/dd(i));
    [R, c] = alignToNearestAxis([rectx,recty]); [~,mx]=max(abs(R(:,1))); if R(mx,1)<0, R=-R; end
    pp_R=pp(:,2:3)*R'+c; pp_R=pp_R-min(pp_R);
    parts = pointBoundariesImage(pp_R(:,1), pp_R(:,2), ...
            'Mode','distance', ...
            'Radius', rr(i), ...         % world units
            'NumPixels', num, ...
            'MorphClose', 2, ...
            'MinRegionArea', 2, ...  % drop tiny CCs (approx via px→world)
            'MinBoundaryArea', 1, ...% drop tiny loops (world-units²) AFTER tracing
            'BoundarySnapTol', [], ...  % auto from raster pixel size; or set e.g. 2e-3
            'KeepHoles', true, ...
            'Plot', false);
    spr=size(parts,1);
    warning 'off'
    if spr==1
        O1=parts.holePolys{1,1};
    else
        O1=parts(1).holePolys{1,1};
        for j=2:spr
            I1=parts(j).outerPoly;
            iix = mstPreorderOrder(I1(:,1:2)); I2=I1(iix,:);
            if Regularize==1
                idx = minCoverSeg_confidence(I2(:,1:2), 0.25, 0.9, ...
                      'Knn', 5, 'Params', params_re, 'Verbose', 1);
            else
                idx=1:length(iix);
            end
            Inner{i,j-1}=I1(sort(iix(idx)),1:2);
            pin{i,j-1}=polyshape(I1(:,1:2),'Simplify',true);
        end
    end
    iix = mstPreorderOrder(O1(:,1:2)); O2=O1(iix,:);
    if Regularize==1
        idx = minCoverSeg_confidence(O2, 0.25, 0.9, ...
              'Knn', 5, 'Params', params_re, 'Verbose', 1);
    else
        idx=1:length(iix);
    end
    Outer{i,1}=O1(sort(iix(idx)),:);
    pg{i,1}=polyshape(O1(:,1:2),'Simplify',true);
end



%% Plot 1: Detected Boundaries
if Plot_1==1
    figure; hold on; axis equal; ct=0; 
    for ll=1:mp+1
        i=ib(ll);
        ptt=pg{i,1}; ptt.Vertices(:,1)=ptt.Vertices(:,1)+ct; plot(ptt);
        scatter(Outer{i,1}(:,1)+ct,Outer{i,1}(:,2),sz_1,rand(1,3));
        if ~isempty(Inner{i,1})
           for j=1:size(Inner,2)
               plot(polyshape(Inner{i,j}(:,1)+ct,Inner{i,j}(:,2),'Simplify',true));
               scatter(Inner{i,j}(:,1)+ct,Inner{i,j}(:,2),sz_1,rand(1,3));
           end
        end
        ct=ct+max(Outer{i,1}(:,1))+off_1; 
    end
end

%% MILP Modularization
startTime = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');   
fprintf('[%s] Start of Rectangular Modularization ...\n', startTime);
Modules=cell(mp+1,1); ppg=cell(mp+1,1); Costs=ppg;
for ll=1:mp+1
    i=ib(ll);
    if isempty(pin{i,1})
        ppg=pg{i,1};
    else
        ptt=pin{i,1};
        for tt=2:length(pin(i,:))
            ptt=union(ptt,pin{i,tt});
        end
        ppg=subtract(pg{i,1},ptt);
        pg{i,1}=ppg;
    end
    xg=0:dx:max(ppg.Vertices(:,1)); yg=0:dy:max(ppg.Vertices(:,2));
    data = makeRectData(ppg, 'Xgrid',xg,'Ygrid',yg);
    ss=zeros(size(data.Sij));
    ss(round(data.Sij,5)./round(data.Aij,5)==1)=1; ss=ss';
    params_mod.MIPFocus= 1;
    params_mod.Heuristics= 0.2;
    params_mod.Cuts= 2;

    % Call the lexicographic solver
    [~, ~, rects, ~, st] = minK_maxCoverage_lexi(ss, floor(wm/dx), floor(wM/dx), floor(hm/dy), floor(hM/dy), params_mod);
    aa=table2array(rects);
    Modules{ll,1}=[aa(:,2).*dx-dx,aa(:,1).*dy-dy,aa(:,4).*dx,aa(:,3).*dy];
    Costs{ll,1}=[area(ppg)*(1-st.coverFracStar),st.coverFracStar];
    startTime = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');   
    fprintf('[%s] Facade %d: Modularization Completed ...\n', startTime,ll);
end

%% Plot 2: Modularization
if Plot_2==1
    figure; hold on; axis equal; axis tight; box on; ct=0; 
    for ll=1:mp+1 
        i=ib(ll);
        pc=pg{i,1};
        pc.Vertices(:,1)=pc.Vertices(:,1)+ct; plot(pc);
        scatter(Outer{i,1}(:,1)+ct,Outer{i,1}(:,2),sz_2,'k','filled'); 
        if ~isempty(Inner{i,1}) 
           for j=1:size(Inner,2)
                scatter(Inner{i,j}(:,1)+ct,Inner{i,j}(:,2),sz_2,'k','filled'); 
           end 
        end
        for kk=1:length(Modules{ll,1})
            drawrectangle('Position',Modules{ll,1}(kk,:)+[ct,0,0,0],'Color',rand(1,3),'MarkerSize',P_sz,'LineWidth',L_sz);
        end
        ct=ct+max(Outer{i,1}(:,1))+off_2;
    end
end
Outer=Outer(ib);
Inner=Inner(ib,:);
Att=Att(ib,:);
ID_P=ID_P(ib);

%% Module Connectivity Mapping
startTime = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');   
fprintf('[%s] Building Facade Connectivity Map...\n', startTime);
n = numel(Modules); xn=zeros(n,1); 
X_p = cell(n,1); M_n = cell(n,1); X_Gp = X_p; ct = 0;
for i=1:n
    X_p{i,1} = Modules{i,1}(:,1:2)+Modules{i,1}(:,3:4)/2+[xn(i)-min(Modules{i,1}(:,1)),0];
    X_Gp{i,1} = Modules{i,1}(:,1:2)+Modules{i,1}(:,3:4)/2+[ct,0];
    M_n{i,1} = Modules{i,1}+[xn(i)-min(Modules{i,1}(:,1)),0,0,0];
    xn(i+1,:) = xn(i)+max(Modules{i,1}(:,1)+Modules{i,1}(:,3)); ct=ct+max(Outer{i,1}(:,1))+0.5;
end
Md = cell2mat(M_n);
E = Build_Connectivity(Md,0.01);
A = graph(E(:,1),E(:,2)); PM=cell2mat(X_p);
G = connectComponentsByNearest(PM, A);
[s, t] = findedge(G);    % s(k), t(k) are the endpoints of edge k
Rb.E = [s t];
Rb.X_n = cell2mat(X_Gp);

%% Plot 3: Connectivity Graph
if Plot_3==1
    sigma = sum(pdist2(Rb.X_n,Rb.X_n),2)/n;
    Rb.s = sigma.*Md(:,3).*Md(:,4);
    ms=min(Rb.s);Ms=max(Rb.s); ss=(Rb.s-ms)/(Ms-ms)*2*sz_3+sz_3;
    aa = adjacency(G);
    gplot(aa, Rb.X_n, ['-','k']); axis equal; axis tight; grid on
    scatter(Rb.X_n(:,1),Rb.X_n(:,2),ss,'k','filled');
end

end

function opts = parseParams(defaults, varargin)
% parseParams  Merge defaults with user-specified Name/Value pairs or struct.
% - Case-insensitive, partial name matching
% - Validates each parameter
%
% Usage:
%   opts = parseParams(defaults);                      % returns defaults
%   opts = parseParams(defaults, 'PlaneIter', 100);    % name/value
%   opts = parseParams(defaults, struct('Plot_1', 0)); % struct

% Allow struct input as a convenience
if isscalar(varargin) && isstruct(varargin{1})
    S = varargin{1};
    fn = fieldnames(S);
    nv = reshape([fn.'; struct2cell(S).'], 1, []);
else
    nv = varargin;
end

p = inputParser;
p.CaseSensitive   = true;
p.PartialMatching = false;
p.KeepUnmatched   = false;

% --- Validators ------------------------------------------------
isPosScalar  = @(x) isnumeric(x) && isscalar(x) && isfinite(x) && (x>0);
isNonNegScal = @(x) isnumeric(x) && isscalar(x) && isfinite(x) && (x>=0);
isBinary     = @(x) (islogical(x) && isscalar(x)) || (isnumeric(x) && isscalar(x) && any(x==[0 1]));
isGap        = @(x) isnumeric(x) && isscalar(x) && isfinite(x) && (x>0) && (x<=1);
isLen7VecPos = @(x) isnumeric(x) && isvector(x) && numel(x)==7 && all(isfinite(x)) && all(x>0);

% --- Register parameters + validation --------------------------
addParameter(p,'dd',       defaults.dd,       isLen7VecPos);
addParameter(p,'rr',       defaults.rr,       isLen7VecPos);
addParameter(p,'PlaneIter',defaults.PlaneIter,isPosScalar);

addParameter(p,'Regularize',defaults.Regularize,isBinary);
addParameter(p,'Tm_Re',     defaults.Tm_Re,     isPosScalar);
addParameter(p,'Out_Re',    defaults.Out_Re,    isBinary);
addParameter(p,'Gap_Re',    defaults.Gap_Re,    isGap);

addParameter(p,'Plot_1', defaults.Plot_1, isBinary);
addParameter(p,'sz_1',   defaults.sz_1,   isPosScalar);
addParameter(p,'off_1',  defaults.off_1,  isNonNegScal);

addParameter(p,'wm', defaults.wm, isPosScalar);
addParameter(p,'wM', defaults.wM, isPosScalar);
addParameter(p,'hm', defaults.hm, isPosScalar);
addParameter(p,'hM', defaults.hM, isPosScalar);
addParameter(p,'dx', defaults.dx, isPosScalar);
addParameter(p,'dy', defaults.dy, isPosScalar);

addParameter(p,'Tm_Mod', defaults.Tm_Mod, isPosScalar);
addParameter(p,'Out_Mod',defaults.Out_Mod,isBinary);
addParameter(p,'Gap_Mod',defaults.Gap_Mod,isGap);

addParameter(p,'Plot_2', defaults.Plot_2, isBinary);
addParameter(p,'sz_2',   defaults.sz_2,   isPosScalar);
addParameter(p,'off_2',  defaults.off_2,  isNonNegScal);
addParameter(p,'L_sz',   defaults.L_sz,   isPosScalar);
addParameter(p,'P_sz',   defaults.P_sz,   isPosScalar);

addParameter(p,'Plot_3', defaults.Plot_3, isBinary);
addParameter(p,'sz_3',   defaults.sz_3,   isPosScalar);

% --- Parse -----------------------------------------------------
parse(p, nv{:});
opts = p.Results;

% --- Extra logical normalization (optional) --------------------
logicalFields = {'Plot_1','Plot_2','Regularize','Out_Re','Out_Mod'};
for k = 1:numel(logicalFields)
    f = logicalFields{k};
    if isfield(opts,f), opts.(f) = logical(opts.(f)); end
end

% --- Optional: simple sanity checks across parameters ----------
% e.g., ensure mins <= maxes
if opts.wm > opts.wM
    error('ParameterError:wmwM','wm (%.3g) must be <= wM (%.3g).', opts.wm, opts.wM);
end
if opts.hm > opts.hM
    error('ParameterError:hmhM','hm (%.3g) must be <= hM (%.3g).', opts.hm, opts.hM);
end
end

% File: struct2vars.m
function struct2vars(S)
%STRUCT2VARS Assign fields of struct S as variables in the caller function.
% WARNING: Overwrites any existing variables with the same names in the caller.
    names = fieldnames(S);
    for k = 1:numel(names)
        assignin('caller', names{k}, S.(names{k}));
    end
end

%% Helper Functions
function E = Build_Connectivity(partsGA1,tol)
E=[]; n=size(partsGA1,1);
for i=1:n
    ix2=find(abs(partsGA1(i,1)-partsGA1(i+1:end,1)-partsGA1(i+1:end,3))<=tol);
    ix4=find(abs(partsGA1(i,1)-partsGA1(i+1:end,1)+partsGA1(i,3))<=tol);
    ixx=union(ix2,ix4);
    ipp=(i+1:n);
    if ~isempty(ixx)
        for j=1:length(ixx)
            y1=partsGA1(i,2); y2=partsGA1(i,2)+partsGA1(i,4);
            y3=partsGA1(ipp(ixx(j)),2); y4=partsGA1(ipp(ixx(j)),2)+partsGA1(ipp(ixx(j)),4);
            if (y2 > y3) && (y1 < y4)
                E=[E;[i,ipp(ixx(j))]];
            end
        end
    end
    ix2=find(abs(partsGA1(i,2)-partsGA1(i+1:end,2)-partsGA1(i+1:end,4))<=tol);
    ix4=find(abs(partsGA1(i,2)-partsGA1(i+1:end,2)+partsGA1(i,4))<=tol);
    ixx=union(ix2,ix4);
    if ~isempty(ixx)
        for j=1:length(ixx)
            y1=partsGA1(i,1); y2=partsGA1(i,1)+partsGA1(i,3);
            y3=partsGA1(ipp(ixx(j)),1); y4=partsGA1(ipp(ixx(j)),1)+partsGA1(ipp(ixx(j)),3);
            if (y2 > y3) && (y1 < y4)
                E=[E;[i,ipp(ixx(j))]];
            end
        end
    end
end
end
function [R, c, ptsAligned] = alignToNearestAxis(pts)
% alignToNearestAxis  Rotate a 2D rectangle so one side aligns to the nearest axis
% Inputs:
%   pts        4×2 array of corner coordinates [x, y]
% Outputs:
%   R          2×2 rotation matrix
%   ptsAligned 4×2 array of rotated, axis-aligned corners

  % 1. Compute centroid
  c = mean(pts, 1);

  % 2. List edges (wrap-around)
  edges = [1 2; 2 3; 3 4; 4 1];

  % 3. Compute orientations and distances
  thetas = zeros(4,1);
  d_x = zeros(4,1);
  d_y = zeros(4,1);
  for k = 1:4
    v = pts(edges(k,2),:) - pts(edges(k,1),:);
    theta = atan2(v(2), v(1));                   % edge angle
    thetas(k) = theta;
    theta0 = atan2(sin(theta), cos(theta));      % normalize to [-π, π]
    d_x(k)  = abs(theta0);                       % horizontal distance
    d_y(k)  = abs(abs(theta0) - pi/2);           % vertical distance
  end

  % 4. Pick minimal-rotation edge
  [~, kmin] = min(min(d_x, d_y));
  theta_sel = thetas(kmin);
  if d_x(kmin) <= d_y(kmin)
    phi = -theta_sel;            % align to x-axis
  else
    phi = pi/2 - theta_sel;      % align to y-axis
  end

  % 5. Build rotation matrix
  R = [cos(phi), -sin(phi);
       sin(phi),  cos(phi)];

  % 6. Rotate about centroid
  pts0       = pts - c;
  ptsRotated = (R * pts0')';
  ptsAligned = ptsRotated + c;
end
function [R, u, s] = bestRotation_flattenZ_svd(N, w)
% Minimizes sum (z after rotation)^2 with SVD (stable).
% N: Nx3 (rows=normals). w: optional Nx1 weights (nonnegative).

    % Normalize directions
    N = normalizeRows(N);

    % Optional weights -> scale rows by sqrt(w)
    if nargin>=2 && ~isempty(w)
        w = w(:); w = w / max(eps, sum(w));
        N = N .* sqrt(w);
    end

    % SVD of data (econ is fine)
    [~, S, V] = svd(N, 'econ');
    [sMin, j] = min(diag(S));  %#ok<ASGLU>
    u = V(:, j);               % minimizing direction in R^3
    [~,r] = max(abs(u));
    if u(r)<0
        u=-u;
    end
    % Rotate u -> e3
    R = Rodrigues(u',[0,0,1]);
    s = diag(S);
end

function X = normalizeRows(X)
    n = sqrt(sum(X.^2,2)); n(n==0)=1; X = X ./ n;
end
function Rot=Rodrigues(N,V)
    teta=acos(dot(N,V)/norm(N)/norm(V));
    k=cross(N,V);
    if norm(k)==0
        Rot=eye(length(N),length(N));
    else
        k=k/norm(k);
        K=[0,-k(3),k(2);k(3),0,-k(1);-k(2),k(1),0];
        R=eye(3)+sin(teta)*K+(1-cos(teta))*K^2;
        Rot=R';
    end
end

function [G2, bridges] = connectComponentsByNearest(P, G)
%CONNECTCOMPONENTSBYNEAREST Repeatedly add the shortest inter-component edge.
% Inputs:
%   P : N-by-d coordinates (double)
%   G : graph with N nodes (undirected)
% Outputs:
%   G2      : connected graph after adding bridges
%   bridges : table of added edges [s t w] (s,t are node ids; w is distance)

    arguments
        P double
        G graph
    end
    N = size(P,1);
    if numnodes(G) ~= N
        error('G must have %d nodes to match P.', N);
    end

    % Pairwise distances; set self to inf
    D = squareform(pdist(P));
    D(1:N+1:end) = inf;

    G2 = G;
    [comp, ~] = conncomp(G2);  % initial components
    K = max(comp);

    % Block out all within-component pairs initially
    u = unique(comp);
    for c = u
        idx = find(comp==c);
        D(idx, idx) = inf;
    end

    sAdd = []; tAdd = [];

    while K > 1
        % Find globally closest cross-component pair
        [wmin, linIdx] = min(D(:));
        if ~isfinite(wmin)
            error('No finite cross-component distances found. Check your P.');
        end
        [i, j] = ind2sub([N N], linIdx);

        % Add the bridge and record it
        G2 = addedge(G2, i, j);
        sAdd(end+1,1) = i; %#ok<AGROW>
        tAdd(end+1,1) = j; %#ok<AGROW>

        % Recompute components (simple & robust)
        comp = conncomp(G2);
        K = max(comp);

        % Update D: disable all now-within-component pairs
        % (simple but safe; for speed you could only update the merged block)
        D(:) = D(:); % no-op; keeps D same size
        for c = unique(comp)
            idx = find(comp==c);
            D(idx, idx) = inf;
        end
    end

    bridges = [sAdd, tAdd];
end
function data = makeRectData(pg, varargin)
% MAKERECTDATA  Build grid and per-cell tables for rectangle MILP.

p = inputParser;
p.addParameter('Xgrid', [], @(v)isvector(v)&&isreal(v));
p.addParameter('Ygrid', [], @(v)isvector(v)&&isreal(v));
p.addParameter('GridStep', [], @(v)isempty(v)||(isscalar(v)&&v>0));
p.addParameter('BBox', [], @(v)isvector(v)&&numel(v)==4);
p.parse(varargin{:});
Xg=p.Results.Xgrid; Yg=p.Results.Ygrid; h=p.Results.GridStep; BBox=p.Results.BBox;

tol=1e-12;

% bbox
if isempty(BBox)
    try
        V = pg.Vertices;
        x1=min(V(:,1)); x2=max(V(:,1)); y1=min(V(:,2)); y2=max(V(:,2));
    catch
        [vx,vy]=boundary(pg); if iscell(vx), vx=vertcat(vx{:}); vy=vertcat(vy{:}); end
        m=~isnan(vx)&~isnan(vy); vx=vx(m); vy=vy(m);
        x1=min(vx); x2=max(vx); y1=min(vy); y2=max(vy);
    end
else
    x1=BBox(1); x2=BBox(2); y1=BBox(3); y2=BBox(4);
end

% grids
if isempty(Xg), if isempty(h), h=(x2-x1)/60; end, X=unique([x1; (x1:h:x2).'; x2]);
else, X=sort(Xg(:)); end
if isempty(Yg), if isempty(h), h=(y2-y1)/60; end, Y=unique([y1; (y1:h:y2).'; y2]);
else, Y=sort(Yg(:)); end
if numel(X)>1, X=X([true; diff(X)>tol]); end
if numel(Y)>1, Y=Y([true; diff(Y)>tol]); end

nX=numel(X); nY=numel(Y); nCx=nX-1; nCy=nY-1;
if nCx<=0||nCy<=0, error('Grid must have at least one cell.'); end

dx=diff(X); dy=diff(Y);
Aij = dx(:)*dy(:)';

% overlaps
Sij=zeros(nCx,nCy);
for i=1:nCx
    xa=X(i); xb=X(i+1);
    for j=1:nCy
        ya=Y(j); yb=Y(j+1);
        R=polyshape([xa xb xb xa],[ya ya yb yb]);
        Sij(i,j)=area(intersect(pg,R));
    end
end

data=struct('X',X,'Y',Y,'nX',nX,'nY',nY,'nCx',nCx,'nCy',nCy, ...
            'dx',dx(:),'dy',dy(:),'Aij',Aij,'Sij',Sij,'area_pg',area(pg));
end
function [rectx,recty,area,perimeter] = minboundrect(x,y,metric)
% minboundrect: Compute the minimal bounding rectangle of points in the plane
% usage: [rectx,recty,area,perimeter] = minboundrect(x,y,metric)
%
% arguments: (input)
%  x,y - vectors of points, describing points in the plane as
%        (x,y) pairs. x and y must be the same lengths.
%
%  metric - (OPTIONAL) - single letter character flag which
%        denotes the use of minimal area or perimeter as the
%        metric to be minimized. metric may be either 'a' or 'p',
%        capitalization is ignored. Any other contraction of 'area'
%        or 'perimeter' is also accepted.
%
%        DEFAULT: 'a'    ('area')
%
% arguments: (output)
%  rectx,recty - 5x1 vectors of points that define the minimal
%        bounding rectangle.
%
%  area - (scalar) area of the minimal rect itself.
%
%  perimeter - (scalar) perimeter of the minimal rect as found
%
%
% Note: For those individuals who would prefer the rect with minimum
% perimeter or area, careful testing convinces me that the minimum area
% rect was generally also the minimum perimeter rect on most problems
% (with one class of exceptions). This same testing appeared to verify my
% assumption that the minimum area rect must always contain at least
% one edge of the convex hull. The exception I refer to above is for
% problems when the convex hull is composed of only a few points,
% most likely exactly 3. Here one may see differences between the
% two metrics. My thanks to Roger Stafford for pointing out this
% class of counter-examples.
%
% Thanks are also due to Roger for pointing out a proof that the
% bounding rect must always contain an edge of the convex hull, in
% both the minimal perimeter and area cases.
%
%
% Example usage:
%  x = rand(50000,1);
%  y = rand(50000,1);
%  tic,[rx,ry,area] = minboundrect(x,y);toc
%
%  Elapsed time is 0.105754 seconds.
%
%  [rx,ry]
%  ans =
%      0.99994  -4.2515e-06
%      0.99998      0.99999
%   2.6441e-05            1
%  -5.1673e-06   2.7356e-05
%      0.99994  -4.2515e-06
%
%  area
%  area =
%      0.99994
%
%
% See also: minboundcircle, minboundtri, minboundsphere
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 3.0
% Release date: 3/7/07
% default for metric
if (nargin<3) || isempty(metric)
  metric = 'a';
elseif ~ischar(metric)
  error 'metric must be a character flag if it is supplied.'
else
  % check for 'a' or 'p'
  metric = lower(metric(:)');
  ind = strmatch(metric,{'area','perimeter'});
  if isempty(ind)
    error 'metric does not match either ''area'' or ''perimeter'''
  end
  
  % just keep the first letter.
  metric = metric(1);
end
% preprocess data
x=x(:);
y=y(:);
% not many error checks to worry about
n = length(x);
if n~=length(y)
  error 'x and y must be the same sizes'
end
% start out with the convex hull of the points to
% reduce the problem dramatically. Note that any
% points in the interior of the convex hull are
% never needed, so we drop them.
if n>3
  edges = convhull(x,y);
  % exclude those points inside the hull as not relevant
  % also sorts the points into their convex hull as a
  % closed polygon
  
  x = x(edges);
  y = y(edges);
  
  % probably fewer points now, unless the points are fully convex
  nedges = length(x) - 1;
elseif n>1
  % n must be 2 or 3
  nedges = n;
  x(end+1) = x(1);
  y(end+1) = y(1);
else
  % n must be 0 or 1
  nedges = n;
end
% now we must find the bounding rectangle of those
% that remain.
% special case small numbers of points. If we trip any
% of these cases, then we are done, so return.
switch nedges
  case 0
    % empty begets empty
    rectx = [];
    recty = [];
    area = [];
    perimeter = [];
    return
  case 1
    % with one point, the rect is simple.
    rectx = repmat(x,1,5);
    recty = repmat(y,1,5);
    area = 0;
    perimeter = 0;
    return
  case 2
    % only two points. also simple.
    rectx = x([1 2 2 1 1]);
    recty = y([1 2 2 1 1]);
    area = 0;
    perimeter = 2*sqrt(diff(x).^2 + diff(y).^2);
    return
end
% 3 or more points.
% will need a 2x2 rotation matrix through an angle theta
Rmat = @(theta) [cos(theta) sin(theta);-sin(theta) cos(theta)];
% get the angle of each edge of the hull polygon.
ind = 1:(length(x)-1);
edgeangles = atan2(y(ind+1) - y(ind),x(ind+1) - x(ind));
% move the angle into the first quadrant.
edgeangles = unique(mod(edgeangles,pi/2));
% now just check each edge of the hull
nang = length(edgeangles);
area = inf;
perimeter = inf;
met = inf;
xy = [x,y];
for i = 1:nang
  % rotate the data through -theta 
  rot = Rmat(-edgeangles(i));
  xyr = xy*rot;
  xymin = min(xyr,[],1);
  xymax = max(xyr,[],1);
  
  % The area is simple, as is the perimeter
  A_i = prod(xymax - xymin);
  P_i = 2*sum(xymax-xymin);
  
  if metric=='a'
    M_i = A_i;
  else
    M_i = P_i;
  end
  
  % new metric value for the current interval. Is it better?
  if M_i<met
    % keep this one
    met = M_i;
    area = A_i;
    perimeter = P_i;
    
    rect = [xymin;[xymax(1),xymin(2)];xymax;[xymin(1),xymax(2)];xymin];
    rect = rect*rot';
    rectx = rect(:,1);
    recty = rect(:,2);
  end
end
% get the final rect
% all done
end % mainline end
function [idx, kOpt, res, info] = minCoverSeg_confidence(X, dmax, covFrac, varargin)
% MINCOVERSEG_CONFIDENCE  Minimal subset size covering ≥ covFrac of points within dmax
%                         to at least one selected segment (endpoints chosen).
% Pure MILP, no Big-M.
%
% [idx, kOpt, res, info] = minCoverSeg_confidence(X, dmax, covFrac, 'Name',Value,...)
%
% Inputs
%   X        : n-by-d matrix of points
%   dmax     : scalar max orthogonal distance
%   covFrac  : required coverage fraction in [0,1]
%
% Name–Value options:
%   'Pairs'   : m-by-2 candidate segments [i j] (i<j). Default: k-NN graph.
%   'Knn'     : neighbors for k-NN graph if Pairs empty                  [8]
%   'Params'  : Gurobi params struct                                     [struct()]
%   'Verbose' : 0/1                                                      [1]
%
% Outputs
%   idx   : indices of selected points
%   kOpt  : |idx|
%   res   : Gurobi result struct
%   info  : struct (coverage masks, distances, active segs, etc.)
%
% -------------------------------------------------------------------------
p = inputParser;
addParameter(p,'Pairs',[],@(z)isnumeric(z)&&size(z,2)==2);
addParameter(p,'Knn',8,@(z)isscalar(z)&&z>=1);
addParameter(p,'Params',struct());
addParameter(p,'Verbose',1,@(z)isscalar(z));
parse(p,varargin{:});
P       = p.Results.Pairs;
Knn     = p.Results.Knn;
Params  = p.Results.Params;
verb    = p.Results.Verbose;

[n,~] = size(X);
if covFrac<0 || covFrac>1
    error('covFrac must be in [0,1]');
end
Umax = floor((1-covFrac)*n);

% ----- Build candidate segments -----
if isempty(P)
    if verb, fprintf('[MILP] Building %d-NN segment set...\\n',Knn); end
    if exist('knnsearch','file')
        [idxNN,~] = knnsearch(X,X,'K',Knn+1);
    else
        D = pdist2(X,X);
        [~,idxNN] = sort(D,2,'ascend');
        idxNN = idxNN(:,1:Knn+1);
    end
    Ptmp = [];
    for i=1:n
        nbr = idxNN(i,2:end);
        Ptmp = [Ptmp; [i*ones(Knn,1) nbr(:)]];
    end
    P = unique(sort(Ptmp,2),'rows');
end
m = size(P,1);
if verb, fprintf('[MILP] n=%d, segments=%d\\n',n,m); end

% ----- Distances and coverage -----
if verb, fprintf('[MILP] Computing distances & coverage...\\n'); end
Delta = zeros(n,m);
for pidx=1:m
    i = P(pidx,1); j = P(pidx,2);
    Delta(:,pidx) = point_to_segment_dist(X, X(i,:), X(j,:));
end
CoverMask = (Delta <= dmax);

% prune useless segments
usefulSeg = any(CoverMask,1);
P      = P(usefulSeg,:);
Delta  = Delta(:,usefulSeg);
CoverMask = CoverMask(:,usefulSeg);
m = size(P,1);
if verb, fprintf('[MILP] Useful segments after prune = %d\\n',m); end

forcedUncov = find(~any(CoverMask,2));
if numel(forcedUncov)>Umax
    warning('Forced uncovered (%d) > allowed (%d). Might be infeasible.',numel(forcedUncov),Umax);
end

% ----- Variables -----
ix = 1:n;                 % x_i
iy = n + (1:m);           % y_p
iu = n + m + (1:n);       % u_r
nv = n + m + n;

obj = zeros(nv,1); obj(ix)=1;
modelsense = 'min';

rows = {}; rhs=[]; sense=[];

% (1) segment activation link
for pidx=1:m
    i = P(pidx,1); j = P(pidx,2);
    rows{end+1} = sparse(1,[iy(pidx) ix(i)],[1 -1],1,nv); rhs(end+1,1)=0; sense(end+1,1)='<';
    rows{end+1} = sparse(1,[iy(pidx) ix(j)],[1 -1],1,nv); rhs(end+1,1)=0; sense(end+1,1)='<';
    rows{end+1} = sparse(1,[iy(pidx) ix(i) ix(j)],[1 -1 -1],1,nv); rhs(end+1,1)=-1; sense(end+1,1)='>';
end

% (2) coverage with slack
for r=1:n
    cols = find(CoverMask(r,:));
    rowIdx = [iy(cols) iu(r)];
    vals   = [ones(1,numel(cols)) 1];
    rows{end+1} = sparse(1,rowIdx,vals,1,nv);
    rhs(end+1,1)=1; sense(end+1,1)='>';
end

% (3) total uncovered <= Umax
rows{end+1} = sparse(1,iu,1,1,nv); rhs(end+1,1)=Umax; sense(end+1,1)='<';

model.A     = cat(1,rows{:});
model.rhs   = double(rhs);
model.sense = char(sense);

vtype = [repmat('B',n,1); repmat('B',m,1); repmat('B',n,1)];
lb = zeros(nv,1); ub = ones(nv,1);

model.obj        = obj;
model.modelsense = modelsense;
model.lb         = lb;
model.ub         = ub;
model.vtype      = vtype(:)';

if ~isfield(Params,'OutputFlag'), Params.OutputFlag = verb; end

res = gurobi(model, Params);

if ~isfield(res,'x')
    error('No solution. Status=%s',fieldOr(res,'status','UNKNOWN'));
end

if verb
    fprintf('Status: %s | Obj=%g | Gap=%.3g | Time=%.1fs\\n',...
        fieldOr(res,'status','?'),fieldOr(res,'objval',NaN),...
        fieldOr(res,'mipgap',NaN),fieldOr(res,'runtime',NaN));
end

xsol = res.x(ix)>0.5;
ysol = res.x(iy)>0.5;
usol = res.x(iu)>0.5;

idx   = find(xsol);
kOpt  = numel(idx);
covered = ~usol;
covRate = mean(covered);

% robust min distance
minDist = inf(n,1);
activeSegMask = ysol(:)';           % 1-by-m logical
for r=1:n
    rowMask = CoverMask(r,:) & activeSegMask;
    if any(rowMask)
        minDist(r) = min( Delta(r,rowMask) );
    end
end

info.P          = P;
info.Delta      = Delta;
info.CoverMask  = CoverMask;
info.covered    = covered;
info.uncovered  = find(~covered);
info.yActive    = find(ysol);
info.Umax       = Umax;
info.covRate    = covRate;
info.minDist    = minDist;
info.runtime    = fieldOr(res,'runtime',NaN);
info.gap        = fieldOr(res,'mipgap',NaN);
end

% ---- helpers ----
function d = point_to_segment_dist(P, A, B)
    AP = bsxfun(@minus,P,A);
    AB = B - A;
    denom = AB*AB';
    if denom==0
        d = sqrt(sum((P-A).^2,2)); return
    end
    t = (AP*AB')/denom;
    t = max(0,min(1,t));
    proj = A + t.*AB;
    d = sqrt(sum((P-proj).^2,2));
end

function v = fieldOr(s,f,def)
    if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = def; end
end

function [res1, res2, selectedRects, assignMap, stats] = minK_maxCoverage_lexi(F, xm, xM, ym, yM, gurobiParams)
% minK_maxCoverage_lexi
% LEXICOGRAPHIC SOLVE:
%   Phase 1: maximize coverage (minimize #uncovered u_p)
%   Phase 2: among all max-coverage solutions, minimize #rectangles
%
% Model (both phases share constraints):
%   Vars: z_r in {0,1}   (choose rectangle r)
%         u_p in {0,1}   (foreground pixel p is UNcovered)
%   (A) Foreground assignment (enforces non-overlap):
%       sum_{r: p in r} z_r + u_p = 1     for all p in F
%   (B) Explicit fullness:
%       sum_{r: q in r} z_r = 0           for all q in B
% Phase 1 OBJ: minimize sum_p u_p            (== maximize coverage)
% Phase 2 OBJ: minimize sum_r z_r  subject to sum_p u_p <= U*
%
% Inputs:
%   F  : HxW logical/double foreground mask (1=foreground)
%   xm,xM, ym,yM : integer size bounds on width/height
%   gurobiParams (optional) : struct of Gurobi params
%
% Outputs:
%   res1        : Gurobi result for Phase 1 (max coverage)
%   res2        : Gurobi result for Phase 2 (min rectangles at max coverage)
%   selectedRects : table of rectangles chosen in Phase 2
%   assignMap     : HxW uint8 map over FOREGROUND: 1=covered, 0=uncovered (Phase 2)
%   stats       : struct with fields H,W,nF,nB,totalR, Ustar, coverFracStar, Kstar

if nargin < 5, error('Usage: F, xm, xM, ym, yM , [gurobiParams]'); end
if nargin < 6, gurobiParams = struct(); end

F = logical(F);
[H,W] = size(F);
if xm < 1 || ym < 1 || xM < xm || yM < ym, error('Invalid size bounds.'); end

idxF = find(F(:));     nF = numel(idxF);
idxB = find(~F(:));    nB = numel(idxB);
nU  = H*W;

% ---------- Enumerate ALL size-bounded rectangles and build membership ----------
totalR = 0; totNZ = 0;
for h = ym:yM
    for w = xm:xM
        nPos = max(H - h + 1, 0) * max(W - w + 1, 0);
        totalR = totalR + nPos;
        totNZ  = totNZ + nPos * (h*w);
    end
end
if totalR == 0, error('No candidate rectangles under given bounds.'); end

topRow = zeros(totalR,1,'uint32');
leftCol= zeros(totalR,1,'uint32');
height = zeros(totalR,1,'uint16');
width  = zeros(totalR,1,'uint16');
area   = zeros(totalR,1,'uint32');

I = zeros(totNZ,1,'uint32');   % pixel indices
J = zeros(totNZ,1,'uint32');   % rectangle indices
V = ones (totNZ,1);
ridx = 0; t = 0;

for h = ym:yM
  for w = xm:xM
    max_i = H - h + 1; max_j = W - w + 1;
    if max_i <= 0 || max_j <= 0, continue; end
    for i = 1:max_i
      for j = 1:max_j
        ridx = ridx + 1;
        topRow(ridx)  = i; leftCol(ridx) = j;
        height(ridx)  = h; width(ridx)   = w;
        area(ridx)    = uint32(h*w);

        rows = i:(i+h-1); cols = j:(j+w-1);
        [RR,CC] = ndgrid(rows, cols);
        pidx = sub2ind([H,W], RR(:), CC(:));

        L = numel(pidx);
        I(t+(1:L)) = uint32(pidx);
        J(t+(1:L)) = uint32(ridx);
        t = t + L;
      end
    end
  end
end
if t ~= totNZ, error('NNZ mismatch.'); end

M = sparse(double(I), double(J), double(V), nU, double(totalR));
M_F = M(idxF, :);   % foreground rows
M_B = M(idxB, :);   % background rows

% ---------- Base model pieces (shared) ----------
% Vars: [ z(1:totalR); u(1:nF) ]
nVars = totalR + nF;

% (A) Foreground assignment: M_F * z + u = 1
A1 = [ M_F, speye(nF) ];
rhs1 = ones(nF,1);
sense1 = repmat('=', nF, 1);

% (B) Explicit fullness on background: M_B * z = 0
A2 = [ M_B, sparse(nB, nF) ];
rhs2 = zeros(nB,1);
sense2 = repmat('=', nB, 1);

% Common bounds/types
lb   = zeros(nVars,1);
ub   = ones (nVars,1);
vtype= [repmat('B', totalR, 1); repmat('B', nF, 1)];

% ---------- Phase 1: maximize coverage -> minimize sum u ----------
model1.A     = [A1; A2];
model1.rhs   = [rhs1; rhs2];
model1.sense = [sense1; sense2];
model1.modelsense = 'min';
model1.obj   = [zeros(totalR,1); ones(nF,1)];  % min sum u
model1.lb    = lb; model1.ub = ub; model1.vtype = vtype;

params1.OutputFlag = 1; params1.MIPGap = 1e-4;
fn = fieldnames(gurobiParams);
for k=1:numel(fn), params1.(fn{k}) = gurobiParams.(fn{k}); end

res1 = gurobi(model1, params1);
if ~isfield(res1,'x'), error('Phase 1 failed to return a solution.'); end
Ustar = round(res1.objval);          % min #uncovered -> maximal coverage
% Ustar = floor(res1.objval + 1e-6);   % min # uncovered from Phase 1 (max coverage)
coverFracStar = (nF - Ustar) / max(1,nF);

% % ---------- Phase 2: minimize #rectangles subject to max coverage ----------
% Add constraint sum u <= Ustar
Acov = sparse(1, nVars); Acov(1, totalR+(1:nF)) = 1;
model2.A     = [A1; A2; Acov];
model2.rhs   = [rhs1; rhs2; Ustar];
model2.sense = [sense1; sense2; '<'];

% allow slack of, say, 1% uncovered beyond Ustar:
% slackFrac = 0.00001;                     % <-- TUNE THIS (e.g., 0.005 .. 0.02)
% Delta     = floor(slackFrac * nF);
% Ucap      = Ustar + Delta;
% 
% Acov = sparse(1, nVars); Acov(1, totalR+(1:nF)) = 1;   % sum u
% model2.A     = [A1; A2; Acov];
% model2.rhs   = [rhs1; rhs2; Ucap];
% model2.sense = [sense1; sense2; '<'];

% Objective: minimize sum z
model2.modelsense = 'min';
model2.obj  = [ones(totalR,1); zeros(nF,1)];
model2.lb   = lb; model2.ub = ub; model2.vtype = vtype;

params2 = params1;
res2 = gurobi(model2, params2);
if ~isfield(res2,'x'), error('Phase 2 failed to return a solution.'); end

% ---------- Extract Phase 2 solution ----------
sol = res2.x; z = sol(1:totalR); u = sol(totalR+1:end);
assignMap = zeros(H,W,'uint8'); assignMap(idxF) = uint8(1 - (u>0.5));

sel = z > 0.5;
selectedRects = table( ...
    double(topRow(sel)), double(leftCol(sel)), ...
    double(height(sel)), double(width(sel)), ...
    double(area(sel)), double(z(sel)), ...
    'VariableNames', {'topRow','leftCol','height','width','area','zVal'});

Kstar = sum(sel);
fprintf('Lexicographic result: max coverage = %.3f (U*=%d of %d), min K = %d\n', coverFracStar, Ustar, nF, Kstar);

% ---------- Stats ----------
stats = struct('H',H, 'W',W, 'nF',nF, 'nB',nB, 'totalR',totalR, ...
               'Ustar',Ustar, 'coverFracStar',coverFracStar, 'Kstar',Kstar);
end
function [idx,root] = mstPreorderOrder(P)
    [rx,ry]=minboundrect(P(:,1),P(:,2));
    [R, c] = alignToNearestAxis([rx,ry]);
    xx=[rx,ry]*R'+c;
    root=knnsearch(P,(min(xx)-c)*R);
    D = squareform(pdist(P));                     % pairwise distances
    G = graph(D, 'upper');                        % complete (small n)
    T = minspantree(G,'Method','dense','Root', root);          % MST
    % DFS preorder from a reasonable root (e.g., farthest from centroid)
    order = dfsearch(T, root);       % nodes in visitation order
    idx = order(:).';                             % make row vector
end
% XBIM XS alfa cen fx fy
function [Att,RMSE,ID,r]=RobustConcentrationPlaneFit(XS,Miter)
cut=sqrt(chi2inv(0.99,3));
id=0;
RMSE=Inf;
iter=1;
N=length(XS);
w=ones(N,1);
X=XS;
while id==0 && iter<=Miter
    [Att,RMSE1,r]=PlaneFitW(X,w);
    if abs(RMSE-RMSE1)<=0.0005 
        id=1;
    else
        md=modect(r);
        if isempty(md)
            md=median(r);
        end
        % md=median(r);
        rj=r-md;
        % s=max(0.001,min(MADN(r),std(r)));
        s=MADN(r);
        u=rj/s/cut;
        w=(1-u.^2).^0.5; %Bisquare
        w(abs(u)>=1)=0;
        iter=iter+1;
        RMSE=RMSE1;
    end
end
r= [XS,ones(N,1)]*Att';
u=(r-md)/s/cut;
ID=find(abs(u)<1);
[Att,RMSE,r]=PlaneFitW(XS(ID,:),ones(length(ID),1));
end
function [Attribute,RMSE,RMS]=PlaneFitW(X,w)
    [N,~]=size(X);
    C=sum(w.*X)/sum(w);
    cv=N/(N-1)*(w.*(X-C))'*(X-C)/sum(w);
    [c,b]=eig(cv);
    [~,r]=min(diag(b));
    aa=c(:,r);
    [~,a]=max(abs(aa));
    if aa(a)<0
        aa=-aa;
    end
    Attribute=[aa',-dot(C,aa')];  
    RMS= [X,ones(N,1)]*Attribute';
    RMSE=sqrt(mean(RMS.^2));
end
%% The Mode detection using Matlab findpeak function
% Only for registered members of the "Digital Technologies for Field
% Information Modeling" Course at KIT. Use of the software for commercial
% purposes is strictly prohibited. Distribution only allowed with written
% consent from Jun.-Prof. Dr. Reza Maalek (reza.maalek@kit.edu). For
% educational use only.
% Inputs: RMS (univariate data)
% Outputs: Mode (Robust mean)
function Mod=modect(RMS)
    [f,xi]=ksdensity(RMS);
    [L,loc]=findpeaks(f);
    [~,r]=max(L);
    Mod=xi(loc(r));
end

%%MADN of a univariate dataset
function [MADN,MAD]=MADN(D)
med=abs(D-repmat(median(D),length(D),1));
MAD=median(med);
MADN=MAD/0.67449;
end
function [parts, meta] = pointBoundariesImage(x, y, varargin)
% pointBoundariesImage
% Extract ordered outer + inner boundaries from 2-D points by rasterizing to an image,
% and map those boundaries back to the ACTUAL original points (ordered along each loop).
%
% Usage:
%   parts = pointBoundariesImage(x, y)
%   parts = pointBoundariesImage([x y], 'Mode','distance', 'Radius',0.05, 'Plot',true)
%   [parts, meta] = pointBoundariesImage(x, y, 'Mode','density', 'Sigma',0.02, 'Thresh',0.2)
%
% Key options:
%   'Mode'            - 'distance' (union of disks) | 'density' (Gaussian splat). Default 'distance'.
%   'NumPixels'       - long-side image size (px). Default 800. (Overrides PixelsPerUnit if set)
%   'PixelsPerUnit'   - pixels per world unit (used when NumPixels is []).
%   'Padding'         - bbox padding (world units). Default 1% of long span.
%   'Radius'          - (distance) disk radius (world units). Default 1% of long span.
%   'Sigma'           - (density) Gaussian sigma (world units). Default 0.5% of long span.
%   'Thresh'          - (density) threshold in [0,1] after max-normalization. Default 0.15.
%   'MorphOpen'       - opening radius (px). Default 0.
%   'MorphClose'      - closing radius (px). Default 2.
%   'MinRegionArea'   - drop tiny connected components before tracing (world-units²; approx via px). Default 0.
%   'MinBoundaryArea' - drop tiny loops (outer/holes) AFTER tracing (world-units²). Default 0.   % NEW
%   'BoundarySnapTol' - snap tolerance (world units) to map original points to each boundary.     % NEW
%                       Default: ~1.5 * average world-pixel size of the raster.
%   'KeepHoles'       - true keeps inner holes (default true). If false, holes are filled in output.
%   'Plot'            - quick visualization (default false).
%
% Output (per connected component i):
%   parts(i).outerPoly          - [N x 2] CCW polygon (closed) in world coordinates
%   parts(i).holePolys{h}       - [M x 2] CW polygons (closed)
%   parts(i).allLoops           - {1xK} all loops (closed) for the component
%   parts(i).areas              - signed areas for allLoops (CCW>0, CW<0; world units²)
%
%   % Mapped "actual" original points ON the boundaries (ordered along the loops):
%   parts(i).outerIdx           - indices into original (x,y) that lie on the outer boundary (ordered)
%   parts(i).outerPts           - coordinates of those points (ordered)
%   parts(i).holeIdx{h}         - indices on hole h (ordered)
%   parts(i).holePts{h}         - coordinates for hole h (ordered)
%
% Meta:
%   meta.mask   - logical raster used for boundary extraction
%   meta.R      - imref2d spatial referencing object
%   meta.params - struct of resolved parameters (for reproducibility)

%% -------- Flexible inputs (Nx2 allowed) --------
if nargin >= 1 && (isempty(y) || (isvector(y) && numel(y) ~= numel(x))) && size(x,2) == 2
    y = x(:,2); x = x(:,1);
end
x = x(:); y = y(:);
assert(numel(x) == numel(y), 'x and y must have the same length.');

%% -------- Parse options --------
p = inputParser;
p.addParameter('Mode','distance', @(s)ischar(s)&&ismember(lower(s),{'distance','density'}));
p.addParameter('NumPixels', 800, @(v) isempty(v) || (isscalar(v) && v>=64));
p.addParameter('PixelsPerUnit', [], @(v) isempty(v) || (isscalar(v) && v>0));
p.addParameter('Padding', [], @(v) isempty(v) || (isscalar(v) && v>=0));
p.addParameter('Radius', [], @(v) isempty(v) || (isscalar(v) && v>=0));
p.addParameter('Sigma', [], @(v) isempty(v) || (isscalar(v) && v>=0));
p.addParameter('Thresh', 0.15, @(v) isscalar(v) && v>=0 && v<=1);
p.addParameter('MorphOpen', 0, @(v) isscalar(v) && v>=0);
p.addParameter('MorphClose', 2, @(v) isscalar(v) && v>=0);
p.addParameter('MinRegionArea', 0, @(v) isscalar(v) && v>=0);
p.addParameter('MinBoundaryArea', 0, @(v) isscalar(v) && v>=0);       % NEW
p.addParameter('BoundarySnapTol', [], @(v) isempty(v) || (isscalar(v) && v>0)); % NEW
p.addParameter('KeepHoles', true, @(v) islogical(v) || ismember(v,[0 1]));
p.addParameter('Plot', false, @(v) islogical(v) || ismember(v,[0 1]));
p.parse(varargin{:});
opt = p.Results;

%% -------- Early outs / cleaning --------
maskFinite = isfinite(x) & isfinite(y);
x = x(maskFinite); y = y(maskFinite);
if numel(x) < 3
    parts = struct('outerPoly', [], 'holePolys', {{}}, 'allLoops', {{}}, 'areas', [], ...
                   'outerIdx', [], 'outerPts', [], 'holeIdx', {{}}, 'holePts', {{}});
    meta  = struct('mask', false(2,2), 'R', imref2d([2 2]), 'params', struct());
    warning('Not enough points to form boundaries.');
    return;
end

%% -------- World bbox and padding --------
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
spanX = max(xmax - xmin, eps);
spanY = max(ymax - ymin, eps);
spanL = max(spanX, spanY);
pad = defaultIfEmpty(opt.Padding, 0.01 * spanL);
xmin = xmin - pad; xmax = xmax + pad;
ymin = ymin - pad; ymax = ymax + pad;

%% -------- Image geometry (resolution) --------
if ~isempty(opt.NumPixels)
    % Target long side, preserve aspect
    aspect = spanY / spanX;
    if aspect >= 1
        H = max(64, round(opt.NumPixels));
        W = max(64, round(H / max(aspect, eps)));
    else
        W = max(64, round(opt.NumPixels));
        H = max(64, round(W * max(aspect, eps)));
    end
    ppuX = (W-1) / max(xmax - xmin, eps);
    ppuY = (H-1) / max(ymax - ymin, eps);
else
    % Use PixelsPerUnit
    ppuX = opt.PixelsPerUnit;  ppuY = opt.PixelsPerUnit;
    W = max(64, 1 + round(ppuX * (xmax - xmin)));
    H = max(64, 1 + round(ppuY * (ymax - ymin)));
end
R = imref2d([H W], [xmin xmax], [ymin ymax]);   % spatial ref
dx = (xmax - xmin) / max(W-1,1);                % world units per pixel (x)
dy = (ymax - ymin) / max(H-1,1);                % world units per pixel (y)
pixWorld = mean([dx, dy]);                      % average world-pixel size

%% -------- Rasterize points (nearest pixel hits) --------
[u, v] = worldToIntrinsic(R, x, y);   % u ~ column (x), v ~ row (y)
col = max(1, min(W, round(u)));
row = max(1, min(H, round(v)));
I = false(H, W);
if ~isempty(row), I(sub2ind([H W], row, col)) = true; end

%% -------- Resolve mode parameters (always init) --------
rpx   = 0;    % distance-mode radius in pixels
sigpx = 0;    % density-mode sigma in pixels
pxScale = mean([ppuX, ppuY]);  % world->px isotropic approx

if strcmpi(opt.Mode,'distance')
    radiusWorld = defaultIfEmpty(opt.Radius, 0.01*spanL);
    rpx = max(0, radiusWorld * pxScale);
else % 'density'
    sigmaWorld = defaultIfEmpty(opt.Sigma, 0.005*spanL);
    sigpx = max(0.1, sigmaWorld * pxScale);
end

%% -------- Build binary mask --------
switch lower(opt.Mode)
    case 'distance'
        % Union of disks of radius rpx around each point
        if nnz(I) == 0
            mask = false(H, W);
        else
            D = bwdist(~I);      % distance (px) to nearest hit
            mask = D <= rpx;     % inside any disk (correct inequality)
        end
    case 'density'
        if nnz(I) == 0
            mask = false(H, W);
        else
            G = gaussianBlur(I, sigpx);     % smooth "splat"
            if max(G(:)) > 0, G = G / max(G(:)); end
            mask = G >= opt.Thresh;
        end
    otherwise
        error('Unknown Mode "%s".', opt.Mode);
end

%% -------- Morphological cleanup (optional) --------
if opt.MorphOpen > 0
    mask = imopen(mask, strel('disk', round(opt.MorphOpen)));
end
if opt.MorphClose > 0
    mask = imclose(mask, strel('disk', round(opt.MorphClose)));
end

%% -------- Remove tiny connected components by world area (approx via pixels) --------
if opt.MinRegionArea > 0 && any(mask(:))
    CC0 = bwconncomp(mask, 8);
    pxArea = dx * dy;                       % world area per pixel
    keep = false(1, CC0.NumObjects);
    for k = 1:CC0.NumObjects
        keep(k) = (numel(CC0.PixelIdxList{k}) * pxArea) >= opt.MinRegionArea;
    end
    mask2 = false(H, W);
    for k = find(keep)
        mask2(CC0.PixelIdxList{k}) = true;
    end
    mask = mask2;
end

%% -------- Extract boundaries per connected component --------
CC = bwconncomp(mask, 8);
parts = repmat(struct('outerPoly', [], 'holePolys', {{}}, 'allLoops', {{}}, 'areas', [], ...
                      'outerIdx', [], 'outerPts', [], 'holeIdx', {{}}, 'holePts', {{}}), ...
               max(CC.NumObjects,1), 1);

if CC.NumObjects == 0
    meta = struct('mask', mask, 'R', R, 'params', paramsStruct());
    return;
end

% Snap tolerance in world units for mapping original points to boundaries
snapTol = defaultIfEmpty(opt.BoundarySnapTol, 1.5 * pixWorld);

for ci = 1:CC.NumObjects
    compMask = false(H, W);
    compMask(CC.PixelIdxList{ci}) = true;

    % Boundaries (outer + holes) for THIS component
    B = bwboundaries(compMask, 8, 'holes');

    % Convert boundaries to world coords and compute signed areas
    compLoops = cell(numel(B),1);
    compAreas = zeros(numel(B),1);
    for bi = 1:numel(B)
        rc = B{bi};                       % [row col]
        if ~isequal(rc(1,:), rc(end,:)), rc = [rc; rc(1,:)]; end
        [xw, yw] = intrinsicToWorld(R, rc(:,2), rc(:,1));  % col->x, row->y
        P = [xw(:), yw(:)];
        compLoops{bi} = P;
        compAreas(bi) = signedArea(P);    % +CCW / -CW
    end

    % Drop tiny loops by world area, if requested
    if opt.MinBoundaryArea > 0 && ~isempty(compLoops)
        keep = abs(compAreas) >= opt.MinBoundaryArea;
        compLoops = compLoops(keep);
        compAreas = compAreas(keep);
    end

    if isempty(compLoops)
        % Entire component suppressed by MinBoundaryArea
        continue;
    end

    % Pick OUTER = loop with max |area|; orient CCW; holes CW
    [~, io] = max(abs(compAreas));
    outer = compLoops{io};
    if signedArea(outer) < 0, outer = flipud(outer); end
    holes = compLoops; holes(io) = [];
    for h = 1:numel(holes)
        if signedArea(holes{h}) > 0, holes{h} = flipud(holes{h}); end
    end
    if ~opt.KeepHoles, holes = {}; end

    % Save polygons (closed)
    parts(ci).outerPoly = ensureClosed(outer);
    parts(ci).holePolys = cellfun(@ensureClosed, holes, 'UniformOutput', false);
    parts(ci).allLoops  = cellfun(@ensureClosed, compLoops, 'UniformOutput', false);
    parts(ci).areas     = compAreas(:);

    % ===== Map boundaries to ACTUAL original points (ordered) =====
    % Outer:
    if ~isempty(parts(ci).outerPoly)
        [idxO, ptsO] = mapPointsToPolyline([x y], parts(ci).outerPoly, snapTol);
        parts(ci).outerIdx = idxO;
        parts(ci).outerPts = ptsO;
    end
    % Holes:
    parts(ci).holeIdx = cell(size(parts(ci).holePolys));
    parts(ci).holePts = cell(size(parts(ci).holePolys));
    for h = 1:numel(parts(ci).holePolys)
        [idxH, ptsH] = mapPointsToPolyline([x y], parts(ci).holePolys{h}, snapTol);
        parts(ci).holeIdx{h} = idxH;
        parts(ci).holePts{h} = ptsH;
    end
end

%% -------- Optional plot --------
if opt.Plot
    figure('Name','pointBoundariesImage');
    subplot(1,2,1);
    plot(x, y, '.'); axis equal; grid on; title('Points');
    subplot(1,2,2);
    imshow(mask, R, 'InitialMagnification', 'fit'); hold on;
    for i = 1:numel(parts)
        if ~isempty(parts(i).outerPoly)
            plot(parts(i).outerPoly(:,1), parts(i).outerPoly(:,2), '-', 'LineWidth', 2);
            % Also draw mapped actual points:
            if ~isempty(parts(i).outerPts)
                plot(parts(i).outerPts(:,1), parts(i).outerPts(:,2), 'o', 'MarkerSize', 4);
            end
        end
        for h = 1:numel(parts(i).holePolys)
            plot(parts(i).holePolys{h}(:,1), parts(i).holePolys{h}(:,2), '--', 'LineWidth', 1.5);
            if ~isempty(parts(i).holePts{h})
                plot(parts(i).holePts{h}(:,1), parts(i).holePts{h}(:,2), 'o', 'MarkerSize', 4);
            end
        end
    end
    title('Mask + Boundaries + mapped points'); hold off;
end

%% -------- Meta --------
meta = struct('mask', mask, 'R', R, 'params', paramsStruct());

%% ======== Helpers ========

function v = defaultIfEmpty(val, def), if isempty(val), v = def; else, v = val; end, end

function P = ensureClosed(P)
    if isempty(P), return; end
    if ~isequal(P(1,:), P(end,:)), P = [P; P(1,:)]; end
end

function A = signedArea(P)
    if ~isequal(P(1,:), P(end,:)), P = [P; P(1,:)]; end
    xP = P(:,1); yP = P(:,2);
    A  = 0.5 * sum( xP(1:end-1).*yP(2:end) - xP(2:end).*yP(1:end-1) );
end

function G = gaussianBlur(B, s)
    if s <= 0
        G = double(B);
        return;
    end
    if exist('imgaussfilt','file')
        G = imgaussfilt(double(B), s);
    else
        w = max(3, 2*ceil(3*s)+1);
        h = fspecial('gaussian', w, s);
        G = imfilter(double(B), h, 'replicate', 'same');
    end
end

function S = paramsStruct()
    S = struct( ...
        'Mode', opt.Mode, ...
        'NumPixels', opt.NumPixels, ...
        'PixelsPerUnit', opt.PixelsPerUnit, ...
        'Padding', pad, ...
        'RadiusPX', rpx, ...
        'SigmaPX', sigpx, ...
        'Thresh', opt.Thresh, ...
        'MorphOpen', opt.MorphOpen, ...
        'MorphClose', opt.MorphClose, ...
        'MinRegionArea', opt.MinRegionArea, ...
        'MinBoundaryArea', opt.MinBoundaryArea, ...
        'BoundarySnapTol', snapTol );
end

function [idxOrdered, ptsOrdered] = mapPointsToPolyline(XY, polyClosed, tol)
% Map original points XY (Nx2) to a closed polyline polyClosed (Mx2).
% Returns indices and coordinates of the points whose shortest distance to the
% polyline is <= tol, ordered by their projection arclength along the polyline.
    if isempty(XY) || isempty(polyClosed)
        idxOrdered = []; ptsOrdered = []; return;
    end
    % Ensure closed
    if ~isequal(polyClosed(1,:), polyClosed(end,:))
        polyClosed = [polyClosed; polyClosed(1,:)];
    end
    % Precompute segment vectors and lengths
    P0 = polyClosed(1:end-1, :);
    P1 = polyClosed(2:end, :);
    seg = P1 - P0;                          % S x 2
    segLen = hypot(seg(:,1), seg(:,2));     % S x 1
    S = size(P0,1);
    cs = [0; cumsum(segLen)];               % cumulative arclength at segment Starts (length S+1)

    % For each point, compute distance to each segment via projection and pick minimum.
    N = size(XY,1);
    bestDist = inf(N,1);
    bestS    = zeros(N,1);                  % arclength position of nearest projection
    % Vectorized by blocks to control memory for large N,S
    blk = max(1, floor(1e7 / max(S,1)));    % heuristic block size
    for st = 1:blk:N
        en = min(N, st+blk-1);
        Q  = XY(st:en, :);                  % B x 2
        % Expand to B x S x 2 via bsxfun-style ops
        % Compute projection parameter t for each segment: t = dot(Q-P0, seg)/|seg|^2 clamped [0,1]
        % (Compute via matrix multiplications block-wise)
        QmP0x = Q(:,1) - P0(:,1).';   QmP0y = Q(:,2) - P0(:,2).';
        segx  = seg(:,1).';           segy  = seg(:,2).';
        seg2  = max(segx.^2 + segy.^2, eps);

        t = (QmP0x .* segx + QmP0y .* segy) ./ seg2;   % B x S
        t = max(0, min(1, t));                         % clamp

        % Projection points and squared distances
        projx = P0(:,1).' + t .* segx;                 % B x S
        projy = P0(:,2).' + t .* segy;                 % B x S
        dx = Q(:,1) - projx;  dy = Q(:,2) - projy;
        d2 = dx.^2 + dy.^2;                             % B x S

        % For each point (row), pick best segment
        [minD2, idxSeg] = min(d2, [], 2);
        d = sqrt(minD2);
        % Arclength parameter at projection:
        tbest = t(sub2ind(size(t), (1:size(t,1))', idxSeg));
        sHere = cs(idxSeg) + tbest .* segLen(idxSeg);

        % Update best
        upd = d < bestDist(st:en);
        bestDist(st:en) = d .* upd + bestDist(st:en) .* ~upd;
        bestS(st:en)    = sHere .* upd + bestS(st:en) .* ~upd;
    end

    % Select points within tolerance and order them by bestS
    in  = bestDist <= tol;
    idx = find(in);
    [~, order] = sort(bestS(in), 'ascend');
    idxOrdered = idx(order);
    ptsOrdered = XY(idxOrdered, :);
end

end
