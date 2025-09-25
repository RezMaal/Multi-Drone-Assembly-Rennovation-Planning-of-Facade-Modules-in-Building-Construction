function out = parse_solution_viz_hor(out1, G, R, sigma, m, E, opt)
%PARSE_SOLUTION_VIZ  Extract routes, compute stats, and visualize in tiles.
% Supports optional XY coordinates, per-route Colors, and custom arrowheads.
%
% out = parse_solution_viz(result, G, R, m, 'Name', value, ...)
%
% Options:
%   'ShowFigure'        : logical (default true)
%   'FigTitle'          : char, figure title
%   'Layout'            : 'auto' or [rows cols] (default 'auto')
%   'MarkerSizeS'       : scalar (default 90)
%   'MarkerSizeD'       : scalar (default 70)
%   'ArrowWidth'        : scalar (default 1.5)
%   'FontSize'          : scalar (default 10)
%   'SXY'               : k-by-2 source coordinates
%   'DXY'               : n-by-2 destination coordinates
%   'Pad'               : axis padding fraction (default 0.08)
%   'ShowNodeIDs'       : logical (default true)
%   'ShowStartEndText'  : logical (default true)
%   'Colors'            : m-by-3 RGB per route (default lines(m))
%   'ArrowHeadAt'       : 'end' (default) or 'middle'
%   'ArrowHeadLenFrac'  : head length as fraction of diagonal span (default 0.04)
%   'ArrowHeadWidFrac'  : head width  as fraction of diagonal span (default 0.03)
%
% OUTPUT fields:
%   ordersD{r}, startSource(r), bridgeSrc{r}, endSource(r),
%   cost_expand(r), T (table), stats (struct), figure, layout, SXY, DXY, Colors
%
% This function assumes the variable ordering used by the provided solvers:
% [ y0 | yjh | yt | v | u | cost | ... ]. If your model uses a different
% ordering, adapt the indices below.
%
% Copyright 2025.

% defaults
opts.ShowFigure       = true;
opts.FigTitle         = '';
opts.Layout           = 'auto';
opts.MarkerSizeS      = 90;
opts.MarkerSizeD      = 50;
opts.ArrowWidth       = 1.5;
opts.FontSize         = 10;
opts.SXY              = [];
opts.DXY              = [];
opts.Pad              = 0.05;
opts.ShowNodeIDs      = true;
opts.ShowStartEndText = true;
opts.Colors           = [];
opts.ArrowHeadAt      = 'end';
opts.ArrowHeadLenFrac = 0.005;
opts.ArrowHeadWidFrac = 0.005;
opts.PG = 0;
opts.sg = 30;

% parse
% if rem(numel(varargin),2)~=0, error('Name-value arguments must be pairs.'); end
% for i=1:2:numel(varargin)
%     name = varargin{i}; val = varargin{i+1};
%     if isfield(opts,name), opts.(name)=val; else, error('Unknown option %s',name); end
% end
if nargin<6
    opt = struct();
end
opts = validateAndFill(opt, opts);
[k,n] = size(G);
if size(R,1)~=n || size(R,2)~=k, error('R must be n-by-k.'); end

% coordinates
useCoords = (~isempty(opts.SXY) && ~isempty(opts.DXY) && size(opts.SXY,2)==2 && size(opts.DXY,2)==2 ...
             && size(opts.SXY,1)==k && size(opts.DXY,1)==n);
if useCoords
    SXY = opts.SXY; DXY = opts.DXY;
else
    xS=0; xD=1; yS=linspace(1,0,k).'; yD=linspace(1,0,n).';
    SXY=[repmat(xS,k,1) yS]; DXY=[repmat(xD,n,1) yD];
end
xmin=min([SXY(:,1);DXY(:,1)]); xmax=max([SXY(:,1);DXY(:,1)]);
ymin=min([SXY(:,2);DXY(:,2)]); ymax=max([SXY(:,2);DXY(:,2)]);
dx=xmax-xmin; dy=ymax-ymin; if dx==0,dx=1;end, if dy==0,dy=1;end
pad=opts.Pad; xlim_all=[xmin-pad*dx, xmax+pad*dx]; ylim_all=[ymin-pad*dy, ymax+pad*dy];
diagSpan = hypot(diff(xlim_all), diff(ylim_all));
headLen  = opts.ArrowHeadLenFrac * diagSpan;
headWid  = opts.ArrowHeadWidFrac * diagSpan;
P_G = opts.PG;
s_g = opts.sg;
% colors
if isempty(opts.Colors), cols = lines(m); else, cols=opts.Colors; end
if size(cols,1)<m || size(cols,2)~=3, error('Colors must be m-by-3.'); end

% index mapping consistent with solvers
numY0=m*n; numYJH=m*n*(n-1); numYT=m*n; numV=m*n; numU=m*n; % followed by cost vars
offset=0;
idx.y0  = @(r,j) offset + (r-1)*n + j;                    offset=offset+numY0;
idx.yjh = @(r,j,h) offset + (r-1)*n*(n-1) + map_jh(j,h,n); offset=offset+numYJH;
idx.yt  = @(r,j) offset + (r-1)*n + j;                    offset=offset+numYT;
idx.v   = @(r,j) offset + (r-1)*n + j;                    offset=offset+numV;
idx.u   = @(r,j) offset + (r-1)*n + j;                    offset=offset+numU;

% extract routes
Mm=zeros(m,1);
if m==1
    Mm=min(DXY(out1.orderD,1));
    [~,il]=sort(Mm);
    out.orderD=out1.orderD(il); out.startSource=out1.startSource(il); out.endSource=out1.endSource(il);
    out.bridgeSrc=out1.bridgeSrc(il); out.cost=out1.cost(il);
    orders={out.orderD}; startSrc=out.startSource; endSrc=out.endSource;
    bridge={out.bridgeSrc}; cost_expand=out.cost;
else
    for i=1:m
        Mm(i)=min(DXY(out1.ordersD{1,i},1));
    end
    [~,il]=sort(Mm);
    out.ordersD=out1.ordersD(il); out.startSource=out1.startSource(il); out.endSource=out1.endSource(il);
    out.bridgeSrc=out1.bridgeSrc(il); out.routeCost=out1.routeCost(il);
    orders=out.ordersD; startSrc=out.startSource; endSrc=out.endSource;
    bridge=out.bridgeSrc; cost_expand=out.routeCost;  
end

% table and stats
NumDest = cellfun(@numel, orders).';
Route=(1:m).'; StartSource=startSrc(:); EndSource=endSrc(:); 
if m==1 
    Cost=out.cost;
else
    Cost=out.routeCost';
end
T = table(Route,NumDest,StartSource,EndSource,Cost);
stats.Max=max(Cost); stats.Min=min(Cost); stats.Mean=mean(Cost); stats.Std=std(Cost,1); stats.Total=sum(Cost);

% visualization
if opts.ShowFigure
    % if ischar(opts.Layout) && strcmpi(opts.Layout,'auto')
    %     rows=ceil(sqrt(m)); cols1=ceil(m/rows);
    % else
    %     rows=opts.Layout(1); cols1=opts.Layout(2);
    % end
    if P_G==1
        cols1=3;
    else
        cols1=2;
    end
    rows=m;
    fig = figure('Name','m-ATSPP Routes','Color','w');
    
    tl = tiledlayout(rows,cols1,'TileSpacing','compact','Padding','compact');
    axx=nexttile(tl,[m 1]);
    if ~isempty(opts.FigTitle)
        title(opts.FigTitle,'FontWeight','bold','FontName','Times New Roman', 'FontSize',30);
    end
    plot_Gantt_Schedule(axx,G,R,sigma,out,m,cols);
    ax = gca;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 10;
    camroll(-90);
    set(gca,'YDir','reverse')

    colS_used=[0.20 0.20 0.70]; colS_unused=[0.75 0.78 0.90];
    colD_used=[0.10 0.55 0.10]; colD_unused=[0.75 0.90 0.75];

    for r=1:m
        nexttile(tl,[1 1]); hold on; axis equal; axis off;
        xlim(xlim_all); ylim(ylim_all);

        % plot sources
        for i=1:k
            used = (i==startSrc(r)) || (i==endSrc(r)) || any(bridge{r}==i);
            c = colS_used; if ~used, c = colS_unused; end
            scatter(SXY(i,1),SXY(i,2),opts.MarkerSizeS,'s','filled','MarkerEdgeColor','k','MarkerFaceColor',c);
            if opts.ShowNodeIDs
                text(SXY(i,1),SXY(i,2)+0.02*dy,sprintf('S%d',i),'HorizontalAlignment','center','FontSize',opts.FontSize);
            end
        end

        % plot destinations
        usedD = orders{r};
        for j=1:n
            c = colD_unused; if ismember(j,usedD), c = colD_used; end
            scatter(DXY(j,1),DXY(j,2),opts.MarkerSizeD,'o','filled','MarkerEdgeColor','k','MarkerFaceColor',c);
            if opts.ShowNodeIDs
                text(DXY(j,1),DXY(j,2)-0.02*dy,sprintf('D%d',j),'HorizontalAlignment','center','FontSize',opts.FontSize);
            end
        end

        % arrows in route color
        colRoute = cols(r,:);
        if ~isempty(usedD)
            i0 = startSrc(r); d1 = usedD(1);
            draw_arrow2(SXY(i0,:), DXY(d1,:), colRoute, opts.ArrowWidth, opts.ArrowHeadAt, headLen, headWid);
            for t=1:numel(usedD)-1
                j = usedD(t); h = usedD(t+1); ib = bridge{r}(t);
                draw_arrow2(DXY(j,:), SXY(ib,:), colRoute, opts.ArrowWidth, opts.ArrowHeadAt, headLen, headWid);
                draw_arrow2(SXY(ib,:), DXY(h,:), colRoute, opts.ArrowWidth, opts.ArrowHeadAt, headLen, headWid);
            end
            il=endSrc(r); jlast=usedD(end);
            draw_arrow2(DXY(jlast,:), SXY(il,:), colRoute, opts.ArrowWidth, opts.ArrowHeadAt, headLen, headWid);

            if opts.ShowStartEndText
                % text(SXY(i0,1),SXY(i0,2)+0.1*dy,'S','HorizontalAlignment','center','FontAngle','italic','FontSize',opts.FontSize);
                % text(SXY(il,1),SXY(il,2)-0.1*dy,'E','HorizontalAlignment','center','FontAngle','italic','FontSize',opts.FontSize);
                scatter(SXY(i0,1),SXY(i0,2),opts.MarkerSizeS,'diamond','filled','MarkerEdgeColor','k','MarkerFaceColor','r');
                scatter(SXY(il,1),SXY(il,2),2*opts.MarkerSizeS,'pentagram','filled','MarkerEdgeColor','k','MarkerFaceColor','m');
            end
        end

        title(sprintf('Route : %d ',r),sprintf('Total Sub-route Time = %.3f',Cost(r)),'Color',colRoute,'FontName','Times New Roman','FontSize',12);
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        hold off; view(2);           % ensure a flat 2D view
        daspect([2 1 1]);
    end
    
else
    fig = [];
    rows=0; cols1=0;
end

out.ordersD=orders; out.startSource=startSrc; out.bridgeSrc=bridge; out.endSource=endSrc;
out.cost_expand=cost_expand; out.T=T; out.stats=stats; out.figure=fig; out.layout=[rows cols1];
out.SXY=SXY; out.DXY=DXY; out.Colors=cols;
%% Separation Graph
if P_G==1
    nexttile(tl,[1 m]);
    G=graph(E(:,1),E(:,2)); ms=min(sigma);Ms=max(sigma); ss=(sigma-ms)/(Ms-ms)*2.5*s_g+s_g;
    aa = adjacency(G);
    gplot(aa, DXY, ['-','k']); axis equal; axis tight; grid on; hold
    scatter(DXY(:,1),DXY(:,2),ss,'k','filled');
    for i=1:m
        is=out.ordersD{1,i};
        scatter(DXY(is,1),DXY(is,2),ss(is),cols(i,:),'filled');
    end
       title('Segmented Routes','FontWeight','bold','FontSize',30,'FontName','Times New Roman');
       ax = gca;                     % get current axes
       ax.FontName = 'Times New Roman';
       ax.FontSize = 20;
end
end
function opts = validateAndFill(user, defaults)

    %--- merge ---
    opts = defaults;
    fields = fieldnames(user);
    for k = 1:numel(fields)
        name = fields{k};
        if ~isfield(defaults, name)
            error('Unknown option "%s". Check spelling.', name);
        end
        opts.(name) = user.(name);
    end

end
function id = map_jh(j,h,n)
    if j==h, error('map_jh: j==h'); end
    if h<j, pos=h; else, pos=h-1; end
    id = (j-1)*(n-1) + pos;
end

function draw_arrow2(p0, p1, col, lw, headAt, headLen, headWid)
    plot([p0(1) p1(1)],[p0(2) p1(2)],'-','Color',col,'LineWidth',lw);
    v = p1 - p0; L = hypot(v(1),v(2)); if L==0, return; end
    u = v / L;
    switch lower(headAt)
        case 'middle'
            pt = p0 + 0.5*v;
        otherwise
            pt = p1;
    end
    base = pt - headLen*u;
    perp = [-u(2), u(1)];
    left  = base + 0.5*headWid*perp;
    right = base - 0.5*headWid*perp;
    patch('XData',[pt(1) left(1) right(1)],'YData',[pt(2) left(2) right(2)],'FaceColor',col,'EdgeColor',col);
end
function plot_Gantt_Schedule(axx,tau_go,tau_back,sigma,outM,m,Col)
    hold(axx,'on');
    
    off = 0;
    if m==1
        sources=[outM.startSource,outM.bridgeSrc,outM.endSource];
        destination=outM.orderD;
        ss=string(destination);
        travel=tau_go(sub2ind(size(tau_go),sources(1:end-1),destination));
        retur=tau_back(sub2ind(size(tau_back),destination,sources(2:end)));
        trips = table( ...
        ss', ...       % dest_id
        travel',  ...       % travel
        sigma(destination),  ...       % service
        retur', ...      % return
        'VariableNames', ...
        {'dest_id','travel','service','return'});
           
        plotTravelGanttBlock(axx, trips, Col, off, m);   % reddish 
    else
        for ll=1:m
            sources=[outM.startSource(ll),outM.bridgeSrc{1,ll},outM.endSource(ll)];
            destination=outM.ordersD{1,ll};
            ss=string(destination);
            travel=tau_go(sub2ind(size(tau_go),sources(1:end-1),destination));
            retur=tau_back(sub2ind(size(tau_back),destination,sources(2:end)));
            trips = table( ...
            ss', ...       % dest_id
            travel',  ...       % travel
            sigma(destination),  ...       % service
            retur', ...      % return
            'VariableNames', ...
            {'dest_id','travel','service','return'});
               
            off = plotTravelGanttBlock(axx, trips, Col(ll,:), off, ll);   % reddish
        end
    end
end   

function nextOffset = plotTravelGanttBlock(ax, trips, lineColor, xOffset, m, sequential)
% Draw ONE coloured vertical-time Gantt block on axes AX.
%
% nextOffset = plotTravelGanttBlock(ax, trips, colour, xOffset, sequential)
%
%   ax         : target axes handle (gca for current axes)
%   trips      : table with dest_id, travel, service, return, [start_time]
%   colour     : MATLAB colour
%   xOffset    : integer column offset where this block starts
%   sequential : true  (default)  → calculate start_time cumulatively
%                false            → missing start_time set to 0
%
% Returns nextOffset = xOffset + #destinations_in_this_block

if nargin < 6, sequential = true; end           % default

%--- cosmetics you may tweak -------------------------------------------
boldLW     = 3;     thinLW = boldLW/2;     markerSize = 6;
%-----------------------------------------------------------------------

% ----- normalise data --------------------------------------------------
if ~istable(trips), trips = struct2table(trips); end

need = {'dest_id','travel','service','return'};
assert(all(ismember(need,trips.Properties.VariableNames)), ...
    'trips must contain dest_id, travel, service, return.');

if ~ismember('start_time',trips.Properties.VariableNames)
    if sequential
        runs = trips.travel + trips.service + trips.return;
        trips.start_time = [0; cumsum(runs(1:end-1))];   % chain jobs
    else
        trips.start_time = zeros(height(trips),1);
    end
end

% ----- map destination → local columns --------------------------------
[destU,~,xIdx] = unique(trips.dest_id,'stable');
x      = xOffset + xIdx;
nDest  = numel(destU);

% ----- plotting --------------------------------------------------------
hold(ax,'on');
for i = 1:height(trips)
    y0 = trips.start_time(i);
    y1 = y0 + trips.travel(i);
    y2 = y1 + trips.service(i);
    y3 = y2 + trips.return(i);

    % travel
    line(ax,[x(i) x(i)],[y0 y1],'Color',lineColor,'LineWidth',thinLW);
    plot(ax,x(i),(y0+y1)/2,'^','MarkerSize',markerSize,'Color',lineColor,...
         'MarkerFaceColor',lineColor);

    % service
    % line(ax,[x(i) x(i)],[y1 y2],'Color',lineColor,'LineWidth',boldLW);
    plot(ax,[x(i) x(i)],[y1 y2],'-','MarkerSize',markerSize,'Color',lineColor,...
         'MarkerFaceColor',lineColor,'LineWidth',boldLW);
    % return
    line(ax,[x(i) x(i)],[y2 y3],'Color',lineColor,'LineStyle','--',...
         'LineWidth',thinLW);
    % plot(ax,x(i),y3,'^','MarkerSize',markerSize,'Color',lineColor,...
    %      'MarkerFaceColor',lineColor);
end

% ----- update X-ticks / labels ----------------------------------------
oldTicks   = get(ax,'XTick');
oldLabels  = cellstr(get(ax,'XTickLabel'));
if numel(oldTicks) ~= numel(oldLabels)
    % brand-new axes or previous mis-match → start fresh
    oldTicks  = [];  oldLabels = {};
end

newTicks   = xOffset + (1:nDest);
newLabels  = cellstr(string(destU));

if m==1
    set(ax,'XTick',newTicks,'XTickLabel',newLabels,...
       'YDir','normal','XDir','normal','Box','on'); xtickangle(30);
else
    allTicks   = [oldTicks(:) ; newTicks(:)];
    allLabels  = [oldLabels   ; newLabels  ];
    [allTicks, idx] = unique(allTicks,'stable');   % keep order, drop dups
    allLabels  = allLabels(idx);

    set(ax,'XTick',allTicks,'XTickLabel',allLabels,...
       'YDir','normal','XDir','normal','Box','on'); xtickangle(90);
end

xlabel(ax,'Model Installation Zone ID','FontSize',20,'FontName','Times New Roman');  ylabel(ax,'Time','FontSize',20,'FontName','Times New Roman');  grid(ax,'on');

% ----- limits ----------------------------------------------------------
xlim(ax,[0.5, xOffset + nDest + 0.5]);

totTime = trips.start_time + trips.travel + trips.service + trips.return;
yMax    = max(totTime);  if yMax<=0, yMax = 1; end
ylim(ax,[-0.05*yMax, 1.05*yMax]);

% ----- next block offset ----------------------------------------------
nextOffset = xOffset + nDest;
end

