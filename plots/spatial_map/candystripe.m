function hg = candystripe(h,varargin)
% CANDYSTRIPE Candystripes patches
%
% CANDYSTRIPE creates a candystripe pattern to fill all patch-like objects
%  in the current axis. Patch like objects are any objects that are plotted
%  as patches in MATLAB. This includes, but is not limited to: bars, fills,
%  patches, etc.
%
% CANDYSTRIPE(h) creates a candystripe pattern to fill the object specified
%  by the handle H.
%
% CANDYSTRIPE(...,'param1',value1,'param2',value2,...) specifies optional
%  parameters to be used in the candystripe creation. The parameters, data
%  types, and default values are listed below.
%
%  Optional Paramters:
%   PARAMETER NAME      DATA TYPE           ORIGINAL DEFAULT VALUE
%
%   'Angle'             scalar              45
%       Specifies the angle in degrees at which the candystripes are to be
%       drawn.
%
%   'Color'         [1x3] OR colorspec      'w'
%       Specifies the color of the candystripes that are to be drawn
%
%   'Units'              string             'pixels'
%       Specifies that units that are to be used when drawing the
%       candystripes. This influences how the angle appears. For example,
%       if the units are in pixels, a candystripe with a 45 degree angle
%       will appear as such, regardless of what the data aspect ratio is.
%       Valid options are {'inches' , 'centimeters' , 'normalized' , 
%       'points' , 'pixels' , 'characters' , and 'native'}, where 'native'
%       corresponds to the plot units
%
%   'Width'              scalar             10
%       Specifies the width of the candystripe to be drawn. The units of
%       this are specified by the optional input parameters UNITS. Note: By
%       default, a cap is placed on how large candystripes may be. This is
%       done so that the candystriping won't end up being a solid color.
%
% NOTE: If the object is already in the legend, CANDYSTRIPE will
%  automatically update the thumbnail provided in the legend. However, if
%  it is not in the legend, then the thumbnail will appear as the patch
%  originally did.
%
%
%   EXAMPLE 1: All patches with legend
%       x = rand(3,1);
%       y = rand(3,1);
%       figure
%       patch(x,y,'r');
%       legend('Triangle')
%       candystripe;
%
%   EXAMPLE 2: Bar plots
%       data = rand(10,2);
%       figure
%       h = bar(data,'stacked');
%       caxis([3 100])
%       legend('Current','Future')
%       candystripe(h(2));
%
%   EXAMPLE 3: Thin candystipes at a 60 degree angle
%       z = peaks;
%       figure
%       [c, h] = contourf(z);
%       candystripe(h,'Units','native','width',1,'Angle',60)
%
%   EXAMPLE 4: Patches with FaceColor set to interp
%       t = 0:pi/5:2*pi;
%       figure
%       axis equal
%       p = patch(sin(t),cos(t),1:length(t),'EdgeColor','k','LineWidth',1.5);
%       candystripe(p)
%
%   EXAMPLE 5
%       data = randn(1000,10);
%       xs = -3:0.1:3;
%       figure;
%       h = hist(data,xs);
%       candystripe(h(5:10));
%
% J Sullivan January 2012


deleteExtras = false;

% All Childeren of the axis
if nargin > 0 && ischar(h)
    varargin = [{h} varargin];
end

% No handle specified. Default to all patches on the axis
if nargin == 0 || ischar(h);
    h = get(gca,'Children');
    deleteExtras = true;
end

% Check if H is a handle
if ~ishandle(h);
    error('JLS:CandyStripes:NotHandle','"H" is not a valid handle');
end

% Recursive call to deal with a vector of handles
n = length(h);
if n > 1;
    hg = NaN(n,1);
    for ii = 1:n
        hg(ii) = candystripe(h(ii),varargin{:});
    end
    if deleteExtras
        hg(isnan(hg)) = [];
    end
    if nargout == 0
        clear hg;
    end
    return;
end

% Parse optional inputs
if length(varargin) == 1 && isstruct(varargin{1});
    S = varargin{1};
else
    S = parseInput(varargin{:});
end

% Extract only the patches
[h, ax, ho] = getPatches(h);
if isempty(h)
    if nargout > 0
        hg = NaN;
    end
    return;
end

if length(h) > 1
    candystripe(h,S);
    return;
end

% Prepare for units transform
hunits = get(ax,'Units');
xl = get(ax,'XLim');
yl = get(ax,'YLim');
pnts = [xl(1) yl(1) diff(xl) diff(yl)];
if ~strcmpi(S.units,'native')
    set(ax,'units',S.units);
    units = get(ax,'Position');
else
    units = pnts;
end
set(ax,'units',hunits);

% Get data
xd = get(h,'XData');
yd = get(h,'YData');

% Change units
[xd, yd] = changeCoord(xd,yd,pnts,units);

% Get Hatch Pattern
if ismatrix(xd)
    xcs = []; ycs = [];
    for ii = 1:size(xd,2)
        [xt, yt] = getHatchPattern(xd(:,ii),yd(:,ii),S);
        xcs = [xcs xt]; ycs = [ycs yt]; %#ok<AGROW>
    end
else
    [xcs, ycs] = getHatchPattern(xd,yd,S);
end

% Plot the hatch pattern
wasHold = ishold;
hold on
hcs = NaN(length(xcs),1);
for ii = 1:length(xcs)
    x1 = xcs{ii}; y1 = ycs{ii};
    if isempty(x1); continue; end
    [x1, y1] = changeCoord(x1,y1,units,pnts);
    if sum(isnan(x1)) == 0
        hcs(ii) = patch(x1,y1,S.color,'EdgeColor','k','FaceColor',S.color);
    else
        ind1 = 1;
        nnans = sum(isnan(x1));
        fs = find(isnan(x1));
        for jj = 1:nnans+1
            if jj > nnans
                f = length(x1)+1;
            else
                f = fs(jj);
            end
            ind2 = f-1;
            hthis = patch(x1(ind1:ind2),y1(ind1:ind2),S.color,'EdgeColor','k','FaceColor',S.color);
            ind1 = f+1;
            if jj == 1
                hcs(ii) = hthis;
            else
                hcs(end+1) = hthis; %#ok<AGROW>
            end
        end
    end
end
hcs(isnan(hcs)) = [];
if ~isempty(hcs)
    set(hcs,'EdgeColor','none')
end

% Make hg group
hg = hggroup;

% Redraw original patch
p = get(h,'Parent');
hl = copyobj(h,p);
set(hl,'FaceColor','none');

% Parent candy stripes to original
ch = get(p,'Children');
ch(ismember(ch,[hcs; hl; hg])) = [];
ch(ch == h) = hg;
set(h,'Parent',hg);
set(hcs,'Parent',hg);
set(hl,'Parent',hg);
set(hg,'Parent',p);
set(hg,'Tag','candystripe');
set(p,'Children',ch);

% Handle the legend
an = get(hg,'Annotation'); an.LegendInformation.IconDisplayStyle = 'on';
[~,hlbl,outh,outm] = legend;
if any(outh == ho | outh == hg)
    candystripe(hlbl(strcmpi(get(hlbl,'tag'),outm(outh == ho | outh == hg))),S);
end

% Reset the hold state
if ~wasHold
    hold off;
end

% Don't force the output
if nargout == 0;
    clear hg;
end

function [x, y] = changeCoord(x,y,from,to)

% Simple coordinate transform
x = (x - from(1)).*to(3)./from(3) + to(1);
y = (y - from(2)).*to(4)./from(4) + to(2);

function [xcs, ycs] = getHatchPattern(xd,yd,S)

% get the bounding box
xlim = [min(xd) max(xd)];
ylim = [min(yd) max(yd)];
xs = linspace(xlim(1),xlim(2),10);
ys = linspace(ylim(1),ylim(2),10);
[xs, ys] = meshgrid(xs,ys);

% Get the rotation matrix
ang = 90 - S.angle;
M = [cosd(ang) -sind(ang); sind(ang) cosd(ang)];

% Rotate the bounding box
if exist('mtimesx','file')
    o = mtimesx(M,[reshape(xs,1,1,[]); reshape(ys,1,1,[])]);
else
    for ii = numel(ys):-1:1
        o(:,:,ii) = M*[xs(ii); ys(ii)];
    end
end

% Get the candy stripes
alim = [min(o(1,:,:),[],3) max(o(1,:,:),[],3)];
blim = [min(o(2,:,:),[],3) max(o(2,:,:),[],3)];

S.width = diff(alim)./max(min(diff(alim)./S.width,100),5);
as = alim(1)-S.width:S.width:alim(2)+S.width;
bs = blim;
n = length(as);
ainds = bsxfun(@plus,[0 1 1 0],(1:2:n-1)');
binds = repmat([1 1 2 2],floor(n/2),1);

% Get the reverse rotation matrix
M = [cosd(-ang) -sind(-ang); sind(-ang) cosd(-ang)];

% Rotate the candy stripes back
if exist('mtimesx','file')
    o2 = mtimesx(M,[reshape(as(ainds),1,1,[]); reshape(bs(binds),1,1,[])]);
else
    x1 = as(ainds); y1 = bs(binds);
    for ii = numel(y1):-1:1
        o2(:,:,ii) = M*[x1(ii); y1(ii)];
    end
end

% Reshape
xs = o2(1,:,:); xs = reshape(xs,[],4)';
ys = o2(2,:,:); ys = reshape(ys,[],4)';

% Find the intersection between the patch and the stripes
cnt = 0;
ycs = {[]}; xcs = {[]};
for ii = 1:size(xs,2)
    [xi,yi] = poly2cw(xs(:,ii),ys(:,ii));
    [xd,yd] = poly2cw(xd,yd);
    [xo, yo] = polybool('intersection',xi,yi,xd,yd);
    if ~isempty(xo)
        cnt = cnt + 1;
        xcs{cnt} = xo;
        ycs{cnt} = yo;
    end
end

function [h, ax, ho] = getPatches(ho)

% Get the original handle, axis, and the patch handle
if ~strcmpi(get(ho,'type'),'patch');
    ax = get(ho,'Parent');
    type = '';
    h = ho;
    while ~strcmpi(type,'patch')
        h = get(h,'Children');
        if isempty(h); break; end
        h(strcmpi(get(h,'tag'),'candystripe')) = [];
        type = get(h,'Type');
    end
else
    h = ho;
    ax = get(ho,'Parent');
end

if ~strcmpi(get(ax,'Type'),'axes')
    type = '';
    while ~strcmpi(type,'axes')
        ax = get(ax,'Parent');
        if isempty(ax); break; end
        type = get(ax,'Type');
    end
end

function S = parseInput(varargin)

% Error checking
isChars = cellfun(@ischar,varargin(1:2:end));
if ~all(isChars)
    error('JLS:CandyStripe:NotValidOptInput','One or more optional inputs are not strings.')
end
if rem(length(varargin),2) ~= 0
    error('JLS:CandyStripe:OptInputNotInPairs','Optional inputs must be in parameter value pairs.')
end

% Defaults
S.width      = 10;
S.units      = 'pixels';
S.color      = 'w';
S.angle      = 45;

% Get valid optional inputs
vInputs = fieldnames(S);
vInputs = strcat({'"'},vInputs,{'", '});
vInputs(end-1) = strcat(vInputs(end-1),{'and '});
vInputs(end)   = strrep(vInputs(end),'", ','"');
vInputs = [vInputs{:}];

% Loop over optional inputs
for ii = 1:2:length(varargin)
    lbl = lower(varargin{ii});
    if ~isfield(S,lbl);
        warning('JLS:CandyStripe:NotValidOptInput','"%s" is not a valid optional input. Valid inputs are %s.',lbl,vInputs)
    end
    S.(lbl) = varargin{ii+1};
end