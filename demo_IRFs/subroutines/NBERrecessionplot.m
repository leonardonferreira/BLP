function varargout = NBERrecessionplot(varargin)
%RECESSIONPLOT Add recession bands to time series plot
%
% Syntax:
%
%   recessionplot
%   recessionplot(param,val,...)
%   hBands = recessionplot(...)
%
% Description:
%
%   Overlay shaded recession bands on a time series plot.
%
% Optional Input Parameter Name/Value Pairs:
%
%   NAME            VALUE
%
%   'axes'          Handle to axes displaying a time series plot with
%                   serial date numbers on the horizontal axis. The default
%                   is gca.
%
%	'recessions'	numRecessions-by-2 matrix of serial date numbers
%                   indicating the beginning (first column) and end (second
%                   column) of historical recessions. The default is the
%                   U.S. recession data in Data_Recessions.mat, reported
%                   by the National Bureau of Economic Research.
%
% Output Arguments:
%
%   hBands - Vector of handles to the recession bands.
%
% Notes:
%
%   o RECESSIONPLOT requires dates on the horizontal axis of a time series
%     plot to be expressed as serial date numbers. To convert other date
%     information to this format before plotting, use the DATENUM command.
%
%   o Use the output handles to change the color and opacity of the
%     recession bands by setting their 'FaceColor' and 'FaceAlpha'
%     properties. This may be necessary to achieve satisfactory displays
%     when working with certain monitors and projectors.
%
% Example:
%
%   load Data_CreditDefaults   
%   X0 = Data(:,1:4);
%   T0 = size(X0,1);
%   
%   % Convert dates to serial date numbers:
% 
%   dates = datenum([dates,ones(T0,2)]);
% 
%   % Create plot:
% 
%   plot(dates,X0,'LineWidth',2)
%   set(gca,'XTick',dates(1:2:end))
%   datetick('x','yyyy','keepticks')
%   xlabel('Year') 
%   ylabel('Level')
%   axis tight
%  
%   % Add recession bands:
% 
%   recessionplot
%
% References:
% 
%   [1] National Bureau of Economic Research. "US Business Cycle Expansions
%       and Contractions." http://www.nber.org/cycles.html.
%
% See also DATENUM.

% Copyright 2012 The MathWorks, Inc.
% $Revision: 1.1.6.2 $ $Date: 2011/11/09 16:44:27 $

% Parse inputs and set defaults:

parseObj = inputParser;
parseObj.addParamValue('axes',get(get(0,'CurrentFigure'),'CurrentAxes'),@axesCheck);
parseObj.addParamValue('recessions',[],@recessionsCheck);

parseObj.addParamValue('dates',[]);
parseObj.addParamValue('ymin',[]);
parseObj.addParamValue('ymax',[]);


parseObj.parse(varargin{:});

hAx = parseObj.Results.axes;
Recessions = parseObj.Results.recessions;

% if isempty(hAx)
%         
% 	error(message('econ:recessionplot:NoCurrentAxes'))
%           
% else
%     
%     axes(hAx)
%     hold on
%     
% end

if isempty(Recessions)
    
    load Data_Recessions
    
end

Dates=parseObj.Results.dates;

Ymin =parseObj.Results.ymin;
Ymax =parseObj.Results.ymax;


% dateRange = get(hAx,'XLim'); % Date range of current axes
% dataRange = get(hAx,'YLim'); % Data range of current axes

dateRange = [Dates(1) Dates(end)]; % Date range of current axes
dataRange = [Ymin Ymax]; % Data range of current axes


inPlot = (Recessions(:,1) < dateRange(2)) & (Recessions(:,2) > dateRange(1));
RPlot = Recessions(inPlot,:); % Recessions to plot

numRecessions = size(RPlot,1);
hBands = zeros(numRecessions,1); % Recession band handles

for n = 1:numRecessions
    
%     hBands(n) = patch([RPlot(n,1),RPlot(n,1),RPlot(n,2),RPlot(n,2)],...
%                       [dataRange,fliplr(dataRange)],'k',...
%                       'FaceAlpha',0.1,'EdgeColor','none');
    timeLine=Dates(Dates>=RPlot(n,1) & Dates<=RPlot(n,2));

    hBands(n) = area(timeLine,repmat(dataRange(end),1,length(timeLine)),...
        'EdgeColor',[.9 .9 .9],'FaceColor',[.9 .9 .9],...
        'basevalue',dataRange(1));
    
end

xlim(dateRange)
% hold off

nargoutchk(0,1);

if nargout > 0
    
    varargout{1} = hBands;
end

%-------------------------------------------------------------------------
% Check value of 'axes' parameter
function OK = axesCheck(hAx)
    
if ~ishghandle(hAx,'axes')

    error(message('econ:recessionplot:AxesHandleInvalid'))

else

    OK = true;

end

%-------------------------------------------------------------------------
% Check value of 'recessions' parameter
function OK = recessionsCheck(Recessions)
    
if ~isnumeric(Recessions) || ~ismatrix(Recessions) || size(Recessions,2) ~= 2

    error(message('econ:recessionplot:RecessionsMatrixInvalid'))

else

    OK = true;

end