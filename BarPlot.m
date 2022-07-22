function H = BarPlot(tt,X)
% =======================================================================
% Creates a bar graph with positive data values stacked on the positive
% quadrant and negative data values stacked on the negative quadrant
% =======================================================================
% H = BarPlot(X)
% -----------------------------------------------------------------------
% INPUT
%   - X: data to plot [nobs x nvars]
%   - tt: dates
% -----------------------------------------------------------------------
% OUTPUT
%   - H: handle to graph
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com
% Miguel Herculano, Note: slightly altered so that the colours of the bars
% on the negative and positive domain would be the same =)
n = size(X,2);
colormap(parula)
H(1,:) = bar(tt,(X).*(X>0),'stacked'); 
hold on;
H(2,:) = bar(tt,(X).*(X<0),'stacked');
% to ensure colors are the same
for k=1:n
H(1,k).FaceColor=H(2,k).FaceColor;
end
