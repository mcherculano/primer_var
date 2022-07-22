% =========================================================================
% =========== A PRIMER ON VARS
% =========================================================================
% This script accompanies the jupyter notebook "A PRIMER ON VARS".
% =========================================================================
% Miguel C. Herculano, UoN September '21
%% 2 Data: First load the data
clear all; clc;
% First load the data
load 'Data_USEconModel'

% Plot the time-series
figure
subplot(3,1,1)
plot(DataTable.Time,DataTable.CPIAUCSL,'b');
title('CPI')
grid on
subplot(3,1,2);
plot(DataTable.Time,DataTable.UNRATE,'b');
title('Unemployment rate')
grid on
subplot(3,1,3);
plot(DataTable.Time,DataTable.TB3MS,'b')
title('3-mo T-bill')
grid on

cpi_inf = price2ret(DataTable.CPIAUCSL); % dif-log
%irate = diff(DataTable.TB3MS); % 1st dif

% notice that by taking dif-log/dif we loose 1 observation and therefore UNRATE which is in levels needs to be adjusted
urate = DataTable.UNRATE(2:end,:); 
irate = DataTable.TB3MS(2:end,:); 

% now put the three variables together in an array called 'Data'
Data = array2timetable([cpi_inf urate irate],...
    'RowTimes',DataTable.Time(2:end),'VariableNames',{'CPI_Inf' 'U_Rate' 'Int_rate'});
    
% standardize variables
 y = standardize_miss(Data{:,1:3});
% store mean and std. for later re-scaling
m_y = nanmean(Data{:,1:3});
s_y = nanstd(Data{:,1:3});

% finally remove all missing values from the begginning of the sample
idx = all(~ismissing(y),2);
y = y(idx,:);

%% 3. What is a Vector Autoregression (VAR) ? 

% construct a model with #series=3 and #lags=2. We will estimate the model with 2 lags. 
model = varm(3,2);
% estimate the model
[var_est,SE_est,logL1,E1] = estimate(model,y);

%% 4.2 Identification Problem 

% recover the Variance-Covariance matrix Sigma from the VAR structure and perform a cholesky decomposition
Sigma = var_est.Covariance;
B_chol = chol(Sigma)';

%% 4.3 Impulse Response Analysis 

% "irf" below is a 3D tensor whose dimensions are (H,#var.shock,#var.response) 
[response,lower,upper] = irf(var_est);
% WARNING: give a different name to an output and a function - otherwise this will cause errors in the notebook.
% Now, for instance, let us examine the IRF of the urate wrt. a MP shock, for a 1 year horizon (1:4 quarters):
response(1:4,3,2)
% Next, plot the IRF for all shocks and all variables neatly
% start by defining some basic quatities
N = 3; % No. variables in the VAR
H = size(response,1); % Horizon of the IRF (by default matlab assumes H=20)
vnames = {'inflation','urate','irate'};
nshock = {'Agg. Supply shock','Agg. Demand shock','Monetary Pol. shock'};
mm = 1;
    for jj = 1:N %nshocks
        for ii = 1:N %nvars
                subplot(N,N,mm)
                plot(1:H,squeeze(lower(:,jj,ii)),'r--')
                hold on
                plot(1:H,squeeze(response(:,jj,ii)),'k-')
                hold on
                plot(1:H,squeeze(upper(:,jj,ii)),'r--')
                if jj==1 %&& jj==1
                    title(vnames(ii), 'FontSize', 7)
                end
                if ii==1 %&& jj==1
                    ylabel(nshock(jj), 'FontSize', 7, 'FontWeight','bold')
                end                
                %title('...')
                xlim([1 H])
                set(gca,'XTick',0:3:H)    
                mm = mm+1;
        end
    end

%% 4.4 Forecast Error Variance Decomposition 

% compute fevd get a 3D tensor
[decomp,decomp_lower,decomp_upper] = fevd(var_est);

% The fevd in element (t,i,j) is the contribution (in % terms) to the variance decomposition 
% of variable j attributable to an innovation shock to variable i at time t.
% Now get the fevd for inflation wrt. a Monetary Policy shock for h={1, 4, 8, 20} quarters.
decomp([1 4 8 20],3,1)

%% 4.5 Historical Decomposition

% define some quantities needed for the calculation
B = cell2mat( var_est.AR);
SIGMA = var_est.Covariance;
hor = 20; % 
% perform HD. - for a peek into the black box, open the function "HD" and see what is going on.
% Long story short: 1st - we write the VAR in MA form; 2nd - we explore the elements in C(L) to desintangle the contribution of 
HD_struct = HD(B,SIGMA,E1,hor,N);

% Now explore the structure and plot HD
% Note: hd(:,:,k) gives the decomposition of y_k due to e1,e2,e3,...
tt = Data.Time(3:end);
figure
H1 = BarPlot(tt,HD_struct.hd_rec(:,:,2)');
hold on 
H2 = plot(tt, HD_struct.y_rec(2,:),'k','LineWidth',1);
legend('Agg. Supply shock','Agg. Demand shock','Monetary Pol. shock')

%% 5. Forecasting with VARs 
% set forecast horizon
h = 16;
% save 16 observations for out-of-sample analysis
Y0 = y(1:end-16,:);
var_f = estimate(model,Y0);
% Forecast based on Y0
[Y_fore, YMSE] = forecast(var_f,h,Y0);

% Calculate confidence bands
extractMSE = @(x)diag(x)';
MSE = cellfun(extractMSE,YMSE,'UniformOutput',false);
SE = sqrt(cell2mat(MSE));
YFI = zeros(h,var_est.NumSeries,2);
YFI(:,:,1) = Y_fore - 2*SE;
YFI(:,:,2) = Y_fore + 2*SE;

% re-scale time-series and forecasts
y_scale = m_y + y.*s_y;
Y_fore_scale = m_y + Y_fore.*s_y;
YFI_scale = m_y + YFI.*s_y;
Y1 = y_scale(end-15:end,:);

% Now plot forecasts 
fh = dateshift(Data.Time(end-20),'end','quarter',1:h);
vnames = {'CPI Inflation' 'URate' 'Int. rate'};
figure;
for k=1:N
    subplot(1,N,k)
    h1 = plot(Data.Time(end-60:end-20),y_scale(end-60:end-20,k),'b');
    hold on;
    h2 = plot(fh,Y_fore_scale(:,k),'k','LineWidth',1);
    hold on 
    h3 =  plot(fh,Y1(:,k),'--b');
    hold on 
    h4 = plot(fh,YFI_scale(:,k,1),'k--');
    hold on
    h5 = plot(fh,YFI_scale(:,k,2),'k--');
    title(vnames{k});
    h = gca;
    fill([Data.Time(end-20) fh([end end]) Data.Time(end-20)],h.YLim([1 1 2 2]),'k',...
        'FaceAlpha',0.1,'EdgeColor','none');
    legend([h1 h2 h3],'Sample','Forecast','observed', 'Location', 'southoutside')
    hold off;
end