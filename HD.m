function HD_struct = HD(B,SIGMA,E1,hor,n);
% --- DESCRIPTION: Historical Decomposition (HD)
% This script performs HD for a VAR described by B and SIGMA (see below).
% Historical decomposition (HD) represent each endogenous variable included
% in the VAR as a sum of all the shocks in the system. HD allows us to assess
% the historical contribution of each structural shock in driving each variable
% away from its long-run mean.
%
%  Main object is hd(:,:,k) gives the decomposition of y1 due to e1,e2,e3,...
% ------------------------ INPUT
% - B : Autoregressive polynomial as given by matlab's default estimate()
% - SIGMA : reduced form VAR COvariance matrix
% - E1 : reduced form VAR innovations
% - hor : horizon for the Wold representation of the VAR
% - n : #variables in the VAR
% - m_y and s_y : are the mean and sdev. of the original time-series for
% reconstruction. Just leave them in blank and the program with perform HD.
% on the standardized series. 
% ------------------------ OUTPUT

% HD_struct.y_rec : reconstruct y (data) based on the wold representation -
% y_t = C(L)u_t
% HD_struct.hd_pc : HD in percentage points
% HD_struct.hd_rec : HD in original units (same as original data)
% -------------------------------------------------------------------------
% Miguel C. Herculano, mcherculano.github.io, miguelcbherculano@gmail.com
% August 2021, University of Nottingham
% -------------------------------------------------------------------------
t = size(E1,1);
C = Wold(B,SIGMA,hor,n);

y_rec = zeros(n,t);
hd = zeros(n,t,n); %  hd(:,:,k) gives the decomposition of y1 due to e1,e2,e3,...
hd_pc = zeros(n,t,n);
hd_rec = zeros(n,t,n);
for h=1:hor
    e_lag = lagmatrix(E1,h-1);
    e_lag(isnan(e_lag)) = 0;
    aux = C(:,:,h)*e_lag';
    for k =1:n
        aux1 = C(k,:,h)'.*e_lag';
        hd(:,:,k) = hd(:,:,k) + aux1; %#serie;#t;#shock
    end
    y_rec = y_rec + aux;
end

% Now convert HD into 0-1 scale
for k =1:n 
    hd_pc(:,:,k) = hd(:,:,k)./y_rec(k,:);
    hd_rec(:,:,k) = [y_rec(k,:)'.*hd_pc(:,:,k)']'; 
end
% rescale data consistent with initial time-series unit  
%y_rec = m_y + s_y.*y_rec'; 


% OUTPUT
HD_struct.y_rec = y_rec;
HD_struct.hd_pc = hd_pc;
HD_struct.hd_rec = hd_rec;

end


% -------- AUXILIARY FUNCTIONS

function C = Wold(B,SIGMA,hor,n,varargin)

% INPUTS:
% [1] B - matrix of coefficients from the VAR 
% [2] shock: 1 - 

%Insert appropriate restrictions:
% if isempty(varargin)==0 & strcmp(varargin{1},'cholimpact_sd') 
%     % Cholesky on impact restriction
%     S_ir=chol(SIGMA,'lower'); 
%     d = diag(diag(S_ir));
%     S_ir = inv(d)*S_ir;
% elseif isempty(varargin)==0 & strcmp(varargin{1},'none')
%     % No restriction
%     S_ir=eye(n);
% elseif isempty(varargin)==0 & strcmp(varargin{1},'cholimpact_unit')
%     S_ir=chol(SIGMA,'lower'); 
%     d = diag(diag(S_ir));
%     S_ir = inv(d)*S_ir;
% else     

        % Cholesky on impact (default)
        S_ir=chol(SIGMA,'lower');
        d = diag(diag(S_ir));
        S_ir = inv(d)*S_ir;
        %warning('Misspecified or missing restriction. Set by default to cholimpact.')
%    end

% i) Find the companion form matrix
A=companion(B);

% ii) Create a 3d array to store the powers of A in each layer
wald=zeros(size(A,1),size(A,2),hor);
for i=1:hor 
wald(:,:,i)=A^(i-1);
end

% iii) Select the n*n upper left elements to get the Cj of the wald representation
C=wald(1:n,1:n,:);

% iv) Multiply times S 
for i=1:hor
    C(:,:,i)=C(:,:,i)*S_ir;   
end

% % STEP 8: Rearrange the ir in a bidimenstional matrix easier to plot
% C_plot=reshape(permute(C, [3 2 1]),[],n*n);
% 
% STEP 9: Compute structural errors:
%e_chol=(S_ir\e_ir')';

end


function [comp]=companion(PI)
% Inputs:
% PI - matrix of coefficents

%Outputs
% comp - n*p suquare matrix in the companion form

s=size(PI,2);
comp = (eye(s));
comp=[PI; comp];
comp=comp(1:s,:);
end