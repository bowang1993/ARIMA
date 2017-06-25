function sys = armax(varargin)
%ARMAX  Estimate ARMAX polynomial model using time domain data.
%
%  M = ARMAX(Z, [na nb nc nk])
%    estimates an ARMAX model, M, represented by:
%       A(q) y(t) = B(q) u(t-nk) +  C(q) e(t)
%    where:
%       na = order of A polynomial     (Ny-by-Ny matrix)
%       nb = order of B polynomial + 1 (Ny-by-Nu matrix)
%       nc = order of C polynomial     (Ny-by-1 matrix)
%       nk = input delay (in number of samples, Ny-by-Nu entries)
%       (Nu = number of inputs; Ny = number of outputs)
%
%    The estimated model, M, is delivered as an @idpoly object. M contains
%    the estimated values for A, B, and C polynomials along with their
%    covariances and structure information.
%
%    Z is the time-domain estimation data given as an IDDATA object. Type
%    "help iddata" for more information. You cannot use frequency-domain
%    data for estimation of ARMAX models. na, nb, nc and nk are the
%    polynomial orders associated with the ARMAX model.
%
%  M = ARMAX(Z, [na nb nc nk], 'Name1', Value1, 'Name2', Value2,...)
%    specifies additional model structure properties as name-value
%    pairs. You can specify as one or more of the following:
%    'InputDelay': Specify input delay as a double vector of length equal
%                  to number of inputs. Entries must be nonnegative
%                  integers denoting the delay as multiples of sample
%                  time.
%      'ioDelay': Input-to-output delay (double matrix). Specify as an
%                 Ny-by-Nu matrix of nonnegative integers denoting the
%                 delays as multiples of sample time. Useful as a
%                 replacement for "nk" order - max(nk-1,0) lags can be
%                 factored out as "ioDelay" value.
%   'IntegrateNoise': Add integrator to noise channel. Logical vector of
%                 length Ny. Default: false(Ny,1). Setting IntegrateNoise
%                 to true (for a particular output) results in models of
%                 structure:
%                                                C(q)      
%                 A(q) y(t) =  B(q) u(t-nk) +  -------- e(t)
%                                              (1-q^-1)
%                 Use this property, for example, to create "ARIMA" models: 
%                 Estimate a 4th order ARIMA model for univariate
%                 time series data.
%                 load iddata9
%                 z9.y = cumsum(z9.y); % integrated data
%                 model = armax(z9, [4 1], 'IntegrateNoise', true); 
%                 compare(z9, model, 10) % 10-step ahead prediction
%
%  M = ARMAX(Z, [na nb nc nk], ..., OPTIONS)
%    specifies estimation options that configure the estimation objective,
%    initial conditions and numerical search method to be used for
%    estimation. Use the "armaxOptions" command to create the option set
%    OPTIONS.
%
%  M = ARMAX(Z, M0)
%  M = ARMAX(Z, M0, OPTIONS)
%    uses the IDPOLY model M0 to configure the initial parameterization of
%    the resulting model M. M0 must be a model of ARMAX structure (only A,
%    B and C polynomials must be active). M0 may be created using the
%    IDPOLY constructor or could be the result of a previous estimation.
%    The initial model argument, M0, may be followed by estimation options
%    to configure estimation options. If OPTIONS is not specified and M0
%    was created by estimation, the options are taken from M0.Report.OptionsUsed.
%
%  Continuous Time Model Estimation: This command cannot be used for
%  estimating continuous-time models. Some alternatives are to estimate a
%  continuous-time transfer function using TFEST command or a state-space
%  model using the SSEST command.
%
%   See also ARMAXOPTIONS, ARX, BJ, OE, POLYEST, SSEST, TFEST, IDPOLY,
%   IDDATA, IDPARAMETRIC/FORECAST.

%   Lennart Ljung 10-10-86
%   Copyright 1986-2011 The MathWorks, Inc.

narginchk(2,Inf)

% Set estimation data name.
I = find(cellfun(@(x)isa(x,'iddata') || isa(x,'frd'),varargin(1:2)));
if ~isempty(I) && isempty(varargin{I(1)}.Name);
   varargin{I(1)}.Name = inputname(I(1));
end

% Validate input arguments and create a template system if required.
try
   [sys, EstimData, Orders] = validatePEMInputs('armax',varargin{:});
catch E
   throw(E)
end

Options = getDefaultOptions(sys);
Disp = ~strcmpi(Options.Display,'off');
if Disp
   W = Options.ProgressWindow;
   Str = ctrlMsgUtils.message('Ident:estimation:msgDispPolyest1','ARMAX');
   idDisplayEstimationInfo('Intro',{Str, ' '},W);
end

%% Perform estimation.
try
   sys = pem_(sys, EstimData, Orders);
catch E
   if Disp
      S{1} = sprintf('<font color="red">%s</font>',E.message);
      S{2} = ctrlMsgUtils.message('Ident:estimation:msgAbortEstimation');
      idDisplayEstimationInfo('Error',S,W);
   end
   throw(E)
end
if Disp, W.STOP = true; end

%% Reconcile metadata between model and data.
sys = copyEstimationDataMetaData(sys, EstimData);
Report = sys.Report; Report.Method = 'ARMAX'; sys = setReport(sys, Report);
