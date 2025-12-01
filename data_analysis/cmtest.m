function [H, pValue, CvMstat, criticalValue] = cmtest(x, varargin)
%CMTEST Single sample Cramer-Von Mises goodness-of-fit hypothesis test.
%   H = CMTEST(X) performs a Cramer-Von Mises (CvM) test to determine if
%   a random sample X could have come from a specified distribution.
%   This is a 2-sided test.
%   H indicates the result of the hypothesis test:
%      H = 0 => Do not reject the null hypothesis at the 5% significance
%      level. 
%      H = 1 => Reject the null hypothesis at the 5% significance
%      level.
%
%   X is a vector representing a random sample from some underlying
%   distribution, with cumulative distribution function F. Missing 
%   observations in X, indicated by NaNs (Not-a-Number), are ignored.
%
%   [H,P] = CMTEST(...) also returns the asymptotic P-value P.
%
%   [H,P,CvMSTAT] = CMTEST(...) also returns the CvM test statistic CvMSTAT
%   defined above.
%
%   [H,P,CvMSTAT,CV] = CMTEST(...) returns the critical value of the test CV.
%
%   [...] = CMTEST(X,'PARAM1',val1,'PARAM2',val2,...) specifies one or
%   more of the following parameter name/value pairs to specify whether to
%   perform a simple or composite hypothesis test, to specify the
%   distribution being tested for, to control for the significance level,
%   and to specify whether to perform the test using Monte-Carlo
%   simulations:
%
%   Parameter       Value
%
%   'alpha'         A value ALPHA between 0 and 1 specifying the
%                   significance level. Default is 0.05 for 5% significance.
%
%   'CDF'           CDF is the c.d.f. under the null hypothesis.  It can
%                   be specified either as a ProbabilityDistribution object
%                   or as a two-column matrix. CDF must be completely
%                   specified. If missing, the default is the standard
%                   normal with sample mean and variance.
%
%   'MCTol'         Monte-Carlo tolerance value. In this case, an
%                   approximation for the P value is computed directly,
%                   using Monte-Carlo simulations. Default = [].
%
%   Let S(X) be the empirical c.d.f. estimated from the sample vector X, F(X)
%   be the corresponding true (but unknown) population c.d.f., and CDF be the
%   known input c.d.f. specified under the null hypothesis.
%   Test Statistic: W2 = integral from -Inf to Inf of n*(S(x)-F(x))^2 dF(x)
%
%   In the matrix version of CDF, column 1 contains the x-axis data and
%   column 2 the corresponding y-axis c.d.f data. Since the CvM test
%   statistic will occur at one of the observations in X, the calculation
%   is most efficient when CDF is only specified at the observations in X.
%   When column 1 of CDF represents x-axis points independent of X, CDF is
%   're-sampled' at the observations found in the vector X via
%   interpolation. In this case, the interval along the x-axis (the column
%   1 spread of CDF) must span the observations in X for successful
%   interpolation.
%
%   The decision to reject the null hypothesis is based on comparing the
%   statistic CVMSTAT with the critical value CV, not by comparing the
%   p-value P with ALPHA. CV is computed separately by interpolation in a
%   table. CV is returned as NaN if ALPHA is outside this range.
%
%   See also UTEST, KSTEST, ADTEST.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Copyright (c) 20 March 2015 by Ahmed Ben Saïda           %
%                 LaREMFiQ Laboratory, IHEC Sousse - Tunisia             %
%                       Email: ahmedbensaida@yahoo.com                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% References:
% S. Csorgo & J. J. Faraway, 1996, "The Exact and Asymptotic
%   Distributions of Cramer-von Mises Statistics", Journal of the Royal
%   Statistical Society: Series B, Vol. 58, No. 1, pp. 221-234.
%
% D'Agostino and Stephens, Goodness-Of-Fit Techniques, Marcel-Dekker, New
%   York, 1986. 
%
% M.A Stephens, 1974, EDF Statistics for Goodness of Fit and Some
%   Comparisons, Journal of the American Statistical Association.
%
% M. A. Stephens, 1970, "Use of the Kolomogorov-Smirnov, Cramer-Von Mises
%   and Related Statistics Without Extensive Tables", Journal of the Royal
%   Statistical Society: Series B, Vol. 32, No. 1, pp. 115-122.
%
% W. F. Scott, 2000, "Tables of the cramér-von mises distributions",
%   Communications in Statistics - Theory and Methods, Vol. 29, No. 1,
%   pp. 227-235.

% Parse arguments and check if parameter/value pairs are valid
paramNames = {'Alpha', 'CDF', 'MCTol'};
dflts  =     {  0.05,   [],     [],  };

[valAlpha, CDF, valMCTol] = parseArgs(paramNames, dflts, varargin{:});

% Ensure the sample data is a real vector.
if ~isvector(x) || ~isreal(x)
    error('Sample data X must be a real vector.');
end

% Ensure the significance level, ALPHA, is a scalar between 0 and 1.
if ~isscalar(valAlpha) || ~(valAlpha > 0 && valAlpha < 1)
    error('Significance level ALPHA must be a scalar between 0 and 1.');
else
    alpha = valAlpha;
end

% Ensure the Monte-Carlo tolerance, MCTOL, is a numeric scalar greater than
% 0.
if ~isempty(valMCTol) && (~isscalar(valMCTol) || valMCTol <=0)
    % Invalid Monte-Carlo tolerance
    error('Monte-Carlo standard error tolerance MCTOL must be a positive scalar value.');
else
    mctol = valMCTol;
end

% Remove missing values.
x = x(~isnan(x));

% Compute sample size n.
n = length(x);

% Check & scrub the hypothesized CDF specified under the null hypothesis.
% If CDF has been specified, remove any rows with NaN's in them and sort
% x-axis data found in column 1 of CDF. If CDF has not been specified, then
% allow the convenience of x ~ N(0,1) under the null hypothesis.
if (isa(CDF,'ProbDist') || isa(CDF,'prob.ProbabilityDistribution'))
   xCDF = x;
   yCDF = cdf(CDF,x);
elseif ~isempty(CDF)
    if size(CDF,2) ~= 2
      error('Hypothesized CDF matrix must have 2 columns.');
    end

    CDF  =  CDF(~isnan(sum(CDF,2)),:);

    if size(CDF,1) == 0
      error('Hypothesized CDF matrix must have at least 1 valid row.');
    end

    [xCDF,i] =  sort(CDF(:,1));    % Sort the theoretical CDF.
    yCDF = CDF(i,2);

    ydiff = diff(yCDF);
    if any(ydiff<0)
      error('CDF must define an increasing function of X.');
    end

    % Remove duplicates, but it's an error if they are not consistent
    dups = find(diff(xCDF) == 0);
    if ~isempty(dups)
      if ~all(ydiff(dups) == 0)
         error('CDF must not have duplicate X values.');
      end
      xCDF(dups) = [];
      yCDF(dups) = [];
    end
    
else
    xCDF = x;
    yCDF = normcdf(x,mean(x),std(x));
end

% If CDF's x-axis values have been specified by the data sample X, then just
% assign column 2 to the null CDF; if not, then we interpolate subject to the
% check that the x-axis interval of CDF must bound the observations in X. Note
% that both X and CDF have been sorted and have had NaN's removed.
if isequal(x,xCDF)
   nullCDF  =  yCDF;       % CDF has been specified at the observations in X.
else
   if (x(1) < xCDF(1)) || (x(end) > xCDF(end))
     error('Hypothesized CDF matrix must span the observations interval in X.');
   end
   nullCDF  =  interp1(xCDF, yCDF, x, 'linear');
end

% Compute the Cramer Von-Mises statistic.
nullCDF = sort(reshape(nullCDF,n,1));
CvMstat = ComputeCvMstat(nullCDF,n);

if isempty(mctol)
    
    % Find critival value.
    [alphas, CVs, sampleSizes] = findCriticalValues;
    
    % Make sure alpha is within the lookup table
    validateAlpha(alpha, alphas);
    
    % The funtion findCriticalValues gives the lower tail, so transform the
    % alphas before computing the critical value.
    alphas = 1 - alphas;
    
    [OneOverSampleSizes, LogAlphas] = meshgrid(1./sampleSizes, log(alphas));
    criticalValue = interp2(OneOverSampleSizes, LogAlphas, CVs', 1/n, log(alpha));
    
    % Compute the P-value based on the Cramer-von Mises statistic
    % P-value = Pr(Wn >= d) = 1-Pr(Wn < d).
    % Both of the following methods works.
%     pValue = CvMpValue(CvMstat, n);
    pValue = inverse_interp2(OneOverSampleSizes, LogAlphas, CVs', CvMstat, 1/n);
    
else
    % Compute the critical value and p-value on the fly using Monte-Carlo simulation.
    [criticalValue, pValue] = cmtestMC(CvMstat, n, alpha, 'unif', mctol);
end

% "H = 0" implies that we "Do not reject the null hypothesis at the
% significance level of alpha," and "H = 1" implies that we "Reject null
% hypothesis at significance level of alpha."
% H  =  (pValue < alpha);
H = (CvMstat > criticalValue); % this is more sure than the asymptotic p-value.

%----------------Subfunctions--------------------------------------------%
function varargout = parseArgs(paramNames, dflts, varargin)

varargout = dflts;
for i = 1:length(dflts)
    varargout{i} = dflts{i};
end

n = length(varargin);

% Check the input argument list:
if rem(n,2) ~= 0
    error('Wrong input variables.')
else
    parameters = varargin(1:2:n-1);
    values = varargin(2:2:n);
end

for j = 1:length(parameters)
    % Ensure that all parameter names are strings:
    if ~ischar(parameters{j})
       error('Parameter name at input argument number %s must be a string.', num2str(2*j - 1))
    end
    % Ensure that all parameter names represent valid field names:
    matches = find(strncmpi(parameters{j},paramNames,length(parameters{j})));
    if isempty(matches)
       % No matches found:
       error('Unrecognized parameter name ''%s''.', parameters{j})
    elseif length(matches) > 1
       % More than one match:
       exacts = find(strcmpi(parameters{j},paramNames));
       if length(exacts) == 1
        % An exact match:
          matches = exacts;
       else
          % The parameter name is ambiguous; additional information needed:
          error('Ambiguous parameter name ''%s''. Additional information needed.', parameters{j})
       end
    end
    varargout{matches} = values{j};
end

%------------------------------------------
function CvMstat = ComputeCvMstat(z,n)
% Compute Cramer Von-Mises Statistic
% Sort the data and compute the statistic
% z = reshape(z,n,1);
% z = sort(z);
w = (2*(1:n)' - 1) ./ (2*n);
CvMstat = 1/(12*n) + sum((w - z).^2);

% ------------------------------------------
function cm = CvMpValue(z,n)
% Evaluate the Cramer-von Mises asymptotic p-value.

% Transform the CvM statistic using Stephens (1970) approximation
% for known F(x).
z = (z - 0.4/n + 0.6/n^2) * (1 + 1/n);

% Compute the upper-tail p-value (Stephen, 1970).
cm = 0.05 * exp(2.79 - 6 * z);

%
% This is the analytical form of the p-value. However, near infinity it
% generates a warning and give a complexe number, (Scott, 2000).
%

% fun = @(k) (-1).^(k-1) * 2/pi .* integral(@(r) exp(-r.^2*z/2)./sqrt(-r.^2.*sin(r.^2)), 2*(k-1)*pi, 2*k*pi);
% k   = 1:1000; % k tends to infinity.
% cm   = sum(arrayfun(fun,k));

% ------------------------------------------
function Yq = inverse_interp2(X, Y, V, Vq, Xq)
% Inverse 2D interpolation to find the p-value for a given critical value.

f  = @(y) abs(Vq - interp2(X, Y, V, Xq, y));

% Since p-value is always between 0 and 1 (0 <= x <= 1), use the logistic
% function to bound it.
y0  = log(0.5); % initial guess, 0 <= x <= 1.
Yq  = fminsearch(f, y0);
Yq  = exp(Yq);

Yq(Yq > 1) = 1; % Avoid constraint violation.

% ------------------------------------------
function validateAlpha(alpha, alphas)
% Make sure alpha is within the lookup table
if alpha< alphas(1) || alpha>alphas(end)
    error(message('Significance level ALPHA is outside range of lookup table.\nUse a value in the interval [%g, %g], or use the MCTOL input argument.', sprintf('%g', alphas(1)),...
        sprintf('%g', alphas(end))));
end

% ------------------------------------------
function [alphas, CVs, sampleSizes] = findCriticalValues
% Find rows of critical values at relevant significance levels from Csorgo
% & Faraway (1996). The following table gives the lower tail P[Wn<=x] = p

alphas = [0.01	0.025	0.05	0.1	0.15	0.20	0.25	0.50	0.75	0.80	0.85	0.90	0.95	0.975	0.99	0.999];

%0.01	0.025	0.05	0.1     0.15	0.20	0.25	0.50	0.75	0.80	0.85	0.90	0.95	0.975	0.99	0.999 %alphas
CVs = [
0.04326	0.04565	0.04962	0.05758	0.06554	0.07350	0.08146	0.12659	0.21522	0.24743	0.28854	0.34343	0.42480	0.48901	0.55058	0.62858;... %n=2
0.03319	0.03774	0.04360	0.05289	0.06091	0.06887	0.07683	0.12542	0.21338	0.24167	0.27960	0.33785	0.43939	0.53318	0.63980	0.82240;... %n=3
0.03002	0.03536	0.04149	0.05093	0.05896	0.06681	0.07493	0.12406	0.21171	0.24260	0.28336	0.34184	0.44206	0.54200	0.67017	0.92970;... %n=4
0.02869	0.03422	0.04036	0.04969	0.05800	0.06610	0.07427	0.12250	0.21164	0.24237	0.28305	0.34238	0.44697	0.55056	0.68352	0.98730;... %n=5
0.02796	0.03344	0.03959	0.04911	0.05747	0.06548	0.07351	0.12200	0.21110	0.24198	0.28331	0.34352	0.44911	0.55572	0.69443	1.02000;... %n=6
0.02741	0.03292	0.03914	0.04869	0.05698	0.06492	0.07297	0.12158	0.21087	0.24197	0.28345	0.34397	0.45100	0.55935	0.70154	1.04250;... %n=7
0.02702	0.03257	0.03880	0.04830	0.05660	0.06460	0.07260	0.12120	0.21070	0.24190	0.28350	0.34450	0.45240	0.56220	0.70720	1.05910;... %n=8
0.02680	0.03230	0.03860	0.04810	0.05630	0.06430	0.07240	0.12100	0.21060	0.24180	0.28360	0.34480	0.45340	0.56430	0.71150	1.07220;... %n=9
0.02650	0.03212	0.03840	0.04780	0.05610	0.06410	0.07210	0.12070	0.21040	0.24170	0.28360	0.34500	0.45420	0.56590	0.71470	1.08220;... %n=10
0.02564	0.03120	0.03742	0.04689	0.05515	0.06312	0.07117	0.11978	0.20989	0.24148	0.28384	0.34617	0.45778	0.57331	0.72895	1.11898;... %n=20
0.02512	0.03068	0.03690	0.04636	0.05462	0.06258	0.07062	0.11924	0.20958	0.24132	0.28396	0.34682	0.45986	0.57754	0.73728	1.14507;... %n=50
0.02488	0.03043	0.03665	0.04610	0.05435	0.06231	0.07035	0.11897	0.20943	0.24125	0.28402	0.34715	0.46091	0.57968	0.74149	1.15783;... %n=200
0.02481	0.03037	0.03658	0.04603	0.05428	0.06224	0.07027	0.11889	0.20938	0.24123	0.28403	0.34724	0.46119	0.58026	0.74262	1.16120;... %n=1000
0.02480	0.03035	0.03656	0.04601	0.05426	0.06222	0.07026	0.11888	0.20939	0.24124	0.28406	0.34730	0.46136	0.58061	0.74346	1.16204];   %n=Inf

sampleSizes = [2 3 4 5 6 7 8 9 10 20 50 200 1000 Inf];

%------------------------------------------
function [crit, p] = cmtestMC(CvMstat, n, alpha, distr, mctol)
%CMTESTMC Simulated critical values and p-values for Cramer-von Mises test.
%   [CRIT,P] = CMTESTMC(CvMSTAT,N,ALPHA,DISTR,MCTOL) returns the critical
%   value CRIT and p-value P for the Cramer-von Mises test of the null
%   hypothesis that data were drawn from a distribution in the family
%   DISTR, for a sample size N and confidence level 100*(1-ALPHA)%.  P is
%   the p-value for the observed value CvMSTAT of the Cramer-von Mises
%   statistic.  DISTR is 'norm', 'exp', 'ev' 'or 'unif'. ALPHA is a scalar
%   or vector.  CMTESTMC uses Monte-Carlo simulation to approximate CRIT
%   and P, and chooses the number of MC replications, MCREPS, large enough
%   to make the standard error for P, SQRT(P*(1-P)/MCREPS), less than
%   MCTOL.

vartol = mctol^2;
crit = 0;
p = 0;
mcRepsTot = 0;
mcRepsMin = 1000;

while true
    mcRepsOld = mcRepsTot;
    mcReps = ceil(mcRepsMin - mcRepsOld);
    CvMstatMC = zeros(mcReps,1);
    
    switch distr
        % Simulate critical values for the normal
        case 'normal'
            mu0 = 0;
            sigma0 = 1;
            for rep = 1:length(CvMstatMC)
                x = normrnd(mu0, sigma0,n,1);
                xCDF = sort(x);
                nullCDF = normcdf(xCDF, mean(x), std(x));
                CvMstatMC(rep) = ComputeCvMstat(nullCDF, n);
%                 w = (2*(1:n)' - 1) ./ (2*n);
%                 CvMstatMC(rep) = 1/(12*n) + sum((w - nullCDF).^2);
            end
        case 'exponential'
            beta0 = 1;
            for rep = 1:length(CvMstatMC)
                x = exprnd(beta0,n,1);
                xCDF = sort(x);
                nullCDF = expcdf(xCDF, mean(x));
                CvMstatMC(rep) = ComputeCvMstat(nullCDF, n);
%                 w = (2*(1:n)' - 1) ./ (2*n);
%                 CvMstatMC(rep) = 1/(12*n) + sum((w - nullCDF).^2);
            end
        case {'ev', 'extreme value'}
            mu0 = 0;
            sigma0 = 1;
            for rep = 1:length(CvMstatMC)
                x = evrnd(mu0,sigma0,n,1);
                pHat = evfit(x); %MLE Estimate
                xCDF = sort(x);
                nullCDF = evcdf(xCDF, pHat(1), pHat(2));
                CvMstatMC(rep) = ComputeCvMstat(nullCDF, n);
%                 w = (2*(1:n)' - 1) ./ (2*n);
%                 CvMstatMC(rep) = 1/(12*n) + sum((w - nullCDF).^2);
            end
        case 'unif'
            for rep = 1:length(CvMstatMC)
                nullCDF = sort(rand(n,1));
                CvMstatMC(rep) = ComputeCvMstat(nullCDF, n);
%                 w = (2*(1:n)' - 1) ./ (2*n);
%                 CvMstatMC(rep) = 1/(12*n) + sum((w - nullCDF).^2);
            end
    end
    
    critMC = prctile(CvMstatMC, 100*(1-alpha));
    pMC = sum(CvMstatMC > CvMstat)./mcReps;
    
    mcRepsTot = mcRepsOld + mcReps;
    crit = (mcRepsOld*crit + mcReps*critMC) / mcRepsTot;
    p = (mcRepsOld*p + mcReps*pMC) / mcRepsTot;
    
    % Compute a std err for p, with lower bound (1/N)*(1-1/N)/N when p==0.
    sepsq = max(p*(1-p)/mcRepsTot, 1/mcRepsTot^2);
    if sepsq < vartol
        break
    end
    % Based on the current estimate, find the number of trials needed to
    % make the MC std err less than the specified tolerance.
    mcRepsMin = 1.2 * (mcRepsTot*sepsq)/vartol;
end