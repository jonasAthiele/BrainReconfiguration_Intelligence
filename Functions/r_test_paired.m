function [p, T2, df, r_jk, r_jh, r_kh]=r_test_paired(j,k,h,tail)
% r_test_paired() - Tests the null hypothesis that two variables (h and k)
%                   are equally correlated with a third variable (j) when 
%                   all three sets of observations were derived from the 
%                   same "individuals" (i.e., a repeated measures/paired 
%                   observations design).
%           
%
% Usage:
%  >> [p, T2, df, r_jk, r_jh, r_kh]=r_test_paired(j,k,h,tail);
%
% Required Inputs:
%   j    - A vector of values of the "dependent variable" 
%   k    - A vector of values of an "independent variable" 
%   h    - A vector of values of another "independent variable" 
%
% Optional Input:
%   tail - The tail of the test.  There are three possible values.
%            0: the alternative hypothesis is that the correlation between
%               variables j and k DOES NOT EQUAL the correlation between
%               variables j and h (i.e., a two tailed test). {default}
%            1: the alternative hypothesis is that the correlation between
%               variables j and k is GREATER THAN the correlation between
%               variables j and h (i.e., an upper tailed test).
%           -1: the alternative hypothesis is that the correlation between
%               variables j and k is LESS THAN the correlation between
%               variables j and h (i.e., a lower tailed test).
%
%
% Outputs:
%   p    - the p-value of the test statistic
%   T2   - the test statisic, a t-score
%   df   - the degrees of freedom of the t distribution (number of
%          "individuals"-3)
%   r_jk - Pearson's linear correlation coefficient between variables j and k
%   r_jh - Pearson's linear correlation coefficient between variables j and h
%   r_kh - Pearson's linear correlation coefficient between variables k and h
%
% 
% Example:
% >> k=rand(1,25);
% >> h=k+rand(1,25)*2;
% >> j=4*k+randn(1,25);
% >> [p, T2, df, r_jk, r_jh, r_kh]=r_test_paired(j,k,h,0);
%
%
% References:
% This test was taken from:
%   Williams, E. J. (1959) The comparison of regression variables. Journal
%   of the Royal Statistical Society: Series B, 21, 396-399.
%
% as recommended by:
%   Steiger, J. H. (1980) Tests for comparing elements of a correlation
%   matrix. Psychological Bulletin. 87(2), 245-251.
%
%
% Author: 
% David M. Groppe (November, 2009)
% Kutaslab
% Department of Cognitive Science
% University of California, San Diego
% La Jolla, CA, USA


%% Error check inputs 
if nargin<3,
   help r_test_paired;
   error('r_test_paired requires at least three arguments.');
elseif nargin==3,
   tail=0; %default: two-tailed test 
elseif (tail~=0) && (tail~=-1) && (tail~=1),
   error('Argument "tail" needs to be set to -1, 0, or 1.'); 
end

% Make sure all data vectors are indeed vectors
if ~isvector(j) || ~isvector(k) || ~isvector(h),
   error('Input variables j, k, and h all need to be vectors.'); 
end


% Make sure all data vectors are column vectors
sj=size(j);
sk=size(k);
sh=size(h);
if sj(1)==1,
   j=j'; 
   n=sj(2); % sample size
else
   n=sj(1); % sample size 
end

if sk(1)==1,
   k=k'; 
end

if sh(1)==1,
   h=h'; 
end

if sum(sj-sk) || sum(sj-sh) || sum(sk-sh)
   error('Vectors j, k, and h need to have the same number of elements.');
end


%% Test

% Compute correlations
r_jk=corr(j,k);
r_jh=corr(j,h);
r_kh=corr(h,k);

% Compute test statistic
detR=1-(r_jk^2)-(r_jh^2)-(r_kh^2)+2*r_jk*r_jh*r_kh; %determinant of correlation matrix
r_tilda=(r_jk+r_jh)/2;

T2=(r_jk-r_jh)*sqrt( (n-1)*(1+r_kh)/(2*(n-1)*detR/(n-3) + r_tilda*((1-r_kh)^3)) );

% Compute p-value
df=n-3;
if tail==0,
    p=2*(1-cdf('t',abs(T2),df)); %two tailed test
elseif tail==1,
    p=1-cdf('t',T2,df); %upper tailed test
else
    p=1-cdf('t',-T2,df); %lower tailed test
end