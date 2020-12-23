function [model,stats,residuals,cint,yhat]=cdom_spectral_slope(lambda,cdom_absorption,lambda_start,lambda_stop,lam0,K)

% cdom_spectral_slope
%
% [model,stats,residuals,cint,yhat]=cdom_spectral_slope(lambda,cdom_absorption,lambda_start,lambda_stop,lam0,K)
%
% Function uses measured CDOM absorption to fit an exponential function 
% and calculate the spectral slope, or shape, of the exponential curve 
% over the specificed wavelength range

% Function uses Matlab's curve fitting toolbox and the fit function for 
% model fitting. Model contains absorption at reference wavelength (lam0)
% and S value for spectral slope over designated range. These parameters as
% accessible using dot notation (e.g., model.s = spectral slope)
%
%
% Inputs:
%
%     lambda             = wavelength values affiliated with absorption 
%                          spectra, format = vector
%     absorption         = CDOM absorption, format = vector
%     lambda_start       = lower wavelength range for spectral slope 
%                          calculation (e.g. 275 nm)
%     lambda_stop        = maximum wavelength range for spectral slope 
%                          calculation (e.g. 275 nm)
%     lam0               = initial wavelength used to calculate exponential
%                          model, used to reference absorption at that
%                          wavelength
%     K                  = offset term, if used (not required)
%
% Returns:
%
%     model              = 
%     stats              = goodness of fit statistics
%     residuals          = residuals, iterations used for model fit, etc.
%     cint               = confidence intervals
%     yhat               = predicted values from model
%
%%%%%%%
%
% Suggested reading
% Twardowski, Michael S., et al. "Modeling the spectral shape of absorption 
% by chromophoric dissolved organic matter." Marine Chemistry 89.1-4 
% (2004): 69-88.
%
% Helms, John R., et al. "Absorption spectral slopes and slope ratios as 
% indicators of molecular weight, source, and photobleaching of 
% chromophoric dissolved organic matter." Limnology and Oceanography 53.3 
% (2008): 955-969.
%
%
% copyright (c) 2020 Brice K. Grunert
% email: bricegrunert@gmail.com
%
%Last modified on 23 December 2020 by BG
%
%%%%%%%


% find wavelengths corresponding to minimum and maximum of spectral range
% being considered

for jj=1:length(lambda)
    if lambda(jj)==lambda_start
        ind=jj;
    end
    if lambda(jj)==lambda_stop
        ind2=jj;
    end
end

% reorient wavelength to column format, if necessary
x=lambda(ind:ind2);
if size(x,2) > 1
    x=x';
end

% reorient cdom absorption to column format, if necessary
y=cdom_absorption(ind:ind2);
if size(y,2) > 1
    y=y';
end

ind=find(lambda==lam0);

% check whether to include K in CDOM model
if exist('K','var')==0
    ft=fittype('cdom_model_noK(x,a,s,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s'},'problem',{'lam0'}); % 'problem',{'variable'} specifices as constant
    fopts=fitoptions(ft);
    fopts.StartPoint=[cdom_absorption(ind) 0.015];
    
    [model,stats,residuals]=fit(x,y,ft,fopts,'problem',{lam0}); % define problem variables as fixed; coefficients should not be fixed
else
    ft=fittype('cdom_model(x,a,s,lam0,K)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s'},'problem',{'lam0','K'});
    fopts=fitoptions(ft);
    fopts.StartPoint=[cdom_absorption(ind) 0.015];
    
    [model,stats,residuals]=fit(x,y,ft,fopts,'problem',{lam0,K});     
end

% Confidence intervals
cint=confint(model);
% Predicted values
yhat=feval(model,x);
