function y=cdom_model_noK(x,a,s,lam0)

% Model for CDOM spectral slope


y=a.*exp(-(x-lam0).*s);
