function y=cdom_model(x,a,s,lam0,K)

% Model for CDOM spectral slope


y=a.*exp(-(x-lam0).*s)+K;