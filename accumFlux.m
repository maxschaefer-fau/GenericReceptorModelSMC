% Function to accumulate fluxes 
% 
% i:    accumulated flux in time step k 
% T:    time step 
% ibar: value of ibar (o or d) in time step k 
% ikm1: accumulated flux (o or d) in time step k-1

function [i] = accumFlux(ikm1, ibar, T) 

    i = ikm1 + T*ibar;

end