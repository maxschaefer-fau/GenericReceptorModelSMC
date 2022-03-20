% Function for saturation of receptors 
%
% csat:     receptor saturation in time step k 
% C:        number of available receptors 
% io:       accumulated flux i_o in time step k 
% id:       accumulated flux i_d in time step k 

function [csat, cs] = saturate(cak, C, io, id) 

    cs = (1 - 1/C*(io + id));
    csat = cs*cak;

end