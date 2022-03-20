% Receptor circuit 
%
%
% iokp1:    flux io_bar in time step k+1
% idkp1:    flux id_bar in time step k+1
% k:        receptor state transition rates
% csatk:    receptor saturation in time step k 
% iok:      flux io in time step k
% iod:      flux id in time step k


function [iokp1, idkp1] = receptorCircuit(csatk, iok, idk, ...
    kco, koc, kod, kdo, kcd, kdc)

    iokp1 = kco*csatk - koc*iok - kod*iok + kdo*idk; 
    idkp1 = kcd*csatk - kdc*idk - kdo*idk + kod*iok;
    
end