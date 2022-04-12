function [kco, koc, kod, kdo, kcd, kdc] = rateVariation(kco_, koc_, kod_, kdo_, ...
    kcd_, kdc_, sc)

    kco = kco_*sc(1); 
    koc = koc_*sc(2); 
    kod = kod_*sc(3); 
    kdo = kdo_*sc(4); 
    kcd = kcd_*sc(5); 
    kdc = kdc_*sc(6); 
end