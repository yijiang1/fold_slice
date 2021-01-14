function [wav,wavA] = electronwavelength(kev)
            %wavelength in a0
            a0 = 0.52917720859;
            wav = 12.3986./sqrt((2*511.0+kev).*kev)/a0;
            wavA = wav*a0;
end 
