function [xdm] = FDMDemux(muxSignal,t,MagSpec,freqshift,PhaseSpec)
    
    xdm = zeros(length(freqshift),length(t));
    
    for i = 1:length(freqshift)
        xdm(i,:) = muxSignal .* cos(2 * pi * freqshift(i) * t + PhaseSpec(i)) * MagSpec(i); 
    end

end