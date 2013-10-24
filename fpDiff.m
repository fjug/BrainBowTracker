function [ fpdiff ] = fpDiff( enc1, enc2 )
%fpDiff get difference between to encoded FPs
muh1 = decodeFPfromUint32(enc1);
muh2 = decodeFPfromUint32(enc2);
fpdiff = muh1-muh2;
end

