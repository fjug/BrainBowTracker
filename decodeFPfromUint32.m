function [ fp ] = decodeFPfromUint32( encodedFP )
%decodeFPfromUint32 gets an encodedFP and decodes it.
%   byte 1 (most significant): volume (number of voxels)
%   byte 2: average red tone
%   byte 3: average green tone
%   byte 4: average blue tone
%   Assumptions: volume always < 255 voxels; colors < 65535 (got rescaled
%   to [0-255]!

muh = bitget(encodedFP,32:-1:25);
volume = sum(bitset(0,8:-1:1,muh,'int8'));

muh = bitget(encodedFP,24:-1:17);
red = sum(bitset(0,8:-1:1,muh,'int8'));

muh = bitget(encodedFP,16:-1:9);
green = sum(bitset(0,8:-1:1,muh,'int8'));

muh = bitget(encodedFP,8:-1:1);
blue = sum(bitset(0,8:-1:1,muh,'int8'));

fp = [ volume, red, green, blue ];

end