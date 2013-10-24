function [ integer ] = encodeFPinUint32( fingerprint )
%encodeFPinUint32 gets an fingerprint and encodes it in an UInt32.
%   byte 1 (most significant): volume (number of voxels)
%   byte 2: average red tone
%   byte 3: average green tone
%   byte 4: average blue tone
%   Assumptions: volume always < 255 voxels; colors < 65535 (get rescaled
%   to [0-255]!
integer = bitshift( uint32(fingerprint(1)),     24 ) +...
          bitshift( uint32(fingerprint(2)/257), 16 ) +...
          bitshift( uint32(fingerprint(3)/257),  8 ) +...
          bitshift( uint32(fingerprint(4)/257),  0 ) ;
end