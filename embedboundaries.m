function Image = embedboundaries(Input, Boundaries, channels)
%Embeds the boundaries in the binary image "Boundaries" into the color image
%"Input" for visualization. "channels" is used to select the boundary color.

%determine datatype
if(isa(Input, 'uint8'))
    Maximum = 255;
elseif(isa(Input,'double') || isa(Input,'single'))
    Maximum = 1.0;
elseif(isa(Input, 'uint16'))
    Maximum = 65535;
end

%extract channels
if(size(Input,3) == 1)
    Red = Input; Green = Input; Blue = Input;
else
    Red = Input(:,:,1); Green = Input(:,:,2); Blue = Input(:,:,3);
end

%embed boundaries
Red(Boundaries) = Maximum * channels(1);
Green(Boundaries) = Maximum * channels(2);
Blue(Boundaries) = Maximum * channels(3);

%generate output
Image = cat(3, Red, Green, Blue);
