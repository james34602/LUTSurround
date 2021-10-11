# LUTSurround
As its name suggests, this is a real time stereo to multichannel upmix algorithm that based on lookup table.

Processing speed is very fast, this is generally good for low power consumption device.

## What's upmix? And what's the motivation to deploy such technique?
Upmix is a technique of surround sound generation that uses the correlation between stereo Left, Right channel.

Most audio material in the world is mainly exist in Stereo format, it's de facto standard, as manufactorer see the advantage that Stereo recording encode spatial information good enough, entropy encoding stereo of the signal can be benefit by Stereo Mid/Side, Stereo also provide good quality of sound reproduction.

However, multichannel configuration exists, convert exist Stereo material to multichannel material have it's value.

## Why most audio material is still recorded in mono/stereo while multichannel format exist?
Reason is simple, transmit large amount of data is still a problem nowaday, most user rather need better video quality than more 3D-ish surround sound.

## How it works?
Short time Fourier transform -> Map Magnitude, Phase difference to cartesian plane with LUT -> Map cartesian data to mask function with LUT -> Filtering -> Inverse short time Fourier transform of N channels

## What sample rate is supported and what's the latency?
Consistent latency over 8000Hz to 384000Hz sample rate with intrinsic delay less than 17.5ms