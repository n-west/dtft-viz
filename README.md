dtft
===

visualization (targeted at undergrad jr/sr level) of DTFT, z-transform, and digital filter design with transform methods

The basic idea is:
a) Provide a transfer function in either numerator/denominator or
pole/zero format. Plot the 2d pole-zero plot
b) Plot a 3d DTFT with a translucent cylinder showing the unit circle
c) Cut out just the unit circle in a 3d plot
d) unravel the unit circle, which becomes the 2-d z-transform. This
shows how poles/zeros affect the look of the transfer function.

e) Show (what I think is called) bilinear transform to see how old
frequencies get mapped to new frequencies when designing a
high/band/low-pass filter from a butterworth/reference design.
f) Plot pzmap, DTFT and z-transform of corresponding filter.
