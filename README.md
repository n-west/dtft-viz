dtft
===

visualization (targeted at undergrad jr/sr level) of DTFT, z-transform, and digital filter design with transform methods

The basic idea is:
1) Provide a transfer function in either numerator/denominator or
pole/zero format. Plot the 2d pole-zero plot
2) Plot a 3d z-transform with a translucent cylinder showing the unit circle
3) Cut out just the unit circle in a 3d plot
4) unravel the unit circle, which becomes the 2-d DTFT. This
shows how poles/zeros affect the look of the transfer function.

5) Show (what I think is called) bilinear transform to see how old
frequencies get mapped to new frequencies when designing a
high/band/low-pass filter from a butterworth/reference design.
6) Plot pzmap, DTFT and z-transform of corresponding filter.
