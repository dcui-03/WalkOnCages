# WalkOnCages
C++ implementation of Stochastic Computation of Barycentric Coordinates [de Goes and Desbrun, 2024]


# Update 08/17

Figured out smoothing issue (smooth u, not alpha!) and implemented MVC and positive MVC. Just need to do some testing with these. I also uniformifed all the different projection structs so they now all use the Utils version and variants.

Cage samplers also live in the cage classes, NOT the original classes.

# Next TODO's

- Find a WoS speedup method. Two problems currently: wasted walks, and slow walks. The second I can try to update by optimizing the code, but the first is a little difficult.