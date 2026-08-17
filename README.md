# StochasticBC
C++ implementation of Stochastic Barycentric Coordinates [de Goes and Desbrun, 2024]


# Update 08/16

Completed stochastic harmonic coordinates and updated the visualizer to handle it. As part of this, I completed the closest point code for curvenets. This is a bit rough, though.

# Next TODO's

- Add a random sampler for curvenets and polynets. 
- Add raycasting function to the mesh class (return ALL positive intersections in order).
- Add mean value coordinates
- Add positive mean value coordinates
- Find a WoS speedup method. Two problems currently: wasted walks, and slow walks. The second I can try to update by optimizing the code, but the first is a little difficult.