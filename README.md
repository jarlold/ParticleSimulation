# Particle Toy
This simulates a lot of particles. They interact with eachother according to Newtonian physics. This is O(N^2), which is pretty tricky, so it makes estimates for far away particles using "spatial partioning". 

Oh except neighboring cells aren't counted so it's possible for two particles to be next to eachother but not interact.
