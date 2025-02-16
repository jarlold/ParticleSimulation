use rand::Rng;
use kiss3d::window::Window;
extern crate nalgebra as na;
use std::collections::VecDeque;

pub struct Particle {
    pub pos: na::Point<f32, 3>,
    pub vel: na::Point<f32, 3>,
    pub inverse_mass: f32,
    pub color: na::Point<f32, 3>,
}

pub struct ParticleSim {
    pub particles: Vec<Particle>,
    pub time_step: f32,
    pub num_cells_per_axis: usize,
    pub partition_edge_len: f32,
}

impl ParticleSim {
    pub fn tick(&mut self) {
        // We'll have to split our particles up into "spatial partitions" or
        // computing the forces will be O(n^2) and my computer he will cry
        let partitions = self.partition_particles();
        let mut centroids : Vec<na::Vector3<f32>> = Vec::with_capacity(partitions.len());
        let mut masses : Vec<f32> = Vec::with_capacity(partitions.len());

        // We'll go through and compute the centroid for each partition
        for i in 0 .. partitions.len() {
            centroids.push( na::Vector3::new(0.0, 0.0, 0.0));
            masses.push(0.0);
            for j in 0 .. partitions[i].len() {
                centroids[i] += self.particles[partitions[i][j]].pos.coords;
                masses[i] += 1.0 / self.particles[partitions[i][j]].inverse_mass;
            }

            if partitions[i].len() <= 0 { continue; }
            centroids[i] /= partitions[i].len() as f32;
        }

        // Then we'll use those centroids to estimate the force of gravity from
        // particles in far away cells- THIS O(# partitions^2 *# particles) BTW
        for i in 0 .. partitions.len() {
            for j in 0 .. partitions[i].len() {
                let pp =  &mut self.particles[partitions[i][j]];
                // start with fnet = 0
                let mut fnet : na::Vector3<f32> = na::Vector3::new(0.0, 0.0, 0.0);
                for k in 0 .. partitions.len() {
                    // Except for those of the current centroid of course
                    if i == k { continue; }

                    // Then add the estimated forces at all the centroids to it
                    fnet += newtonian_gravity(
                        1.0/pp.inverse_mass,
                        masses[k],
                        pp.pos.coords,
                        centroids[k]
                    );
                }

                // Finally we can multiply it by dt and add it to the velocity
                let fi = fnet * pp.inverse_mass * self.time_step;
                pp.vel.coords -= fi;
            }
        }

        // Then we can actually update their positions...
        for i in 0 .. self.particles.len() {
            let local_so_rust_shuts_up = self.particles[i].vel.coords;
            self.particles[i].pos.coords += local_so_rust_shuts_up*self.time_step;
        }
    }

    fn partition_particles(&mut self) -> Vec<Vec<usize>> {
        // We'll just store the indicies instead of a reference/pointer because Rust is being iky about it 
        let mut partitions : Vec<Vec<usize>> = vec![Vec::new(); self.num_cells_per_axis*self.num_cells_per_axis*self.num_cells_per_axis];

        // We'll store the indices of all the points in this stack, so that we don't iterate through
        // the same points multiple times (once they're allocated to a partition we'll stop checking them).
        //let mut stack: Vec<usize> = (0..= self.particles.len()-1).rev().collect();
        let mut queue: VecDeque<usize> = (0..= self.particles.len()-1).rev().collect();

        // It's a little esoteric but we're using modulos operator to iterate throuh 3 dimensions in our
        // 1-dimensional vector because Rust's multi-dimensional vectors are pretty slow.
        let mut escape_flag = false;
        for i in 0..partitions.len() {
            // If we run out of particles, this will let us leave early
            if escape_flag { break; }

            // We'll disinteangle our X, Y, Z from the index
            let (x, y, z) = deserialize_partition_index(i, self.num_cells_per_axis);
           
            // To prevent re-counting the length each time we'll
            // store it and subtract it ourselves.
            let mut l = queue.len();
            let mut j = 0;

            while j < l {
                // Since usize can't be negative we'll use number of particles + 1 as our
                // flag for "no more particles in the queue"
                let index: usize = queue.pop_front().unwrap_or(self.particles.len()+1);

                // If we see it, we'll set this other flag that'll let us break the for 
                // loop early.
                if index > self.particles.len() {
                    escape_flag = true;
                }

                // Check if the point is in the cells hitbox
                let inside = point_in_cube(
                    self.partition_edge_len,
                    (x as f32)*self.partition_edge_len,
                    (y as f32)*self.partition_edge_len,
                    (z as f32)*self.partition_edge_len,
                    self.particles[index].pos.coords.x,
                    self.particles[index].pos.coords.y,
                    self.particles[index].pos.coords.z,
                );

                // If it's inside the bounding box, add it to the partition
                if inside {
                    partitions[i].push(index);
                    l -= 1;
                } else {
                    // Otherwise drop him back into the queue
                    j += 1;
                    queue.push_back(index);
                }

            }
        }

        // Then we'll find any particles that left the cells, and bounce them back into
        // the cells.
        for i in 0 .. queue.len() {
            // Invert their velocity
            self.particles[queue[i]].vel *= -0.2;

            // Then since they're out-of-bounds we'll give them 5 free tick to get back into
            // bounds. (That's 5 because we're bouncing him back with 1/5th the velocity).
            let d = self.particles[queue[i]].vel.coords;
            self.particles[queue[i]].pos += d * self.time_step*5.0;
        }

        return partitions;
    }

    pub fn new(num_particles: i32) -> ParticleSim {
        // Make a new particle simulation
        let mut p = ParticleSim {
            particles: Vec::new(),
            time_step: 0.01,
            partition_edge_len: 10.0,
            num_cells_per_axis: 10,
        };

        // Fill it with n random particles
        for _i in 0..num_particles {
            p.particles.push( random_particle() );
        }

        // Return it
        return p;
    }

    pub fn draw_partitions(&mut self, window: &mut Window) {
        // Draw the blue "cage" that surrounds the simulation
        let edge_len = (self.num_cells_per_axis as f32) * self.partition_edge_len;
        window.draw_line(
            &na::Point3::new(0.0, 0.0, 0.0),
            &na::Point3::new(0.0, 0.0, edge_len),
            &na::Point3::new(0.0, 0.0, 0.5),
        );

        window.draw_line(
            &na::Point3::new(0.0, 0.0, edge_len),
            &na::Point3::new(edge_len, 0.0, edge_len),
            &na::Point3::new(0.0, 0.0, 0.5),
        );

        window.draw_line(
            &na::Point3::new(edge_len, 0.0, edge_len),
            &na::Point3::new(edge_len, 0.0, 0.0),
            &na::Point3::new(0.0, 0.0, 0.5),
        );
        window.draw_line(
            &na::Point3::new(edge_len, 0.0, 0.0),
            &na::Point3::new(0.0, 0.0, 0.0),
            &na::Point3::new(0.0, 0.0, 0.5),
        );

        // Draw the cell divisions on the floor 
        // We're iterating self.num_cells_per_axis times more than we need to here but that's ok
        for i in 0 .. self.num_cells_per_axis*self.num_cells_per_axis*self.num_cells_per_axis {
            let (x, _y, z) = deserialize_partition_index(i, self.num_cells_per_axis);
            window.draw_line(
                &na::Point3::new(x as f32*self.partition_edge_len, 0.0, z as f32*self.partition_edge_len),
                &na::Point3::new((x+1) as f32*self.partition_edge_len, 0.0, (z+0) as f32*self.partition_edge_len),
                &na::Point3::new(0.0, 0.5, 0.0),
            );

            window.draw_line(
                &na::Point3::new(x as f32*self.partition_edge_len, 0.0, z as f32*self.partition_edge_len),
                &na::Point3::new((x+0) as f32*self.partition_edge_len, 0.0, (z+1) as f32*self.partition_edge_len),
                &na::Point3::new(0.0, 0.5, 0.0),
            );
        }
    }
}

fn deserialize_partition_index(i: usize, num_cells_per_axis: usize) -> (usize, usize, usize) {
    let x: usize = i / (num_cells_per_axis*num_cells_per_axis);
    let y: usize = (i / num_cells_per_axis) % num_cells_per_axis;
    let z: usize = i % num_cells_per_axis;
    return (x, y, z);
}

pub fn random_particle() -> Particle {
    let mut rng = rand::rng();

    const POS_RANGE : std::ops::Range<f32> = 0.0 .. 99.0;

    let m = 1.0/rng.random_range(1.0 .. 10.0);
    let p1 = na::Point3::new(rng.random_range(POS_RANGE), rng.random_range(POS_RANGE), rng.random_range(POS_RANGE));
    let p2 = na::Point3::new(rng.random(), rng.random(), rng.random());
    let p3 = na::Point3::new(m, 0.0, 0.0);

    Particle {
        pos: p1,
        vel: p2*100.0,
        inverse_mass: m,
        color: p3
    }
}

fn newtonian_gravity(m1: f32, m2: f32, p1: na::Vector3<f32>, p2: na::Vector3<f32>) -> na::Vector3<f32> {
    return (0.0001)*(m1*m2)/( (p1 - p2).magnitude().powf(2.0) ) * (p1 - p2);
}

fn point_in_cube(edge_len: f32, offx:f32, offy:f32, offz:f32, px:f32, py:f32, pz:f32) -> bool {
    // If it's below the floor, or the two negative walls, it's not in our cube
    if px < offx || py < offy || pz < offz {
        return false;
    }

    // If it's above the cieling, or past any of the two positive walls, it's not in our cube
    if px > offx + edge_len || py > offy + edge_len || pz > offz + edge_len {
        return false;
    }

    // otherwise I guess it's in our cube
    return true;

}
