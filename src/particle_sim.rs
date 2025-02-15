use rand::Rng;
use kiss3d::window::Window;
extern crate nalgebra as na;

pub struct Particle {
    pub pos: na::Point<f32, 3>,
    pub vel: na::Point<f32, 3>,
    pub inverse_mass: f32,
    pub color: na::Point<f32, 3>,
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

pub struct ParticleSim {
    pub particles: Vec<Particle>,
    pub time_step: f32,
    num_cells_per_axis: usize,
    partition_edge_len: f32,
}

impl ParticleSim {
    pub fn draw_partitions(&mut self, window: &mut Window) {

        let edge_len = (self.num_cells_per_axis as f32) * self.partition_edge_len;

        // 1
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
    }

    fn partition_particles(&mut self) -> Vec<Vec<usize>> {
        // We'll just store the indicies instead of a reference/pointer because Rust is being iky about it 
        let mut partitions : Vec<Vec<usize>> = vec![Vec::new(); self.num_cells_per_axis*self.num_cells_per_axis*self.num_cells_per_axis];

        // We'll store the indices of all the points in this stack, so that we don't iterate through
        // the same points multiple times (once they're allocated to a partition we'll stop checking them).
        let mut stack: Vec<usize> = (0..= self.particles.len()-1).rev().collect();

        // It's a little esoteric but we're using modulos operator to iterate throuh 3 dimensions in our
        // 1-dimensional vector because Rust's multi-dimensional vectors are pretty slow.
        for i in 0..partitions.len() {
            // We'll disinteangle our X, Y, Z from the index
            let x: usize = i / (self.num_cells_per_axis*self.num_cells_per_axis);
            let y: usize = (i / self.num_cells_per_axis) % self.num_cells_per_axis;
            let z: usize = i % self.num_cells_per_axis;
           
            // To prevent re-counting the length each time we'll
            // store it and subtract it ourselves.
            let mut l = stack.len();
            let mut j = 0;

            while j < l {
                let index = stack[j];

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

                // If it is, add it to the cell's partition,
                // then take it out of the stack of particle ID's
                // we're sorting.
                if inside {
                    partitions[i].push(index);
                    stack.remove(j);
                    l -= 1;
                } else {
                    j += 1;
                }

            }
        }

        // Then we'll find any particles that left the cells, and bounce them back into
        // the cells.
        for i in 0 .. stack.len() {
            // Invert their velocity
            self.particles[stack[i]].vel *= -0.2;

            // Then since they're out-of-bounds we'll give them 5 free tick to get back into
            // bounds. (That's 5 because we're bouncing him back with 1/5th the velocity).
            let d = self.particles[stack[i]].vel.coords;
            self.particles[stack[i]].pos += d * self.time_step*5.0;
        }

        return partitions;
    }

    pub fn tick(&mut self) {
        // These are computer-science vectors of mathematical vectors...
        let mut displacements : Vec<na::Vector3<f32>> = Vec::new();
        let mut accelerations : Vec<na::Vector3<f32>> = Vec::new();

        // We'll have to split our particles up into "spatial partitions" or
        // computing the forces will be O(n^2) and my computer he will cry
        let partitions = self.partition_particles();

        // Iterate once just to compute the changes. We'll have to iterate
        // again after to actually modify the values. Otherwise we'll be messing
        // with the values of some particles before computing the values of others.
        for i in 0..partitions.len() {
            for j in 0..partitions[i].len() {
                displacements.push(
                    self.particles[ partitions[i][j] ].vel.coords*self.time_step
                );

                let mut accl: na::Vector3<f32> = self.local_gravity(
                    i,
                    j,
                    &partitions,
                );

                accl += global_gravity(
                    &self.particles[ partitions[i][j]],
                    na::Vector3::new(50.0, 50.0, 50.0), 2000.0, self.time_step
                );

                accelerations.push(accl);
            }

        }

        // Then we'll iterate another time to add all the changes
        let mut count = 0;
        for i in 0..partitions.len() {
            for j in 0..partitions[i].len() {
                let p_index: usize = partitions[i][j];
                self.particles[p_index].pos += displacements[count];
                self.particles[p_index].vel += accelerations[count];
                count += 1; 
            }
        }

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

    fn local_gravity(&mut self, self_partition: usize, self_particle:usize, partitions: &Vec<Vec<usize>>) -> na::Vector3<f32>{
        // Compute actual real gravity forces actually between all particles (in the cells). Essentially making this an
        // N-body problem where N is the number of particles in the cell.
        let mut force = na::Vector3::new(0.0, 0.0, 0.0);

        // We'll go dig up the particle we're actually working with here
        let p = &self.particles[partitions[self_partition][self_particle]];

        //let (x, y, z) = deserialize_partition_index(self_partition, self.num_cells_per_axis);

        for other_index in 0 .. partitions[self_partition].len() {
            let l = partitions[self_partition][other_index];
            let other = &self.particles[l];

            if std::ptr::eq(other, p) { continue };

            let r2 : f32 = (p.pos.coords - other.pos.coords).magnitude().powf(2.0);
            let f = 1.0 / (r2 * other.inverse_mass);
            force += (p.pos.coords - other.pos.coords) * f;
        }

        return -force*(1.0/p.inverse_mass)*(11.0)*self.time_step;
    }

}

fn global_gravity(p: &Particle, centroid_pos : na::Vector3<f32>, centroid_mass: f32, dt:f32) -> na::Vector3<f32> {
    // Based off the centre-of-mass FOR ALL THE PARTICLES we'll simulate long-distance interactions.
    let r2: f32 = (p.pos.coords - centroid_pos).magnitude().powf(3.0);
    let f = 1.0 / (r2 * (1.0/centroid_mass));
    let force = (p.pos.coords - centroid_pos) * f;
    return -force*(1.0/p.inverse_mass)*(10.0)*dt;
}

fn deserialize_partition_index(i: usize, num_cells_per_axis: usize) -> (usize, usize, usize) {
    let x: usize = i / (num_cells_per_axis*num_cells_per_axis);
    let y: usize = (i / num_cells_per_axis) % num_cells_per_axis;
    let z: usize = i % num_cells_per_axis;
    return (x, y, z);
}


pub fn random_particle() -> Particle{
    let mut rng = rand::rng();

    const POS_RANGE : std::ops::Range<f32> = 0.0 .. 100.0;

    let p1 = na::Point3::new(rng.random_range(POS_RANGE), rng.random_range(POS_RANGE), rng.random_range(POS_RANGE));
    let p2 = na::Point3::new(rng.random(), rng.random(), rng.random());
    let p3 = na::Point3::new(rng.random(), rng.random(), rng.random());

    Particle {
        pos: p1,
        vel: p2*10.0,
        inverse_mass: 1.0/rng.random_range(1.0 .. 10.0),
        color: p3
    }
}


