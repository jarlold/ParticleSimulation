extern crate kiss3d;
extern crate nalgebra as na;

use kiss3d::camera::ArcBall;
use kiss3d::light::Light;
use kiss3d::window::Window;
use na::Point3;
mod particle_sim;

fn main() {
    // Setup the window
    let mut window = Window::new("Kiss3d: camera");

    // Setup the camera
    let eye = Point3::new(200.0f32, 200.0, 200.0);
    let mut arc_ball = ArcBall::new(eye, Point3::origin());

    // Add in a light
    window.set_light(Light::StickToCamera);

    // Then we'll make a particle simulation
    let mut simulation = particle_sim::ParticleSim::new(2000);

    window.set_background_color(1.0, 1.0, 1.0);

    window.set_point_size(100.0);

    // Main draw loop
    while !window.should_close() {
        // Draw all the particles
        for p in &simulation.particles {
            window.draw_point(
                &p.pos,
                &p.color
            );
        }

        // Draw a compass so I don't get lost
        // Draw x-axis
        window.draw_line(
            &Point3::new(0.0, 0.0, 0.0),
            &Point3::new(1.0, 0.0, 0.0),
            &Point3::new(1.0, 0.0, 0.0),
        );

        // Draw y-axis
        window.draw_line(
            &Point3::new(0.0, 0.0, 0.0),
            &Point3::new(0.0, 1.0, 0.0),
            &Point3::new(0.0, 1.0, 0.0),
        );

        // Draw z-axis
        window.draw_line(
            &na::Point3::new(0.0, 0.0, 0.0),
            &na::Point3::new(0.0, 0.0, 1.0),
            &na::Point3::new(0.0, 0.0, 1.0),
        );

        // Draw lines where the partition cells are
        simulation.draw_partitions(&mut window);

        // Then we'll do a physics tick in the draw thread because
        // that's life.
        simulation.tick();

        // And keep our orbit camera
        window.render_with_camera(&mut arc_ball);
    }
}
