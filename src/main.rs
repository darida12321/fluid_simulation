extern crate sdl2;

pub mod input_handler;
pub mod fluid_field;

use input_handler::Input;
use fluid_field::*;
use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use std::time::Duration;


const SCREEN_WIDTH: u32 = 800;
const SCREEN_HEIGHT: u32 = 800;

fn main() -> Result<(), String> {
    let sdl = sdl2::init()?;
    let video_subsystem = sdl.video()?;

    let window = video_subsystem
        .window("Thrust", SCREEN_WIDTH, SCREEN_HEIGHT)
        .position_centered()
        .opengl()
        .build()
        .map_err(|e| e.to_string())?;

    let mut canvas = window.into_canvas().build().map_err(|e| e.to_string())?;
    let mut input = Input::new();

    let mut d = create_field();
    let mut u = create_field();
    let mut v = create_field();

    let mut event_pump = sdl.event_pump()?;
    'running: loop {
        input.update();
        for event in event_pump.poll_iter() {
            input.handle_event(&event);
            match event {
                Event::Quit { .. }
                | Event::KeyDown {
                    keycode: Some(Keycode::Escape),
                    ..
                } => break 'running,
                _ => {}
            }
        }

        canvas.set_draw_color(Color::RGB(20, 20, 20));
        canvas.clear();

        let diff = 0.001;
        let visc = 0.001;

        let dt = 1.0/30.0;
        add_source(&input, &mut d, &mut u, &mut v, dt);
        let mut d0 = d.clone();
        let mut u0 = u.clone();
        let mut v0 = v.clone();
        vel_step(&mut u, &mut v, &mut u0, &mut v0, visc, dt);
        dens_step(&mut d, &mut d0, &u, &v, diff, dt);


        display_density(&canvas, &d);
        //display_velocity(&canvas, &u, &v);

        canvas.present();

        std::thread::sleep(Duration::from_millis(1_000u64 / 60));
    }

    Ok(())
}
