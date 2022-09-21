extern crate sdl2;

pub mod input_handler;
pub mod fluid_field;
pub mod color_density;
pub mod to_color;

use color_density::ColorDensity;
use input_handler::Input;
use sdl2::mouse::MouseButton;
use fluid_field::{ FluidField, FIELD_WIDTH, FIELD_HEIGHT };
use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use std::time::Duration;


const SCREEN_WIDTH: u32 = 800;
const SCREEN_HEIGHT: u32 = 800;

fn get_source(
    input: &Input,
    dt: f64,
) -> (
    Vec<(usize, usize, ColorDensity)>,
    Vec<(usize, usize, f64)>,
    Vec<(usize, usize, f64)>,
) {
    let mut d_add = vec![];
    let mut u_add = vec![];
    let mut v_add = vec![];

    let tile_w = (SCREEN_WIDTH as f64 / (FIELD_WIDTH-1) as f64) as i16;
    let tile_h = (SCREEN_HEIGHT as f64 / (FIELD_HEIGHT-1) as f64) as i16;
    let tx = (input.mouse_position().x / tile_w as i32) as usize + 1;
    let ty = (input.mouse_position().y / tile_h as i32) as usize + 1;

    let d_amount = 9.0;
    let v_amount = 5.0;
    if input.is_mouse_down(&MouseButton::Left) {
        d_add.push((tx, ty, ColorDensity::new(0.0, d_amount, 0.0)));
    }
    if input.is_mouse_down(&MouseButton::Right) {
        d_add.push((tx, ty, ColorDensity::new(d_amount, 0.0, d_amount*0.5)));
    }
    u_add.push((tx, ty, input.mouse_movement().x as f64 * v_amount * dt));
    v_add.push((tx, ty, input.mouse_movement().y as f64 * v_amount * dt));

    (d_add, u_add, v_add)
}

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

    // Create fluid field
    let mut field = FluidField::new();

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
        
        // Add sources to the fluid field
        let dt = 1.0/30.0;
        let (d_add, u_add, v_add) = get_source(&input, dt);
        field.add_source(d_add);
        field.add_velocity(u_add, v_add);

        // Update fluid field
        field.update(dt);
        field.display_density(&canvas);

        canvas.present();

        std::thread::sleep(Duration::from_millis(1_000u64 / 60));
    }

    Ok(())
}
