use crate::input_handler::Input;
use crate::{SCREEN_HEIGHT, SCREEN_WIDTH};
use sdl2::gfx::primitives::DrawRenderer;
use sdl2::mouse::MouseButton;
use sdl2::pixels::Color;
use sdl2::render::Canvas;
use sdl2::video::Window;

const FIELD_WIDTH: usize = 50;
const FIELD_HEIGHT: usize = 50;

type Field = Vec<Vec<f64>>;

pub fn create_field() -> Field {
    return vec![vec![0.0; FIELD_HEIGHT+2]; FIELD_WIDTH+2]
}

pub fn display_density(canvas: &Canvas<Window>, d: &Field) {
    let tile_w = (SCREEN_WIDTH as f64 / (FIELD_WIDTH-1) as f64) as i16;
    let tile_h = (SCREEN_HEIGHT as f64 / (FIELD_HEIGHT-1) as f64) as i16;
    for i in 1..=FIELD_WIDTH {
        for j in 1..=FIELD_HEIGHT {
            let x = (i - 1) as i16 * tile_w;
            let y = (j - 1) as i16 * tile_h;
            canvas
                .filled_polygon(
                    &[x, x, (x + tile_w), (x + tile_w)],
                    &[y, (y + tile_h), (y + tile_h), y],
                    Color::RGB(0, 0, (d[i][j] * 256.0) as u8),
                )
                .unwrap();
        }
    }
}

pub fn display_velocity(canvas: &Canvas<Window>, u: &Field, v: &Field) {
    let tile_w = (SCREEN_WIDTH as f64 / (FIELD_WIDTH-1) as f64) as i16;
    let tile_h = (SCREEN_HEIGHT as f64 / (FIELD_HEIGHT-1) as f64) as i16;
    for i in 1..=FIELD_WIDTH {
        for j in 1..=FIELD_HEIGHT {
            let x = (i - 1) as i16 * tile_w + tile_w / 2;
            let y = (j - 1) as i16 * tile_h + tile_h / 2;
            let vel_x = (u[i][j] * tile_w as f64) as i16;
            let vel_y = (v[i][j] * tile_h as f64) as i16;
            canvas.line(x, y, x + vel_x, y + vel_y, Color::RED).unwrap();
        }
    }
}

pub fn add_source(input: &Input, d: &mut Field, u: &mut Field, v: &mut Field, dt: f64) {
    let d_amount = 9.0;
    let v_amount = 5.0;
    let tile_w = (SCREEN_WIDTH as f64 / (FIELD_WIDTH-1) as f64) as i16;
    let tile_h = (SCREEN_HEIGHT as f64 / (FIELD_HEIGHT-1) as f64) as i16;
    let tx = (input.mouse_position().x / tile_w as i32) as usize + 1;
    let ty = (input.mouse_position().y / tile_h as i32) as usize + 1;
    // Add sources
    if input.is_mouse_down(&MouseButton::Left) {
        d[tx][ty] = d_amount;
    }
    if input.is_mouse_down(&MouseButton::Right) {
        d[tx][ty] = -d_amount;
    }
    u[tx][ty] += input.mouse_movement().x as f64 * v_amount * dt;
    v[tx][ty] += input.mouse_movement().y as f64 * v_amount * dt;
}

pub fn diffuse(b: i32, x: &mut Field, x0: &Field, diff: f64, dt: f64) {
    let a = dt * diff * FIELD_WIDTH as f64 * FIELD_HEIGHT as f64;
    for _ in 0..20 {
        for i in 1..=FIELD_WIDTH {
            for j in 1..=FIELD_HEIGHT {
                x[i][j] = x0[i][j] + a * (x[i - 1][j] + x[i + 1][j] + x[i][j - 1] + x[i][j + 1]);
                x[i][j] /= 1.0 + 4.0 * a;
            }
        }
        set_bnd(b, x);
    }
}

pub fn advect(b: i32, d: &mut Field, d0: &Field, u: &Field, v: &Field, dt: f64) {
    let dt0 = dt * FIELD_WIDTH as f64;
    for i in 1..=FIELD_WIDTH {
        for j in 1..=FIELD_HEIGHT {
            let mut x = i as f64 - dt0 * u[i][j];
            let mut y = j as f64 - dt0 * v[i][j];

            x = x.max(0.5).min(FIELD_WIDTH as f64 + 0.5);
            y = y.max(0.5).min(FIELD_HEIGHT as f64 + 0.5);

            let i0 = x.floor() as usize;
            let i1 = i0 + 1;
            let j0 = y.floor() as usize;
            let j1 = j0 + 1;

            let s1 = x - i0 as f64;
            let s0 = 1.0 - s1;
            let t1 = y - j0 as f64;
            let t0 = 1.0 - t1;

            d[i][j] =
                s0 * (t0 * d0[i0][j0] + t1 * d0[i0][j1]) + s1 * (t0 * d0[i1][j0] + t1 * d0[i1][j1]);
        }
    }
    set_bnd(b, d);
}

pub fn dens_step(x: &mut Field, x0: &mut Field, u: &Field, v: &Field, diff: f64, dt: f64) {
    std::mem::swap(x0, x);
    diffuse(0, x, x0, diff, dt);
    std::mem::swap(x0, x);
    advect(0, x, x0, u, v, dt);
}

pub fn vel_step(u: &mut Field, v: &mut Field, u0: &mut Field, v0: &mut Field, visc: f64, dt: f64) {
    std::mem::swap(u0, u);
    diffuse(1, u, u0, visc, dt);
    std::mem::swap(v0, v);
    diffuse(2, v, v0, visc, dt);
    project(u, v, u0, v0);
    std::mem::swap(u0, u);
    std::mem::swap(v0, v);
    advect(1, u, u0, u0, v0, dt);
    advect(2, v, v0, u0, v0, dt);
    project(u, v, u0, v0);
}

fn project(u: &mut Field, v: &mut Field, p: &mut Field, div: &mut Field) {
    let h = 1.0 / FIELD_WIDTH as f64;
    for i in 1..=FIELD_WIDTH {
        for j in 1..=FIELD_HEIGHT {
            div[i][j] = -0.5 * h * (u[i+1][j] - u[i-1][j] + v[i][j+1] - v[i][j-1]);
            p[i][j] = 0.0;
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);

    for _ in 0..20 {
        for i in 1..=FIELD_WIDTH {
            for j in 1..=FIELD_HEIGHT {
                p[i][j] = (div[i][j] + p[i-1][j] + p[i+1][j] + p[i][j-1] + p[i][j+1])/4.0;
            }
        }
        set_bnd(0, p);
    }

    for i in 1..=FIELD_WIDTH {
        for j in 1..=FIELD_HEIGHT {
            u[i][j] -= 0.5*(p[i+1][j] - p[i-1][j])/h;
            v[i][j] -= 0.5*(p[i][j+1] - p[i][j-1])/h;
        }
    }
    set_bnd(1, u);
    set_bnd(2, v);
}

pub fn set_bnd(b: i32, x: &mut Field) {
    const W: usize = FIELD_WIDTH;
    const H: usize = FIELD_HEIGHT;
    for i in 1..=W {
        x[i][0] = if b==2 {-x[i][1]} else {x[i][1]};
        x[i][H+1] = if b==2 {-x[i][H]} else {x[i][H]};
    }
    for i in 1..=H {
        x[0][i] = if b==1 {-x[1][i]} else {x[1][i]};
        x[W+1][i] = if b==1 {-x[W][i]} else {x[W][i]};
    }
    x[0][0] = 0.5*(x[1][0] + x[0][1]);
    x[0][H+1] = 0.5*(x[1][H+1] + x[0][H]);
    x[W+1][0] = 0.5*(x[W][0] + x[W+1][1]);
    x[W+1][H+1] = 0.5*(x[W][H+1] + x[W+1][H]);
}
