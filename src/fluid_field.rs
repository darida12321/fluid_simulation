
use crate::{SCREEN_HEIGHT, SCREEN_WIDTH};
use sdl2::gfx::primitives::DrawRenderer;
use sdl2::render::Canvas;
use sdl2::video::Window;
use crate::to_color::ToColor;

pub const FIELD_WIDTH: usize = 50;
pub const FIELD_HEIGHT: usize = 50;



type Field<T> = Vec<Vec<T>>;
pub struct FluidField<T> {
    d: Field<T>,
    u: Field<f64>,
    v: Field<f64>,
    d0: Field<T>,
    u0: Field<f64>,
    v0: Field<f64>,
}

impl<T> FluidField<T> where
    T: std::ops::Add<Output = T>,
    T: std::ops::AddAssign<T>,
    T: std::ops::Mul<f64, Output = T>,
    T: std::ops::DivAssign<f64>,
    T: std::ops::Neg<Output = T>,
    T: Copy,
    T: Default,
    T: ToColor, // TODO Use the correct ToColor method from sdl2
{
    pub fn new() -> Self {
        FluidField {
            d: vec![vec![T::default(); FIELD_HEIGHT+2]; FIELD_WIDTH+2],
            u: vec![vec![0.0; FIELD_HEIGHT+2]; FIELD_WIDTH+2],
            v: vec![vec![0.0; FIELD_HEIGHT+2]; FIELD_WIDTH+2],
            d0: vec![vec![T::default(); FIELD_HEIGHT+2]; FIELD_WIDTH+2],
            u0: vec![vec![0.0; FIELD_HEIGHT+2]; FIELD_WIDTH+2],
            v0: vec![vec![0.0; FIELD_HEIGHT+2]; FIELD_WIDTH+2],
        }
    }

    pub fn add_source(&mut self, d_add: Vec<(usize, usize, T)>) {
        for (x, y, d) in d_add.iter() {
            self.d[*x][*y] += *d;
        }
    }

    pub fn add_velocity(
        &mut self,
        u_add: Vec<(usize, usize, f64)>,
        v_add: Vec<(usize, usize, f64)>,
    ) {
        for (x, y, u) in u_add.iter() {
            self.u[*x][*y] += u;
        }
        for (x, y, v) in v_add.iter() {
            self.v[*x][*y] += v;
        }
    }

    pub fn update(&mut self, dt: f64) {
        // Set diffusion and viscosity
        let diff = 0.001;
        let visc = 0.000001;
        
        // Update fluid field
        self.vel_step(visc, dt);
        self.dens_step(diff, dt);
    }

    pub fn display_density(&self, canvas: &Canvas<Window>) {
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
                        self.d[i][j].to_color()
                    )
                    .unwrap();
            }
        }
    }

    // Diffuse a field with its neighbors
    fn diffuse<U>(
        b: i32, x: &mut Field<U>, x0: &Field<U>, diff: f64, dt: f64
    ) where
        U: std::ops::Add<Output = U>,
        U: std::ops::Mul<f64, Output = U>,
        U: std::ops::DivAssign<f64>,
        U: std::ops::Neg<Output = U>,
        U: Copy,
    {
        let a = dt * diff * FIELD_WIDTH as f64 * FIELD_HEIGHT as f64;
        for _ in 0..20 {
            for i in 1..=FIELD_WIDTH {
                for j in 1..=FIELD_HEIGHT {
                    x[i][j] = x0[i][j] + (x[i - 1][j] + x[i + 1][j] + x[i][j - 1] + x[i][j + 1]) * a;
                    x[i][j] /= 1.0 + 4.0 * a;
                }
            }
            Self::set_bnd(b, x);
        }
    }

    fn advect<U>(b: i32, d: &mut Field<U>, d0: &Field<U>, u: &Field<f64>, v: &Field<f64>, dt: f64
    ) where
        U: std::ops::Add<Output = U>,
        U: std::ops::Mul<f64, Output = U>,
        U: std::ops::DivAssign<f64>,
        U: std::ops::Neg<Output = U>,
        U: Copy,
    {
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
                    (d0[i0][j0]*t0 + d0[i0][j1]*t1)*s0 + (d0[i1][j0]*t0 + d0[i1][j1]*t1)*s1;
            }
        }
        Self::set_bnd(b, d);
    }

    fn dens_step(&mut self, diff: f64, dt: f64) {
        std::mem::swap(&mut self.d0, &mut self.d);
        Self::diffuse(0, &mut self.d, &self.d0, diff, dt);
        std::mem::swap(&mut self.d0, &mut self.d);
        Self::advect(0, &mut self.d, &mut self.d0, &mut self.u, &mut self.v, dt);
    }

    fn vel_step(&mut self, visc: f64, dt: f64) {
        std::mem::swap(&mut self.u0, &mut self.u);
        Self::diffuse(1, &mut self.u, &self.u0, visc, dt);
        std::mem::swap(&mut self.v0, &mut self.v);
        Self::diffuse(2, &mut self.v, &self.v0, visc, dt);
        Self::project(&mut self.u, &mut self.v, &mut self.u0, &mut self.v0);
        std::mem::swap(&mut self.u0, &mut self.u);
        std::mem::swap(&mut self.v0, &mut self.v);
        Self::advect(1, &mut self.u, &self.u0, &self.u0, &self.v0, dt);
        Self::advect(2, &mut self.v, &self.v0, &self.u0, &self.v0, dt);
        Self::project(&mut self.u, &mut self.v, &mut self.u0, &mut self.v0);
    }

    fn project(u: &mut Field<f64>, v: &mut Field<f64>, p: &mut Field<f64>, div: &mut Field<f64>) {
        let h = 1.0 / FIELD_WIDTH as f64;
        for i in 1..=FIELD_WIDTH {
            for j in 1..=FIELD_HEIGHT {
                div[i][j] = -0.5 * h * (u[i+1][j] - u[i-1][j] + v[i][j+1] - v[i][j-1]);
                p[i][j] = 0.0;
            }
        }
        Self::set_bnd(0, div);
        Self::set_bnd(0, p);

        for _ in 0..20 {
            for i in 1..=FIELD_WIDTH {
                for j in 1..=FIELD_HEIGHT {
                    p[i][j] = (div[i][j] + p[i-1][j] + p[i+1][j] + p[i][j-1] + p[i][j+1])/4.0;
                }
            }
            Self::set_bnd(0, p);
        }

        for i in 1..=FIELD_WIDTH {
            for j in 1..=FIELD_HEIGHT {
                u[i][j] -= 0.5*(p[i+1][j] - p[i-1][j])/h;
                v[i][j] -= 0.5*(p[i][j+1] - p[i][j-1])/h;
            }
        }
        Self::set_bnd(1, u);
        Self::set_bnd(2, v);
    }

    fn set_bnd<U>(b: i32, x: &mut Field<U>)
    where
        U: std::ops::Add<Output = U>,
        U: std::ops::Mul<f64, Output = U>,
        U: std::ops::Neg<Output = U>,
        U: Copy,
    {
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
        x[0][0] = (x[1][0] + x[0][1]) * 0.5;
        x[0][H+1] = (x[1][H+1] + x[0][H]) * 0.5;
        x[W+1][0] = (x[W][0] + x[W+1][1]) * 0.5;
        x[W+1][H+1] = (x[W][H+1] + x[W+1][H]) * 0.5;
    }
}
