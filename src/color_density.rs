use crate::to_color::ToColor;


#[derive(Clone, Copy, Default)]
pub struct ColorDensity {
    pub r: f64,
    pub g: f64,
    pub b: f64,
}

impl ColorDensity {
    pub fn new(r: f64, g: f64, b: f64) -> Self {
        ColorDensity { r, g, b }
    }
}

impl std::ops::Add for ColorDensity {
    type Output = ColorDensity;
    fn add(self, rhs: ColorDensity) -> ColorDensity {
        ColorDensity {
            r: self.r + rhs.r,
            g: self.g + rhs.g,
            b: self.b + rhs.b,
        }
    }
}

impl std::ops::AddAssign for ColorDensity {
    fn add_assign(&mut self, rhs: Self) {
        self.r += rhs.r;
        self.g += rhs.g;
        self.b += rhs.b;
    }
}

impl std::ops::Mul<f64> for ColorDensity {
    type Output = ColorDensity;
    fn mul(self, rhs: f64) -> ColorDensity {
        ColorDensity {
            r: self.r * rhs,
            g: self.g * rhs,
            b: self.b * rhs,
        }
    }
}

impl std::ops::DivAssign<f64> for ColorDensity {
    fn div_assign(&mut self, rhs: f64) {
        self.r /= rhs;
        self.g /= rhs;
        self.b /= rhs;
    }
}

impl std::ops::Neg for ColorDensity {
    type Output = ColorDensity;
    fn neg(self) -> ColorDensity {
        ColorDensity {
            r: -self.r,
            g: -self.g,
            b: -self.b,
        }
    }
}

impl ToColor for ColorDensity {
    fn to_color(&self) -> sdl2::pixels::Color {
        sdl2::pixels::Color::RGB(
            (self.r * 256.0) as u8,
            (self.g * 256.0) as u8,
            (self.b * 256.0) as u8,
        )
    }
}







