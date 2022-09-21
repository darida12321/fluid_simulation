
use sdl2::pixels::Color;

pub trait ToColor {
    fn to_color(&self) -> Color;
}
impl ToColor for f64 {
    fn to_color(&self) -> Color {
        Color::RGB(0, 0, (*self*256.0) as u8)
    }
}

