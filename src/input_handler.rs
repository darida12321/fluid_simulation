use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::mouse::MouseButton;
use std::collections::HashSet;
use vector2d::Vector2D;

pub struct Input {
    keyboard_state: HashSet<Keycode>,
    keyboard_state_pressed: HashSet<Keycode>,
    keyboard_state_released: HashSet<Keycode>,
    mouse_pos: Vector2D<i32>,
    mouse_pos_diff: Vector2D<i32>,
    mouse_state: HashSet<MouseButton>,
    mouse_state_pressed: HashSet<MouseButton>,
    mouse_state_released: HashSet<MouseButton>,
}
impl Input {
    pub fn new() -> Self {
        Input {
            keyboard_state: HashSet::new(),
            keyboard_state_pressed: HashSet::new(),
            keyboard_state_released: HashSet::new(),
            mouse_pos: Vector2D::new(0, 0),
            mouse_pos_diff: Vector2D::new(0, 0),
            mouse_state: HashSet::new(),
            mouse_state_pressed: HashSet::new(),
            mouse_state_released: HashSet::new(),
        }
    }

    pub fn update(&mut self) {
        self.keyboard_state_pressed.clear();
        self.keyboard_state_released.clear();
        self.mouse_state_pressed.clear();
        self.mouse_state_released.clear();
        self.mouse_pos_diff = Vector2D::new(0, 0);
    }

    pub fn handle_event(&mut self, event: &Event) {
        match event {
            Event::KeyDown {
                keycode: Some(code),
                ..
            } => {
                self.keyboard_state_pressed.insert(*code);
                self.keyboard_state.insert(*code);
            }
            Event::KeyUp {
                keycode: Some(code),
                ..
            } => {
                self.keyboard_state_released.insert(*code);
                self.keyboard_state.remove(&code);
            }
            Event::MouseMotion { x, y, xrel, yrel, .. } => {
                self.mouse_pos = Vector2D::new(*x, *y);
                self.mouse_pos_diff += Vector2D::new(*xrel, *yrel);
            }
            Event::MouseButtonDown { mouse_btn, .. } => {
                self.mouse_state_pressed.insert(*mouse_btn);
                self.mouse_state.insert(*mouse_btn);
            }
            Event::MouseButtonUp { mouse_btn, .. } => {
                self.mouse_state_released.insert(*mouse_btn);
                self.mouse_state.remove(&mouse_btn);
            }
            _ => (),
        }
    }

    pub fn is_key_down(&self, key: &Keycode) -> bool {
        self.keyboard_state.contains(key)
    }
    pub fn is_key_pressed(&self, key: &Keycode) -> bool {
        self.keyboard_state_pressed.contains(key)
    }
    pub fn is_key_released(&self, key: &Keycode) -> bool {
        self.keyboard_state_released.contains(key)
    }

    pub fn is_mouse_down(&self, key: &MouseButton) -> bool {
        self.mouse_state.contains(key)
    }
    pub fn is_mouse_pressed(&self, key: &MouseButton) -> bool {
        self.mouse_state_pressed.contains(key)
    }
    pub fn is_mouse_released(&self, key: &MouseButton) -> bool {
        self.mouse_state_released.contains(key)
    }

    pub fn mouse_position(&self) -> Vector2D<i32> {
        self.mouse_pos.clone()
    }
    pub fn mouse_movement(&self) -> Vector2D<i32> {
        self.mouse_pos_diff.clone()
    }
}
