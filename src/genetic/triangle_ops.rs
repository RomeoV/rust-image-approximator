extern crate rand;
// use imageproc::drawing::Point;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct Point<T: Copy + PartialEq + Eq> {
    pub x: T,
    pub y: T,
}

impl<T: Copy + PartialEq + Eq> Point<T> {
    // <T: Copy + PartialEq + Eq>
    pub fn to_imageproc_point(&self) -> imageproc::drawing::Point<T> {
        imageproc::drawing::Point::new(self.x, self.y)
    }
}

#[derive(Clone)]
pub struct Triangle(pub [Point<i32>; 3], pub image::Rgba<u8>);

pub mod triangle_ops {
    use super::Triangle;
    use rand::Rng;
    use std::f64;

    impl Triangle {
        pub fn resize(&mut self) {
            let center_x = self.0.iter().map(|p| p.x).sum::<i32>() as f64 / self.0.len() as f64;
            let center_y = self.0.iter().map(|p| p.y).sum::<i32>() as f64 / self.0.len() as f64;
            let factor: f64 = rand::thread_rng().gen_range(-0.5, 0.3);
            for p in &mut self.0 {
                p.x += (factor*(p.x as f64 - center_x)) as i32;
                p.y += (factor*(p.y as f64 - center_y)) as i32;
            }
        }

        pub fn shift(&mut self, (w, h): (u32, u32)) {
            let shift_x = rand::thread_rng().gen_range(-1.*w as f64/3., w as f64/3.) as i32;
            let shift_y = rand::thread_rng().gen_range(-1.*h as f64/3., h as f64/3.) as i32;
            for p in &mut self.0 {
                p.x += shift_x as i32;
                p.y += shift_y as i32;
            }
        }

        pub fn rotate(&mut self) {
            let center_x = self.0.iter().map(|p| p.x).sum::<i32>() as f64 / self.0.len() as f64;
            let center_y = self.0.iter().map(|p| p.y).sum::<i32>() as f64 / self.0.len() as f64;
            let angle: f64 = rand::thread_rng().gen_range(-f64::consts::PI/4., f64::consts::PI/4.);
            for p in &mut self.0 {
                let dx: f64 = p.x as f64 - center_x;
                let dy: f64 = p.y as f64 - center_y;
                p.x = (angle.cos()*dx - angle.sin()*dy) as i32;
                p.y = (angle.sin()*dx + angle.cos()*dy) as i32;
            }

        }
    }
}
