extern crate image;
extern crate imageproc;
extern crate rand;

pub trait GeneticSample<'a, D> {
    fn new_from_rng(data: &'a D) -> Self;

    #[allow(non_snake_case)]
    fn combine_DNA(lhs: &Self, rhs: &Self) -> Self;

    fn maybe_mutate(&mut self, p_mutate: f64);

    fn calc_fitness(&self) -> u32;
}

pub mod rgba {
    use super::GeneticSample;
    use image::RgbaImage;
    use rand::distributions::Distribution;

    pub struct RgbaApproximator<'a> {
        reference_img: &'a RgbaImage,
        triangles: Vec<[imageproc::drawing::Point<i32>; 3]>,
    }

    impl<'a> GeneticSample<'a, RgbaImage> for RgbaApproximator<'a> {
        fn new_from_rng(data: &'a RgbaImage) -> RgbaApproximator<'a> {
            RgbaApproximator {
                reference_img: data,
                triangles: [].to_vec(),
            }
        }

        fn maybe_mutate(&mut self, p_mutate: f64) {
            let decider = rand::distributions::Bernoulli::new(p_mutate).unwrap();
            if decider.sample(&mut rand::thread_rng()) {
                self.triangles
                    .push([imageproc::drawing::Point::new(0, 0); 3]);
            }
        }

        fn combine_DNA<'b>(
            lhs: &RgbaApproximator<'b>,
            rhs: &RgbaApproximator,
        ) -> RgbaApproximator<'b> {
            RgbaApproximator {
                reference_img: &lhs.reference_img,
                triangles: [
                    lhs.triangles.clone().as_slice(),
                    rhs.triangles.clone().as_slice(),
                ]
                .concat(),
            }
        }

        fn calc_fitness(&self) -> u32 {
            self.triangles.len() as u32
        }
    }
}

pub mod placeholder {
    use super::GeneticSample;
    use rand::distributions::Distribution;
    use rand::Rng;

    #[derive(Clone)]
    struct Point {
        x: u32,
        y: u32,
    }

    pub struct ApproximatorPlaceholder {
        i: usize,
        triangles: Vec<[Point; 3]>,
    }

    impl GeneticSample<'_, u32> for ApproximatorPlaceholder {
        fn new_from_rng<'a>(_data: &'a u32) -> ApproximatorPlaceholder {
            ApproximatorPlaceholder {
                i: rand::thread_rng().gen_range(0, 10),
                triangles: [].to_vec(),
            }
        }

        fn maybe_mutate(&mut self, p_mutate: f64) {
            let decider = rand::distributions::Bernoulli::new(p_mutate).unwrap();
            if decider.sample(&mut rand::thread_rng()) {
                self.i += 1;
            }
        }

        fn combine_DNA(
            lhs: &ApproximatorPlaceholder,
            rhs: &ApproximatorPlaceholder,
        ) -> ApproximatorPlaceholder {
            ApproximatorPlaceholder {
                i: lhs.i + rhs.i,
                triangles: [
                    lhs.triangles.clone().as_slice(),
                    rhs.triangles.clone().as_slice(),
                ]
                .concat(),
            }
        }

        fn calc_fitness(&self) -> u32 {
            self.i as u32
        }
    }
}
