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
    use image::Pixel;
    use rand::distributions::Distribution;
    use rand::Rng;

    pub struct RgbaApproximator<'a> {
        reference_img: &'a RgbaImage,
        triangles: Vec<[imageproc::drawing::Point<i32>; 3]>,
    }

    impl RgbaApproximator<'_> {
        fn random_triangle(img: &RgbaImage) -> [imageproc::drawing::Point<i32>; 3] {
            let (w, h) = img.dimensions();
            let mut rng = rand::thread_rng();
            let mut ret = [imageproc::drawing::Point::<i32>::new(0,0); 3];
            for p in &mut ret {
                *p = imageproc::drawing::Point::<i32>::new(
                        rng.gen_range(0, w-1) as i32,
                        rng.gen_range(0, h-1) as i32);
            }
            while ret[0] == ret[1] || ret[1] == ret[2] || ret[0] == ret[2] {
                ret = RgbaApproximator::random_triangle(&img);
            }
            return ret
        }
    }

    impl<'a> GeneticSample<'a, RgbaImage> for RgbaApproximator<'a> {
        fn new_from_rng(img: &'a RgbaImage) -> RgbaApproximator<'a> {
            const N_INITIAL_TRIANGLES: usize = 50;
            RgbaApproximator {
                reference_img: img,
                triangles: (0..N_INITIAL_TRIANGLES).map(|_| RgbaApproximator::random_triangle(img)).collect(),
            }
        }

        fn maybe_mutate(&mut self, p_mutate: f64) {
            let decider = rand::distributions::Bernoulli::new(p_mutate).unwrap();
            if decider.sample(&mut rand::thread_rng()) {
                self.triangles
                    .push(RgbaApproximator::random_triangle(&self.reference_img));
            }
        }

        fn combine_DNA<'b>(
            lhs: &RgbaApproximator<'b>,
            rhs: &RgbaApproximator,
        ) -> RgbaApproximator<'b> {
            const P_USE_TRIANGLE: f64 = 0.5;
            let decider = rand::distributions::Bernoulli::new(P_USE_TRIANGLE).unwrap();
            let lhs_triangles = lhs
                .triangles
                .clone()
                .into_iter()
                .filter(|_| decider.sample(&mut rand::thread_rng()));

            let rhs_triangles = rhs
                .triangles
                .clone()
                .into_iter()
                .filter(|_| decider.sample(&mut rand::thread_rng()));

            RgbaApproximator {
                reference_img: &lhs.reference_img,
                triangles: lhs_triangles.chain(rhs_triangles).collect(),
            }
        }

        fn calc_fitness(&self) -> u32 {
            /* use imageproc::integral_image::{integral_image, sum_image_pixels};
             * let integral = integral_squared_image(&image);
             * sum_image_pixels(&integral, 0, 0, w-1, h-1);
             */
            let (w, h) = self.reference_img.dimensions();
            let mut triangle_image = image::RgbaImage::new(w,h);
            for t in &self.triangles {
                triangle_image = imageproc::drawing::draw_convex_polygon(&mut triangle_image, t, image::Pixel::from_channels(120, 120, 120, 120));
            }
            let lhs_iter = triangle_image.pixels();
            let rhs_iter = self.reference_img.pixels();
            let fitness: f64 = lhs_iter.zip(rhs_iter).map(|(lhs_pixel, rhs_pixel)| {
                let lhs = lhs_pixel.channels();
                let rhs = rhs_pixel.channels();
                (lhs[0] as f64 * lhs[3] as f64/256. - rhs[0] as f64 * rhs[3] as f64/256.).powf(2.)
            }).sum();

            //self.triangles.len() as u32
            return fitness as u32;
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
