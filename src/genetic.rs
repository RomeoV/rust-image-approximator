extern crate image;
extern crate imageproc;
extern crate rand;

pub mod triangle_ops;

pub trait GeneticSample<'a, D> {
    fn new_from_rng(data: &'a D) -> Self;

    #[allow(non_snake_case)]
    fn combine_DNA(lhs: &Self, rhs: &Self) -> Self;

    fn maybe_mutate(&mut self, p_mutate: f64);

    fn calc_fitness(&self) -> u32;
}

pub mod rgba {
    use super::GeneticSample;
    use image::{GenericImage, ImageBuffer, Pixel, RgbaImage};
    use imageproc::definitions::Image;
    use imageproc::drawing::Blend;
    use rand::distributions::Bernoulli;
    use rand::distributions::Distribution;
    use rand::Rng;

    use super::triangle_ops::Point;
    use super::triangle_ops::Triangle;

    // Reimplements method from imageproc::drawing but with blend
    pub fn draw_convex_polygon_with_blend<I>(
        image: &I,
        poly: &[Point<i32>],
        color: I::Pixel,
    ) -> Image<I::Pixel>
    where
        I: GenericImage,
        I::Pixel: 'static,
    {
        let mut out = ImageBuffer::new(image.width(), image.height());
        out.copy_from(image, 0, 0);
        let mut out_blend = Blend(out);
        imageproc::drawing::draw_convex_polygon_mut(
            &mut out_blend,
            poly.iter()
                .map(|p| p.to_imageproc_point())
                .collect::<Vec<_>>()
                .as_slice(),
            color,
        );
        out_blend.0
    }

    pub struct RgbaApproximator<'a> {
        reference_img: &'a RgbaImage,
        pub triangles: Vec<Triangle>,
    }

    impl RgbaApproximator<'_> {
        fn random_triangle(img: &RgbaImage) -> Triangle {
            let (w, h) = img.dimensions();
            let mut rng = rand::thread_rng();
            let mut ret = [Point { x: 0, y: 0 }; 3];

            // no two points are allowed to overlap
            while ret[0] == ret[1] || ret[1] == ret[2] || ret[0] == ret[2] {
                let x = rng.gen_range(0, w - 1) as i32;
                let y = rng.gen_range(0, h - 1) as i32;
                ret[0] = Point { x: x, y: y }; //imageproc::drawing::Point::<i32>::new(x, y);
                ret[1] = //imageproc::drawing::Point::<i32>::new(
                    Point {
                    x: x + rng.gen_range(1, w / 2) as i32,
                    y: y + rng.gen_range(1, h / 2) as i32,
                    };
                ret[2] = Point {
                    //imageproc::drawing::Point::<i32>::new(
                    x: x + rng.gen_range(1, w / 2) as i32,
                    y: y + rng.gen_range(1, h / 2) as i32,
                };
                /* for p in &mut ret {
                *p = imageproc::drawing::Point::<i32>::new(
                        rng.gen_range(0, w-1) as i32,
                        rng.gen_range(0, h-1) as i32);
                }
                */
            }

            let color = image::Pixel::from_channels(
                rng.gen_range(0, 255),
                rng.gen_range(0, 255),
                rng.gen_range(0, 255),
                rng.gen_range(0, 255),
            );

            return Triangle(ret, color);
        }
    }

    impl<'a> GeneticSample<'a, RgbaImage> for RgbaApproximator<'a> {
        fn new_from_rng(img: &'a RgbaImage) -> RgbaApproximator<'a> {
            const N_INITIAL_TRIANGLES: usize = 50;
            RgbaApproximator {
                reference_img: img,
                triangles: (0..N_INITIAL_TRIANGLES)
                    .map(|_| RgbaApproximator::random_triangle(img))
                    .collect(),
            }
        }

        fn maybe_mutate(&mut self, p_mutate: f64) {
            let p_new_triangle = (2. * 0. + p_mutate) / 3.;
            let p_duplicate = (0. + p_mutate) / 2.;
            let p_shift = p_mutate;
            let p_resize = (p_mutate + 1.) / 2.;
            let p_rotate = p_mutate;
            let p_swap = p_mutate;
            let p_change_r = (0. + 2. * p_mutate) / 3.;
            let p_change_g = (0. + 2. * p_mutate) / 3.;
            let p_change_b = (0. + 2. * p_mutate) / 3.;
            let p_change_a = (0. + 2. * p_mutate) / 3.;

            let decide_new_triangle = Bernoulli::new(p_new_triangle).unwrap();
            let _decide_duplicate = Bernoulli::new(p_duplicate).unwrap();
            let decide_shift = Bernoulli::new(p_shift).unwrap();
            let decide_resize = Bernoulli::new(p_resize).unwrap();
            let decide_rotate = Bernoulli::new(p_rotate).unwrap();
            let decide_swap = Bernoulli::new(p_swap).unwrap();
            let _decide_change_r = Bernoulli::new(p_change_r).unwrap();
            let _decide_change_g = Bernoulli::new(p_change_g).unwrap();
            let _decide_change_b = Bernoulli::new(p_change_b).unwrap();
            let _decide_change_a = Bernoulli::new(p_change_a).unwrap();

            let mut rng = rand::thread_rng();
            let img_dims = self.reference_img.dimensions();

            if decide_new_triangle.sample(&mut rng) {
                self.triangles
                    .push(RgbaApproximator::random_triangle(&self.reference_img));
            }

            // resize
            self.triangles
                .iter_mut()
                .filter(|_| decide_resize.sample(&mut rng))
                .for_each(|t| t.resize());

            // shift
            self.triangles
                .iter_mut()
                .filter(|_| decide_shift.sample(&mut rng))
                .for_each(|t| t.shift(img_dims));
            
            // rotate
            self.triangles
                .iter_mut()
                .filter(|_| decide_rotate.sample(&mut rng))
                .for_each(|t| t.rotate());

            // swap triangle z layer
            for i in 0..self.triangles.len()-2 {
                if decide_swap.sample(&mut rng) {
                    let j = rng.gen_range(i+1, i+6);
                    if j < self.triangles.len() {
                        self.triangles.swap(i, j);
                    }
                }
            }

            for Triangle(p, _) in &mut self.triangles {
                if p[0] == p[1] || p[1] == p[2] || p[0] == p[2] {
                    p[1].x += 1;
                    p[2].y += 1;
                }
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
                triangles: lhs_triangles.zip(rhs_triangles).map(|(a, b)| vec![a, b].into_iter()).flatten().collect(),
            }
        }

        fn calc_fitness(&self) -> u32 {
            /* use imageproc::integral_image::{integral_image, sum_image_pixels};
             * let integral = integral_squared_image(&image);
             * sum_image_pixels(&integral, 0, 0, w-1, h-1);
             */
            let (w, h) = self.reference_img.dimensions();
            let mut triangle_image = image::RgbaImage::new(w, h);
            for Triangle(t, p) in &self.triangles {
                triangle_image = draw_convex_polygon_with_blend(&mut triangle_image, t, *p);
            }
            let lhs_iter = triangle_image.pixels();
            let rhs_iter = self.reference_img.pixels();
            let fitness: f64 = lhs_iter
                .zip(rhs_iter)
                .map(|(lhs_pixel, rhs_pixel)| {
                    let lhs = lhs_pixel.channels();
                    let rhs = rhs_pixel.channels();
                    (lhs[0] as f64 * lhs[3] as f64 / 256. - rhs[0] as f64 * rhs[3] as f64 / 256.)
                        .powf(2.)
                })
                .sum();

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
