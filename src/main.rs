/* Overall plan:
 * - [x] Create project
 * - [x] Create basic trait (mutability) and implementation; Log it
 * - [x] Initialize population, log it
 * - [ ] Implement fitness trait (image manipulation)
 * - [x] Sort population
 * - [x] Generate parents
 * - [ ] Create and implement DNA combination trait
 * - [x] Generate children
 * - [x] Close optimization loop
 * - [ ] DONE
 * Extras:
 * - [ ] Write docstrings
 * - [ ] File output
 * - [ ] Extensive testing and logging
 * - [ ] Benchmark
 */

extern crate image;
extern crate imageproc;
extern crate rand;

use rand::distributions::Distribution;
use rand::seq::SliceRandom;
use rand::Rng;

#[derive(Clone)]
struct Point {
    x: u32,
    y: u32,
}

trait GeneticSample {
    fn new_from_rng() -> Self;

    #[allow(non_snake_case)]
    fn combine_DNA(lhs: &Self, rhs: &Self) -> Self;

    fn maybe_mutate(&mut self, p_mutate: f64);
}

struct ImageApproximator {
    i: usize,
    triangles: Vec<[Point; 3]>,
}

impl GeneticSample for ImageApproximator {
    fn new_from_rng() -> ImageApproximator {
        ImageApproximator {
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

    fn combine_DNA(lhs: &ImageApproximator, rhs: &ImageApproximator) -> ImageApproximator {
        ImageApproximator {
            i: lhs.i + rhs.i,
            triangles: [
                lhs.triangles.clone().as_slice(),
                rhs.triangles.clone().as_slice(),
            ]
            .concat(),
        }
    }
}

fn main() {
    const N_SAMPLES: usize = 100;
    const N_BEST: usize = 20;
    const N_EPOCHS: usize = 10;
    const P_MUTATE: f64 = 0.1;

    let mut samples: Vec<ImageApproximator> = (0..N_SAMPLES)
        .map(|_| GeneticSample::new_from_rng())
        .collect();
    println!("Created {} initial samples.", N_SAMPLES);

    for epoch in 1..N_EPOCHS+1 {
        samples.sort_unstable_by_key(|s| s.i);
        samples.reverse();

        let best_samples = &samples[0..N_BEST];

        let parents: Vec<(&ImageApproximator, &ImageApproximator)> = (0..N_SAMPLES)
            .map(|_| {
                (
                    best_samples.choose(&mut rand::thread_rng()).unwrap(),
                    best_samples.choose(&mut rand::thread_rng()).unwrap(),
                )
            })
            .collect();

        samples = parents
            .iter()
            .map(|(lhs, rhs)| ImageApproximator::combine_DNA(lhs, rhs))
            .collect();

        samples.iter_mut().for_each(|s| s.maybe_mutate(P_MUTATE));

        println!(
            "In epoch {}, after recombination, a maximum of {} was reached.",
            epoch, samples[0].i
        );
    }
}
