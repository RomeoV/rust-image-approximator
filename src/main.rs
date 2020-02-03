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

use rand::seq::SliceRandom;

mod genetic;
use crate::genetic::GeneticSample;

fn main() {
    const N_SAMPLES: usize = 100;
    const N_BEST: usize = 20;
    const N_EPOCHS: usize = 10;
    const P_MUTATE: f64 = 0.1;

    let reference_img = image::open("bliss.png").unwrap().to_rgba();

    let mut samples: Vec<genetic::rgba::RgbaApproximator> = (0..N_SAMPLES)
        .map(|_| GeneticSample::new_from_rng(&reference_img))
        .collect();
    println!("Created {} initial samples.", N_SAMPLES);

    for epoch in 1..N_EPOCHS + 1 {
        samples.sort_unstable_by_key(|s| s.calc_fitness());
        samples.reverse();

        let best_samples = &samples[0..N_BEST];

        let parents = (0..N_SAMPLES).map(|_| {
            (
                best_samples.choose(&mut rand::thread_rng()).unwrap(),
                best_samples.choose(&mut rand::thread_rng()).unwrap(),
            )
        });

        samples = parents
            .map(|(lhs, rhs)| GeneticSample::combine_DNA(lhs, rhs))
            .collect();

        samples.iter_mut().for_each(|s| s.maybe_mutate(P_MUTATE));

        println!(
            "In epoch {}, after recombination, a maximum of {} was reached.",
            epoch,
            samples[0].calc_fitness()
        );
    }
}
