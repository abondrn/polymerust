use polymerust::polymerust::*;

fn main() {

    let mut fc = FoldContext::new(
        "ACCGGUAGU",
        &EnergyParams::<i32> {
            stacking: -3,
            pair: -2,
            bulge_loop: 0,
            internal_loop: 0,
        },
    );

    fc.zuker();
    let secondary_structure = fc.traceback(0, fc.sequence.len() - 1);
    
    println!("RNA Sequence: {}", fc.sequence);
    println!("Predicted Secondary Structure: {}", secondary_structure);
}