use std::collections::HashMap;
use ndarray::{prelude::*, Zip};

/******************************************************************************

This file defines structs needed to contain information about RNA free energy
params, and the types needed to specify which set of energy parameter values
to parse.

For more information on parsing the energy params, please see `parse.go`.

For more information on scaling the energy params, please see `scale.go`.

******************************************************************************/

// NUM_BASE_PAIRS is the number of distinguishable base pairs:
// CG, GC, GU, UG, AU, UA, & non-standard
const NUM_BASE_PAIRS: usize = 7;
// NUM_NT is the number of distinguishable nucleotides
const NUM_NT: usize = 4;
// MAX_LOOP_LEN is the maximum length of a loop (hairpin, or multi-loop)
const MAX_LOOP_LEN: usize = 30;
// ZeroCelsiusInKelvin is 0 deg Celsius in Kelvin
const ZeroCelsiusInKelvin: f64 = 273.15;

// EnergyParams contains all the energy parameters needed for the free energy
// calculations.
//
// The order of entries to access theses matrices always uses the closing pair or
// pairs as the first indices followed by the unpaired bases in 5' to 3' direction.
// For example, if we have a 2x2 interior loop:
// ```
// 		      5'-GAUA-3'
// 		      3'-CGCU-5'
// ```
// The closing pairs for the loops are GC and UA (not AU!), and the unpaired bases
// are (in 5' to 3' direction, starting at the first pair) A U C G.
// Thus, the energy for this sequence is:
// ```
// 	pairs:                    GC UA A  U  C  G
// 						interior2x2Loop[1][5][1][4][2][3]
// ```
// Note that this sequence is symmetric so the sequence is equivalent to:
// ```
// 					5'-UCGC-3'
// 					3'-AUAG-5'
// ```
// which means the energy of the sequence is equivalent to:
// ```
// 	pairs:                    UA GC C  G  A  U
// 						interior2x2Loop[5][1][1][4][2][3]
// ```
#[derive(Default)]
pub struct EnergyParams<T: Default> {

	// The matrix of free energies for stacked pairs, indexed by the two encoded closing
	// pairs. The list should be formatted as symmetric a `7*7`
	// (`NUM_BASE_PAIRS`*`NUM_BASE_PAIRS`) matrix,
	// conforming to the order explained above. As an example the stacked pair
	// ```
	// 					5'-GU-3'
	// 					3'-CA-5'
	// ```
	// corresponds to the entry StackingPair[1][4] (GC=1, AU=4) which should be
	// identical to StackingPair[4][1] (AU=4, GC=1).
	// size: [NUM_BASE_PAIRS][NUM_BASE_PAIRS]int
	StackingPair: Array2<T>,
	// Free energies of hairpin loops as a function of size. The list should
	// contain 31 (MAX_LOOP_LEN + 1) entries. Since the minimum size of a hairpin loop
	// is 3 and we start counting with 0, the first three values should be `inf` to
	// indicate a forbidden value.
	// size: [MAX_LOOP_LEN + 1]int
	HairpinLoop: Array1<T>,
	// Free energies of Bulge loops. Should contain 31 (MAX_LOOP_LEN + 1) entries,
	// the first one being `inf`.
	// size: [MAX_LOOP_LEN + 1]int
	Bulge: Array1<T>,

	// Free energies of interior loops. Should contain 31 (MAX_LOOP_LEN + 1) entries,
	// the first 4 being `inf` (since smaller loops are tabulated).
	//
	// This field was previous called internal_loop, but has been renamed to
	// interior loop to remain consistent with the names of other interior loops
	// size: [MAX_LOOP_LEN + 1]int
	InteriorLoop: Array1<T>,

	// Free energies for the interaction between the closing pair of an interior
	// loop and the two unpaired bases adjacent to the helix. This is a three
	// dimensional array indexed by the type of the closing pair and the two
	// unpaired bases. Since we distinguish 4 bases (A, C, G, & U) the list
	// contains `7*5*5` entries. The order is such that for example the mismatch
	//
	// ```
	// 							5'-CU-3'
	// 							3'-GC-5'
	// ```
	// corresponds to entry MismatchInteriorLoop[0][3][1] (CG=0, U=3, C=1).
	// Note that the matrix uses dims of length 5 (instead of 4) for the bases
	// since indexing the energy params uses the `NucleotideEncodedIntMap` map
	// which starts at 1 instead of 0.
	//
	// More information about mismatch energy: https://rna.urmc.rochester.edu/NNDB/turner04/tm.html
	// size: [NUM_BASE_PAIRS][NUM_NT + 1][NUM_NT + 1]int
	MismatchInteriorLoop: Array3<T>,
	// size: [NUM_BASE_PAIRS][NUM_NT + 1][NUM_NT + 1]int
	Mismatch1xnInteriorLoop: Array3<T>,
	// size: [NUM_BASE_PAIRS][NUM_NT + 1][NUM_NT + 1]int
	Mismatch2x3InteriorLoop: Array3<T>,
	// size: [NUM_BASE_PAIRS][NUM_NT + 1][NUM_NT + 1]int
	MismatchExteriorLoop: Array3<T>,
	// size: [NUM_BASE_PAIRS][NUM_NT + 1][NUM_NT + 1]int
	MismatchHairpinLoop: Array3<T>,
	// size: [NUM_BASE_PAIRS][NUM_NT + 1][NUM_NT + 1]int
	MismatchMultiLoop: Array3<T>,

	// Energies for the interaction of an unpaired base on the 5' side and
	// adjacent to a helix in multiloops and free ends (the equivalent of mismatch
	// energies in interior and hairpin loops). The array is indexed by the type
	// of pair closing the helix and the unpaired base and, therefore, forms a `7*5`
	// matrix. For example the dangling base in
	// ```
	// 						5'-GC-3'
	// 						3'- G-5'
	// ```
	// corresponds to entry DanglingEndsFivePrime[0][3] (CG=0, G=3).
	//
	// More information about dangling ends: https://rna.urmc.rochester.edu/NNDB/turner04/de.html
	// size: [NUM_BASE_PAIRS][NUM_NT + 1]int
	DanglingEndsFivePrime: Array2<T>,

	// Same as above for bases on the 3' side of a helix.
	// ```
	// 			       5'- A-3'
	// 			       3'-AU-5'
	// ```
	// corresponds to entry DanglingEndsThreePrime[4][1] (AU=4, A=1).
	// size: [NUM_BASE_PAIRS][NUM_NT + 1]int
	DanglingEndsThreePrime: Array2<T>,

	// Free energies for symmetric size 2 interior loops. `7*7*5*5` entries.
	// Example:
	// ```
	// 						5'-CUU-3'
	// 						3'-GCA-5'
	// ```
	// corresponds to entry Interior1x1Loop[0][4][4][2] (CG=0, AU=4, U=4, C=2),
	// which should be identical to Interior1x1Loop[4][0][2][4] (AU=4, CG=0, C=2, U=4).
	// size: [NUM_BASE_PAIRS][NUM_BASE_PAIRS][NUM_NT + 1][NUM_NT + 1]int
	Interior1x1Loop: Array4<T>,


	// Free energies for 2x1 interior loops, where 2 is the number of unpaired
	// nucleotides on the larger 'side' of the interior loop and 1 is the number of
	// unpaired nucleotides on the smaller 'side' of the interior loop.
	// `7*7*5*5*5` entries.
	// Example:
	// ```
	// 						5'-CUUU-3'
	// 						3'-GC A-5'
	// ```
	// corresponds to entry Interior2x1Loop[0][4][4][4][2] (CG=0, AU=4, U=4, U=4, C=2).
	// Note that this matrix is always accessed in the 5' to 3' direction with the
	// larger number of unpaired nucleotides first.
	// size: [NUM_BASE_PAIRS][NUM_BASE_PAIRS][NUM_NT + 1][NUM_NT + 1][NUM_NT + 1]int
	Interior2x1Loop: Array5<T>,

	// Free energies for symmetric size 4 interior loops. `7*7*5*5*5*5` entries.
	// Example:
	// ```
	// 						5'-CUAU-3'
	// 						3'-GCCA-5'
	// ```
	// corresponds to entry Interior2x2Loop[0][4][4][1][2][2] (CG=0, AU=4, U=4, A=1, C=2, C=2),
	// which should be identical to Interior2x2Loop[4][0][2][2][1][4] (AU=4, CG=0, C=2, C=2, A=1, U=4).
	// size: [NUM_BASE_PAIRS][NUM_BASE_PAIRS][NUM_NT + 1][NUM_NT + 1][NUM_NT + 1][NUM_NT + 1]int
	Interior2x2Loop: Array6<T>,

	// LogExtrapolationConstant is used to scale energy parameters for hairpin,
	// bulge and interior loops when the length of the loop is greather than
	// `MAX_LOOP_LEN`.
	LogExtrapolationConstant: f64,
	MultiLoopUnpairedNucleotideBonus: T,
    MultiLoopClosingPenalty: T,
    TerminalAUPenalty: T,
	MultiLoopIntern: T,

	// Some tetraloops particularly stable tetraloops are assigned an energy
	// bonus. For example:
	// ```
	// 	GAAA    -200
	// ```
	// assigns a bonus energy of -2 kcal/mol to TetraLoop containing
	// the sequence GAAA.
	TetraLoop: HashMap<String, T>,
	TriLoop:   HashMap<String, T>,
	HexaLoop:  HashMap<String, T>,

    // TODO: what are these?
    Ninio: i32,
	MaxNinio: i32,
}

// BasePairType is a type to hold information of the type of a base pair.
// The chosen numbers denote where the energy paramater values can be found
// for the base pair type in the `EnergyParams` energy parameter matrices.
pub enum BasePairType {
	CG = 0,
	GC = 1,
	GU = 2,
	UG = 3,
	AU = 4,
	UA = 5,
	// NoPair denotes that two bases don't pair
	NoPair = -1,
}

// encode_base_pair returns the type of a base pair encoded as an `BasePairType`,
// which is used to access loop energy paramater values in the `EnergyParams` struct.
// by the base pair closing the loop.
pub fn encode_base_pair(five_prime_base: &char, three_prime_base: &char) -> BasePairType {
	match (five_prime_base, three_prime_base) {
        ('C', 'G') => BasePairType::CG,
        ('G', 'C') => BasePairType::GC,
        ('G', 'U') => BasePairType::GU,
        ('U', 'G') => BasePairType::UG,
        ('A', 'U') => BasePairType::AU,
        ('U', 'A') => BasePairType::UA,
        _ => BasePairType::NoPair,
    }
}

pub fn encode_nt(nt: char) -> usize {
    match nt {
        'A' => 1,
        'C' => 2,
        'G' => 3,
        'U' => 4,
	    _ => panic!("Invalid nucleotide {nt}"),
    }
}

use constcat;

// EnergyParamsSet is used to specify the set of RNA free energy paramaters to
// parse.
pub mod energy_params_set {
	const ROOT: &str = "src/energy_params/param_files/";

	// Langdon2018 specifies the set of RNA energy parameters obtained from
	// Grow and Graft Genetic Programming (GGGP) as published
	// in Langdon et al. 2018, "Evolving Better RNAfold
	// Structure Prediction", EuroGP-2018, M. Castelli,
	// L. Sekanina, M. Zhang Eds., Parma. 4-6 April 2018
    pub const Langdon2018: &str = constcat::concat!(ROOT, "rna_langdon2018.par");

	// Andronescu2007 specifies the set of RNA energy parameters obtained from
	// Andronescu M, Condon A, Hoos HH, Mathews DH, Murphy KP. Efficient
	// parameter estimation for RNA secondary structure prediction.
	// Bioinformatics. 2007 Jul 1;23(13):i19-28
	pub const Andronescu2007: &str = constcat::concat!(ROOT, "rna_andronescu2007.par");

	// Turner2004 specifies the set of RNA energy parameters obtained from
	// Mathews DH, Disney MD, Childs JL, Schroeder SJ, Zuker M, Turner DH.
	// Incorporating chemical modification constraints into a dynamic
	// programming algorithm for prediction of RNA secondary structure.
	// Proc Natl Acad Sci U S A. 2004;101(19):7287-7292.
	pub const Turner2004: &str = constcat::concat!(ROOT, "rna_turner2004.par");

	// Turner1999 specifies the set of RNA energy parameters obtained from
	// Mathews DH, Sabina J, Zuker M, Turner DH. Expanded sequence dependence
	// of thermodynamic parameters improves prediction of RNA secondary
	// structure. J Mol Biol. 1999 May 21;288(5):911-40.
	pub const Turner1999: &str = constcat::concat!(ROOT, "rna_turner1999.par");
}

// rawEnergyParams is an un-exported intermediate struct used to store parsed
// energy param values. It is converted into usable energy params by the
// function `scaleByTemperature()`. For more information on the fields of
// this struct, read the documentation of the fields of the `EnergyParams`
// struct.
#[derive(Default)]
pub struct RawEnergyParams {
    for37C:                   EnergyParams<i32>,
    forEnthalpy:              EnergyParams<i32>,

	maxNinio:                 i32,
    logExtrapolationConstant: f64,
}

// The constant used to extrapolate energy values when length of loop > `MAX_LOOP_LEN`
const defaultLogExtrapolationConstantAt37C: f64 = 107.856;

impl RawEnergyParams {
    fn set_multiloop_params(&mut self, parsed_multi_loop_params: Vec<i32>) {
        self.for37C.MultiLoopUnpairedNucleotideBonus = parsed_multi_loop_params[0];
        self.forEnthalpy.MultiLoopUnpairedNucleotideBonus = parsed_multi_loop_params[1];
        self.for37C.MultiLoopClosingPenalty = parsed_multi_loop_params[2];
        self.forEnthalpy.MultiLoopClosingPenalty = parsed_multi_loop_params[3];
        self.for37C.MultiLoopIntern = parsed_multi_loop_params[4];
        self.forEnthalpy.MultiLoopIntern = parsed_multi_loop_params[5];
    }
    
    fn set_ninio_params(&mut self, parsed_ninio_params: Vec<i32>) {
        self.for37C.Ninio = parsed_ninio_params[0];
        self.forEnthalpy.Ninio = parsed_ninio_params[1];
        self.maxNinio = parsed_ninio_params[2];
    }
    
    fn set_misc_params(&mut self, parsed_misc_params: Vec<f64>) {
        self.for37C.TerminalAUPenalty = parsed_misc_params[2] as i32;
        self.forEnthalpy.TerminalAUPenalty = parsed_misc_params[3] as i32;
        // parameter files may or may not include the log extrapolation constant
        if parsed_misc_params.len() > 4 {
            self.logExtrapolationConstant = parsed_misc_params[5];
        } else {
            // no log extrapolation constant so set to default
            self.logExtrapolationConstant = defaultLogExtrapolationConstantAt37C;
        }
    }
}


/******************************************************************************

This file defines the functions needed to scale free energy params by a
specified temperature.

******************************************************************************/

// The temperature at which the energy parameters have been measured at
const measurementTemperatureInCelsius: f64 = 37.0;

impl RawEnergyParams {

    // scaleByTemperature scales energy paramaters according to the specificed temperatue.
    // See `rescaleDg` for more information of how energy values are rescaled.

    pub fn scale_by_temperature(&self, temperatureInCelsius: f64) -> EnergyParams<i32> {
        let for37C = &self.for37C;
        let forEnthalpy = &self.forEnthalpy;
        
        let clip_le_zero = |x| if x > 0 { 0 } else { x };

        EnergyParams{
            LogExtrapolationConstant:         rescale_dG(self.logExtrapolationConstant, 0.0, temperatureInCelsius),
            TerminalAUPenalty:                rescale_dG_round(for37C.TerminalAUPenalty, forEnthalpy.TerminalAUPenalty, temperatureInCelsius),
            MultiLoopUnpairedNucleotideBonus: rescale_dG_round(for37C.MultiLoopUnpairedNucleotideBonus, forEnthalpy.MultiLoopUnpairedNucleotideBonus, temperatureInCelsius),
            MultiLoopClosingPenalty:          rescale_dG_round(for37C.MultiLoopClosingPenalty, forEnthalpy.MultiLoopClosingPenalty, temperatureInCelsius),
            Ninio:                            rescale_dG_round(for37C.Ninio, forEnthalpy.Ninio, temperatureInCelsius),
            MaxNinio:                         self.maxNinio,

            HairpinLoop: rescale_dG_array(&for37C.HairpinLoop, &forEnthalpy.HairpinLoop, temperatureInCelsius),
            Bulge: rescale_dG_array(&for37C.Bulge, &forEnthalpy.Bulge, temperatureInCelsius),
            InteriorLoop: rescale_dG_array(&for37C.InteriorLoop, &forEnthalpy.InteriorLoop, temperatureInCelsius),
    
            MultiLoopIntern: rescale_dG_round(self.for37C.MultiLoopIntern, self.forEnthalpy.MultiLoopIntern, temperatureInCelsius),
            
            TetraLoop: rescale_dG_loop(&for37C.TetraLoop, &forEnthalpy.TetraLoop, temperatureInCelsius),
            TriLoop: rescale_dG_loop(&for37C.TriLoop, &forEnthalpy.TriLoop, temperatureInCelsius),
            HexaLoop: rescale_dG_loop(&for37C.HexaLoop, &forEnthalpy.HexaLoop, temperatureInCelsius),
        
            StackingPair: rescale_dG_array(&for37C.StackingPair, &forEnthalpy.StackingPair, temperatureInCelsius),
        
            /* mismatches */
            MismatchInteriorLoop: rescale_dG_array(&for37C.MismatchInteriorLoop, &forEnthalpy.MismatchInteriorLoop, temperatureInCelsius),
            MismatchHairpinLoop: rescale_dG_array(&for37C.MismatchHairpinLoop, &forEnthalpy.MismatchHairpinLoop, temperatureInCelsius),
            Mismatch1xnInteriorLoop: rescale_dG_array(&for37C.Mismatch1xnInteriorLoop, &forEnthalpy.Mismatch1xnInteriorLoop, temperatureInCelsius),
            Mismatch2x3InteriorLoop: rescale_dG_array(&for37C.Mismatch2x3InteriorLoop, &forEnthalpy.Mismatch2x3InteriorLoop, temperatureInCelsius),
        
            MismatchMultiLoop: rescale_dG_array(&for37C.MismatchMultiLoop, &forEnthalpy.MismatchMultiLoop, temperatureInCelsius).mapv_into(clip_le_zero),
            MismatchExteriorLoop: rescale_dG_array(&for37C.MismatchExteriorLoop, &forEnthalpy.MismatchExteriorLoop, temperatureInCelsius).mapv_into(clip_le_zero),
        
            /* dangling ends energies */
            DanglingEndsFivePrime: rescale_dG_array(&for37C.DanglingEndsFivePrime, &forEnthalpy.DanglingEndsFivePrime, temperatureInCelsius).mapv_into(clip_le_zero),
            DanglingEndsThreePrime: rescale_dG_array(&for37C.DanglingEndsThreePrime, &forEnthalpy.DanglingEndsThreePrime, temperatureInCelsius).mapv_into(clip_le_zero),
        
            /* interior 1x1 loops */
            Interior1x1Loop: rescale_dG_array(&for37C.Interior1x1Loop, &forEnthalpy.Interior1x1Loop, temperatureInCelsius),
        
            /* interior 2x1 loops */
            Interior2x1Loop: rescale_dG_array(&for37C.Interior2x1Loop, &forEnthalpy.Interior2x1Loop, temperatureInCelsius),
        
            /* interior 2x2 loops */
            Interior2x2Loop: rescale_dG_array(&for37C.Interior2x2Loop, &forEnthalpy.Interior2x2Loop, temperatureInCelsius),
        }
    }

}


// Rescale Gibbs free energy according to the equation dG = dH - T * dS
// where dG is the change in Gibbs free energy
// 			dH is the change in enthalpy
// 			dS is the change in entropy
// 			T is the temperature
// more information: https://chemed.chem.purdue.edu/genchem/topicreview/bp/ch21/gibbs.php
fn rescale_dG(dG: f64, dH: f64, temperatureInCelsius: f64) -> f64 {
	// if temperate == measurementTemperatureInCelsius then below calculation will
	// always return dG. So we save some computation with this check.
	if temperatureInCelsius == measurementTemperatureInCelsius {
		dG
	} else {
        let measurementTemperatureInKelvin = measurementTemperatureInCelsius + ZeroCelsiusInKelvin;
        let temperatureInKelvin = temperatureInCelsius + ZeroCelsiusInKelvin;
        let temp: f64 = temperatureInKelvin / measurementTemperatureInKelvin;
        let dS = dH - dG;

        dH - dS*temp
    }
}

fn rescale_dG_round(dG: i32, dH: i32, temperatureInCelsius: f64) -> i32 {
	rescale_dG(dG as f64, dH as f64, temperatureInCelsius) as i32
}

fn rescale_dG_array<D>(dG: &Array<i32, D>, dH: &Array<i32, D>, temperatureInCelsius: f64) -> Array<i32, D> 
where D: Dimension {
	let mut result = dG.clone();

    // Apply a function element-wise on the two arrays using Zip
    Zip::from(&mut result).and(dH)
        .for_each(|g, &h| *g = rescale_dG_round(*g, h, temperatureInCelsius));

	result
}

fn rescale_dG_loop(loops1: &HashMap<String, i32>, loops2: &HashMap<String, i32>, temperatureInCelsius: f64) -> HashMap<String, i32> {
    let mut loops3 = HashMap::new();
    for l in loops1.keys() {
        loops3.insert(l.to_string(), rescale_dG_round(loops1[l], loops2[l], temperatureInCelsius));
    }
    loops3
}


/******************************************************************************

This file defines the structs and functions needed to parse RNA energy
parameters as specified by the `RNAfold parameter file v2.0` file format from
[ViennaRNA](https://github.com/ViennaRNA/ViennaRNA#energy-parameters).

******************************************************************************/

use std::io::{self, BufRead};

type AnyErr = Box<dyn std::error::Error>;

// The main function where the parsing of energy param files starts from.
//
// Note that the `RNAfold parameter file v2.0` file format has a flaw.
// Since the nucleotide encoded int map (`NucleotideEncodedIntMap`) starts at
// 1 (instead of 0), all parameters specified for nucleotides in the energy
// params file have an extra first dimension even though the values can never
// be accessed.
// For example, for `mismatch_exterior`, we parse parameters of the dimesions
// `NUM_BASE_PAIRS` x `NUM_NT+1` x
// `NUM_NT+1`.
// The extra `+ 1` comes from the fact that `NucleotideEncodedIntMap` starts
// at 1 instead of 0. The first dimension that we parse doesn't have any
// useful information since it can't be accessed, but we still have to parse
// it since the parameter file includes the dimension. Thus, in the future,
// if we'd like to make `NucleotideEncodedIntMap` start at 0, we have to
// include a function that removes the first dimension whenever a
// `NUM_NT+1` dimension is parsed.
//
// Note that for some reason int22 and int22_enthalpies parameters have the
// right dimension sizes for the nucleotide dimensions, so we have to add
// a pre offset to those nucleotide dimensions.
pub fn parse_param_file(filename: &str) -> Result<RawEnergyParams, AnyErr> {
	use std::fs::File;
	use regex::Regex;

    let file = File::open(filename)?;
    let reader = io::BufReader::new(file);

	let mut current_header = "".to_owned();
	let mut current_matrix = Vec::new();
	let mut params: RawEnergyParams = Default::default();

	let mut inTriTetraHexaLoop = false;
	let mut energyLoop: HashMap<String, i32> = HashMap::new();
	let mut enthalpyLoop: HashMap<String, i32> = HashMap::new();

	let comment_pattern = Regex::new(r"/\*.*?\*/").unwrap();

	//headerLine, lineAvailable := readLine(scanner)
	//if lineAvailable && headerLine != "## RNAfold parameter file v2.0" {
	//	panic("Missing header line in file.\nThis file is not of the RNAfold parameter file v2.0 format.\n")
	//}

    for untrimmed in reader.lines() {
		let line = untrimmed?.trim().to_string();

        if line.starts_with("#") {
            // Line starts with '#', indicating a new matrix
			if inTriTetraHexaLoop {
				match current_header.as_str() {
					"Triloops" => {
						params.for37C.TriLoop = energyLoop;
						params.forEnthalpy.TriLoop = enthalpyLoop;
					},
					"Tetraloops" => {
						params.for37C.TetraLoop = energyLoop;
						params.forEnthalpy.TetraLoop = enthalpyLoop;
					},
					"Hexaloops" => {
						params.for37C.HexaLoop = energyLoop;
						params.forEnthalpy.HexaLoop = enthalpyLoop;
					},
					_ => {},
				};

				inTriTetraHexaLoop = false;
				energyLoop = HashMap::new();
				enthalpyLoop = HashMap::new();
			} else if current_header != "" && !current_matrix.is_empty() {
				store_matrix(&mut params, &current_header, current_matrix);
				current_matrix = Vec::new();
			}

			current_header = line[2..].to_owned();
			inTriTetraHexaLoop = match current_header.as_str() {
				"Triloops" => true,
				"Tetraloops" => true,
				"Hexaloops" => true,
				"END" => { break; }
				_ => false,
			};
			continue;
		// parse tri-, tetra-, hexa- loop values
        }
		
		let trimmed = comment_pattern.replace_all(&line, "").trim().to_owned();

		if  trimmed.len() == 0 {
			// encountered a line with only comments so continue parsing
			continue;
		}
		
		if inTriTetraHexaLoop {
			let values: Vec<&str> = trimmed.split_whitespace().collect();
			if values.len() != 3 {
				let error_message = format!("encountered incorrect number of values. expected 3, got {}", values.len());
				return Err(error_message.into());
			}
			let (loopKmer, energy, enthalpy) = (values[0], parse_int(values[1]), parse_int(values[2]));
			energyLoop.insert(loopKmer.to_string(), energy?);
			enthalpyLoop.insert(loopKmer.to_string(), enthalpy?);
		} else {
            let row = trimmed.split_whitespace();
			match current_header.as_str() {
				"Misc" => {
					params.set_misc_params(row.filter_map(|tok| parse_float(tok).ok()).collect())
				},
				_ => {
					current_matrix.extend(row.filter_map(|tok| parse_int(tok).ok()));
				}
			}
        }
    }

    Ok(params)
}

const MISMATCH_SHAPE: [usize; 3] = [NUM_BASE_PAIRS, NUM_NT+1, NUM_NT+1];
const LOOP_SHAPE: [usize; 1] = [MAX_LOOP_LEN+1];
const BP_SHAPE: [usize; 2] = [NUM_BASE_PAIRS, NUM_BASE_PAIRS];
const DANGLE_SHAPE: [usize; 2] = [NUM_BASE_PAIRS, NUM_NT+1];

use ndarray::{OwnedRepr, StrideShape, ShapeError};
fn reshape<A, D, Sh>(shape: Sh, v: Vec<A>) -> Result<ArrayBase<OwnedRepr<A>, D>, ShapeError>
where Sh: Into<StrideShape<D>>, D: Dimension {
	Array::from_shape_vec(shape, v)
}


fn store_field<T>(params: &mut EnergyParams<T>, header: &str, vec: Vec<T>)
where T: Default {
	match header {
		"stack" | "stack_enthalpies" => {
			params.StackingPair = reshape(BP_SHAPE, vec).unwrap();
		},
		"hairpin" | "hairpin_enthalpies" => {
			params.HairpinLoop = reshape(LOOP_SHAPE, vec).unwrap();
		},
		"bulge" | "bulge_enthalpies" => {
			params.Bulge = reshape(LOOP_SHAPE, vec).unwrap();
		},
		"interior" | "interior_enthalpies" => {
			params.InteriorLoop = reshape(LOOP_SHAPE, vec).unwrap();
		},

		"mismatch_exterior" | "mismatch_exterior_enthalpies" => {
			params.MismatchExteriorLoop = reshape(MISMATCH_SHAPE, vec).unwrap();
		},
		"mismatch_hairpin" | "mismatch_hairpin_enthalpies" => {
			params.MismatchHairpinLoop = reshape(MISMATCH_SHAPE, vec).unwrap();
		},
		"mismatch_interior" | "mismatch_interior_enthalpies" => {
			params.MismatchInteriorLoop = reshape(MISMATCH_SHAPE, vec).unwrap();
		},
		"mismatch_interior_1n" | "mismatch_interior_1n_enthalpies" => {
			params.Mismatch1xnInteriorLoop = reshape(MISMATCH_SHAPE, vec).unwrap();
		},
		"mismatch_interior_23" | "mismatch_interior_23_enthalpies" => {
			params.Mismatch2x3InteriorLoop = reshape(MISMATCH_SHAPE, vec).unwrap();
		},
		"mismatch_multi" | "mismatch_multi_enthalpies" => {
			params.MismatchMultiLoop = reshape(MISMATCH_SHAPE, vec).unwrap();
		},
		
		"int11" | "int11_enthalpies" => {
			params.Interior1x1Loop = reshape((NUM_BASE_PAIRS, NUM_BASE_PAIRS, NUM_NT+1, NUM_NT+1), vec).unwrap();
		},
		"int21" | "int21_enthalpies" => {
			params.Interior2x1Loop = reshape((NUM_BASE_PAIRS, NUM_BASE_PAIRS, NUM_NT+1, NUM_NT+1, NUM_NT+1), vec).unwrap();
		},
		
		"int22" | "int22_enthalpies" => {
			params.Interior2x2Loop = reshape((NUM_BASE_PAIRS-1, NUM_BASE_PAIRS-1, NUM_NT, NUM_NT, NUM_NT, NUM_NT), vec).unwrap();
		},

		"dangle5" | "dangle5_enthalpies" => {
			params.DanglingEndsFivePrime = reshape(DANGLE_SHAPE, vec).unwrap();
		},
		"dangle3" | "dangle3_enthalpies" => {
			params.DanglingEndsThreePrime = reshape(DANGLE_SHAPE, vec).unwrap();
		},
		
		_ => {},
	}
}


fn store_matrix(params: &mut RawEnergyParams, current_header: &str, current_matrix: Vec<i32>) {
	match current_header {

		"NINIO" => {
			// len = 3
			params.set_ninio_params(current_matrix);
		},

		"ML_params" => {
			// len = 6
			params.set_multiloop_params(current_matrix);
		},

		_ => {
			let side = if current_header.ends_with("_enthalpies") {
				&mut params.forEnthalpy
			} else {
				&mut params.for37C
			};

			store_field(side, current_header, current_matrix);
		},
	}
}

use std::num::{ParseIntError, ParseFloatError};

fn parse_int(token: &str) -> Result<i32, ParseIntError> {
	if token == "INF" {
		Ok(i32::MAX)
	} else {
		token.parse()
	}
}

fn parse_float(token: &str) -> Result<f64, ParseFloatError> {
	if token == "INF" {
		Ok(f64::MAX)
	} else {
		token.parse()
	}
}


/*****************************************************************************
End Section: Miscellaneous Parsing Functions
*****************************************************************************/


/*
func addOffset1Dim(values []int, dim1Offset int, offsetType offsetType) (ret []int) {
	switch offsetType {
	case preOffset:
		ret = prependInfsToSlice(values, dim1Offset)
	case postOffset:
		ret = appendInfsToSlice(values, dim1Offset)
	}
	return
}

// adds `length` `inf`s to the front of a slice
func prependInfsToSlice(slice []int, length int) []int {
	return append(newIntSlice(inf, length), slice...)
}

func appendInfsToSlice(slice []int, length int) []int {
	return append(slice, newIntSlice(inf, length)...)
}
*/


#[cfg(test)]
mod tests {
    use super::*;

	#[test]
	fn test_parse_param_files() {
		println!("{}", energy_params_set::Langdon2018);
		parse_param_file(energy_params_set::Langdon2018).unwrap();
		parse_param_file(energy_params_set::Andronescu2007).unwrap();
		parse_param_file(energy_params_set::Turner2004).unwrap();
		parse_param_file(energy_params_set::Turner1999).unwrap();
	}

	#[test]
	fn test_load() {
		parse_param_file(energy_params_set::Langdon2018).unwrap().scale_by_temperature(40.0);
	}
}