use seq_io::parallel::ParallelRecordsets;

#[derive(PartialEq, Debug)]
struct SecondaryStructure {
    structures: Vec<Structure>,
    length: usize,
    exterior_loops_energy: i32,
    energy: i32,
    annotated: Vec<u8>,
    pair_table: Vec<i32>,
}

// five_prime_idx, three_prime_idx
#[derive(PartialEq, Debug, Clone, Copy)]
struct Region(usize, usize);

impl Region {
    fn five_prime_index(&self) -> usize {
        return self.0
    }

    fn three_prime_index(&self) -> usize {
        return self.1
    }
}

#[derive(PartialEq, Debug)]
struct Stem {
    top: StemStructure, // the stem that encompasses all the substructures
	substructures: Vec<StemStructure>,
}

#[derive(PartialEq, Debug)]
struct StemStructure {
    closing: Region,
    enclosed: Region,
	energy: i32, // free energy of the `StemStructure` in dcal / mol
}

#[derive(PartialEq, Debug)]
enum Structure {
    MultiLoop { stem: Stem, substructures: Vec<Structure>, energy: i32, substructures_energy: i32 },
    Hairpin { stem: Stem, single_stranded: Region, energy: i32, },
    SingleStranded(Region),
}


mod annotated_nt {
    pub const EXTERIOR_LOOP: u8 = b'e';
    pub const INTERIOR_LOOP: u8 = b'i';
    pub const HAIRPIN_LOOP: u8 = b'h';
    pub const MULTI_LOOP: u8 = b'm';
}

impl SecondaryStructure {
    pub fn annotated_str(&self) -> &str {
        unsafe {
            std::str::from_utf8_unchecked(&self.annotated)
        }
    }

    pub fn from_dot_bracket(structure: &str) -> Self {
        let bytes = structure.as_bytes();
        let mut ss = SecondaryStructure{
            structures: Vec::new(),
            length: structure.len(),
            exterior_loops_energy: 0,
            energy: 0,
            annotated: bytes.to_vec(),
            pair_table: pair_table_from_dot_bracket(bytes),
        };

        ss.parse_exterior();

        ss
    }

    pub fn has_pair(&self, idx: usize) -> bool {
        return self.pair_table[idx] != -1;
    }

    fn are_paired(&self, r: Region) -> bool {
        return self.has_pair(r.0) &&  self.pair_table[r.0] == r.1 as i32;
    }

    fn parse_exterior(&mut self) {
        let mut len_exterior_loop = 0;
        let mut i = 0;
        while i < self.length {
            if self.pair_table[i] == -1 {
                self.annotated[i] = annotated_nt::EXTERIOR_LOOP;
                len_exterior_loop += 1;
            } else {
                if len_exterior_loop != 0 {
                    let ssr = Region(i - len_exterior_loop, i - 1);
    
                    self.structures.push(Structure::SingleStranded(ssr));
                }

                let sub = self.parse_loop(i);
                self.structures.push(sub);
                i = self.pair_table[i] as usize;
            }

            i += 1;
        }

        if len_exterior_loop != 0 {
            let ssr = Region(self.length - len_exterior_loop, self.length - 1);
    
            self.structures.push(Structure::SingleStranded(ssr));
        }
    }

    fn parse_loop(&mut self, closing_5p_idx: usize) -> Structure {
        assert!(self.has_pair(closing_5p_idx));
        let closing_3p_idx = self.pair_table[closing_5p_idx] as usize;
        let mut closing = Region(closing_5p_idx, closing_3p_idx);
        
        let mut enclosed = Region(closing_5p_idx, closing_3p_idx);
        let mut stem = Stem {
            top: self.new_stem_structure(closing.clone(), enclosed.clone()),
            substructures: Vec::new(),
        };

        while enclosed.0 < enclosed.1 {
            
            loop {
                enclosed.0 += 1;
                if self.pair_table[enclosed.0] != -1 {
                    break;
                }
            }

            loop {
                enclosed.1 -= 1;
                if self.pair_table[enclosed.1] != -1 {
                    break;
                }
            }

            if self.pair_table[enclosed.1] as usize != enclosed.0 || enclosed.0 > enclosed.1 {
                break;
            }

            stem.substructures.push(self.new_stem_structure(closing.clone(), enclosed.clone()));
            closing = enclosed;
        }

        if closing.0 == stem.top.closing.0 {
            assert!(stem.substructures.len() == 0);
        } else {
            stem.top.enclosed = enclosed.clone();
        }

        if enclosed.0 > enclosed.1 {
            self.parse_hairpin(closing, stem)
        } else {
            self.parse_multi_loop(closing_5p_idx, stem)
        }
    }

    fn new_stem_structure(&mut self, closing: Region, enclosed: Region) -> StemStructure {
        assert!(self.are_paired(closing));
        //self.annotated[enclosed.0] = b'(';
        assert!(self.are_paired(enclosed));
        //self.annotated[enclosed.1] = b')';
        
        for i in closing.0+1 .. enclosed.0 {
            self.annotated[i] = annotated_nt::INTERIOR_LOOP;
        }

        for i in enclosed.1+1 .. closing.1 {
            self.annotated[i] = annotated_nt::INTERIOR_LOOP;
        }
        
        StemStructure {
            closing: closing,
            enclosed: enclosed,
            energy: 0,
        }
    }

    fn parse_hairpin(&mut self, closing: Region, stem: Stem) -> Structure {
        assert!(self.are_paired(closing));
        for i in closing.0+1 .. closing.1 {
            assert!(!self.has_pair(i));
            self.annotated[i] = annotated_nt::HAIRPIN_LOOP;
        }
        
        Structure::Hairpin { stem: stem, single_stranded: closing, energy: 0 }
    }

    fn skip_unpaired_in_multi_loop(&mut self, start: usize, closing_3p_idx: usize) -> usize {
        for i in start .. closing_3p_idx {
            if self.pair_table[i] != -1 {
                return i
            }
            self.annotated[i] = annotated_nt::MULTI_LOOP;
        }
        return closing_3p_idx
    }

    fn parse_multi_loop(&mut self, closing_5p_idx: usize, stem: Stem) -> Structure {
        let closing_3p_idx = self.pair_table[closing_5p_idx] as usize;
        let mut substructures: Vec<Structure> = Vec::new();

        let mut enclosed_5p_idx = self.skip_unpaired_in_multi_loop(closing_5p_idx+1, closing_3p_idx);
        let len_single_stranded_region = enclosed_5p_idx - closing_5p_idx - 1;

        //println!("{:?} {}", self.pair_table[closing_5p_idx ..= enclosed_5p_idx].to_vec(), len_single_stranded_region);

        if len_single_stranded_region != 0 {
            let ssr = Structure::SingleStranded(Region(closing_5p_idx + 1, enclosed_5p_idx - 1));
            substructures.push(ssr);
        }

        while enclosed_5p_idx < closing_3p_idx {
            let sub = self.parse_loop(enclosed_5p_idx);
            substructures.push(sub);
            enclosed_5p_idx = (self.pair_table[enclosed_5p_idx]+1) as usize;
            
            let next_enclosed_5p_idx = self.skip_unpaired_in_multi_loop(enclosed_5p_idx, closing_3p_idx);
            let len_single_stranded_region = next_enclosed_5p_idx - enclosed_5p_idx;
            if len_single_stranded_region != 0 {
                let ssr = Structure::SingleStranded(Region(enclosed_5p_idx, next_enclosed_5p_idx - 1));
                substructures.push(ssr);
            }

            enclosed_5p_idx = next_enclosed_5p_idx;
        }
        
        Structure::MultiLoop { stem: stem, substructures: substructures, energy: 0, substructures_energy: 0 }
    }
}

pub fn pair_table_from_dot_bracket(structure: &[u8]) -> Vec<i32> {
    let mut stack: Vec<usize> = Vec::new();
    let mut pt: Vec<i32> = vec![-1; structure.len()];

    for (i, &c) in structure.iter().enumerate() {
        match c {
            b'(' => stack.push(i),
            b')' => {
                if let Some(open) = stack.pop() {
                    pt[open] = i as i32;
                    pt[i] = open as i32;
                }
            },
            _ => {},
        }
    }

    pt
}


mod tests {
    use super::*;

	#[test]
	fn test_ss_from_db() {
		let ss = SecondaryStructure::from_dot_bracket("..((((...))))...((........))..");
        assert_eq!(ss.annotated_str(), "ee((((hhh))))eee((hhhhhhhh))ee");
	}

    #[test]
    fn test_ss_from_db_ml() {
        let ss = SecondaryStructure::from_dot_bracket("..(..((((...))))...((........))..).");
        assert_eq!(ss.annotated_str(), "ee(mm((((hhh))))mmm((hhhhhhhh))mm)e");
    }

    #[test]
    fn test_ss_from_db_hairpin_wo_ssr() {
        let ss = SecondaryStructure::from_dot_bracket("..((()))..");
        assert_eq!(ss.annotated_str(), "ee((()))ee");
    }

    #[test]
    fn test_ss_from_db_interior_bulge() {
        let ss = SecondaryStructure::from_dot_bracket("(((.).....).)");
        assert_eq!(ss.annotated_str(), "(((h)iiiii)i)");
    }

    #[test]
    fn test_ss_from_db_interior_1xn_loops() {
        let ss = SecondaryStructure::from_dot_bracket("(.(.(...(..).)..).)");
        assert_eq!(ss.annotated_str(), "(i(i(iii(hh)i)ii)i)");
    }

    #[test]
    fn test_ss_from_db_interior_2xn_loops() {
        let ss = SecondaryStructure::from_dot_bracket("(..(...(..(..)....)..)..)");
        assert_eq!(ss.annotated_str(), "(ii(iii(ii(hh)iiii)ii)ii)");
    }

    #[test]
    fn test_ss_from_db_interior_generic_loops() {
        let ss = SecondaryStructure::from_dot_bracket("(..(...(....().....)...)....)");
        assert_eq!(ss.annotated_str(), "(ii(iii(iiii()iiiii)iii)iiii)");
    }

    #[test]
    fn test_ss_from_db_interior_stem_wo_enclosed_pair() {
        let ss = SecondaryStructure::from_dot_bracket(".(.).");
        assert_eq!(ss.annotated_str(), "e(h)e");
    }
}