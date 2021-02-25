extern crate clap;
extern crate flate2;


use std::fs::File;
use std::io::{BufReader, BufRead, Error};
use flate2::read::MultiGzDecoder;
use std::collections::{HashSet, HashMap};
use std::str::Chars;
use clap::{Arg, App};

#[derive (Default)]
struct DeMultiSeq {
    data : HashMap<String, HashMap<String, HashSet<String>>>,
    tolerance : u8,
    barcode_size: u32
}
impl DeMultiSeq {

    fn read_record(&self, reader:&mut BufReader<MultiGzDecoder<File>>) -> Result<String, i32> {

        /*
        Reads an expected 4-line record from a fastq file,
        exits if line lengths are malformed.
        */

        let mut header = String::new();
        let mut sequence = String::new();
        let mut plus = String::new();
        let mut qual = String::new();

        let eof = reader.read_line(&mut header);
        if eof.unwrap() == 0 {return Err(-1);}

        let eof = reader.read_line(&mut sequence);
        if eof.unwrap() == 0 {return Err(-2);}

        let eof = reader.read_line(&mut plus);
        if eof.unwrap() == 0 {return Err(-3);}

        let eof = reader.read_line(&mut qual);
        if eof.unwrap() == 0 {return Err(-4);}

        // println!("{}", sequence);
        Ok(sequence.trim().to_string())
    }

    fn char_slice(&self, chars:&mut Chars, size:u32) -> String {

        /*
        Pop the first <size> characters from a list of chars and return a string
        */

        let mut vec = Vec::<char>::new();

        let mut count = 0;
        while count < size {
            vec.push(chars.next().unwrap());
            count += 1
        }

        vec.iter().collect()
    }

    fn parse_r1(&self, seq:String) -> (String, String) {

        /*
        Parses the cell-barcode and UMI from R1

        expects :
            barcode : 16 basepairs
            UMI     : 10 basepairs

        */

        let mut chars = seq.chars();
        let barcode:String = self.char_slice(&mut chars, 16);
        let umi:String = self.char_slice(&mut chars, 10);
        return (barcode, umi)
    }

    fn parse_r2(&self, seq:String) -> String {

        /*
        Parses the Multiseq barcode from R2

        expects :
            Multiseq : 8 basepairs

        */

        let mut chars = seq.chars();
        let multiseq:String = self.char_slice(&mut chars, self.barcode_size);
        return multiseq
    }

    fn add_record(&mut self, cbs:&BarcodeSet, mbs:&BarcodeSet, bc:String, umi:String, multiseq:String) {

        /*
        Adds a barcode,umi,multiseq record to the existing dataset

        Validates that the barcode/multiseq exist OR are within a given hamming distance tolerance
        of the whitelist sets
        */

        let (cb_exists, bc) = cbs.exists(&bc, self.tolerance);
        let (ms_exists, multiseq) = mbs.exists(&multiseq, self.tolerance);

        if cb_exists & ms_exists {

            if !self.data.contains_key(&bc) {
                self.data.insert(bc.clone(), HashMap::new());
            }

            if !self.data.get(&bc).unwrap().contains_key(&multiseq) {
                self.data.get_mut(&bc)
                .unwrap()
                .insert(multiseq.clone(), HashSet::new());
            }

            if !self.data.get(&bc).unwrap().get(&multiseq).unwrap().contains(&umi) {
                self.data.get_mut(&bc)
                .unwrap()
                .get_mut(&multiseq)
                .unwrap()
                .insert(umi);
            }


        }

    }

    fn run(&mut self, r1:&str, r2:&str, cbs:&BarcodeSet, mbs:&BarcodeSet) -> Result<(), Error>{

        /*
        Iterates through both fastq records simulatenously and populates dataset.
        */

        let f1 = File::open(r1).expect("Error : File not found");
        let gz1 = MultiGzDecoder::new(f1);
        let mut buf1 = BufReader::new(gz1);

        let f2 = File::open(r2).expect("Error : File not found");
        let gz2 = MultiGzDecoder::new(f2);
        let mut buf2 = BufReader::new(gz2);

        let mut count = 0;
        loop {
            let seq1 = self.read_record(&mut buf1);
            let seq2 = self.read_record(&mut buf2);

            if seq1.is_err() | seq2.is_err() {
                break;
            }

            let (barcode, umi) = self.parse_r1(seq1.unwrap());
            let multiseq = self.parse_r2(seq2.unwrap());


            self.add_record(cbs, mbs, barcode, umi, multiseq);

            count += 1;
            if count % 10000 == 0 {
                eprintln!("Read pairs processed : {}", count);
            }
        }

        Ok(())
    }

    fn pprint(&self) {

        /*
        Format prints dataset into tab file.
        */

        println!("{}\t{}\t{}", "Barcode", "Multiseq", "nUMI");

        for barcode in self.data.keys() {

            for multiseq in self.data.get(barcode).unwrap().keys() {

                let num_umi = self.data.get(barcode)
                    .unwrap()
                    .get(multiseq)
                    .unwrap()
                    .len();

                println!("{}\t{}\t{}", barcode, multiseq, num_umi);
            }
        }
    }
}

#[derive (Default)]
struct BarcodeSet {
    barcodes : HashSet<String>
}
impl BarcodeSet {

    fn parse(&mut self, filename:&str) -> Result<(), Error> {

        let f = File::open(filename).expect("Error : Filename not found");
        let mut reader = BufReader::new(f);

        loop {
            let mut line = String::new();
            let is_eof = reader.read_line(&mut line)?;
            if is_eof == 0 {break;}

            // strip new line
            line = line.trim().to_string();

            self.barcodes.insert(line);
        }

        Ok(())
    }

    fn hamming_distance(&self, s1:&str, s2:&str) -> u8 {

        let mut distance = 0;

        let c1 = s1.chars();
        let c2 = s2.chars();

        for (i,j) in c1.zip(c2) {
            if i != j {
                distance += 1
            }
        }

        return distance;
    }

    fn pdist(&self, barcode:&str, tolerance:u8) -> (bool, String) {

        for seq in self.barcodes.iter() {
            let distance = self.hamming_distance(seq, barcode);
            if distance <= tolerance {
                return (true, seq.to_string())
            }
        }
        return (false, barcode.to_string())
    }

    fn exists(&self, barcode:&str, tolerance:u8) -> (bool, String) {

        // hash match
        if self.barcodes.contains(barcode) {
            return (true, barcode.to_string())
        }

        // pairwise hamming distance across all barcodes
        else if tolerance > 0 {

            // true if within tolerance else false
            let (within_tolerance, tolerant_seq) = self.pdist(barcode, tolerance);
            return (within_tolerance, tolerant_seq)

        }

        else {
            return (false, barcode.to_string())
        }
    }
}

fn main() {

    let matches = App::new("DeMultiSeq")
        .version("0.1")
        .author("Noam Teyssier")
        .about("Demultiplexes MultiSeq paired-end fastq sequences (expected barcode[16]-umi[10] & multiseq[8])")
        .arg(
            Arg::with_name("r1")
                .short("i")
                .long("read_1")
                .value_name("Read 1")
                .help("Read 1 of paired-end sequencing (expected fastq.gz format)")
                .takes_value(true)
                .required(true)
            )
        .arg(
            Arg::with_name("r2")
                .short("I")
                .long("read_2")
                .value_name("Read 2")
                .help("Read 2 of paired-end sequencing (expected fastq.gz format)")
                .takes_value(true)
                .required(true)
            )
        .arg(
            Arg::with_name("cell_barcodes")
                .short("c")
                .long("cell_barcodes")
                .value_name("Cell Barcodes")
                .help("White list of cell barcodes to match against")
                .takes_value(true)
                .required(true)
            )
        .arg(
            Arg::with_name("multiseq_barcodes")
                .short("m")
                .long("multiseq_barcodes")
                .value_name("Multiseq Barcodes")
                .help("White list of multiseq barcodes to match against")
                .takes_value(true)
                .required(true)
            )
        .arg(
            Arg::with_name("tolerance")
                .short("t")
                .long("tol")
                .value_name("tolerance")
                .help("Tolerance of hamming distance to implement on imperfect sequences (default = 0)")
                .takes_value(true)
                .default_value("0")
            )
        .arg(
            Arg::with_name("multiseq_barcode_size")
                .short("s")
                .long("size")
                .value_name("multiseq_barcode_size")
                .help("Size of multiseq barcode to extract from R2 (default = 8)")
                .takes_value(true)
                .default_value("8")
        )
        .get_matches();

    let r1 = matches.value_of("r1").unwrap();
    let r2 = matches.value_of("r2").unwrap();
    let fn_cell_barcodes = matches.value_of("cell_barcodes").unwrap();
    let fn_ms_barcodes = matches.value_of("multiseq_barcodes").unwrap();
    let tolerance = matches.value_of("tolerance").unwrap().parse::<u8>().unwrap();
    let barcode_size = matches.value_of("multiseq_barcode_size").unwrap().parse::<u32>().unwrap();

    let mut cbs = BarcodeSet::default();
    let _cbs_parse_result = cbs.parse(fn_cell_barcodes);

    let mut mbs = BarcodeSet::default();
    let _mbs_parse_result = mbs.parse(fn_ms_barcodes);

    let mut dms = DeMultiSeq {
        tolerance,
        barcode_size,
        ..Default::default()};
    let _dms_result = dms.run(r1, r2, &cbs, &mbs);
    dms.pprint();
}
