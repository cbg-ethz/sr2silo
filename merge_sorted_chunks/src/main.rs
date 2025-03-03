use clap::Command;
use serde_json::Value;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::env;
use std::fs::File;
use std::io::{stdin, stdout, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use zstd::stream::Decoder;
use zstd::Encoder;

// Wrapper struct to allow sorting JSON values in a min-heap
#[derive(Eq, PartialEq)]
struct HeapEntry {
    sort_field: i64,
    value: Value,
    index: usize,
}

impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        other.sort_field.cmp(&self.sort_field) // Reverse order to make BinaryHeap a min-heap
    }
}

impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// Merging function that reads from readers and writes to any object implementing `Write`
fn merge_sorted_readers<R: Read, W: Write>(
    output: &mut W,
    sorted_readers: &mut [BufReader<R>],
    sort_field_name: &String,
) -> std::io::Result<()> {
    let mut heap = BinaryHeap::new();

    // Store an iterator for each reader
    let mut reader_iters: Vec<_> = sorted_readers.iter_mut().map(|r| r.lines()).collect();

    // Initialize heap with the first line from each reader
    for (index, iter) in reader_iters.iter_mut().enumerate() {
        if let Some(Ok(line)) = iter.next() {
            let json: Value = serde_json::from_str(&line)?;
            heap.push(HeapEntry {
                sort_field: json[sort_field_name]
                    .as_i64()
                    .expect("the specified sort_column is not of type i64"),
                value: json,
                index,
            });
        }
    }

    let mut writer = BufWriter::new(output);
    while let Some(HeapEntry {
        sort_field: _sort_field,
        value,
        index,
    }) = heap.pop()
    {
        writeln!(writer, "{}", value)?;
        if let Some(Ok(line)) = reader_iters[index].next() {
            let json: Value = serde_json::from_str(&line)?;
            heap.push(HeapEntry {
                sort_field: json[sort_field_name]
                    .as_i64()
                    .expect("the specified sort_column is not of type i64"),
                value: json,
                index,
            });
        }
    }

    Ok(())
}

fn main() -> std::io::Result<()> {
    let matches = Command::new("Merge Sorted Chunks")
        .version("1.0")
        .author("Alexander Taepper")
        .arg(
            clap::Arg::new("sort_field")
                .long("sort-field")
                .value_name("FIELD")
                .help("Specifies the field in the json to sort by")
                .required(true),
        )
        .arg(
            clap::Arg::new("tmp_dir")
                .long("tmp-directory")
                .value_name("PATH")
                .help("Specifies a directory for placing intermediate files"),
        )
        .get_matches();

    let tmp_dir = matches
        .get_one::<String>("tmp_dir")
        .map(PathBuf::from)
        .unwrap_or(env::temp_dir());

    let sort_field = matches.get_one("sort_field").unwrap();

    let reader = stdin();
    let reader = reader.lock();

    let mut input_files: Vec<PathBuf> = Vec::new();

    let mut merge_iteration = 0;
    let mut chunk_id = 0;

    let mut next_input_files = Vec::new();

    for file in BufReader::new(reader).lines() {
        let file = file?;
        if input_files.len() == 512 {
            let mut sorted_readers = open_zst_files(&mut input_files);

            let file_name = Path::join(
                &tmp_dir,
                format!("merged_chunk_{}_{}.ndjson.zst", merge_iteration, chunk_id),
            );

            let file = File::create(file_name.clone())?;
            let mut encoder = Encoder::new(file, 3)?;

            merge_sorted_readers(&mut encoder, &mut sorted_readers, sort_field)?;

            encoder.finish()?;

            next_input_files.push(file_name);
            chunk_id += 1;
            input_files = Vec::new();
        }

        input_files.push(PathBuf::from(file));
    }
    input_files = next_input_files;
    merge_iteration += 1;

    if input_files.len() == 0 {
        panic!("No input files received");
    }

    while input_files.len() > 1 {
        let mut next_input_files: Vec<PathBuf> = Vec::new();
        next_input_files.reserve(input_files.len() / 512);
        for (chunk_id, chunk) in input_files.chunks_mut(512).enumerate() {
            let mut sorted_readers = open_zst_files(chunk);

            let file_name = Path::join(
                &tmp_dir,
                format!("merged_chunk_{}_{}.ndjson.zst", merge_iteration, chunk_id),
            );
            let file = File::create(file_name.clone())?;
            let mut encoder = Encoder::new(file, 3)?;

            merge_sorted_readers(&mut encoder, &mut sorted_readers, sort_field)?;

            encoder.finish()?;

            next_input_files.push(file_name);
        }
        merge_iteration += 1;
        input_files = next_input_files;
    }

    println!("Merged chunks to files \"{:?}\"", input_files);

    let mut sorted_readers = open_zst_files(&mut input_files);

    merge_sorted_readers(&mut stdout(), &mut sorted_readers, sort_field)?;

    Ok(())
}

fn open_zst_files(filenames: &mut [PathBuf]) -> Vec<BufReader<Decoder<'static, BufReader<File>>>> {
    filenames
        .iter_mut()
        .map(|f| BufReader::new(Decoder::new(File::open(f).unwrap()).unwrap()))
        .collect()
}
