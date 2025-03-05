use clap::Command;
use itertools::Itertools;
use serde_json::Value;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{stdin, stdout, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::{env, fs};
use zstd::stream::Decoder;
use zstd::Encoder;

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
        .arg(
            clap::Arg::new("parallel_files")
                .long("parallel-files")
                .value_name("u32")
                .help("Specifies the number of files that should be merged per iteration")
                .default_value("512"),
        )
        .get_matches();

    let tmp_dir = if let Some(given_tmp_dir) = matches.get_one::<String>("tmp_dir") {
        if fs::exists(given_tmp_dir)? {
            assert_eq!(
                fs::read_dir(given_tmp_dir)?.count(),
                0,
                "The given tmp directory is not empty"
            );
        } else {
            fs::create_dir(&given_tmp_dir)?
        };
        PathBuf::from(given_tmp_dir)
    } else {
        env::temp_dir()
    };

    let parallel_files: usize = matches
        .get_one::<String>("parallel_files")
        .unwrap()
        .parse()
        .unwrap();
    assert!(
        parallel_files > 1,
        "We need to work on at least 2 files in parallel."
    );

    let sort_field = matches.get_one("sort_field").unwrap();

    let reader = stdin();
    let reader = reader.lock();

    let mut merge_iteration = 0;

    let input_files_stdin = BufReader::new(reader)
        .lines()
        .map(Result::unwrap)
        .map(PathBuf::from);

    let mut input_files = merge_files_in_batches(
        input_files_stdin,
        &tmp_dir,
        sort_field,
        parallel_files,
        merge_iteration,
    )?;

    merge_iteration += 1;

    if input_files.len() == 0 {
        panic!("No input files received");
    }

    while input_files.len() > parallel_files {
        input_files = merge_files_in_batches(
            input_files,
            &tmp_dir,
            sort_field,
            parallel_files,
            merge_iteration,
        )?;
        merge_iteration += 1;
    }

    merge_files(input_files, &mut stdout().lock(), sort_field)?;

    Ok(())
}

fn merge_files_in_batches<I>(
    input_files: I,
    tmp_dir: &PathBuf,
    sort_field: &String,
    batch_size: usize,
    merge_iteration: usize,
) -> std::io::Result<Vec<PathBuf>>
where
    I: IntoIterator<Item = PathBuf>,
{
    let mut next_input_files: Vec<PathBuf> = Vec::new();
    for (batch_id, batch) in input_files
        .into_iter()
        .chunks(batch_size)
        .into_iter()
        .enumerate()
    {
        let file_name = Path::join(
            tmp_dir,
            format!("merged_chunks_{}_{}.ndjson.zst", merge_iteration, batch_id),
        );
        let file = File::create(file_name.clone())?;
        let mut encoder = Encoder::new(file, 3)?;

        merge_files(batch, &mut encoder, sort_field)?;

        encoder.finish()?;

        next_input_files.push(file_name);
    }
    Ok(next_input_files)
}

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
fn merge_files<I, W: Write>(
    files: I,
    output: &mut W,
    sort_field_name: &String,
) -> std::io::Result<()>
where
    I: IntoIterator<Item = PathBuf>,
{
    let mut heap = BinaryHeap::new();

    let sorted_readers = files
        .into_iter()
        .map(|f| BufReader::new(Decoder::new(File::open(f).unwrap()).unwrap()));

    // Store an iterator for each reader
    let mut reader_iters: Vec<_> = sorted_readers.into_iter().map(|r| r.lines()).collect();

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
