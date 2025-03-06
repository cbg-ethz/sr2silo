use clap::Parser;
use itertools::Itertools;
use serde_json::Value;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, VecDeque};
use std::fs::File;
use std::io::{stdin, stdout, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::thread::JoinHandle;
use std::{env, fs, thread};
use zstd::stream::Decoder;
use zstd::Encoder;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(long)]
    sort_field: String,

    #[arg(long)]
    tmp_directory: Option<String>,

    #[arg(long, default_value_t = 64)]
    parallel_files: usize,

    #[arg(long, default_value_t = 16)]
    num_threads: usize,
}

fn main() -> std::io::Result<()> {
    let args = Args::parse();

    let tmp_dir = if let Some(given_tmp_dir) = args.tmp_directory {
        if Path::new(&given_tmp_dir).exists() {
            assert_eq!(
                fs::read_dir(&given_tmp_dir)?.count(),
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

    assert!(
        args.parallel_files > 1,
        "We need to work on at least 2 files in parallel."
    );

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
        &args.sort_field,
        args.parallel_files,
        merge_iteration,
        args.num_threads,
    )?;

    merge_iteration += 1;

    if input_files.len() == 0 {
        panic!("No input files received");
    }

    while input_files.len() > args.parallel_files {
        input_files = merge_files_in_batches(
            input_files,
            &tmp_dir,
            &args.sort_field,
            args.parallel_files,
            merge_iteration,
            args.num_threads,
        )?;
        merge_iteration += 1;
    }

    merge_files(input_files, &mut stdout().lock(), &args.sort_field)?;

    Ok(())
}

fn merge_files_in_batches<I>(
    input_files: I,
    tmp_dir: &PathBuf,
    sort_field: &String,
    batch_size: usize,
    merge_iteration: usize,
    max_threads: usize,
) -> std::io::Result<Vec<PathBuf>>
where
    I: IntoIterator<Item = PathBuf>,
{
    let tmp_dir = Arc::new(tmp_dir.clone());
    let sort_field = Arc::new(sort_field.clone());
    let mut next_input_files = Vec::new();
    let mut handles: VecDeque<JoinHandle<std::io::Result<()>>> = VecDeque::new();

    for (batch_id, batch) in input_files
        .into_iter()
        .chunks(batch_size)
        .into_iter()
        .enumerate()
    {
        let tmp_dir = Arc::clone(&tmp_dir);
        let sort_field = Arc::clone(&sort_field);
        let batch: Vec<PathBuf> = batch.collect();
        let file_name = tmp_dir.join(format!(
            "merged_chunks_{}_{}.ndjson.zst",
            merge_iteration, batch_id
        ));

        next_input_files.push(file_name.clone());

        // Limit parallelism by waiting when we hit max_threads
        while handles.len() >= max_threads {
            if let Some(handle) = handles.pop_front() {
                handle.join().unwrap().expect("IO error in worker thread"); // Wait for one to finish
            }
        }

        let handle = thread::spawn(move || -> std::io::Result<()> {
            let file = File::create(&file_name)?;
            let mut encoder = Encoder::new(file, 3)?;
            merge_files(batch, &mut encoder, &sort_field)?;
            encoder.finish()?;
            Ok(())
        });

        handles.push_back(handle);
    }

    // Wait for all remaining threads to complete
    for handle in handles {
        handle.join().unwrap()?;
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
