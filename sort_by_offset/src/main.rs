use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write, Read, stdout, stdin};
use serde_json::Value;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

// Function to read ndjson lines from a reader
fn read_ndjson_lines<R: Read>(reader: R) -> std::io::Result<Vec<Value>> {
    let reader = BufReader::new(reader);
    let mut lines = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let json: Value = serde_json::from_str(&line)?;
        lines.push(json);
    }
    Ok(lines)
}

// Function to write ndjson lines to a writer
fn write_ndjson_lines<W: Write>(writer: &mut W, lines: &[Value]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(writer);
    for line in lines {
        writeln!(writer, "{}", line)?;
    }
    Ok(())
}

// Function to sort ndjson lines by "N_length"
fn sort_by_n_length(mut lines: Vec<Value>) -> Vec<Value> {
    lines.sort_by(|a, b| {
        let a_length = a["N_length"].as_i64().unwrap_or(0);
        let b_length = b["N_length"].as_i64().unwrap_or(0);
        a_length.cmp(&b_length)
    });
    lines
}

// Merging function that reads from readers and writes to any object implementing `Write`
fn merge_sorted_readers<R: Read, W: Write>(output: &mut W, sorted_readers: &mut [BufReader<R>]) -> std::io::Result<()> {
    let mut heap = BinaryHeap::new();

    // Initialize heap with the first line from each reader
    for (index, reader) in sorted_readers.iter_mut().enumerate() {
        if let Some(Ok(line)) = reader.lines().next() {
            let json: Value = serde_json::from_str(&line)?;
            heap.push((json, index));
        }
    }

    let mut writer = BufWriter::new(output);

    while let Some((json, index)) = heap.pop() {
        writeln!(writer, "{}", json)?;
        if let Some(Ok(line)) = sorted_readers[index].lines().next() {
            let json: Value = serde_json::from_str(&line)?;
            heap.push((json, index));
        }
    }

    Ok(())
}

fn main() -> std::io::Result<()> {
    let chunk_size = 1000; // Adjust based on available memory
    let mut chunk_counter = 0;
    let mut sorted_files = Vec::new();

    // Read and process input from stdin in chunks
    let reader = stdin();
    let reader = reader.lock();
    let mut lines = Vec::new();

    for line in BufReader::new(reader).lines() {
        let line = line?;
        let json: Value = serde_json::from_str(&line)?;
        lines.push(json);

        if lines.len() >= chunk_size {
            let sorted_lines = sort_by_n_length(lines);
            let chunk_file = format!("chunk_{}.ndjson", chunk_counter);
            let mut file = File::create(&chunk_file)?;
            write_ndjson_lines(&mut file, &sorted_lines)?;
            sorted_files.push(File::open(&chunk_file)?);
            lines = Vec::new();
            chunk_counter += 1;
        }
    }

    // Process any remaining lines
    if !lines.is_empty() {
        let sorted_lines = sort_by_n_length(lines);
        let chunk_file = format!("chunk_{}.ndjson", chunk_counter);
        let mut file = File::create(&chunk_file)?;
        write_ndjson_lines(&mut file, &sorted_lines)?;
        sorted_files.push(File::open(&chunk_file)?);
    }

    // Convert sorted files to BufReaders
    let mut sorted_readers: Vec<_> = sorted_files.iter_mut().map(|f| BufReader::new(f)).collect();

    // Merge all sorted chunk files and write to stdout
    merge_sorted_readers(&mut stdout(), &mut sorted_readers)?;

    // Optionally, clean up temporary chunk files
    for file in sorted_files {
        drop(file);
    }

    Ok(())
}

