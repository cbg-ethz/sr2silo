# Snakemake Workflow Configuration

## New Dynamic Configuration Approach

The workflow now supports dynamic sample selection based on date ranges and location codes instead of hardcoded `SAMPLE_BATCH_IDS`.

### Configuration Parameters

Replace the old `SAMPLE_BATCH_IDS` with these new parameters:

```yaml
LOCATIONS:
  - 5    # Lugano (TI)
  - 10   # Zürich (ZH)
  - 15   # Basel (BS)
  - 17   # Chur (GR)
  - 25   # Additional location code

START_DATE: "2025-02-26"
END_DATE: "2025-03-08"
TIMELINE_FILE: "/path/to/timeline.tsv"
```

### Timeline File Format

The timeline file should be a TSV with the following columns (no header):

```
sample	batch	reads	proto	location_code	date	location
A1_05_2025_05_21	20250613_2427498204	250	v532	5	2025-05-21	Lugano (TI)
A2_15_2025_05_22	20250613_2427498204	250	v532	15	2025-05-22	Basel (BS)
```

### Location Code Mapping

Common location codes used in the examples:
- 5: Lugano (TI)
- 10: Zürich (ZH)
- 15: Basel (BS)
- 17: Chur (GR)
- 25: Additional location

### Usage

1. Update your config file with the desired location codes, date range, and timeline file path
2. Run the workflow as usual:

```bash
snakemake --use-conda -c1
```

The workflow will automatically:
1. Parse the timeline file
2. Filter samples by the specified date range (START_DATE to END_DATE)
3. Filter samples by the specified location codes
4. Generate the sample-batch mapping dynamically
5. Process all matching samples

### Migration from Old Config

**Old approach:**
```yaml
SAMPLE_BATCH_IDS:
  A1_05_2025_02_26 : "20250321_2429695737"
  A2_15_2025_02_27 : "20250321_2429695737"
  # ... many more entries
```

**New approach:**
```yaml
LOCATIONS: [5, 15]
START_DATE: "2025-02-26"
END_DATE: "2025-02-27"
TIMELINE_FILE: "/path/to/timeline.tsv"
```

This approach is much more flexible and allows easy selection of samples based on biological/temporal criteria rather than maintaining large hardcoded lists.
