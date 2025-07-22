#!/usr/bin/env python3
"""
Simple script to analyze and visualize read count statistics from subsampling.

Usage:
    python analyze_read_counts.py logs/subsampling/read_counts.tsv

This script provides:
- Summary statistics of original vs subsampled reads
- Weekly read ingestion trends
- Basic plotting (if matplotlib is available)
"""

from __future__ import annotations

import csv
import sys

try:
    import pandas as pd

    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False


def analyze_read_counts(tsv_file):
    """Analyze read count statistics from the subsampling log."""

    # Read the data
    try:
        if HAS_PANDAS:
            df = pd.read_csv(tsv_file, sep="\t")  # type: ignore
            data = None
        else:
            # Simple CSV fallback
            with open(tsv_file, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")
                data = list(reader)
            df = None
    except FileNotFoundError:
        print(f"Error: File {tsv_file} not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    print(f"📊 Read Count Analysis Report")
    print(f"{'='*50}")
    print(f"Data source: {tsv_file}")

    if HAS_PANDAS and df is not None:
        _analyze_with_pandas(df, tsv_file)
    elif data is not None:
        _analyze_with_csv(data, tsv_file)


def _analyze_with_pandas(df, tsv_file):
    """Analysis using pandas (full-featured)."""
    print(f"Total samples processed: {len(df)}")
    print()

    # Summary statistics
    print("📈 Summary Statistics:")
    print(f"Original reads - Total: {df['original_reads'].sum():,}")
    print(f"Original reads - Mean: {df['original_reads'].mean():,.0f}")
    print(f"Original reads - Median: {df['original_reads'].median():,.0f}")
    print(
        f"Original reads - Range: {df['original_reads'].min():,} - {df['original_reads'].max():,}"
    )
    print()

    print(f"Subsampled reads - Total: {df['subsampled_reads'].sum():,}")
    print(f"Subsampled reads - Mean: {df['subsampled_reads'].mean():,.0f}")
    print(f"Subsampled reads - Median: {df['subsampled_reads'].median():,.0f}")
    print()

    # Reduction statistics
    df["reduction_factor"] = df["original_reads"] / df["subsampled_reads"]
    df["reads_removed"] = df["original_reads"] - df["subsampled_reads"]
    df["percent_kept"] = (df["subsampled_reads"] / df["original_reads"]) * 100

    print(f"💾 Data Reduction:")
    print(f"Average reduction factor: {df['reduction_factor'].mean():.1f}x")
    print(f"Total reads removed: {df['reads_removed'].sum():,}")
    print(f"Average % of reads kept: {df['percent_kept'].mean():.1f}%")
    print()

    # Weekly trends
    if "week_of_year" in df.columns:
        weekly_stats = (
            df.groupby("week_of_year")
            .agg({"original_reads": ["sum", "count"], "subsampled_reads": "sum"})
            .round(0)
        )

        print("📅 Weekly Ingestion Summary:")
        print(weekly_stats)
        print()

    _create_plots_if_available(df, tsv_file)


def _analyze_with_csv(data, tsv_file):
    """Simple analysis using basic CSV (no pandas)."""
    print(f"Total samples processed: {len(data)}")
    print("(Using basic CSV mode - install pandas for full analysis)")
    print()

    if not data:
        print("No data to analyze.")
        return

    # Calculate basic statistics
    original_reads = [int(row["original_reads"]) for row in data]
    subsampled_reads = [int(row["subsampled_reads"]) for row in data]

    print("📈 Summary Statistics:")
    print(f"Original reads - Total: {sum(original_reads):,}")
    print(f"Original reads - Mean: {sum(original_reads) / len(original_reads):,.0f}")
    print()

    print(f"Subsampled reads - Total: {sum(subsampled_reads):,}")
    print(
        f"Subsampled reads - Mean: {sum(subsampled_reads) / len(subsampled_reads):,.0f}"
    )
    print()

    # Basic reduction statistics
    total_removed = sum(
        orig - sub for orig, sub in zip(original_reads, subsampled_reads)
    )
    avg_reduction = sum(
        orig / sub for orig, sub in zip(original_reads, subsampled_reads)
    ) / len(data)

    print(f"💾 Data Reduction:")
    print(f"Average reduction factor: {avg_reduction:.1f}x")
    print(f"Total reads removed: {total_removed:,}")
    print()

    # Weekly trends (simple)
    if data and "week_of_year" in data[0]:
        weeks = {}
        for row in data:
            week = row["week_of_year"]
            if week not in weeks:
                weeks[week] = {"count": 0, "original_total": 0}
            weeks[week]["count"] += 1
            weeks[week]["original_total"] += int(row["original_reads"])

        print("📅 Weekly Ingestion Summary:")
        for week, stats in sorted(weeks.items()):
            print(
                f"  {week}: {stats['count']} samples, {stats['original_total']:,} total reads"
            )
        print()


def _create_plots_if_available(df, tsv_file):
    try:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle("Read Count Analysis", fontsize=16)

        # Original vs subsampled reads scatter
        axes[0, 0].scatter(df["original_reads"], df["subsampled_reads"], alpha=0.6)
        axes[0, 0].set_xlabel("Original Reads")
        axes[0, 0].set_ylabel("Subsampled Reads")
        axes[0, 0].set_title("Original vs Subsampled Reads")
        axes[0, 0].grid(True, alpha=0.3)

        # Read count distributions
        axes[0, 1].hist(
            df["original_reads"], bins=20, alpha=0.7, label="Original", color="blue"
        )
        axes[0, 1].hist(
            df["subsampled_reads"], bins=20, alpha=0.7, label="Subsampled", color="red"
        )
        axes[0, 1].set_xlabel("Read Count")
        axes[0, 1].set_ylabel("Frequency")
        axes[0, 1].set_title("Read Count Distributions")
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)

        # Reduction factors
        axes[1, 0].hist(df["reduction_factor"], bins=20, alpha=0.7, color="green")
        axes[1, 0].set_xlabel("Reduction Factor")
        axes[1, 0].set_ylabel("Frequency")
        axes[1, 0].set_title("Data Reduction Factors")
        axes[1, 0].grid(True, alpha=0.3)

        # Weekly trends if available
        if "week_of_year" in df.columns:
            weekly_totals = df.groupby("week_of_year")["original_reads"].sum()
            axes[1, 1].bar(
                weekly_totals.index, weekly_totals.values, alpha=0.7, color="orange"
            )
            axes[1, 1].set_xlabel("Week of Year")
            axes[1, 1].set_ylabel("Total Original Reads")
            axes[1, 1].set_title("Weekly Read Ingestion")
            axes[1, 1].grid(True, alpha=0.3)
        else:
            axes[1, 1].text(
                0.5,
                0.5,
                "Weekly data\nnot available",
                ha="center",
                va="center",
                transform=axes[1, 1].transAxes,
            )
            axes[1, 1].set_title("Weekly Trends")

        plt.tight_layout()

        # Save plot
        output_plot = tsv_file.replace(".tsv", "_analysis.png")
        plt.savefig(output_plot, dpi=300, bbox_inches="tight")
        print(f"📊 Visualization saved to: {output_plot}")

        # Optionally show plot
        # plt.show()

    except ImportError:
        print("📊 Install matplotlib for visualizations: conda install matplotlib")
    except Exception as e:
        print(f"Warning: Could not create plots: {e}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python analyze_read_counts.py <read_counts.tsv>")
        sys.exit(1)

    tsv_file = sys.argv[1]
    analyze_read_counts(tsv_file)
