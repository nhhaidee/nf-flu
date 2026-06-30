#!/usr/bin/env python3
import argparse
import json
import logging
import os
import re
import sys
import tempfile
import numpy as np
import pandas as pd
import torch

# Set up module-level logger
logger = logging.getLogger("waveseekernet_analyze")
logger.setLevel(logging.INFO)

if not logger.handlers:
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter("%(levelname)s | %(name)s | %(message)s"))
    logger.addHandler(handler)

# Fix for headless container crashes (e.g. Docker, HPC)
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patches as patches

from Bio import SeqIO
from WaveSeekerNet import WaveSeekerClassifier, IsotonicCalibrator, compute_calibration_metrics
from WaveSeekerNet.utils import fasta_to_one_hot, fasta_to_fcgr


def detect_segment(header):
    """Detect segment type from FASTA sequence header."""
    header_upper = header.upper()

    # Split header into tokens using non-alphanumeric characters as delimiters
    tokens = [t for t in re.split(r'[^A-Z0-9]+', header_upper) if t]
    if not tokens:
        return None

    # 1. Search tokens from right to left for exact segment name match
    # This ensures that if a sample name contains a segment name (e.g. "MP"),
    # we still match the actual segment identifier near the end.
    for token in reversed(tokens):
        if token in ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]:
            return token

    # 2. Search tokens from right to left for segment number matches (e.g. SEGMENT_1, SEGMENT1, SEG1)
    num_map = {
        "1": "PB2", "2": "PB1", "3": "PA", "4": "HA",
        "5": "NP", "6": "NA", "7": "MP", "8": "NS"
    }
    for i in reversed(range(len(tokens))):
        token = tokens[i]
        # Match SEGMENTX or SEG_X or SEGX
        match = re.match(r'^(?:SEGMENT|SEG)_?([1-8])$', token)
        if match:
            return num_map[match.group(1)]
        # Match separate tokens: e.g. ["SEGMENT", "8"] or ["SEG", "8"]
        if token in num_map and i > 0:
            prev_token = tokens[i - 1]
            if prev_token in ["SEGMENT", "SEG"]:
                return num_map[token]

    # 3. Fallback: search for any occurrence of segment patterns as substrings with boundary regexes
    for seg in ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]:
        pattern = rf"(?:^|[^A-Z0-9]){seg}(?:[^A-Z0-9]|$)"
        if re.search(pattern, header_upper):
            return seg

    # 4. Fallback for segment patterns (SEGMENT_1)
    num_map_full = {
        "GENE_1": "PB2", "GENE_2": "PB1", "GENE_3": "PA", "GENE_4": "HA",
        "GENE_5": "NP", "GENE_6": "NA", "GENE_7": "MP", "GENE_8": "NS",
        "SEGMENT_1": "PB2", "SEGMENT_2": "PB1", "SEGMENT_3": "PA", "SEGMENT_4": "HA",
        "SEGMENT_5": "NP", "SEGMENT_6": "NA", "SEGMENT_7": "MP", "SEGMENT_8": "NS"
    }
    for pattern_key, seg in num_map_full.items():
        pattern = rf"(?:^|[^A-Z0-9]){pattern_key}(?:[^A-Z0-9]|$)"
        if re.search(pattern, header_upper):
            return seg

    return None


def is_at_rich(codon):
    """Returns True if the codon contains 2 or more A/T bases."""
    return sum(1 for b in codon if b in 'AT') >= 2


def calculate_rscu(codon_counts):
    """
    Calculate Relative Synonymous Codon Usage (RSCU) for 61 sense codons.
    """
    CODON_TO_AA = {
        'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S', 'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
        'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'TAC': 'Y', 'TAT': 'Y', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TGC': 'C', 'TGT': 'C',
        'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F', 'TGG': 'W'
    }

    # Group codons by amino acid
    aa_to_codons = {}
    for codon, aa in CODON_TO_AA.items():
        aa_to_codons.setdefault(aa, []).append(codon)

    rscu = {}
    for aa, synonymous_codons in aa_to_codons.items():
        total_count = sum(codon_counts.get(c, 0) for c in synonymous_codons)
        n_synonymous = len(synonymous_codons)

        for c in synonymous_codons:
            if total_count == 0:
                rscu[c] = 0.0
            else:
                expected = total_count / n_synonymous
                rscu[c] = codon_counts.get(c, 0) / expected

    return rscu


def parse_cds_coordinates(header):
    """Extract inclusive, 1-indexed (start, end) coordinate tuples from the FASTA header."""
    fields = header.split('|')
    if len(fields) >= 4:
        loc_str = fields[3].split()[0]
    else:
        match = re.search(r'\[location=([^\]]+)\]', header)
        loc_str = match.group(1) if match else header

    ranges = []
    for r_match in re.finditer(r'(\d+)\.\.(\d+)', loc_str):
        start = int(r_match.group(1))
        end = int(r_match.group(2))
        ranges.append((start, end))

    is_complement = 'complement' in loc_str.lower()
    return ranges, is_complement


def find_matching_consensus(cds_id, consensus_dict):
    """Find matching consensus sequence by substring ID comparison."""
    parts = cds_id.split('|')
    con_id_candidate = parts[0]
    if con_id_candidate in consensus_dict:
        return consensus_dict[con_id_candidate]

    if cds_id in consensus_dict:
        return consensus_dict[cds_id]
    for con_id, record in consensus_dict.items():
        if con_id in cds_id or cds_id in con_id:
            return record
    return None


def plot_shap_given_ax(ax, df, text='', is_at_rich_fn=None, show_legend=True, rotate_xticks=45):
    """Bar plot for ranked SHAP values."""
    s = pd.Series(df).sort_values()

    if is_at_rich_fn is None:
        nuc_color_map = {'A': 'green', 'C': 'blue', 'G': 'orange', 'T': 'red'}
        tick_colors = [next((col for nuc, col in nuc_color_map.items() if f"({nuc})" in c), 'black') for c in s.index]
        ax.bar(s.index.astype(str), s.values)
        legend_patches = None
    else:
        tick_colors = ['orange' if is_at_rich_fn(c) else 'blue' for c in s.index]
        legend_patches = [
            mpatches.Patch(color='blue', label='G/C-rich'),
            mpatches.Patch(color='orange', label='A/T-rich')
        ]
        bar_colors = ['red' if x < 0 else 'green' for x in s.values]
        ax.bar(s.index.astype(str), s.values, color=bar_colors)

    ax.set_ylabel('SHAP Values', fontsize=14)
    ax.tick_params(axis='y', labelsize=12)
    ax.tick_params(axis='x', labelsize=16, labelrotation=rotate_xticks)
    ax.set_title(text, loc='left', fontsize=15)

    for label, color in zip(ax.get_xticklabels(), tick_colors):
        label.set_color(color)
        label.set_ha('right')
        label.set_rotation_mode('anchor')

    ax.grid(axis='y', linestyle='--', alpha=0.7)

    if show_legend and legend_patches is not None:
        ax.legend(handles=legend_patches, loc='upper left', fontsize='small')


# ==========================================
# GLYPH HELPER FUNCTIONS (For Seq Logo)
# ==========================================
def plot_a(ax, base, left_edge, height, color):
    a_polygon_coords = [
        np.array([[0.0, 0.0], [0.5, 1.0], [0.5, 0.8], [0.2, 0.0]]),
        np.array([[1.0, 0.0], [0.5, 1.0], [0.5, 0.8], [0.8, 0.0]]),
        np.array([[0.225, 0.45], [0.775, 0.45], [0.85, 0.3], [0.15, 0.3]])
    ]
    for polygon_coords in a_polygon_coords:
        ax.add_patch(patches.Polygon(
            (np.array([1, height])[None, :] * polygon_coords + np.array([left_edge, base])[None, :]),
            facecolor=color, edgecolor=color))

def plot_c(ax, base, left_edge, height, color):
    ax.add_patch(patches.Ellipse(xy=[left_edge + 0.65, base + 0.5 * height], width=1.3, height=height,
                                 facecolor=color, edgecolor=color))
    ax.add_patch(patches.Ellipse(xy=[left_edge + 0.65, base + 0.5 * height], width=0.7 * 1.3, height=0.7 * height,
                                 facecolor='white', edgecolor='white'))
    ax.add_patch(patches.Rectangle(xy=[left_edge + 1, base], width=1.0, height=height,
                                   facecolor='white', edgecolor='white', fill=True))

def plot_g(ax, base, left_edge, height, color):
    ax.add_patch(patches.Ellipse(xy=[left_edge + 0.65, base + 0.5 * height], width=1.3, height=height,
                                 facecolor=color, edgecolor=color))
    ax.add_patch(patches.Ellipse(xy=[left_edge + 0.65, base + 0.5 * height], width=0.7 * 1.3, height=0.7 * height,
                                 facecolor='white', edgecolor='white'))
    ax.add_patch(patches.Rectangle(xy=[left_edge + 1, base], width=1.0, height=height,
                                   facecolor='white', edgecolor='white', fill=True))
    ax.add_patch(patches.Rectangle(xy=[left_edge + 0.825, base + 0.085 * height], width=0.174, height=0.415 * height,
                                   facecolor=color, edgecolor=color, fill=True))
    ax.add_patch(patches.Rectangle(xy=[left_edge + 0.625, base + 0.35 * height], width=0.374, height=0.15 * height,
                                   facecolor=color, edgecolor=color, fill=True))

def plot_t(ax, base, left_edge, height, color):
    ax.add_patch(patches.Rectangle(xy=[left_edge + 0.4, base],
                                   width=0.2, height=height, facecolor=color, edgecolor=color, fill=True))
    ax.add_patch(patches.Rectangle(xy=[left_edge, base + 0.8 * height],
                                   width=1.0, height=0.2 * height, facecolor=color, edgecolor=color, fill=True))

default_colors = {0: 'green', 1: 'blue', 2: 'orange', 3: 'red'}
default_plot_funcs = {0: plot_a, 1: plot_c, 2: plot_g, 3: plot_t}


# ==========================================
# CORE PLOTTING FUNCTIONS
# ==========================================
def plot_weights_given_ax(ax, array,
                          height_padding_factor,
                          length_padding,
                          subticks_frequency,
                          highlight,
                          start_index=0,
                          colors=default_colors,
                          plot_funcs=default_plot_funcs):
    if len(array.shape) == 3:
        array = np.squeeze(array)
    if array.shape[0] == 4 and array.shape[1] != 4:
        array = array.transpose(1, 0)

    seq_len = array.shape[0]

    max_pos_height = 0.0
    min_neg_height = 0.0
    heights_at_positions = []
    depths_at_positions = []

    for i in range(seq_len):
        acgt_vals = sorted(enumerate(array[i, :]), key=lambda x: abs(x[1]))
        positive_height_so_far = 0.0
        negative_height_so_far = 0.0

        for letter_idx, score in acgt_vals:
            plot_func = plot_funcs[letter_idx]
            color = colors[letter_idx]
            if score > 0:
                h = positive_height_so_far
                positive_height_so_far += score
            else:
                h = negative_height_so_far
                negative_height_so_far += score

            plot_func(ax=ax, base=h, left_edge=i, height=score, color=color)

        max_pos_height = max(max_pos_height, positive_height_so_far)
        min_neg_height = min(min_neg_height, negative_height_so_far)
        heights_at_positions.append(positive_height_so_far)
        depths_at_positions.append(negative_height_so_far)

    ax.set_xlim(-length_padding, seq_len + length_padding)
    tick_indices = np.arange(0, seq_len, subticks_frequency)
    tick_locs = tick_indices + 0.5
    tick_labels = (tick_indices + start_index + 1).astype(int)

    ax.set_xticks(tick_locs)
    ax.set_xticklabels(tick_labels, fontsize=14)
    ax.set_ylabel('SHAP Values', fontsize=14)
    ax.tick_params(axis='y', labelsize=12)

    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())

    first_codon = (start_index // 3) + 1
    last_codon = ((start_index + seq_len) // 3) + 1

    codon_ticks = []
    codon_labels = []

    for c in range(first_codon, last_codon + 2):
        global_center_coord = (c - 1) * 3 + 1.5
        local_center = global_center_coord - start_index

        if 0 <= local_center <= seq_len:
            codon_ticks.append(local_center)
            codon_labels.append(c)

    ax2.set_xticks(codon_ticks)
    ax2.set_xticklabels(codon_labels, fontsize=14)
    ax2.grid(False)

    for color, regions in highlight.items():
        for item in regions:
            label = None
            if len(item) == 3:
                h_start, h_end, label = item
            else:
                h_start, h_end = item

            local_start = h_start - start_index
            local_end = h_end - start_index

            if local_end < 0 or local_start > seq_len:
                continue

            draw_start = max(0, local_start)
            draw_end = min(seq_len, local_end)
            slice_idx_start = int(max(0, local_start))
            slice_idx_end = int(min(seq_len, local_end))

            if slice_idx_end > slice_idx_start:
                valid_depths = depths_at_positions[slice_idx_start:slice_idx_end]
                valid_heights = heights_at_positions[slice_idx_start:slice_idx_end]
                if valid_depths and valid_heights:
                    min_depth = np.min(valid_depths)
                    max_height = np.max(valid_heights)
                    rect = patches.Rectangle(
                        xy=[draw_start, min_depth],
                        width=draw_end - draw_start,
                        height=max_height - min_depth,
                        edgecolor=color, facecolor='none', linewidth=2
                    )
                    ax.add_patch(rect)
                    if label:
                        ax.text(x=draw_start + (draw_end - draw_start) / 2,
                                y=max_height + (max_pos_height * 0.05),
                                s=label, color=color, ha='center', va='bottom', fontweight='bold')

    height_pad = max(abs(min_neg_height), abs(max_pos_height)) * height_padding_factor
    ax.set_ylim(min_neg_height - height_pad, max_pos_height + height_pad)


def plot_scatter_given_ax(ax, shap_dict, title='', markers=None, highlight_region=None):
    """Genome-wide SHAP scatter/line trace."""
    positions = sorted(shap_dict.keys())
    values = [shap_dict[p] for p in positions]
    
    ax.plot(positions, values, color='pink', lw=2, zorder=1)
    ax.scatter(positions, values, color='crimson', s=50, zorder=2)
    
    if markers:
        for m in markers:
            pos = m['pos']
            label = m.get('label', str(pos))
            color = m.get('color', 'dodgerblue')
            
            # Fetch and display the exact SHAP score at this codon position
            score_val = shap_dict.get(pos, 0.0)
            label_with_score = f"{label}\nSHAP: {score_val:.4f}"
            
            ax.axvline(x=pos, color=color, linestyle='--', lw=2, zorder=0)
            
            y_max = max(values) if values else 1.0
            ax.text(pos, y_max, label_with_score, color='crimson',
                    fontsize=12, va='top', ha='center', fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'))

    if highlight_region:
        start, end = highlight_region
        ax.axvspan(start, end, color='gray', alpha=0.15, zorder=0)

    ax.set_xlabel('Codon Position in CDS', fontsize=14)
    ax.set_ylabel('SHAP Values', fontsize=14)
    ax.set_title(title, loc='left', fontsize=15)
    ax.tick_params(axis='y', labelsize=12)
    ax.tick_params(axis='x', labelsize=12)
    ax.grid(True, linestyle=':', alpha=0.6)


def plot_dashboard_grid(
    scatter_data,      # Dict: {position_index: shap_score}
    logo_array,        # Array: (L, 4)
    bar_data,          # Dict: {codon_name: shap_score} (Ranked Codons)
    top_pos_data,      # List of tuples: [(pos, score, nuc), ...]
    codon=['Lysine K', 'AAG', 'AAA'],
    fig_title='',
    save_path=None,
    figsize=(20, 14),
    scatter_markers=None,
    logo_region=None,
    logo_highlight=None,
    subticks_frequency=1,
    is_at_rich_fn=None
):
    fig = plt.figure(figsize=figsize, constrained_layout=True)

    gs = fig.add_gridspec(
        3, 2,
        height_ratios=[0.9, 1.1, 1.0],
        width_ratios=[1, 1]
    )

    ax_ranked_codons = fig.add_subplot(gs[0, :])
    ax_logo = fig.add_subplot(gs[1, :])
    ax_scatter = fig.add_subplot(gs[2, 0])
    ax_top_pos = fig.add_subplot(gs[2, 1])

    if fig_title:
        fig.suptitle(fig_title, fontsize=16)

    plot_shap_given_ax(
        ax_ranked_codons,
        bar_data,
        text="A - Ranked Codons (Global Feature Importance)",
        is_at_rich_fn=is_at_rich_fn,
        show_legend=True
    )

    if logo_region:
        s_start, s_end = logo_region
        s_start, s_end = max(0, s_start), min(logo_array.shape[0], s_end)
        sliced_array, start_idx = logo_array[s_start:s_end], s_start
        logo_title = f"B - Sequence Logo (Nucleotides {s_start}-{s_end})"
    else:
        sliced_array, start_idx = logo_array, 0
        logo_title = "B - Sequence Logo (Full Sequence)"

    plot_weights_given_ax(
        ax=ax_logo,
        array=sliced_array,
        height_padding_factor=0.2,
        length_padding=1.0,
        subticks_frequency=subticks_frequency,
        highlight=logo_highlight if logo_highlight else {},
        start_index=start_idx
    )
    ax_logo.set_title(logo_title, loc='left', fontsize=15)

    scatter_highlight = (logo_region[0] // 3, logo_region[1] // 3) if logo_region else None
    plot_scatter_given_ax(
        ax_scatter,
        scatter_data,
        title=f"C - Segment-wide Codon SHAP Values (Highlighting Pos {codon[0]}: {codon[1]}/{codon[2]})",
        markers=scatter_markers,
        highlight_region=scatter_highlight
    )

    top_pos_dict = {f"{p}({n})": s for p, s, n in top_pos_data}
    if not top_pos_dict:
        top_pos_dict = {'None': 0}

    plot_shap_given_ax(
        ax_top_pos,
        top_pos_dict,
        text=f"D - Top {len(top_pos_dict)} Positive Nucleotide Positions",
        is_at_rich_fn=None,
        show_legend=False
    )

    if save_path:
        import os
        directory = os.path.dirname(save_path)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Figure saved to: {save_path}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Unified prediction and SHAP analysis for WaveSeekerNet")
    parser.add_argument("--consensus-fasta", required=True, help="Input full consensus FASTA file containing segments")
    parser.add_argument("--cds-fasta", required=True,
                        help="Input CDS nucleotide FASTA file (.ffn) containing coordinates")
    parser.add_argument("--config", required=True, help="Path to segments JSON config file")
    parser.add_argument("--weights-dir", required=True, help="Directory containing the weights subfolders")
    parser.add_argument("--encoding", choices=["one-hot", "fcgr"], default="one-hot", help="Encoding format")
    parser.add_argument("--k", type=int, default=6, help="k-mer size for FCGR")

    # Defaults
    parser.add_argument("--patch_size", type=int, nargs=2, default=[12, 5], help="Fallback patch size")
    parser.add_argument("--emb_dim", type=int, default=64, help="Fallback emb_dim")
    parser.add_argument("--res_len", type=int, default=5, help="Fallback res_len")
    parser.add_argument("--final_hidden_size", type=int, default=24, help="Fallback final hidden size")
    parser.add_argument("--batch_size", type=int, default=256, help="Fallback batch size")

    parser.add_argument("--output-predictions", required=True, help="Output TSV file path for ensemble predictions")
    parser.add_argument("--output-html", help="Output HTML report file path")
    parser.add_argument("--output-shap-prefix", required=True, help="Prefix path for SHAP output files")
    parser.add_argument("--device", default=None, help="Device ('cpu' or 'cuda')")

    args = parser.parse_args()

    SENSE_CODONS = [
        'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT',
        'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
        'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
        'TAC', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGC', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT', 'TGG'
    ]

    with open(args.config, 'r') as f:
        seg_configs = json.load(f)

    classes = ["Human", "Avian", "Non-human Mammals"]
    n_out = len(classes)

    # 1. Store predictions and SHAP matrices per consensus segment
    pred_results = []
    consensus_shap_data = {}
    consensus_dict = {rec.id: rec for rec in SeqIO.parse(args.consensus_fasta, "fasta")}

    logger.info("Running WaveSeekerNet predictions and SHAP calculations on consensus segments...")
    for seq_id, seq_record in consensus_dict.items():
        segment = detect_segment(seq_id)
        logger.info("Processing segment: %s", segment)
        if not segment or segment not in seg_configs:
            print(f"Skipping prediction and SHAP for {seq_id}: Segment not configured.")
            continue

        cfg = seg_configs[segment]
        seq_len = cfg["seq_len"]
        weights_list = cfg["weights"]
        backgrounds_list = cfg.get("backgrounds", [])

        if seq_len < len(seq_record):
            logger.warning(f"Padding one-hot encoding length is less than actual sequence length: {seq_len} vs {len(seq_record)}. Skipping segment {seq_id}. ")
            continue
        logger.info(f"Padding one-hot encoding length vs actual sequence length: {seq_len} vs {len(seq_record)}")

        res_len = cfg.get("res_len", args.res_len)
        patch_size = tuple(cfg.get("patch_size", args.patch_size))
        emb_dim = cfg.get("emb_dim", args.emb_dim)
        batch_size = cfg.get("batch_size", args.batch_size)
        final_hidden_size = cfg.get("final_hidden_size", args.final_hidden_size)

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f_temp:
            SeqIO.write(seq_record, f_temp, "fasta")
            temp_fasta_path = f_temp.name

        try:
            X, _ = fasta_to_one_hot(temp_fasta_path, seq_len, res_len, True, out_filename=None)

            clf = WaveSeekerClassifier(
                n_channels=1, seq_L=seq_len, res_L=res_len, patch_size=patch_size,
                n_out=n_out, batch_size=batch_size, emb_dim=emb_dim,
                final_hidden_size=final_hidden_size, patch_mode="patch",
                wavelet_names=["sym4"], n_blocks=1, lr=0.0025, device=args.device
            )

            # Ensemble Prediction
            accumulated_probs = np.zeros((X.shape[0], n_out))
            individual_preds = {}
            valid_models = 0

            for idx, weight_rel_path in enumerate(weights_list):
                weight_path = os.path.join(args.weights_dir, weight_rel_path)
                if not os.path.exists(weight_path):
                    logger.warning("Model weight file not found at: %s", weight_path)
                    continue
                clf.load_weights(weight_path)
                pred, logits = clf.predict(X, return_logits=True)
                pred_idx = pred[0]
                pred_label = classes[pred_idx] if pred_idx < len(classes) else "Unknown"
                individual_preds[f"model_{idx}_class"] = pred_label

                exp_logits = np.exp(logits - np.max(logits, axis=-1, keepdims=True))
                probs = exp_logits / np.sum(exp_logits, axis=-1, keepdims=True)
                for c_idx, class_name in enumerate(classes):
                    col_class_name = class_name.replace(" ", "_").replace("-", "_")
                    individual_preds[f"model_{idx}_{col_class_name}_prob"] = f"{probs[0][c_idx]:.4f}"

                accumulated_probs += probs
                valid_models += 1

            if valid_models == 0:
                continue

            avg_probs = accumulated_probs / valid_models
            
            # Calibrate predicted probabilities using IsotonicCalibrator if validation data is available
            calibrated_probs = avg_probs.copy()
            calibration_applied = False
            
            # Read calibration paths from config, fallback to directory of weights
            logits_rel_path = cfg.get("logits")
            y_true_rel_path = cfg.get("y_true")
            
            if logits_rel_path and y_true_rel_path:
                logits_path = os.path.join(args.weights_dir, logits_rel_path)
                y_true_path = os.path.join(args.weights_dir, y_true_rel_path)
            elif weights_list:
                first_weight_path = os.path.join(args.weights_dir, weights_list[0])
                segment_dir = os.path.dirname(first_weight_path)
                logits_path = os.path.join(segment_dir, "pre2020_logits.npy")
                y_true_path = os.path.join(segment_dir, "pre2020_y_true.npy")
            else:
                logits_path = None
                y_true_path = None
                
            if logits_path and y_true_path and os.path.exists(logits_path) and os.path.exists(y_true_path):
                try:
                    logger.info("Training IsotonicCalibrator using validation logits (%s) and true labels (%s)...", logits_path, y_true_path)
                    calib_logits = np.load(logits_path)
                    calib_y = np.load(y_true_path)
                    
                    # Convert validation logits to probabilities using softmax
                    exp_calib_logits = np.exp(calib_logits - np.max(calib_logits, axis=-1, keepdims=True))
                    calib_probs = exp_calib_logits / np.sum(exp_calib_logits, axis=-1, keepdims=True)
                    
                    calibrator = IsotonicCalibrator(n_classes=n_out)
                    calibrator.fit(calib_probs, calib_y)
                    
                    # Compute validation calibration metrics before and after calibration
                    ece_pre, mce_pre, ace_pre, mcs_pre, _ = compute_calibration_metrics(calib_probs, calib_y, n_classes=n_out)
                    val_calib_probs = calibrator.predict_proba(calib_probs)
                    ece_post, mce_post, ace_post, mcs_post, _ = compute_calibration_metrics(val_calib_probs, calib_y, n_classes=n_out)
                    
                    calibrated_probs = calibrator.predict_proba(avg_probs)
                    calibration_applied = True
                    logger.info("Probability calibration successfully applied.")
                    
                    val_metrics = {
                        "val_ece_pre": ece_pre,
                        "val_ece_post": ece_post,
                        "val_ace_pre": ace_pre,
                        "val_ace_post": ace_post,
                        "val_mce_pre": mce_pre,
                        "val_mce_post": mce_post,
                        "val_mcs_pre": mcs_pre,
                        "val_mcs_post": mcs_post
                    }
                except Exception as e:
                    logger.error("Failed to perform probability calibration: %s", str(e))
                    val_metrics = {}
            else:
                val_metrics = {}
            
            # Uncalibrated metrics for reference
            uncal_pred_idx = np.argmax(avg_probs, axis=-1)[0]
            uncal_confidence = np.max(avg_probs, axis=-1)[0]
            uncal_pred_class = classes[uncal_pred_idx]

            # Final predicted metrics (calibrated if available)
            ensemble_pred_idx = np.argmax(calibrated_probs, axis=-1)[0]
            ensemble_confidence = np.max(calibrated_probs, axis=-1)[0]
            ensemble_pred_class = classes[ensemble_pred_idx]

            # Save ensemble outputs
            record = {
                "sequence_id": seq_record.id,
                "segment": segment,
                "ensemble_class": ensemble_pred_class,
                "ensemble_confidence": f"{ensemble_confidence:.4f}",
                "uncalibrated_class": uncal_pred_class,
                "uncalibrated_confidence": f"{uncal_confidence:.4f}",
                "calibration_applied": calibration_applied
            }
            
            # Save validation metrics to record
            for k, v in val_metrics.items():
                record[k] = f"{v:.4f}" if isinstance(v, float) else str(v)
                
            record.update(individual_preds)
            for c_idx, class_name in enumerate(classes):
                col_name = class_name.replace(" ", "_").replace("-", "_") + "_prob"
                uncal_col_name = "uncalibrated_" + col_name
                record[col_name] = f"{calibrated_probs[0][c_idx]:.4f}"
                record[uncal_col_name] = f"{avg_probs[0][c_idx]:.4f}"
            record["models_evaluated"] = valid_models
            pred_results.append(record)

            # Gradient SHAP computation with fold-specific backgrounds
            accumulated_shap = np.zeros((res_len, seq_len, n_out))
            for idx, weight_rel_path in enumerate(weights_list):
                weight_path = os.path.join(args.weights_dir, weight_rel_path)
                if not os.path.exists(weight_path):
                    logger.warning("Model weight file not found at: %s during SHAP calculation.", weight_path)
                    continue
                clf.load_weights(weight_path)

                # Fetch corresponding background
                background = None
                if idx < len(backgrounds_list):
                    bg_rel_path = backgrounds_list[idx]
                    bg_path = os.path.join(args.weights_dir, bg_rel_path)
                    if os.path.exists(bg_path):
                        background = np.load(bg_path)

                        max_bg_samples = 1000
                        if background.shape[0] > max_bg_samples:
                            indices = np.random.choice(background.shape[0], max_bg_samples, replace=False)
                            background = background[indices]
                        logger.info("  Background data shape: %s", str(background.shape))
                    else:
                        logger.warning("Background file not found at: %s", bg_path)

                if background is None:
                    logger.info("Falling back to zero-initialized baseline background for SHAP.")
                    background = np.zeros_like(X)

                shap_vals = clf.explain(
                    X_explain=X, background_data=background,
                    explainer_type="gradient", output_type="logits", batch_size=batch_size
                )
                # Multiply SHAP values by one-hot encoded input sequence (X) to get actual base attributions
                # X has shape (1, 5, seq_len), shap_vals[0] has shape (5, seq_len, n_out)
                x_one_hot = X[0]  # shape (5, seq_len)
                shap_vals_actual = shap_vals[0] * np.expand_dims(x_one_hot, axis=-1)
                accumulated_shap += shap_vals_actual

            avg_shap = accumulated_shap / valid_models

            # Cache SHAP matrix and predicted class
            consensus_shap_data[seq_id] = {
                "avg_shap": avg_shap,
                "ensemble_pred_class": ensemble_pred_class,
                "segment": segment,
                "X": X
            }
        finally:
            if os.path.exists(temp_fasta_path):
                os.remove(temp_fasta_path)

    # Save consolidated TSV file
    pred_dir = os.path.dirname(args.output_predictions)
    if pred_dir:
        os.makedirs(pred_dir, exist_ok=True)
    df_pred = pd.DataFrame(pred_results)
    df_pred.to_csv(args.output_predictions, sep='\t', index=False)
    logger.info("Ensemble predictions saved to: %s", args.output_predictions)

    # 2. Slice SHAP for each CDS sequence in .ffn
    logger.info("Slicing CDS coordinates and generating SHAP outputs...")
    for cds_record in SeqIO.parse(args.cds_fasta, "fasta"):
        fields = cds_record.id.split('|')
        if len(fields) >= 2:
            if fields[1] != "CDS":
                continue
        else:
            if "CDS" not in cds_record.id.upper():
                continue

        consensus_record = find_matching_consensus(cds_record.id, consensus_dict)
        logger.info("Processing CDS: %s in Consensus: %s", cds_record.id, consensus_record.id if consensus_record else "None")
        if not consensus_record or consensus_record.id not in consensus_shap_data:
            logger.warning("No precomputed SHAP data found for %s. Skipping.", cds_record.id)
            continue

        shap_info = consensus_shap_data[consensus_record.id]
        avg_shap = shap_info["avg_shap"]
        ensemble_pred_class = shap_info["ensemble_pred_class"]
        segment = shap_info["segment"]

        cfg = seg_configs[segment]
        seq_len = cfg["seq_len"]

        shap_dir = os.path.dirname(args.output_shap_prefix)
        shap_base = os.path.basename(args.output_shap_prefix)
        # Organize outputs into {sample_name}/{segment}/ directories
        segment_dir = os.path.join(shap_dir, shap_base, segment) if shap_dir else os.path.join(shap_base, segment)
        os.makedirs(segment_dir, exist_ok=True)
        clean_cds_id = cds_record.id.replace('|', '_')

        ranges, is_complement = parse_cds_coordinates(cds_record.description)
        logger.info("CDS ranges: %s, is_complement: %s", str(ranges), str(is_complement))
        if not ranges:
            ranges, is_complement = parse_cds_coordinates(cds_record.id)
        if not ranges:
            continue

        # Output per class
        for c_idx, class_name in enumerate(classes):
            con_seq_str = str(consensus_record.seq).upper()
            n_bases = len(con_seq_str)
            seq_shap = np.zeros(n_bases)
            base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'U': 3}

            for p in range(min(n_bases, seq_len)):
                base = con_seq_str[p]
                if base in base_to_idx:
                    idx = base_to_idx[base]
                    seq_shap[p] = avg_shap[idx, p, c_idx]

            cds_shap_list = []
            cds_seq_chars = []

            if is_complement:
                for start, end in reversed(ranges):
                    slice_start = max(0, start - 1)
                    slice_end = min(n_bases, end)
                    cds_shap_list.extend(seq_shap[slice_start: slice_end][::-1])
                    for base in reversed(con_seq_str[slice_start: slice_end]):
                        comp = {'A': 'T', 'T': 'A', 'U': 'A', 'G': 'C', 'C': 'G', 'N': 'N', '-': '-'}.get(base, base)
                        cds_seq_chars.append(comp)
            else:
                for start, end in ranges:
                    slice_start = max(0, start - 1)
                    slice_end = min(n_bases, end)
                    cds_shap_list.extend(seq_shap[slice_start: slice_end])
                    cds_seq_chars.extend(list(con_seq_str[slice_start: slice_end]))

            cds_shap = np.array(cds_shap_list)
            cds_seq = "".join(cds_seq_chars)

            # Validate coordinates matches ffn sequence
            cds_fasta_str = str(cds_record.seq).upper()
            if cds_seq != cds_fasta_str:
                logger.warning("Warning: Extracted coordinate seq (%s bp) does not match ffn seq (%s bp)!", len(cds_seq), len(cds_fasta_str))

            n_codons = len(cds_seq) // 3
            codon_shap_sums = {codon: 0.0 for codon in SENSE_CODONS}
            codon_counts = {codon: 0 for codon in SENSE_CODONS}

            for c_pos in range(n_codons):
                codon = cds_seq[c_pos * 3: (c_pos + 1) * 3]
                if codon in codon_shap_sums:
                    val = cds_shap[c_pos * 3] + cds_shap[c_pos * 3 + 1] + cds_shap[c_pos * 3 + 2]
                    codon_shap_sums[codon] += val
                    codon_counts[codon] += 1

            # Save TSV
            shap_series = pd.Series(codon_shap_sums)
            clean_class_name = class_name.replace(' ', '_')
            output_tsv = os.path.join(segment_dir, f"{shap_base}_{clean_cds_id}_shap_{clean_class_name}.tsv")
            rscu_dict = calculate_rscu(codon_counts)
            df_out = pd.DataFrame({
                "codon": shap_series.index,
                "shap_score": shap_series.values,
                "frequency": [codon_counts[c] for c in shap_series.index],
                "rscu": [rscu_dict[c] for c in shap_series.index]
            })
            df_out.to_csv(output_tsv, sep='\t', index=False)

            # Calculate nucleotide-level SHAP scores
            nuc_shap_sums = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}
            nuc_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            nuc_vals = {b: [] for b in ['A', 'C', 'G', 'T']}

            for char, val in zip(cds_seq, cds_shap):
                if char in nuc_shap_sums:
                    nuc_shap_sums[char] += val
                    nuc_counts[char] += 1
                    nuc_vals[char].append(val)

            nuc_shap_means = {
                b: (nuc_shap_sums[b] / nuc_counts[b] if nuc_counts[b] > 0 else 0.0)
                for b in ['A', 'C', 'G', 'T']
            }

            # Save Nucleotide TSV
            output_nuc_tsv = os.path.join(segment_dir, f"{shap_base}_{clean_cds_id}_shap_nuc_{clean_class_name}.tsv")
            df_nuc_out = pd.DataFrame({
                "nucleotide": ['A', 'C', 'G', 'T'],
                "shap_sum": [nuc_shap_sums[b] for b in ['A', 'C', 'G', 'T']],
                "shap_mean": [nuc_shap_means[b] for b in ['A', 'C', 'G', 'T']],
                "frequency": [nuc_counts[b] for b in ['A', 'C', 'G', 'T']]
            })
            df_nuc_out.to_csv(output_nuc_tsv, sep='\t', index=False)

            # Save Nucleotide Plot
            fig_nuc, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
            bases = ['A', 'C', 'G', 'T']
            colors = ['#4CAF50', '#2196F3', '#FF9800', '#F44336'] # A: Green, C: Blue, G: Orange, T: Red

            # Panel 1: Sum of SHAP values
            shap_sums = [nuc_shap_sums[b] for b in bases]
            bars = ax1.bar(bases, shap_sums, color=colors, edgecolor='#333333', linewidth=1.2, alpha=0.85)
            ax1.set_ylabel('Total (Sum) SHAP Value', fontsize=12)
            ax1.set_title('Sum of SHAP Scores per Nucleotide', fontsize=14, fontweight='bold', pad=10)
            ax1.grid(axis='y', linestyle='--', alpha=0.5)
            ax1.axhline(0, color='black', linewidth=0.8, linestyle='-')

            # Annotate bar values
            for bar in bars:
                height = bar.get_height()
                va = 'bottom' if height >= 0 else 'top'
                xytext = (0, 3) if height >= 0 else (0, -12)
                ax1.annotate(f'{height:.2f}',
                             xy=(bar.get_x() + bar.get_width() / 2, height),
                             xytext=xytext, textcoords="offset points",
                             ha='center', va=va, fontsize=10, fontweight='semibold')

            # Panel 2: Boxplot of SHAP distributions
            data_to_plot = [nuc_vals[b] if len(nuc_vals[b]) > 0 else [0.0] for b in bases]
            bp = ax2.boxplot(data_to_plot, patch_artist=True, showfliers=False,
                             boxprops=dict(linewidth=1.2, edgecolor='#333333'),
                             whiskerprops=dict(linewidth=1.2, color='#333333'),
                             capprops=dict(linewidth=1.2, color='#333333'),
                             medianprops=dict(linewidth=1.5, color='#111111'))
            ax2.set_xticks(range(1, len(bases) + 1))
            ax2.set_xticklabels(bases)

            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)

            ax2.set_ylabel('SHAP Value per Position', fontsize=12)
            ax2.set_title('Distribution of SHAP Scores per Position', fontsize=14, fontweight='bold', pad=10)
            ax2.grid(axis='y', linestyle='--', alpha=0.5)
            ax2.axhline(0, color='black', linewidth=0.8, linestyle='-')

            fig_nuc.suptitle(f"{cds_record.id} ({segment}) - Nucleotide SHAP Analysis\nClass: {class_name} (Predicted: {ensemble_pred_class})",
                             fontsize=14, fontweight='bold', y=0.98)
            plt.tight_layout(rect=[0, 0, 1, 0.95])

            output_nuc_png = os.path.join(segment_dir, f"{shap_base}_{clean_cds_id}_shap_nuc_{clean_class_name}.png")
            plt.savefig(output_nuc_png, dpi=150)
            plt.close(fig_nuc)

            # Save Plot
            fig, ax = plt.subplots(figsize=(15, 6))
            plot_shap_given_ax(
                ax=ax, df=codon_shap_sums,
                text=f"{cds_record.id} ({segment}) - SHAP: {class_name} (Predicted: {ensemble_pred_class})",
                is_at_rich_fn=is_at_rich, show_legend=True, rotate_xticks=90
            )
            plt.tight_layout()
            output_png = os.path.join(segment_dir, f"{shap_base}_{clean_cds_id}_shap_{clean_class_name}.png")
            plt.savefig(output_png, dpi=150)
            plt.close()

            if segment == "PB2":
                # Build logo_array of shape (len(cds_seq), 4)
                cds_logo_list = []
                if is_complement:
                    for start, end in reversed(ranges):
                        slice_start = max(0, start - 1)
                        slice_end = min(n_bases, end)
                        for p in reversed(range(slice_start, slice_end)):
                            logo_vals = [avg_shap[3, p, c_idx], avg_shap[2, p, c_idx], avg_shap[1, p, c_idx], avg_shap[0, p, c_idx]]
                            cds_logo_list.append(logo_vals)
                else:
                    for start, end in ranges:
                        slice_start = max(0, start - 1)
                        slice_end = min(n_bases, end)
                        for p in range(slice_start, slice_end):
                            logo_vals = [avg_shap[0, p, c_idx], avg_shap[1, p, c_idx], avg_shap[2, p, c_idx], avg_shap[3, p, c_idx]]
                            cds_logo_list.append(logo_vals)
                logo_array = np.array(cds_logo_list)

                # Calculate codon-level position-specific SHAP scores
                scatter_data = {}
                for c_pos in range(n_codons):
                    val = cds_shap[c_pos * 3] + cds_shap[c_pos * 3 + 1] + cds_shap[c_pos * 3 + 2]
                    scatter_data[c_pos + 1] = val
                # Get top positive positions
                pos_scores_nucs = []
                for i, (nuc, score) in enumerate(zip(cds_seq, cds_shap)):
                    pos_scores_nucs.append((i + 1, score, nuc))
                top_pos_data = sorted(pos_scores_nucs, key=lambda x: x[1], reverse=True)[:20]

                # 1. Genotype Screening for E627K and D701N
                codon_627 = cds_seq[1878:1881] if len(cds_seq) >= 1881 else "N/A"
                codon_701 = cds_seq[2100:2103] if len(cds_seq) >= 2103 else "N/A"

                # Map 627
                if codon_627 in ['GAA', 'GAG']:
                    aa_627 = 'E (Glutamic Acid - Avian/WT)'
                elif codon_627 in ['AAA', 'AAG']:
                    aa_627 = 'K (Lysine - Mammalian/Mutant)'
                elif codon_627 == "N/A":
                    aa_627 = 'Sequence too short'
                else:
                    aa_627 = f'Other ({codon_627})'

                # Map 701
                if codon_701 in ['GAC', 'GAT']:
                    aa_701 = 'D (Aspartic Acid - Avian/WT)'
                elif codon_701 in ['AAC', 'AAT']:
                    aa_701 = 'N (Asparagine - Mammalian/Mutant)'
                elif codon_701 == "N/A":
                    aa_701 = 'Sequence too short'
                else:
                    aa_701 = f'Other ({codon_701})'

                # 2. E627K Dashboard
                plot_dashboard_grid(
                    scatter_data=scatter_data,
                    logo_array=logo_array,
                    bar_data=codon_shap_sums,
                    top_pos_data=top_pos_data,
                    codon=['627', 'E', 'K'],
                    fig_title=f"{cds_record.id} ({segment}) - E627K Mutation Screening Dashboard ({class_name})\n[Genotype at position 627: {codon_627} -> {aa_627}]",
                    save_path=os.path.join(segment_dir, f"{shap_base}_{clean_cds_id}_shap_dashboard_E627K_{clean_class_name}.png"),
                    scatter_markers=[{'pos': 627, 'label': '627', 'color': 'red'}],
                    logo_region=(1863, 1896), # codons 622 to 632
                    logo_highlight={'red': [(1878, 1881, 'Codon 627 (E/K)')]},
                    is_at_rich_fn=is_at_rich
                )

                # 3. D701N Dashboard
                plot_dashboard_grid(
                    scatter_data=scatter_data,
                    logo_array=logo_array,
                    bar_data=codon_shap_sums,
                    top_pos_data=top_pos_data,
                    codon=['701', 'D', 'N'],
                    fig_title=f"{cds_record.id} ({segment}) - D701N Mutation Screening Dashboard ({class_name})\n[Genotype at position 701: {codon_701} -> {aa_701}]",
                    save_path=os.path.join(segment_dir, f"{shap_base}_{clean_cds_id}_shap_dashboard_D701N_{clean_class_name}.png"),
                    scatter_markers=[{'pos': 701, 'label': '701', 'color': 'red'}],
                    logo_region=(2085, 2118), # codons 696 to 706
                    logo_highlight={'red': [(2100, 2103, 'Codon 701 (D/N)')]},
                    is_at_rich_fn=is_at_rich
                )

    # 3. Generate HTML Report
    output_html = getattr(args, "output_html", None)
    if not output_html:
        output_html = args.output_predictions.replace(".tsv", ".html")
    
    try:
        sample_name = os.path.basename(args.output_predictions).split(".")[0]
        generate_html_report(pred_results, output_html, sample_name, args.output_shap_prefix)
    except Exception as e:
        logger.error("Failed to generate HTML report: %s", str(e))

    print("WaveSeekerNet analysis complete.")


def generate_html_report(pred_results, output_html_path, sample_name, shap_prefix):
    """
    Generates a beautiful, standalone HTML report summarizing the WaveSeekerNet predictions
    and embedding the generated SHAP and mutation screening figures as Base64 data URIs.
    """
    import glob
    import base64
    
    def get_base64_image(image_path):
        try:
            with open(image_path, "rb") as img_file:
                encoded = base64.b64encode(img_file.read()).decode('utf-8')
                return f"data:image/png;base64,{encoded}"
        except Exception as e:
            logger.error("Failed to base64 encode image %s: %s", image_path, str(e))
            return ""
            
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>WaveSeekerNet Analysis Report - {sample_name}</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background-color: #f5f7fa;
            color: #333;
            margin: 0;
            padding: 20px;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 12px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.05);
        }}
        h1 {{
            color: #1e3a8a;
            border-bottom: 2px solid #e2e8f0;
            padding-bottom: 10px;
            margin-top: 0;
        }}
        h2 {{
            color: #2b6cb0;
            margin-top: 30px;
            border-bottom: 1px solid #e2e8f0;
            padding-bottom: 8px;
        }}
        .summary-box {{
            background-color: #ebf8ff;
            border-left: 4px solid #3182ce;
            padding: 15px;
            margin-bottom: 25px;
            border-radius: 4px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
            margin-bottom: 25px;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #e2e8f0;
        }}
        th {{
            background-color: #f7fafc;
            color: #4a5568;
            font-weight: 600;
        }}
        tr:hover {{
            background-color: #f8fafc;
        }}
        .badge {{
            display: inline-block;
            padding: 4px 8px;
            border-radius: 12px;
            font-size: 12px;
            font-weight: bold;
            text-transform: uppercase;
        }}
        .badge-avian {{ background-color: #feebc8; color: #c05621; }}
        .badge-human {{ background-color: #c6f6d5; color: #22543d; }}
        .badge-mammal {{ background-color: #e2e8f0; color: #4a5568; }}
        .badge-yes {{ background-color: #c6f6d5; color: #22543d; }}
        .badge-no {{ background-color: #fed7d7; color: #9b2c2c; }}
        
        .grid {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
            margin-top: 20px;
        }}
        .card {{
            background: #fff;
            border: 1px solid #e2e8f0;
            border-radius: 8px;
            padding: 15px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.02);
            text-align: center;
        }}
        .card img {{
            max-width: 100%;
            height: auto;
            border-radius: 4px;
            margin-top: 10px;
            border: 1px solid #edf2f7;
        }}
        .footer {{
            margin-top: 50px;
            text-align: center;
            font-size: 12px;
            color: #a0aec0;
            border-top: 1px solid #e2e8f0;
            padding-top: 15px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>WaveSeekerNet Analysis Report</h1>
        <div class="summary-box">
            <strong>Sample ID:</strong> {sample_name}<br>
            <strong>Date:</strong> {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}<br>
            <strong>Pipeline Version:</strong> 1.0 (with Probability Calibration)
        </div>

        <h2>Ensemble Predictions Summary</h2>
        <table>
            <thead>
                <tr>
                    <th>Segment</th>
                    <th>Consensus Sequence ID</th>
                    <th>Calibrated Prediction</th>
                    <th>Confidence</th>
                    <th>Uncalibrated Prediction</th>
                    <th>Uncalibrated Conf.</th>
                    <th>Calibration Applied</th>
                </tr>
            </thead>
            <tbody>
"""
    
    for r in pred_results:
        seg = r.get("segment", "N/A")
        seq_id = r.get("sequence_id", "N/A")
        ens_class = r.get("ensemble_class", "N/A")
        ens_conf = r.get("ensemble_confidence", "0.0000")
        uncal_class = r.get("uncalibrated_class", "N/A")
        uncal_conf = r.get("uncalibrated_confidence", "0.0000")
        cal_applied = "Yes" if r.get("calibration_applied", False) else "No"
        
        def get_badge_class(cls_name):
            if "avian" in cls_name.lower(): return "badge-avian"
            if "human" in cls_name.lower(): return "badge-human"
            return "badge-mammal"
            
        badge_cal = get_badge_class(ens_class)
        badge_uncal = get_badge_class(uncal_class)
        badge_applied = "badge-yes" if cal_applied == "Yes" else "badge-no"
        
        html_content += f"""
                <tr>
                    <td><strong>{seg}</strong></td>
                    <td>{seq_id}</td>
                    <td><span class="badge {badge_cal}">{ens_class}</span></td>
                    <td><strong>{float(ens_conf)*100:.2f}%</strong></td>
                    <td><span class="badge {badge_uncal}">{uncal_class}</span></td>
                    <td>{float(uncal_conf)*100:.2f}%</td>
                    <td><span class="badge {badge_applied}">{cal_applied}</span></td>
                </tr>
"""
        
    html_content += """
            </tbody>
        </table>
"""

    # Add Calibration Quality Section
    html_content += """
        <h2>Calibration Quality & Validation Performance</h2>
        <p>This section displays the calibration quality metrics evaluated on the validation cohort for each segment before and after applying Isotonic Regression. The Expected Calibration Error (ECE) and Adaptive Calibration Error (ACE) quantify the reliability of the confidence scores (lower is better, perfect is 0.0000).</p>
        <table>
            <thead>
                <tr>
                    <th>Segment</th>
                    <th>ECE (Before &rarr; After)</th>
                    <th>ACE (Before &rarr; After)</th>
                    <th>MCE (Before &rarr; After)</th>
                    <th>MCS (Before &rarr; After)</th>
                    <th>Calibration Status</th>
                </tr>
            </thead>
            <tbody>
"""
    for r in pred_results:
        seg = r.get("segment", "N/A")
        if r.get("calibration_applied", False) and "val_ece_pre" in r:
            ece_pre = r.get("val_ece_pre", "N/A")
            ece_post = r.get("val_ece_post", "N/A")
            ace_pre = r.get("val_ace_pre", "N/A")
            ace_post = r.get("val_ace_post", "N/A")
            mce_pre = r.get("val_mce_pre", "N/A")
            mce_post = r.get("val_mce_post", "N/A")
            mcs_pre = r.get("val_mcs_pre", "N/A")
            mcs_post = r.get("val_mcs_post", "N/A")
            
            try:
                ece_val = float(ece_post)
                if ece_val <= 0.05:
                    status_badge = '<span class="badge badge-yes">🟢 GOOD</span>'
                elif ece_val <= 0.10:
                    status_badge = '<span class="badge badge-mammal">🟡 FAIR</span>'
                else:
                    status_badge = '<span class="badge badge-no">🔴 POOR</span>'
            except ValueError:
                status_badge = '<span class="badge badge-no">N/A</span>'
                
            html_content += f"""
                <tr>
                    <td><strong>{seg}</strong></td>
                    <td>{ece_pre} &rarr; <strong>{ece_post}</strong></td>
                    <td>{ace_pre} &rarr; <strong>{ace_post}</strong></td>
                    <td>{mce_pre} &rarr; <strong>{mce_post}</strong></td>
                    <td>{mcs_pre} &rarr; <strong>{mcs_post}</strong></td>
                    <td>{status_badge}</td>
                </tr>
"""
        else:
            html_content += f"""
                <tr>
                    <td><strong>{seg}</strong></td>
                    <td colspan="4" style="color: #a0aec0; font-style: italic;">No calibration data available (skipped or missing validation files)</td>
                    <td><span class="badge badge-no">Skipped</span></td>
                </tr>
"""
            
    html_content += """
            </tbody>
        </table>
"""

    has_pb2 = any(r.get("segment") == "PB2" for r in pred_results)
    if has_pb2:
        html_content += """
        <h2>PB2 Mutation Screening (Host Adaptation Markers)</h2>
        <p>The PB2 segment of Influenza A virus contains key markers associated with mammalian adaptation and increased virulence, specifically at positions 627 and 701.</p>
        <div class="grid">
        """
        
        shap_dir = os.path.dirname(output_html_path)
        e627k_files = glob.glob(os.path.join(shap_dir, shap_prefix, "PB2", "*_shap_dashboard_E627K_*.png"))
        d701n_files = glob.glob(os.path.join(shap_dir, shap_prefix, "PB2", "*_shap_dashboard_D701N_*.png"))
        
        for f in e627k_files:
            img_b64 = get_base64_image(f)
            if img_b64:
                html_content += f"""
                <div class="card">
                    <h3>E627K Adaptation Marker Dashboard</h3>
                    <img src="{img_b64}" alt="E627K Dashboard">
                </div>
                """
        for f in d701n_files:
            img_b64 = get_base64_image(f)
            if img_b64:
                html_content += f"""
                <div class="card">
                    <h3>D701N Adaptation Marker Dashboard</h3>
                    <img src="{img_b64}" alt="D701N Dashboard">
                </div>
                """
            
        html_content += """
        </div>
        """
        
    html_content += """
        <h2>Segment-specific SHAP Feature Attributions</h2>
        <p>The following figures show the codon and nucleotide-level feature importance (SHAP values) for the predicted classes.</p>
        <div class="grid">
    """
    
    shap_dir = os.path.dirname(output_html_path)
    all_pngs = glob.glob(os.path.join(shap_dir, shap_prefix, "**", "*.png"), recursive=True)
    shap_pngs = [f for f in all_pngs if "_shap_" in f and "dashboard" not in f and "nuc" not in f]
    
    for f in sorted(shap_pngs):
        img_b64 = get_base64_image(f)
        if img_b64:
            filename = os.path.basename(f)
            parts = filename.split("_")
            seg_name = os.path.basename(os.path.dirname(f))
            
            html_content += f"""
            <div class="card">
                <h3>{seg_name} Codon SHAP ({parts[-1].replace('.png', '')})</h3>
                <img src="{img_b64}" alt="{filename}">
            </div>
            """
        
    html_content += """
        </div>
        <div class="footer">
            Report generated by WaveSeekerNet Deep Learning Suite. All rights reserved.
        </div>
    </div>
</body>
</html>
"""
    with open(output_html_path, "w") as f_html:
        f_html.write(html_content)
    logger.info("HTML report saved to: %s", output_html_path)


if __name__ == "__main__":
    main()
