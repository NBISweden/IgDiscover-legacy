"""
Augment AIRR-formatted IgBLAST output with extra IgDiscover-specific columns

In particular, detect the CDR3 sequence
"""
import json
import sys
import os
import time
import logging

import numpy as np
import pandas as pd
import dnaio
from tinyalign import edit_distance, hamming_distance

from .igblast import igblast
from ..utils import nt_to_aa
from ..parse import parse_header
from ..species import cdr3_start, cdr3_end


logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    arg(
        "--sequence-type",
        default="Ig",
        choices=("Ig", "TCR"),
        help="Sequence type. Default: %(default)s",
    )
    arg(
        "--rename",
        default=None,
        metavar="PREFIX",
        help="Rename reads to PREFIXseqN (where N is a number starting at 1)",
    )
    arg("--stats", metavar="FILE", help="Write statistics in JSON format to FILE")

    arg("database", help="Database directory with V.fasta, D.fasta, J.fasta.")
    arg("table", help="AIRR rearrangement table")


class Database:
    def __init__(self, path, sequence_type):
        """path -- path to database directory with V.fasta, D.fasta, J.fasta"""
        self.path = path
        self.sequence_type = sequence_type
        self._v_records = self._read_fasta(os.path.join(path, "V.fasta"))
        self.v = self._records_to_dict(self._v_records)
        self._d_records = self._read_fasta(os.path.join(path, "D.fasta"))
        self.d = self._records_to_dict(self._d_records)
        self._j_records = self._read_fasta(os.path.join(path, "J.fasta"))
        self.j = self._records_to_dict(self._j_records)
        self._cdr3_starts = dict()
        self._cdr3_ends = dict()
        for chain in ("heavy", "kappa", "lambda", "alpha", "beta", "gamma", "delta"):
            self._cdr3_starts[chain] = {
                name: cdr3_start(s, chain) for name, s in self.v.items()
            }
            self._cdr3_ends[chain] = {
                name: cdr3_end(s, chain) for name, s in self.j.items()
            }
        self.v_regions_nt, self.v_regions_aa = self._find_v_regions()

    @staticmethod
    def _read_fasta(path):
        records = []
        with dnaio.open(path) as sr:
            for record in sr:
                record.name = record.name.split(maxsplit=1)[0]
                records.append(record)
        return records

    @staticmethod
    def _records_to_dict(records):
        return {record.name: record.sequence.upper() for record in records}

    def v_cdr3_start(self, gene, chain):
        return self._cdr3_starts[chain][gene]

    def j_cdr3_end(self, gene, chain):
        return self._cdr3_ends[chain][gene]

    def _find_v_regions(self):
        """
        Run IgBLAST on the V sequences to determine the nucleotide and amino-acid sequences of the
        FR1, CDR1, FR2, CDR2 and FR3 regions
        """
        v_regions_nt = dict()
        v_regions_aa = dict()
        for record in igblast(
            self.path, self._v_records, self.sequence_type, threads=1
        ):
            nt_regions = dict()
            aa_regions = dict()
            for region in ("FR1", "CDR1", "FR2", "CDR2", "FR3"):
                nt_seq = record.region_sequence(region)
                if nt_seq is None:
                    break
                if len(nt_seq) % 3 != 0:
                    logger.warning(
                        "Length %s of %s region in %r is not divisible by three; region "
                        "info for this gene will not be available",
                        len(nt_seq),
                        region,
                        record.query_name,
                    )
                    # not codon-aligned, skip entire record
                    break
                nt_regions[region] = nt_seq
                try:
                    aa_seq = nt_to_aa(nt_seq)
                except ValueError as e:
                    logger.warning(
                        "The %s region could not be converted to amino acids: %s",
                        region,
                        str(e),
                    )
                    break
                if "*" in aa_seq:
                    logger.warning(
                        "The %s region in %r contains a stop codon (%r); region info "
                        "for this gene will not be available",
                        region,
                        record.query_name,
                        aa_seq,
                    )
                    break
                aa_regions[region] = aa_seq
            else:
                v_regions_nt[record.query_name] = nt_regions
                v_regions_aa[record.query_name] = aa_regions

        return v_regions_nt, v_regions_aa


def count_errors(s, t):
    if pd.isna(s) or pd.isna(t):
        return np.nan
    if len(s) != len(t):
        raise ValueError("alignment sequences must be of the same length")
    return sum(1 for c, d in zip(s, t) if c != d)


COLUMN_TYPES = {
    "v_alignment_start": "Int64",
    "v_alignment_end": "Int64",
    "d_alignment_start": "Int64",
    "d_alignment_end": "Int64",
    "j_alignment_start": "Int64",
    "j_alignment_end": "Int64",
    "v_sequence_start": "Int64",
    "v_sequence_end": "Int64",
    "v_germline_start": "Int64",
    "v_germline_end": "Int64",
    "d_sequence_start": "Int64",
    "d_sequence_end": "Int64",
    "d_germline_start": "Int64",
    "d_germline_end": "Int64",
    "j_sequence_start": "Int64",
    "j_sequence_end": "Int64",
    "j_germline_start": "Int64",
    "j_germline_end": "Int64",
    "fwr1_start": "Int64",
    "fwr1_end": "Int64",
    "cdr1_start": "Int64",
    "cdr1_end": "Int64",
    "fwr2_start": "Int64",
    "fwr2_end": "Int64",
    "cdr2_start": "Int64",
    "cdr2_end": "Int64",
    "fwr3_start": "Int64",
    "fwr3_end": "Int64",
    "fwr4_start": "Int64",
    "fwr4_end": "Int64",
    "cdr3_start": "Int64",
    "cdr3_end": "Int64",
}


def main(args):
    logger.info("Running IgBLAST on database sequences to find CDR/FR region locations")
    database = Database(args.database, args.sequence_type)
    start_time = time.time()
    last_status_update = 0
    n = 0
    detected_cdr3s = 0
    first = True
    for table in pd.read_table(args.table, sep="\t", chunksize=1000, dtype=COLUMN_TYPES):
        table = augment_table(table, database)
        if args.rename is not None:
            table["sequence_id"] = [f"seq{i+1}" for i in table.index]
        try:
            table.to_csv(sys.stdout, index=False, header=first, sep="\t")
        except BrokenPipeError:
            sys.exit(1)
        n += len(table)
        detected_cdr3s += sum(table.cdr3.notna())
        first = False

        elapsed = time.time() - start_time
        if elapsed >= last_status_update + 3:  # 60:
            logger.info(
                "Processed {:10,d} sequences at {:.3f} ms/sequence".format(
                    n, elapsed / n * 1e3
                )
            )
            last_status_update = elapsed

    elapsed = time.time() - start_time
    logger.info(
        "Processed {:10,d} sequences at {:.1f} ms/sequence".format(n, elapsed / n * 1e3)
    )

    logger.info("%d IgBLAST assignments parsed and written", n)
    logger.info("CDR3s detected in %.1f%% of all sequences", detected_cdr3s / n * 100)
    if args.stats:
        stats = {"total": n, "detected_cdr3s": detected_cdr3s}
        with open(args.stats, "w") as f:
            json.dump(stats, f)
            print(file=f)


def augment_table(table, database):
    # Fix database sequence identifiers that had to be added to prevent IgBLAST from mangling
    # some identifiers
    for col in "v_call", "d_call", "j_call":
        table[col] = table[col].str.lstrip("%")

    check_table(table, database)

    # Minor transformation needed
    table["V_nt"] = table["v_sequence_alignment"].str.replace("-", "")
    table["V_aa"] = table["v_sequence_alignment_aa"]
    table["D_region"] = table["d_sequence_alignment"].str.replace("-", "")
    table["J_nt"] = table["j_sequence_alignment"].str.replace("-", "")
    table["J_aa"] = table["j_sequence_alignment_aa"]

    table["V_SHM"] = 100.0 - table["v_identity"]
    table["J_SHM"] = 100.0 - table["j_identity"]

    locus_to_chain = {
        "IGH": "VH",
        "IGK": "VK",
        "IGL": "VL",
        "TRA": "VA",
        "TRB": "VB",
        "TRG": "VG",
        "TRD": "VD",
    }
    table["chain"] = table["locus"].map(locus_to_chain, na_action="ignore")

    parsed_headers = table["sequence_id"].apply(lambda s: pd.Series(parse_header(s)))
    table["sequence_id"] = parsed_headers[0]
    table["count"] = parsed_headers[1]  # TODO consensus_count/duplicate_count
    table["barcode"] = parsed_headers[2]
    table["stop"] = table["stop_codon"].map({"F": "no", "T": "yes"})

    table["V_errors"] = table.apply(
        lambda row: count_errors(row.v_germline_alignment, row.v_sequence_alignment),
        axis=1,
    )
    table["D_errors"] = table.apply(
        lambda row: count_errors(row.d_germline_alignment, row.d_sequence_alignment),
        axis=1,
    )
    table["J_errors"] = table.apply(
        lambda row: count_errors(row.j_germline_alignment, row.j_sequence_alignment),
        axis=1,
    )

    def covered(alignment, db, call):
        if pd.isna(alignment) or pd.isna(call):
            return np.nan
        return len(alignment.replace("-", "")) / len(db[call]) * 100.0

    table["V_covered"] = table.apply(
        lambda row: covered(row.v_germline_alignment, database.v, row.v_call), axis=1
    )
    table["D_covered"] = table.apply(
        lambda row: covered(row.d_germline_alignment, database.d, row.d_call), axis=1
    )
    table["J_covered"] = table.apply(
        lambda row: covered(row.j_germline_alignment, database.j, row.j_call), axis=1
    )

    def vdj_nt(row):
        if pd.isna(row.v_call) or pd.isna(row.j_call):
            return np.nan

        return row.sequence[row.v_sequence_start - 1 : row.j_sequence_end]

    table["VDJ_nt"] = table.apply(vdj_nt, axis=1)
    table["VDJ_aa"] = table["VDJ_nt"].map(nt_to_aa, na_action="ignore")

    def v_regions_shm(row):
        """
        Compute SHM (actually mutation rate on nucleotide level) for
        all regions on V
        """
        result = {}
        for airr_col, region in (
            ("fwr1", "FR1"),
            ("cdr1", "CDR1"),
            ("fwr2", "FR2"),
            ("cdr2", "CDR2"),
            ("fwr3", "FR3"),
        ):
            start = getattr(row, airr_col + "_start")
            end = getattr(row, airr_col + "_end")
            if pd.isna(start) or pd.isna(end):
                result[region + "_SHM"] = np.nan
                continue
            sequence = row.sequence[start - 1 : end]
            germline = database.v_regions_nt[row.v_call].get(region)
            dist = edit_distance(germline, sequence)
            result[region + "_SHM"] = 100.0 * dist / len(germline)

        return pd.Series(result)

    def v_regions_aa_mut(row):
        """
        Compute amino acid mutation rate for all regions on V and also for V
        itself as the sum of the regions (that is, excluding the CDR3)
        """
        total_length = 0
        total_dist = 0
        n_regions = 0
        result = {}
        for airr_col, region in (
            ("fwr1", "FR1"),
            ("cdr1", "CDR1"),
            ("fwr2", "FR2"),
            ("cdr2", "CDR2"),
            ("fwr3", "FR3"),
        ):
            result[region + "_aa_mut"] = np.nan
            start = getattr(row, airr_col + "_start")
            end = getattr(row, airr_col + "_end")
            if pd.isna(start) or pd.isna(end):
                continue
            sequence_aa = nt_to_aa(row.sequence[start - 1 : end])
            germline_aa = database.v_regions_aa[row.v_call].get(region)
            if germline_aa is None:
                continue
            # Some FR1 alignments are reported with a frameshift by IgBLAST. By requiring that
            # reference and query lengths are identical, we can filter out these cases (and use
            # Hamming distance to get some speedup)
            if len(germline_aa) != len(sequence_aa):
                continue
            dist = hamming_distance(germline_aa, sequence_aa)
            mut_aa = dist / len(germline_aa)
            if mut_aa >= 0.8:
                # assume something went wrong
                continue
            total_dist += dist
            n_regions += 1
            total_length += len(germline_aa)
            result[region + "_aa_mut"] = 100.0 * mut_aa
        if n_regions == 5:
            result["V_aa_mut"] = 100.0 * total_dist / total_length
        else:
            result["V_aa_mut"] = np.nan
        return pd.Series(result)

    table = pd.concat(
        [
            table,
            table.apply(v_regions_shm, axis=1),
            table.apply(v_regions_aa_mut, axis=1),
        ],
        axis=1,
    )
    locus_to_long_chain = {
        "IGH": "heavy",
        "IGK": "kappa",
        "IGL": "lambda",
        "TRA": "alpha",
        "TRB": "beta",
        "TRG": "gamma",
        "TRD": "delta",
    }

    def find_cdr3(row):
        no_result = pd.Series({
            "cdr3_start": np.nan,
            "cdr3_end": np.nan,
            "cdr3": np.nan,
            "cdr3_aa": np.nan,
        })
        if pd.isna(row.v_call) or pd.isna(row.j_call) or row.locus not in locus_to_chain:
            return no_result
        # CDR3 start
        cdr3_ref_start = database.v_cdr3_start(
            row.v_call, locus_to_long_chain[row.locus]
        )
        if cdr3_ref_start is None:
            return no_result
        cdr3_query_start = query_position(row, "v", reference_position=cdr3_ref_start)
        if cdr3_query_start is None:
            # Alignment is not long enough to cover CDR3 start position; try to rescue it
            # by assuming that the alignment would continue without indels.
            cdr3_query_start = row.v_sequence_end + (
                cdr3_ref_start - row.v_germline_end
            )

        # CDR3 end
        cdr3_ref_end = database.j_cdr3_end(row.j_call, locus_to_long_chain[row.locus])
        if cdr3_ref_end is None:
            return no_result

        cdr3_query_end = query_position(row, "j", reference_position=cdr3_ref_end)
        if cdr3_query_end is None:
            return no_result
        cdr3_nt = row.sequence[cdr3_query_start:cdr3_query_end]
        return pd.Series({
            "cdr3_start": cdr3_query_start + 1,
            "cdr3_end": cdr3_query_end,
            "cdr3": cdr3_nt,
            "cdr3_aa": nt_to_aa(cdr3_nt),
        })

    cdr3 = table.apply(find_cdr3, axis=1)
    table["cdr3_start"] = cdr3["cdr3_start"].astype("Int64")
    table["cdr3_end"] = cdr3["cdr3_end"].astype("Int64")
    table["cdr3"] = cdr3["cdr3"]
    table["cdr3_aa"] = cdr3["cdr3_aa"]

    def find_v_cdr3_start(row):
        """Start of CDR3 within V"""
        if pd.isna(row.v_call) or pd.isna(row.cdr3_start):
            return 0

        # TODO 0-based position because that’s what 'igdiscover discover' wants
        return row.cdr3_start - row.v_sequence_start

    table["V_CDR3_start"] = table.apply(find_v_cdr3_start, axis=1)

    def make_fwr4(row):
        """
        Compute fwr4
        """
        no_result = pd.Series({
            "fwr4_start": np.nan,
            "fwr4_end": np.nan,
            "fwr4": np.nan,
            "fwr4_aa": np.nan,
        })
        j_call = row.j_call
        if pd.isna(j_call) or row.locus not in locus_to_chain:
            return no_result

        cdr3_ref_end = database.j_cdr3_end(row.j_call, locus_to_long_chain[row.locus])
        cdr3_query_end = row.cdr3_end
        if cdr3_ref_end is None or pd.isna(cdr3_query_end):
            return no_result

        fwr4_nt = row.sequence[cdr3_query_end:row.j_sequence_end]
        return pd.Series({
            "fwr4_start": cdr3_query_end + 1,
            "fwr4_end": row.j_sequence_end,
            "fwr4": fwr4_nt,
            "fwr4_aa": nt_to_aa(fwr4_nt),
        })

    # Overwrite some existing columns
    fwr4 = table.apply(make_fwr4, axis=1)
    table["fwr4_start"] = fwr4["fwr4_start"].astype("Int64")
    table["fwr4_end"] = fwr4["fwr4_end"].astype("Int64")
    table["fwr4"] = fwr4["fwr4"]
    table["fwr4_aa"] = fwr4["fwr4_aa"]

    def fr4_shm(row):
        """
        Compute FR4 mutation rate on nucleotide level
        """
        if pd.isna(row.cdr3_end) or pd.isna(row.j_sequence_end):
            return np.nan
        sequence = row.sequence[row.cdr3_end:row.j_sequence_end]
        germline = database.j[row.j_call][row.j_germline_start-1:row.j_germline_end]
        dist = edit_distance(germline, sequence)
        return 100. * dist / len(germline)

    def fr4_aa_mut(row):
        """
        Compute FR4 amino acid mutation rate
        """
        if pd.isna(row.cdr3_end) or pd.isna(row.j_sequence_end):
            return np.nan
        sequence = row.sequence[row.cdr3_end:row.j_sequence_end]
        germline = database.j[row.j_call][row.j_germline_start-1:row.j_germline_end]
        sequence_aa = nt_to_aa(sequence)
        germline_aa = nt_to_aa(germline)
        dist = edit_distance(germline_aa, sequence_aa)
        return 100. * dist / len(germline_aa)

    table["FR4_SHM"] = table.apply(fr4_shm, axis=1)
    table["J_aa_mut"] = table.apply(fr4_aa_mut, axis=1)

    for colname in (
        "FR1_SHM",
        "CDR1_SHM",
        "FR2_SHM",
        "CDR2_SHM",
        "FR3_SHM",
        "FR4_SHM",
        "FR1_aa_mut",
        "CDR1_aa_mut",
        "FR2_aa_mut",
        "CDR2_aa_mut",
        "FR3_aa_mut",
        "J_aa_mut",
        "V_aa_mut",
    ):
        assert colname in table.columns

    # The following columns existed in the non-AIRR version of IgDiscover, but
    # are omitted for now:
    #
    # - leader
    # - UTR
    # - V_end

    # TODO
    # - extend_left_ungapped for V sequences
    # - write with {:.1f} (see TableWriter.write())
    return table


def check_table(table, database):
    for row in table.itertuples(index=False):
        if pd.isna(row.v_call) or pd.isna(row.j_call):
            continue

        if pd.notna(row.complete_vdj) and pd.notna(row.j_call):
            is_complete_vdj = row.v_germline_start == 1 and row.j_germline_end == len(
                database.j[row.j_call]
            )
            assert bool(row.complete_vdj == "T") == is_complete_vdj

        # V
        if pd.notna(row.v_call):
            assert len(row.v_germline_alignment) == len(row.v_sequence_alignment)
            assert row.sequence[
                row.v_sequence_start - 1 : row.v_sequence_end
            ] == row.v_sequence_alignment.replace("-", "")
            assert "-" not in row.v_sequence_alignment_aa

            v_ref = database.v[row.v_call]
            assert v_ref[
                row.v_germline_start - 1 : row.v_germline_end
            ] == row.v_germline_alignment.replace("-", "")

        # D
        if not pd.isna(row.d_call):
            assert len(row.d_germline_alignment) == len(row.d_sequence_alignment)
            assert row.sequence[
                int(row.d_sequence_start) - 1 : int(row.d_sequence_end)
            ] == row.d_sequence_alignment.replace("-", "")
            assert "-" not in row.d_sequence_alignment_aa

        # J
        if not pd.isna(row.j_call):
            assert len(row.j_germline_alignment) == len(row.j_sequence_alignment)
            assert row.sequence[
                row.j_sequence_start - 1 : row.j_sequence_end
            ] == row.j_sequence_alignment.replace("-", "")

            j_ref = database.j[row.j_call]
            assert (
                j_ref[row.j_germline_start - 1 : row.j_germline_end]
                == row.j_germline_alignment.replace("-", "")
            )


def query_position(record, gene: str, reference_position: int):
    """
    Given a position on the reference, return the same position but relative to
    the full query sequence.
    """
    ref_pos = getattr(record, gene + "_germline_start") - 1
    query_pos = getattr(record, gene + "_sequence_start") - 1

    # Iterate over alignment columns
    if ref_pos == reference_position:
        return query_pos

    sequence_alignment = getattr(record, gene + "_sequence_alignment")
    germline_alignment = getattr(record, gene + "_germline_alignment")

    for ref_c, query_c in zip(germline_alignment, sequence_alignment):
        if ref_c != "-":
            ref_pos += 1
        if query_c != "-":
            query_pos += 1
        if ref_pos == reference_position:
            return query_pos
    return None
