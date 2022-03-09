"""
Augment AIRR-formatted IgBLAST output with extra IgDiscover-specific columns

Also, fill in the CDR3 columns
"""
import csv
import json
import sys
import logging
import time

from tinyalign import edit_distance, hamming_distance
from xopen import xopen

from ..igblast import Database
from ..utils import nt_to_aa

EXTRA_COLUMNS = [
    "V_SHM",
    "J_SHM",
    "count",
    "barcode",
    "V_errors",
    "D_errors",
    "J_errors",
    "V_covered",
    "D_covered",
    "J_covered",
    "FR1_SHM",
    "CDR1_SHM",
    "FR2_SHM",
    "CDR2_SHM",
    "FR3_SHM",
    "FR1_aa_mut",
    "CDR1_aa_mut",
    "FR2_aa_mut",
    "CDR2_aa_mut",
    "FR3_aa_mut",
    "V_aa_mut",
    "V_CDR3_start",
    "FR4_SHM",
    "J_aa_mut",
]

COLUMN_TYPES = {
    "v_identity": float,
    "j_identity": float,
    "v_sequence_start": int,
    "v_sequence_end": int,
    "v_germline_start": int,
    "v_germline_end": int,
    "j_sequence_start": int,
    "j_sequence_end": int,
    "j_germline_start": int,
    "j_germline_end": int,
    "fwr1_start": int,
    "fwr1_end": int,
    "fwr2_start": int,
    "fwr2_end": int,
    "fwr3_start": int,
    "fwr3_end": int,
    "cdr1_start": int,
    "cdr1_end": int,
    "cdr2_start": int,
    "cdr2_end": int,
}

ALLOWED_LOCI = {
    "IGH",
    "IGK",
    "IGL",
    "TRA",
    "TRB",
    "TRG",
    "TRD",
}

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


def main(args):
    logger.info("Running IgBLAST on database sequences to find CDR/FR region locations")
    database = Database(args.database, args.sequence_type)
    start_time = time.time()
    last_status_update = 0
    detected_cdr3s = 0

    csv.register_dialect(
        "airr",
        delimiter="\t",
        lineterminator="\n",
        strict=True,
    )
    with xopen(args.table) as f:
        reader = csv.DictReader(f, dialect="airr")
        writer = csv.DictWriter(
            sys.stdout,
            fieldnames=list(reader.fieldnames) + EXTRA_COLUMNS,
            dialect="airr",
        )
        writer.writeheader()
        n = 0
        for i, record in enumerate(reader):
            n += 1  # TODO
            record = parse_record(record)
            if args.rename is not None:
                record["sequence_id"] = f"seq{i + 1}"
            record = augment_record(record, database)
            record = format_float_columns(record)
            try:
                writer.writerow(record)
            except BrokenPipeError:
                sys.exit(1)

            if record["cdr3"]:
                detected_cdr3s += 1
            if n % 10000 == 0:
                elapsed = time.time() - start_time
                if elapsed >= last_status_update + 60:
                    logger.info(
                        "Processed {:10,d} sequences at {:.2f} ms/sequence".format(
                            n, elapsed / n * 1e3
                        )
                    )
                    last_status_update = elapsed

    elapsed = time.time() - start_time
    logger.info(
        "Processed {:10,d} sequences at {:.2f} ms/sequence".format(n, elapsed / n * 1e3)
    )
    logger.info("CDR3s detected in %.1f%% of all sequences", detected_cdr3s / n * 100)
    if args.stats:
        stats = {"total": n, "detected_cdr3s": detected_cdr3s}
        with open(args.stats, "w") as f:
            json.dump(stats, f)
            print(file=f)


def parse_record(record):
    for col, typ in COLUMN_TYPES.items():
        record[col] = typ(record[col]) if record[col] else None
    return record


def format_float_columns(record):
    """
    Round values in some float columns and convert to string

    Unnecessarily precise float values just take up space
    """
    for name in (
        "V_covered",
        "D_covered",
        "J_covered",
        "FR1_SHM",
        "CDR1_SHM",
        "FR2_SHM",
        "CDR2_SHM",
        "FR3_SHM",
        "FR4_SHM",
        "V_SHM",
        "J_SHM",
        "V_aa_mut",
        "J_aa_mut",
        "FR1_aa_mut",
        "CDR1_aa_mut",
        "FR2_aa_mut",
        "CDR2_aa_mut",
        "FR3_aa_mut",
    ):
        record[name] = f"{record[name]:.1f}" if record.get(name) is not None else None
    # for name in ("v_support", "d_support", "j_support"):
    #     record[name] = f"{record[name]:G}" if record[name] is not None else None
    return record


def augment_record(record, database):
    for col in "v_call", "d_call", "j_call":
        record[col] = record[col].lstrip("%")

    record["V_SHM"] = (
        (100.0 - record["v_identity"]) if record["v_identity"] is not None else None
    )
    record["J_SHM"] = (
        (100.0 - record["j_identity"]) if record["j_identity"] is not None else None
    )

    sequence_id, size, barcode = parse_header(record["sequence_id"])
    record["sequence_id"] = sequence_id
    record["count"] = size  # TODO consensus_count/duplicate_count
    record["barcode"] = barcode

    record["V_errors"] = count_errors(
        record["v_germline_alignment"], record["v_sequence_alignment"]
    )
    record["D_errors"] = count_errors(
        record["d_germline_alignment"], record["d_sequence_alignment"]
    )
    record["J_errors"] = count_errors(
        record["j_germline_alignment"], record["j_sequence_alignment"]
    )

    def covered(alignment, db, call):
        if not alignment or not call:
            return None
        return len(alignment.replace("-", "")) / len(db[call]) * 100.0

    record["V_covered"] = covered(
        record["v_germline_alignment"], database.v, record["v_call"]
    )
    record["D_covered"] = covered(
        record["d_germline_alignment"], database.d, record["d_call"]
    )
    record["J_covered"] = covered(
        record["j_germline_alignment"], database.j, record["j_call"]
    )

    set_shm_columns(record, database)
    set_aa_mut_columns(record, database)
    set_cdr3_columns(record, database)

    if record["v_call"] and record["cdr3_start"]:
        # Start of CDR3 within V
        # TODO 0-based position because that’s what 'igdiscover discover' wants
        record["V_CDR3_start"] = record["cdr3_start"] - record["v_sequence_start"]
    else:
        record["V_CDR3_start"] = 0

    set_fwr4_columns(record, database)

    # The following columns existed in the non-AIRR version of IgDiscover, but
    # are omitted for now:
    # - leader
    # - UTR
    # - V_end

    # TODO
    # - extend_left_ungapped for V sequences
    return record


def set_shm_columns(record, database):
    """
    Compute SHM (actually mutation rate on nucleotide level) for
    all regions on V
    """
    for airr_col, region in (
        ("fwr1", "FR1"),
        ("cdr1", "CDR1"),
        ("fwr2", "FR2"),
        ("cdr2", "CDR2"),
        ("fwr3", "FR3"),
    ):
        start = record[airr_col + "_start"]
        end = record[airr_col + "_end"]
        if start is None or end is None:
            record[region + "_SHM"] = None
            continue
        sequence = record["sequence"][start - 1 : end]
        regions = database.v_regions_nt.get(record["v_call"])
        if not regions:
            continue
        germline = regions.get(region)
        dist = edit_distance(germline, sequence)
        record[region + "_SHM"] = 100.0 * dist / len(germline)


def set_aa_mut_columns(record, database):
    """
    Compute amino acid mutation rate for all regions on V and also for V
    itself as the sum of the regions (that is, excluding the CDR3)
    """
    total_length = 0
    total_dist = 0
    n_regions = 0
    for airr_col, region in (
        ("fwr1", "FR1"),
        ("cdr1", "CDR1"),
        ("fwr2", "FR2"),
        ("cdr2", "CDR2"),
        ("fwr3", "FR3"),
    ):
        record[region + "_aa_mut"] = None
        start = record[airr_col + "_start"]
        end = record[airr_col + "_end"]
        if start is None or end is None:
            continue
        sequence_aa = nt_to_aa(record["sequence"][start - 1 : end])
        regions = database.v_regions_aa.get(record["v_call"])
        if not regions:
            continue
        germline_aa = regions.get(region)
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
        record[region + "_aa_mut"] = 100.0 * mut_aa
    if n_regions == 5:
        record["V_aa_mut"] = 100.0 * total_dist / total_length
    else:
        record["V_aa_mut"] = None


def count_errors(s, t):
    if s == "" or t == "":
        return None
    if len(s) != len(t):
        raise ValueError("alignment sequences must be of the same length")
    return sum(1 for c, d in zip(s, t) if c != d)


def set_cdr3_columns(record, database):
    if (
        not record["v_call"]
        or not record["j_call"]
        or record["locus"] not in ALLOWED_LOCI
    ):
        return

    # CDR3 start
    cdr3_ref_start = database.v_cdr3_start(
        record["v_call"], record["locus"]
    )
    if cdr3_ref_start is None:
        return
    cdr3_query_start = query_position(record, "v", reference_position=cdr3_ref_start)
    if cdr3_query_start is None:
        # Alignment is not long enough to cover CDR3 start position; try to rescue it
        # by assuming that the alignment would continue without indels.
        cdr3_query_start = record["v_sequence_end"] + (
            cdr3_ref_start - record["v_germline_end"]
        )

    # CDR3 end
    cdr3_ref_end = database.j_cdr3_end(record["j_call"], record["locus"])
    if cdr3_ref_end is None:
        return

    cdr3_query_end = query_position(record, "j", reference_position=cdr3_ref_end)
    if cdr3_query_end is None:
        return
    cdr3_nt = record["sequence"][cdr3_query_start:cdr3_query_end]

    record["cdr3_start"] = cdr3_query_start + 1
    record["cdr3_end"] = cdr3_query_end
    record["cdr3"] = cdr3_nt
    record["cdr3_aa"] = nt_to_aa(cdr3_nt)


def set_fwr4_columns(record, database):
    j_call = record["j_call"]
    if not j_call or record["locus"] not in ALLOWED_LOCI:
        return

    cdr3_ref_end = database.j_cdr3_end(record["j_call"], record["locus"])
    cdr3_query_end = record["cdr3_end"]
    if cdr3_ref_end is None or not cdr3_query_end:
        return

    fwr4_nt = record["sequence"][cdr3_query_end : record["j_sequence_end"]]

    # This overwrites some existing columns
    record["fwr4_start"] = cdr3_query_end + 1
    record["fwr4_end"] = record["j_sequence_end"]
    record["fwr4"] = fwr4_nt
    record["fwr4_aa"] = nt_to_aa(fwr4_nt)

    # Compute FR4 mutation rate on nucleotide level
    germline = database.j[record["j_call"]][
        record["j_germline_start"] - 1 : record["j_germline_end"]
    ]
    dist = edit_distance(germline, fwr4_nt)
    record["FR4_SHM"] = 100.0 * dist / len(germline)

    # Compute FR4 amino acid mutation rate
    sequence_aa = record["fwr4_aa"]
    germline_aa = nt_to_aa(germline)
    dist = edit_distance(germline_aa, sequence_aa)
    record["J_aa_mut"] = 100.0 * dist / len(germline_aa)


def query_position(record, gene: str, reference_position: int):
    """
    Given a position on the reference, return the same position but relative to
    the full query sequence.
    """
    ref_pos = record[gene + "_germline_start"] - 1
    query_pos = record[gene + "_sequence_start"] - 1

    # Iterate over alignment columns
    if ref_pos == reference_position:
        return query_pos

    sequence_alignment = record[gene + "_sequence_alignment"]
    germline_alignment = record[gene + "_germline_alignment"]

    for ref_c, query_c in zip(germline_alignment, sequence_alignment):
        if ref_c != "-":
            ref_pos += 1
        if query_c != "-":
            query_pos += 1
        if ref_pos == reference_position:
            return query_pos
    return None


def parse_header(header):
    """
    Extract size= and barcode= fields from the FASTA/FASTQ header line

    >>> parse_header("name;size=12;barcode=ACG;")
    ('name', 12, 'ACG')
    >>> parse_header("another name;size=200;foo=bar;")
    ('another name', 200, None)
    """
    fields = header.split(';')
    query_name = fields[0]
    size = barcode = None
    for field in fields[1:]:
        if field == '':
            continue
        if '=' in field:
            key, value = field.split('=', maxsplit=1)
            if key == 'size':
                size = int(value)
            elif key == 'barcode':
                barcode = value
    return query_name, size, barcode
