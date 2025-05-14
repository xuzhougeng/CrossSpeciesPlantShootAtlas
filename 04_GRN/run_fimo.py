#!/usr/bin/env python

r"""
Use fimo to scan genome for motif hits
"""

import os
import re
import tempfile
from argparse import ArgumentParser
from typing import List

import joblib
import numpy as np
import pandas as pd
from Bio.SeqIO import SeqRecord, parse, write
import subprocess
from typing import Any, List, Mapping, Optional

TMP_PREFIX = "GLUETMP"

def run_command(
        command: str,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        log_command: bool = True, print_output: bool = True,
        err_message: Optional[Mapping[int, str]] = None, **kwargs
) -> Optional[List[str]]:
    r"""
    Run an external command and get realtime output

    Parameters
    ----------
    command
        A string containing the command to be executed
    stdout
        Where to redirect stdout
    stderr
        Where to redirect stderr
    echo_command
        Whether to log the command being printed (log level is INFO)
    print_output
        Whether to print stdout of the command.
        If ``stdout`` is PIPE and ``print_output`` is set to False,
        the output will be returned as a list of output lines.
    err_message
        Look up dict of error message (indexed by error code)
    **kwargs
        Other keyword arguments to be passed to :class:`subprocess.Popen`

    Returns
    -------
    output_lines
        A list of output lines (only returned if ``stdout`` is PIPE
        and ``print_output`` is False)
    """
    if log_command:
        run_command.logger.info("Executing external command: %s", command)
    executable = command.split(" ")[0]
    with subprocess.Popen(command, stdout=stdout, stderr=stderr,
                          shell=True, **kwargs) as p:
        if stdout == subprocess.PIPE:
            prompt = f"{executable} ({p.pid}): "
            output_lines = []

            def _handle(line):
                line = line.strip().decode()
                if print_output:
                    print(prompt + line)
                else:
                    output_lines.append(line)

            while True:
                _handle(p.stdout.readline())
                ret = p.poll()
                if ret is not None:
                    # Handle output between last readlines and successful poll
                    for line in p.stdout.readlines():
                        _handle(line)
                    break
        else:
            output_lines = None
            ret = p.wait()
    if ret != 0:
        err_message = err_message or {}
        if ret in err_message:
            err_message = " " + err_message[ret]
        elif "__default__" in err_message:
            err_message = " " + err_message["__default__"]
        else:
            err_message = ""
        raise RuntimeError(
            f"{executable} exited with error code: {ret}.{err_message}")
    if stdout == subprocess.PIPE and not print_output:
        return output_lines


def fimo(
        seqs: List[SeqRecord], motif_file: os.PathLike,
        fimo_args: str = None, n_jobs: int = 1, random_state: int = 0
) -> pd.DataFrame:
    r"""
    Scan for motif in the sequences using ``fimo``.

    Parameters
    ----------
    seqs
        A list of sequences to be scanned
    motif_file
        Path to the motif file (must be in meme format)
    fimo_args
        Additional command line arguments passed to ``fimo``
    n_jobs
        Split peaks into ``n_jobs`` portions and run ``fimo`` in parallel
    random_state
        Random seed that determines peak shuffling before splitting

    Returns
    -------
    fimo_result
        ``fimo`` scanning result stored in a data frame.
        If not hit is found, the return value is ``None``.
    """
    seq_array = np.empty(len(seqs), dtype=object)
    seq_array[:] = seqs  # Otherwise a 2D char array is produced
    batch_size = round(seq_array.size / n_jobs)
    shuffle_idx = np.random.RandomState(random_state).permutation(
        seq_array.size
    )  # Shuffle peaks to prevent bias over genomic region
    tmp_fasta_list, tmp_fimo_list = [], []
    for i in range(n_jobs):
        tmp_fasta_list.append(tempfile.mkstemp(
            prefix=TMP_PREFIX, suffix=".fasta"
        ))
        tmp_fimo_list.append(tempfile.mkstemp(
            prefix=TMP_PREFIX, suffix=".fimo"
        ))
        fimo.logger.info("Writing temporary file: %s", tmp_fasta_list[-1][1])
        write(seq_array[shuffle_idx[
            (i * batch_size):min((i + 1) * batch_size, shuffle_idx.size)
        ]], tmp_fasta_list[-1][1], "fasta")

    def _fimo_job(motif: os.PathLike, fasta: os.PathLike, fimo_args: str, output: os.PathLike) -> None:
        log = fasta.replace(".fasta", ".log")
        run_command(
            f"fimo --text --skip-matched-sequence {fimo_args} "
            f"'{motif}' '{fasta}' > '{output}' 2> '{log}'",
            err_message={"__default__": f"Check {log} for details."}
        )
        os.remove(log)

    fimo_args = "" if fimo_args is None else fimo_args
    joblib.Parallel(
        n_jobs=n_jobs, backend="threading"
    )(joblib.delayed(_fimo_job)(
        motif_file, tmp_fasta, fimo_args, tmp_fimo
    ) for (_, tmp_fasta), (_, tmp_fimo) in zip(
        tmp_fasta_list, tmp_fimo_list
    ))

    fimo.logger.info("Parsing fimo results...")
    fimo_result = []
    for _, tmp_fimo in tmp_fimo_list:
        try:
            fimo_result.append(pd.read_table(tmp_fimo))
        except pd.errors.EmptyDataError:
            pass  # No significant match in this split
    if fimo_result:
        fimo_result = pd.concat(fimo_result, axis=0)
        fimo_result.columns = np.vectorize(
            lambda x: re.sub(r"^# ", "", x)
        )(fimo_result.columns)  # index are immutable, cannot replace
        fimo_result.index = pd.RangeIndex(stop=len(fimo_result))
    else:
        fimo.logger.warning("No motif match found!")
        fimo_result = None

    fimo.logger.info("Removing temporary files...")
    for tmp_fimo_fd, tmp_fimo in tmp_fimo_list:
        os.close(tmp_fimo_fd)
        os.remove(tmp_fimo)
    for tmp_fasta_fd, tmp_fasta in tmp_fasta_list:
        os.close(tmp_fasta_fd)
        os.remove(tmp_fasta)
    return fimo_result


def parse_args():
    r"""
    Parse command line arguments
    """
    parser = ArgumentParser()
    parser.add_argument("-f", "--fasta", dest="fasta", type=str, required=True)
    parser.add_argument("-m", "--motif", dest="motif", type=str, required=True)
    parser.add_argument("-o", "--output", dest="output", type=str, required=True)
    parser.add_argument("-a", "--fimo-args", dest="fimo_args", type=str, default=None)
    parser.add_argument("-j", "--n-jobs", dest="n_jobs", type=int, default=1)
    return parser.parse_args()


def main(args):
    r"""
    Main function
    """
    seqs = list(parse(args.fasta, "fasta"))
    fimo_result = fimo(seqs, args.motif, args.fimo_args, args.n_jobs)
    fimo_result.to_csv(args.output, index=False)


if __name__ == "__main__":
    main(parse_args())