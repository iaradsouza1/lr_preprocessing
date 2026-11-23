"""IsoQuant (3.13) entry point for per-sample re-quantification.

Before writing counts, IsoQuant zeroes every feature without a "confirming" read (a unique read
that is spliced or matches a single-exon transcript). That check runs over ALL read groups of a run,
so per-sample runs would zero transcripts confirmed only in other samples. Here counts are written
without that step and each counter's confirmed transcripts go to $IQ_CONFIRMED_DIR;
merge_requant.py zeroes transcripts confirmed in no sample, which reproduces a joint run.
Counting runs in forked worker processes (one per chromosome), which inherit the patch.
Usage: same arguments as `isoquant`.
"""
import os
import sys

from isoquant_lib import long_read_counter as lrc

CONFIRMED_DIR = os.environ["IQ_CONFIRMED_DIR"]


def dump_unfiltered(self):
    all_features = sorted(f for f in self.all_features if f is not None)
    # the extractor is passed as a class (create_transcript_counter), not an instance
    if self.assignment_extractor is lrc.TranscriptAssignmentExtractor:
        name = os.path.join(CONFIRMED_DIR, f"confirmed.{os.getpid()}.{id(self)}.txt")
        with open(name, "w") as fh:
            fh.writelines(f"{f}\n" for f in sorted(self.confirmed_features) if f is not None)
    if self.ignore_read_groups:
        self.dump_ungrouped(all_features)
    else:
        self.dump_grouped(all_features)
    self.dump_usable()


lrc.AssignedFeatureCounter.dump = dump_unfiltered

from isoquant import main_entry  # noqa: E402

if __name__ == "__main__":
    sys.exit(main_entry())
