"""
# Primer3InputTag Class and Methods

This module contains a class, `Primer3InputTag`, to standardize user input to Primer3. Settings that
are controlled here include general parameters that apply to a singular Primer3 session
as well as query-specific settings that can be changed in between primer queries.

"""

from enum import auto
from enum import unique

from strenum import UppercaseStrEnum


@unique
class Primer3InputTag(UppercaseStrEnum):
    """
    Enumeration of Primer3 input tags.

    Please see the Primer3 manual for additional details:
     https://primer3.org/manual.html#commandLineTags

    This class represents two categories of Primer3 input tags:
     * `SEQUENCE_` tags are those that control sequence-specific attributes of Primer3 jobs.
       These can be modified after each query submitted to Primer3.
     * `PRIMER_` tags are those that describe more general parameters of Primer3 jobs.
       These attributes persist across queries to Primer3 unless they are explicitly reset.
       Errors in these "global" input tags are fatal.
    """

    # Developer note: sequence-specific tags must be specified prior to global tags

    @property
    def is_sequence_arg(self) -> bool:
        """True if this is a sequence input tag (query-specific)"""
        return self.name.startswith("SEQUENCE")

    @property
    def is_global_arg(self) -> bool:
        """True if this is a global input tags (will persist across primer3 queries)"""
        return not self.is_sequence_arg

    # Sequence input tags; query-specific
    SEQUENCE_EXCLUDED_REGION = auto()
    SEQUENCE_INCLUDED_REGION = auto()
    SEQUENCE_PRIMER_REVCOMP = auto()
    SEQUENCE_FORCE_LEFT_END = auto()
    SEQUENCE_INTERNAL_EXCLUDED_REGION = auto()
    SEQUENCE_QUALITY = auto()
    SEQUENCE_FORCE_LEFT_START = auto()
    SEQUENCE_INTERNAL_OLIGO = auto()
    SEQUENCE_START_CODON_POSITION = auto()
    SEQUENCE_FORCE_RIGHT_END = auto()
    SEQUENCE_OVERLAP_JUNCTION_LIST = auto()
    SEQUENCE_TARGET = auto()
    SEQUENCE_FORCE_RIGHT_START = auto()
    SEQUENCE_PRIMER = auto()
    SEQUENCE_TEMPLATE = auto()
    SEQUENCE_ID = auto()
    SEQUENCE_PRIMER_PAIR_OK_REGION_LIST = auto()

    # Global input tags; will persist across primer3 queries
    PRIMER_DNA_CONC = auto()
    PRIMER_MAX_END_GC = auto()
    PRIMER_PAIR_WT_PRODUCT_SIZE_LT = auto()
    PRIMER_DNTP_CONC = auto()
    PRIMER_MAX_END_STABILITY = auto()
    PRIMER_PAIR_WT_PRODUCT_TM_GT = auto()
    PRIMER_EXPLAIN_FLAG = auto()
    PRIMER_MAX_GC = auto()
    PRIMER_PAIR_WT_PRODUCT_TM_LT = auto()
    PRIMER_FIRST_BASE_INDEX = auto()
    PRIMER_MAX_HAIRPIN_TH = auto()
    PRIMER_PAIR_WT_PR_PENALTY = auto()
    PRIMER_GC_CLAMP = auto()
    PRIMER_MAX_LIBRARY_MISPRIMING = auto()
    PRIMER_PAIR_WT_TEMPLATE_MISPRIMING = auto()
    PRIMER_INSIDE_PENALTY = auto()
    PRIMER_MAX_NS_ACCEPTED = auto()
    PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH = auto()
    PRIMER_INTERNAL_DNA_CONC = auto()
    PRIMER_MAX_POLY_X = auto()
    PRIMER_PICK_ANYWAY = auto()
    PRIMER_INTERNAL_DNTP_CONC = auto()
    PRIMER_MAX_SELF_ANY = auto()
    PRIMER_PICK_INTERNAL_OLIGO = auto()
    PRIMER_INTERNAL_MAX_GC = auto()
    PRIMER_MAX_SELF_ANY_TH = auto()
    PRIMER_PICK_LEFT_PRIMER = auto()
    PRIMER_INTERNAL_MAX_HAIRPIN_TH = auto()
    PRIMER_MAX_SELF_END = auto()
    PRIMER_PICK_RIGHT_PRIMER = auto()
    PRIMER_INTERNAL_MAX_LIBRARY_MISHYB = auto()
    PRIMER_MAX_SELF_END_TH = auto()
    PRIMER_PRODUCT_MAX_TM = auto()
    PRIMER_INTERNAL_MAX_NS_ACCEPTED = auto()
    PRIMER_MAX_SIZE = auto()
    PRIMER_PRODUCT_MIN_TM = auto()
    PRIMER_INTERNAL_MAX_POLY_X = auto()
    PRIMER_MAX_TEMPLATE_MISPRIMING = auto()
    PRIMER_PRODUCT_OPT_SIZE = auto()
    PRIMER_INTERNAL_MAX_SELF_ANY = auto()
    PRIMER_MAX_TEMPLATE_MISPRIMING_TH = auto()
    PRIMER_PRODUCT_OPT_TM = auto()
    PRIMER_INTERNAL_MAX_SELF_ANY_TH = auto()
    PRIMER_MAX_TM = auto()
    PRIMER_PRODUCT_SIZE_RANGE = auto()
    PRIMER_INTERNAL_MAX_SELF_END = auto()
    PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION = auto()
    PRIMER_QUALITY_RANGE_MAX = auto()
    PRIMER_INTERNAL_MAX_SELF_END_TH = auto()
    PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION = auto()
    PRIMER_QUALITY_RANGE_MIN = auto()
    PRIMER_INTERNAL_MAX_SIZE = auto()
    PRIMER_MIN_END_QUALITY = auto()
    PRIMER_SALT_CORRECTIONS = auto()
    PRIMER_INTERNAL_MAX_TM = auto()
    PRIMER_MIN_GC = auto()
    PRIMER_SALT_DIVALENT = auto()
    PRIMER_INTERNAL_MIN_GC = auto()
    PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE = auto()
    PRIMER_SALT_MONOVALENT = auto()
    PRIMER_INTERNAL_MIN_QUALITY = auto()
    PRIMER_MIN_QUALITY = auto()
    PRIMER_SEQUENCING_ACCURACY = auto()
    PRIMER_INTERNAL_MIN_SIZE = auto()
    PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE = auto()
    PRIMER_SEQUENCING_INTERVAL = auto()
    PRIMER_INTERNAL_MIN_TM = auto()
    PRIMER_MIN_SIZE = auto()
    PRIMER_SEQUENCING_LEAD = auto()
    PRIMER_INTERNAL_MISHYB_LIBRARY = auto()
    PRIMER_MIN_THREE_PRIME_DISTANCE = auto()
    PRIMER_SEQUENCING_SPACING = auto()
    PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME = auto()
    PRIMER_MIN_TM = auto()
    PRIMER_TASK = auto()
    PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME = auto()
    PRIMER_MISPRIMING_LIBRARY = auto()
    PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = auto()
    PRIMER_INTERNAL_OPT_GC_PERCENT = auto()
    PRIMER_MUST_MATCH_FIVE_PRIME = auto()
    PRIMER_THERMODYNAMIC_PARAMETERS_PATH = auto()
    PRIMER_INTERNAL_OPT_SIZE = auto()
    PRIMER_MUST_MATCH_THREE_PRIME = auto()
    PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT = auto()
    PRIMER_INTERNAL_OPT_TM = auto()
    PRIMER_NUM_RETURN = auto()
    PRIMER_TM_FORMULA = auto()
    PRIMER_INTERNAL_SALT_DIVALENT = auto()
    PRIMER_OPT_GC_PERCENT = auto()
    PRIMER_WT_END_QUAL = auto()
    PRIMER_INTERNAL_SALT_MONOVALENT = auto()
    PRIMER_OPT_SIZE = auto()
    PRIMER_WT_END_STABILITY = auto()
    PRIMER_INTERNAL_WT_END_QUAL = auto()
    PRIMER_OPT_TM = auto()
    PRIMER_WT_GC_PERCENT_GT = auto()
    PRIMER_INTERNAL_WT_GC_PERCENT_GT = auto()
    PRIMER_OUTSIDE_PENALTY = auto()
    PRIMER_WT_GC_PERCENT_LT = auto()
    PRIMER_INTERNAL_WT_GC_PERCENT_LT = auto()
    PRIMER_PAIR_MAX_COMPL_ANY = auto()
    PRIMER_WT_HAIRPIN_TH = auto()
    PRIMER_INTERNAL_WT_HAIRPIN_TH = auto()
    PRIMER_PAIR_MAX_COMPL_ANY_TH = auto()
    PRIMER_WT_LIBRARY_MISPRIMING = auto()
    PRIMER_INTERNAL_WT_LIBRARY_MISHYB = auto()
    PRIMER_PAIR_MAX_COMPL_END = auto()
    PRIMER_WT_NUM_NS = auto()
    PRIMER_INTERNAL_WT_NUM_NS = auto()
    PRIMER_PAIR_MAX_COMPL_END_TH = auto()
    PRIMER_WT_POS_PENALTY = auto()
    PRIMER_INTERNAL_WT_SELF_ANY = auto()
    PRIMER_PAIR_MAX_DIFF_TM = auto()
    PRIMER_WT_SELF_ANY = auto()
    PRIMER_INTERNAL_WT_SELF_ANY_TH = auto()
    PRIMER_PAIR_MAX_LIBRARY_MISPRIMING = auto()
    PRIMER_WT_SELF_ANY_TH = auto()
    PRIMER_INTERNAL_WT_SELF_END = auto()
    PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING = auto()
    PRIMER_WT_SELF_END = auto()
    PRIMER_INTERNAL_WT_SELF_END_TH = auto()
    PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH = auto()
    PRIMER_WT_SELF_END_TH = auto()
    PRIMER_INTERNAL_WT_SEQ_QUAL = auto()
    PRIMER_PAIR_WT_COMPL_ANY = auto()
    PRIMER_WT_SEQ_QUAL = auto()
    PRIMER_INTERNAL_WT_SIZE_GT = auto()
    PRIMER_PAIR_WT_COMPL_ANY_TH = auto()
    PRIMER_WT_SIZE_GT = auto()
    PRIMER_INTERNAL_WT_SIZE_LT = auto()
    PRIMER_PAIR_WT_COMPL_END = auto()
    PRIMER_WT_SIZE_LT = auto()
    PRIMER_INTERNAL_WT_TM_GT = auto()
    PRIMER_PAIR_WT_COMPL_END_TH = auto()
    PRIMER_WT_TEMPLATE_MISPRIMING = auto()
    PRIMER_INTERNAL_WT_TM_LT = auto()
    PRIMER_PAIR_WT_DIFF_TM = auto()
    PRIMER_WT_TEMPLATE_MISPRIMING_TH = auto()
    PRIMER_LIBERAL_BASE = auto()
    PRIMER_PAIR_WT_IO_PENALTY = auto()
    PRIMER_WT_TM_GT = auto()
    PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS = auto()
    PRIMER_PAIR_WT_LIBRARY_MISPRIMING = auto()
    PRIMER_WT_TM_LT = auto()
    PRIMER_LOWERCASE_MASKING = auto()
    PRIMER_PAIR_WT_PRODUCT_SIZE_GT = auto()
