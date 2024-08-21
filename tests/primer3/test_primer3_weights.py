from prymer.primer3.primer3_input_tag import Primer3InputTag
from prymer.primer3.primer3_weights import Primer3Weights


def test_primer_weights_valid() -> None:
    """Test instantiation of Primer3Weights object with valid input"""
    test_weights = Primer3Weights()
    test_dict = test_weights.to_input_tags()
    assert test_dict[Primer3InputTag.PRIMER_PAIR_WT_PRODUCT_SIZE_LT] == 1
    assert test_dict[Primer3InputTag.PRIMER_PAIR_WT_PRODUCT_SIZE_GT] == 1
    assert test_dict[Primer3InputTag.PRIMER_PAIR_WT_PRODUCT_TM_LT] == 0.0
    assert test_dict[Primer3InputTag.PRIMER_PAIR_WT_PRODUCT_TM_GT] == 0.0
    assert test_dict[Primer3InputTag.PRIMER_WT_END_STABILITY] == 0.25
    assert test_dict[Primer3InputTag.PRIMER_WT_GC_PERCENT_LT] == 0.25
    assert test_dict[Primer3InputTag.PRIMER_WT_GC_PERCENT_GT] == 0.25
    assert test_dict[Primer3InputTag.PRIMER_WT_SELF_ANY] == 0.1
    assert test_dict[Primer3InputTag.PRIMER_WT_SELF_END] == 0.1
    assert test_dict[Primer3InputTag.PRIMER_WT_SIZE_LT] == 0.5
    assert test_dict[Primer3InputTag.PRIMER_WT_SIZE_GT] == 0.1
    assert test_dict[Primer3InputTag.PRIMER_WT_TM_LT] == 1.0
    assert test_dict[Primer3InputTag.PRIMER_WT_TM_GT] == 1.0
    assert len((test_dict.values())) == 13


def test_primer_weights_to_input_tags() -> None:
    """Test results from to_input_tags() with and without default values"""
    default_map = Primer3Weights().to_input_tags()
    assert default_map[Primer3InputTag.PRIMER_PAIR_WT_PRODUCT_SIZE_LT] == 1
    customized_map = Primer3Weights(product_size_lt=5).to_input_tags()
    assert customized_map[Primer3InputTag.PRIMER_PAIR_WT_PRODUCT_SIZE_LT] == 5
