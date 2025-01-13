from prymer.primer3 import Primer3InputTag


def test_primer3_input_tag() -> None:
    # 17
    for tag in Primer3InputTag:
        assert tag.is_sequence_arg == tag.startswith("S")
        assert tag.is_global_arg == tag.startswith("P")
        assert tag.is_global_arg != tag.is_sequence_arg
