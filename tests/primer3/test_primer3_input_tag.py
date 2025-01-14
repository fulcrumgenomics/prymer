from prymer.primer3 import Primer3InputTag


def test_primer3_input_tag() -> None:
    for tag in Primer3InputTag:
        assert tag.is_sequence_arg == tag.startswith("S")
        assert tag.is_global_arg == tag.startswith("P")
        assert tag.is_global_arg != tag.is_sequence_arg
        assert isinstance(Primer3InputTag[tag], Primer3InputTag)
        # This is supported in python 3.12+, but not in 3.11
        # assert f"{tag}" in Primer3InputTag
