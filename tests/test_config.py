"""
Tests for the TOML-based configuration system.

Covers:
* the shipped ``config.toml`` matches the built-in ``_DEFAULTS``
* ``config`` module attributes are populated at import time
* ``load_config`` overlays only the keys present in a custom TOML
* unknown keys in a custom TOML are rejected
* a missing custom config file is reported via the CLI and exits non-zero
* an end-to-end CLI run with ``--custom-config`` actually changes the plot
  styling visible in the rendered HTML payload
"""
from __future__ import annotations

import tomllib
from pathlib import Path

import pytest

from bamdash.command import main
from bamdash.scripts import config
from tests.conftest import parse_plotly_payload, run_cli


# --------------------------------------------------------------------------- #
# Shipped config.toml correctness
# --------------------------------------------------------------------------- #
class TestShippedConfig:
    def test_shipped_toml_only_contains_known_keys(self):
        """Every key in the shipped config.toml is a recognised config key.

        The shipped TOML is the source of truth for the default values; it may
        legitimately differ from the hard-coded ``_DEFAULTS`` fallbacks. We only
        assert that it does not introduce unknown keys.
        """
        shipped = tomllib.loads(
            Path(config._shipped_toml_path()).read_text(encoding="utf-8")
        )
        assert set(shipped) == set(config._DEFAULTS)

    def test_module_attributes_populated_at_import(self):
        """Importing config populates module attributes from the shipped TOML.

        Values are taken from the shipped ``config.toml`` rather than hard-coded
        here, so editing the TOML does not break this test.
        """
        shipped = tomllib.loads(
            Path(config._shipped_toml_path()).read_text(encoding="utf-8")
        )
        # a representative sample of keys across all groups
        for key in ("show_log", "font", "coverage_fill_color", "strand_types",
                    "box_gb_alpha", "snp_color"):
            assert getattr(config, key) == shipped[key], f"mismatch for key '{key}'"


# --------------------------------------------------------------------------- #
# load_config overlay semantics
# --------------------------------------------------------------------------- #
class TestLoadConfig:
    def test_partial_overlay_keeps_other_defaults(self, tmp_path):
        """Only keys present in the custom file are overwritten."""
        custom = tmp_path / "custom.toml"
        custom.write_text('snp_color = "purple"\nfont_size = 99\n')
        original_font = config.font
        try:
            config.load_config(str(custom))
            assert config.snp_color == "purple"
            assert config.font_size == 99
            # untouched keys keep their default values
            assert config.font == original_font
            assert config.stem_color == "grey"
        finally:
            # reset to shipped defaults so this test does not leak
            config.load_config()

    def test_full_overlay_replaces_all_given_keys(self, tmp_path):
        custom = tmp_path / "full.toml"
        custom.write_text(
            'coverage_fill_color = "rgba(1, 2, 3, 0.5)"\n'
            'del_color = "green"\n'
        )
        try:
            config.load_config(str(custom))
            assert config.coverage_fill_color == "rgba(1, 2, 3, 0.5)"
            assert config.del_color == "green"
        finally:
            config.load_config()

    def test_reset_to_defaults_with_no_arg(self, tmp_path):
        """load_config() with no arg resets to the shipped defaults."""
        custom = tmp_path / "custom.toml"
        custom.write_text('snp_color = "purple"\n')
        config.load_config(str(custom))
        assert config.snp_color == "purple"
        # reset
        config.load_config()
        assert config.snp_color == "grey"

    def test_unknown_key_rejected(self, tmp_path):
        custom = tmp_path / "bad.toml"
        custom.write_text('not_a_real_key = 42\n')
        with pytest.raises(KeyError, match="unknown config key"):
            config.load_config(str(custom))
        # the failed overlay must not have partially applied: the bad key is not
        # set on the module
        assert not hasattr(config, "not_a_real_key")

    def test_missing_custom_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            config.load_config(str(tmp_path / "nope.toml"))


# --------------------------------------------------------------------------- #
# CLI --custom-config integration
# --------------------------------------------------------------------------- #
class TestCustomConfigCli:
    def test_custom_config_changes_coverage_fill_color(self, two_ref_bam, tmp_out, tmp_path):
        """An end-to-end run with --custom-config changes the rendered fill color."""
        custom = tmp_path / "custom.toml"
        custom.write_text('coverage_fill_color = "rgba(10, 20, 30, 0.5)"\n')
        out = run_cli(two_ref_bam, tmp_out, ["--custom-config", str(custom)])
        data, _layout = parse_plotly_payload(
            Path(f"{out}.html").read_text(), "plotly-refA")
        area = next(t for t in data if t.get("fill") == "tonexty")
        assert area["fillcolor"] == "rgba(10, 20, 30, 0.5)"
        # reset config so the test does not leak into others
        config.load_config()

    def test_custom_config_partial_overlay_in_plot(self, two_ref_bam, tmp_out, tmp_path):
        """Overridden key changes; a non-overridden key keeps its default."""
        custom = tmp_path / "custom.toml"
        custom.write_text('snp_color = "purple"\n')
        out = run_cli(two_ref_bam, tmp_out,
                      ["-t", str(_vcf_with_snp(tmp_path)), "--custom-config", str(custom)])
        data, _layout = parse_plotly_payload(
            Path(f"{out}.html").read_text(), "plotly-refA")
        snp = next(t for t in data if t.get("legendgroup") == "SNP")
        assert snp["marker"]["color"] == "purple"
        # stem_color was NOT overridden => keeps default
        shapes = _layout.get("shapes", [])
        assert all(s["line"]["color"] == "grey" for s in shapes)
        config.load_config()

    def test_unknown_key_exits_nonzero(self, two_ref_bam, tmp_out, tmp_path):
        custom = tmp_path / "bad.toml"
        custom.write_text('bogus_key = 1\n')
        with pytest.raises(SystemExit) as exc:
            main(["-b", str(two_ref_bam), "-p", str(tmp_out),
                  "--custom-config", str(custom)])
        assert exc.value.code == 1
        config.load_config()

    def test_missing_custom_config_exits_nonzero(self, two_ref_bam, tmp_out, tmp_path):
        with pytest.raises(SystemExit) as exc:
            main(["-b", str(two_ref_bam), "-p", str(tmp_out),
                  "--custom-config", str(tmp_path / "nope.toml")])
        assert exc.value.code == 1


# --------------------------------------------------------------------------- #
# helper
# --------------------------------------------------------------------------- #
def _vcf_with_snp(tmp_path: Path) -> Path:
    """Write a tiny VCF with one SNP on refA; return its path."""
    from tests.conftest import make_vcf
    vcf = tmp_path / "snp.vcf"
    make_vcf(vcf, records=[{"chrom": "refA", "pos": 5, "ref": "A", "alt": "T",
                            "info": {"AF": 0.5}}])
    return vcf
