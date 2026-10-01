from pathlib import Path


def test_public_web_ui_has_no_developer_explanatory_captions():
    source = Path("app.py").read_text(encoding="utf-8")

    forbidden = [
        "One PLIP analysis is reused for the original figures and the interactive editor",
        "the approved publication renderer is preserved; imported molecular structure",
        "The main table above is the authoritative record of interactions represented in the figure",
        "protected publication renderer + one reusable PLIP analysis",
    ]
    for text in forbidden:
        assert text not in source


def test_editor_does_not_expose_serialized_annotations_to_system_clipboard():
    source = Path("editor.html").read_text(encoding="utf-8")

    assert "navigator.clipboard?.writeText" not in source
    assert "navigator.clipboard.writeText" not in source
    assert "navigator.clipboard?.readText" not in source
    assert "navigator.clipboard.readText" not in source


def test_editor_status_starts_with_simple_user_facing_copy():
    source = Path("editor.html").read_text(encoding="utf-8")
    assert '<span id="pv-status">Ready</span>' in source
