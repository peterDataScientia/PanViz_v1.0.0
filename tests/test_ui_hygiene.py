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


def test_editor_is_the_single_interaction_table_surface():
    app = Path("app.py").read_text(encoding="utf-8")
    editor = Path("editor.html").read_text(encoding="utf-8")

    assert "Figure interaction records" not in app
    assert "All PLIP detections" not in app
    assert "st.dataframe(" not in app
    assert "Publication PNG" not in app
    assert "Publication SVG" not in app

    assert 'id="pv-interaction-records"' in editor
    assert 'id="pv-interaction-table-body"' in editor
    assert "state.scientificData?.interactions" in editor
    assert "state.customInteractions" in editor
    assert "Manually added" in editor
    assert "Detected" in editor


def test_live_table_is_refreshed_with_editor_render_cycle():
    editor = Path("editor.html").read_text(encoding="utf-8")
    assert "function renderInteractionTable()" in editor
    assert "legendLayer();renderInteractionTable();" in editor


def test_editor_remains_scrollable_below_figure():
    app = Path("app.py").read_text(encoding="utf-8")
    editor = Path("editor.html").read_text(encoding="utf-8")

    assert "scrolling=True" in app
    assert "MutationObserver" in app
    assert "function requestHostResize()" in editor
    assert "requestHostResize();" in editor
    assert "streamlit:setFrameHeight" in editor


def test_live_interaction_table_has_vertical_scroll_fallback():
    editor = Path("editor.html").read_text(encoding="utf-8")
    assert "#pv-interaction-table-wrap{overflow:auto;max-height:min(420px,45vh)}" in editor
    assert "#pv-interaction-table th{position:sticky;top:0;" in editor


def test_interaction_table_is_collapsed_by_default():
    editor = Path("editor.html").read_text(encoding="utf-8")
    assert '#pv-interaction-table-wrap{display:none;' in editor
    assert '#pv-interaction-records.open #pv-interaction-table-wrap{display:block}' in editor
    assert 'id="pv-interaction-toggle"' in editor
    assert 'aria-expanded="false"' in editor
    assert 'Show interaction records' in editor
    assert "Hide interaction records" in editor
