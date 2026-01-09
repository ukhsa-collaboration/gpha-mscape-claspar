from claspar import handle_tables


def test_create_bacterial_analysis_fields():
    handle_tables.create_bacterial_analysis_fields(
        domain="bacteria",
        classifier="sylph",
        record_id="ID_123456",
        thresholds={'stuff': 1},
        headline_result='Found some stuff here.',
        results={}

