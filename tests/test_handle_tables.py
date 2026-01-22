from datetime import datetime

from claspar import handle_tables


def test_create_bacterial_analysis_fields():
    today = datetime.today().strftime("%Y-%m-%d")
    expected_analysis_table_dict = {
        "identifiers": [],
        "name": "bacteria-classifier-parser",
        "description": "This is an analysis to parse and filter the bacteria classifications from sylph",
        "analysis_date": today,
        "pipeline_name": "ClasPar",
        "pipeline_version": "1.0.0",
        "pipeline_url": "https://github.com/ukhsa-collaboration/gpha-mscape-claspar",
        "methods": '{"stuff": 1}',
        "result": "Found some stuff here.",
        "result_metrics": '{"0": {"thing": 10, "type": "little"}, "1": {"thing": 10, "type": "big"}}',
        "mscape_records": ["ID_123456"],
    }
    onyx_analysis_table, exitcode = handle_tables.create_analysis_fields(
        domain="bacteria",
        classifier="sylph",
        record_id="ID_123456",
        thresholds={"stuff": 1},
        headline_result="Found some stuff here.",
        results={0: {"thing": 10, "type": "little"}, 1: {"thing": 10, "type": "big"}},
        server="mscape",
    )
    assert onyx_analysis_table.__dict__ == expected_analysis_table_dict

    print(dir(onyx_analysis_table))
