import pytest

@pytest.mark.parametrize("test_num", [1, 2, 3, 4, 5, 6,7, 8, 9, 10, 11, 12, 13, 14, 15])
def test_cli_simulation_regressions(run_simulation_testset, expected_output_data, test_num):
    _, observed_final_lines = run_simulation_testset(test_num)

    for source_filename, expected_by_test in expected_output_data.items():
        if test_num not in expected_by_test:
            continue

        assert source_filename in observed_final_lines, (
            f"Expected {source_filename} for test_{test_num}, but file was not produced"
        )
        assert observed_final_lines[source_filename] == expected_by_test[test_num]

