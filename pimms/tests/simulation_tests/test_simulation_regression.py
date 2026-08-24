import pytest

@pytest.mark.parametrize("test_num", [1, 2, 3, 4, 5, 6,7, 8, 9, 10, 11, 12, 13, 14, 15])
def test_cli_simulation_regressions(run_simulation_testset, expected_output_data,
                                    expected_line_counts, test_num):
    _, observed_final_lines, observed_line_counts = run_simulation_testset(test_num)

    for source_filename, expected_by_test in expected_output_data.items():
        if test_num not in expected_by_test:
            continue

        assert source_filename in observed_final_lines, (
            f"Expected {source_filename} for test_{test_num}, but file was not produced"
        )
        assert observed_final_lines[source_filename] == expected_by_test[test_num]

        # write-cadence check: the final line alone is blind to dropped or
        # duplicated intermediate writes, so also pin the non-empty line count
        counts_by_test = expected_line_counts.get(source_filename, {})
        if test_num in counts_by_test:
            assert observed_line_counts[source_filename] == counts_by_test[test_num], (
                f"{source_filename} for test_{test_num}: wrote "
                f"{observed_line_counts[source_filename]} lines, expected "
                f"{counts_by_test[test_num]} (write frequency changed?)"
            )
