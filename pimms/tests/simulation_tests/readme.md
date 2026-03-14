# Adding new simulation tests:

The guide below explains the steps needed to add a new simulation test!



1. FIRST: RUN all the tests! This is SUPER important! In adding a new simulation test we will re-generate all the expected outcomes from the extant simulation tests, so if the code is broken we will break the test suite! So, MAKE SURE you only add a new test where the code is in a *known* non-broken state!
   ```bash
   # run in the simulations directory
   pytest -v
   ```

   

2. Add a directory called `test_n ` where `n` is the next number

3. Set up the keyfile and parameter files you want (and other input files) you want for for the simulation, run it, and make sure the simulation is behaving as you expect. This will be used to create the "expected output" data for that simulation. 

4. Add a `.gitignore` file into the new directory (or copy from any of the other `test_n` directories) with the following files listed:

   ```
   *dat
   *pdb
   restart.pimms
   *xtc
   parameters_used.prm
   absolute_energies_of_angles.txt
   log.txt
   ```

   This avoids git tracking the simulation output files! 

4. ENSURING you have previously run all the tests, including the new test, we now generate the expected output. This will regenerate the files in `expected_output/` by cycling through each `test_` directory and parse in the output which gets used as the expected tests results. 

5. NOW CHECK EVERYTHING LOOKS GOOD! Use `git status` to check all the files in`expected_output/` and ensure that the ONLY change is the new test! If rebuilding these files changes anything other than new lines THIS IS PROBABLY A PROBLEM!

6. Next, update `test_simulation_regression.py` file and explicitly add the next test to the lists of tests:
   ``` python
   @pytest.mark.parametrize("test_num", [1, 2, 3, 4, 5, 6,7, 8, 9, 10, 11, 12, 13, 14, 15])
   ```

7. Finally, we run the tests, and hope everything passes!
   ```bash
   pytest -v
   ```

   