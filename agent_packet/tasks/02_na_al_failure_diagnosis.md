# Task 02: Public Failure Diagnosis

Goal: compare the public Kohn-Sham baseline against public ARPES references for
Na and Al.

Use:

- `data/public/Na/arpes_reference.csv`
- `data/public/Al/arpes_reference.csv`
- your `ks_bands.csv` from Task 01

For each ARPES point, compare it with the nearest occupied Kohn-Sham band at
the same `point_id`.

Report:

- RMSE for Na and Al;
- mean signed error (`prediction - reference`);
- deepest occupied public-reference point and the nearest Kohn-Sham prediction;
- a short diagnosis of why a static Kohn-Sham band is not sufficient.

Expected qualitative pattern:

- Na should show substantial overbinding / excessive occupied bandwidth.
- Al should be much closer, with a small correction.

Do not tune the Kohn-Sham bands to match ARPES in this task.

