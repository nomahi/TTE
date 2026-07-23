# Contributing to TTE

Contributions are welcome.

Before submitting a pull request:

1. Open an issue for substantial changes to the public interface or statistical
   behavior.
2. Add tests for new behavior and regression tests for bug fixes.
3. Update the relevant help file and public tutorial.
4. Run:

```r
devtools::test()
devtools::check()
```

Bug reports should include a minimal reproducible example and `sessionInfo()`.
Use synthetic or otherwise non-confidential data. Never upload identifiable
patient-level data.
