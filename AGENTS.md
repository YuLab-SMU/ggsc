# AGENTS.md

## Known false-positive `R CMD check` NOTE — do not "fix" it in the package

**Symptom** (with `devtools::check()` / `R CMD check --as-cran`):

```
* checking compilation flags used ... NOTE
Compilation used the following non-portable flag(s):
  '-mno-omit-leaf-frame-pointer'
```

Same class of NOTE can show other flags (`-fno-omit-frame-pointer`, `-march=...`) depending on the distro.

**Cause is the local R install, not this package.**

- The flag lives in the system R's `/etc/R/Makeconf` (Ubuntu/Debian `r-base` hardening flags) and is prepended to
  `CFLAGS`/`CXXFLAGS`/`FCFLAGS`. `src/Makevars` and `src/Makevars.win` in this repo set no such flag.
- `R CMD check` reports every token matching `^[-]m` (and every `-W...`) that is not in
  `_R_CHECK_COMPILATION_FLAGS_KNOWN_`. Ubuntu's `/usr/lib/R/etc/Renviron.site` whitelists only
  `-Wformat -Werror=format-security -Wdate-time`, so newer hardening flags always trip this NOTE.
- It is checked only under `--as-cran`; a plain `R CMD check` skips the step.
- Verified by reproducing the identical NOTE with an empty throwaway package containing no `Makevars`.

**Not fixable in the package:** `src/Makevars` is read *before* `/etc/R/Makeconf`, so `CXXFLAGS` is still empty
there and `$(filter-out ...)` is a no-op. Overriding the `%.o` recipe or shipping a full `Makefile` to drop a
harmless system flag is not worth it.

**Not a Bioconductor problem:** Bioc builders compile R from source, so nebbiolo1/2 (Ubuntu 24.04) compile with
`-g -O2 -Wall -Werror=format-security`; `ggsc` checks `OK` there. BiocCheck and submission are unaffected.

**Optional, local only (never commit):** to silence it, add to `~/.Renviron` — all four flags are needed because
the user file replaces the site value entirely. A plain `export` does not work, since `Renviron.site` overrides it.

```
_R_CHECK_COMPILATION_FLAGS_KNOWN_='-Wformat -Werror=format-security -Wdate-time -mno-omit-leaf-frame-pointer'
```
