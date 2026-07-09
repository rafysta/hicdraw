# hicdraw

Drawing and handling of Hi-C contact maps and domain boundaries in R.

`hicdraw` reads, slices and plots Hi-C contact maps together with
boundary-strength and domain-annotation tracks. It supports:

- the **sparse genome-wide bin matrix** produced by
  [hic200-cpp](https://github.com/mbyamaguchi/hic200-cpp)
  (`bin1 bin2 score`, one genome-wide 0-based bin index per resolution step), and
- the **`.hic`** format, read on demand through
  [`strawr`](https://github.com/aidenlab/straw).

A configuration-file driven interface (`draw_from_config()`) generates many
figures — for example a Hi-C triangle stacked with its boundary-strength track
— in a single call.

## Installation

```r
# install.packages("remotes")
remotes::install_github("rafysta/hicdraw")
```

Runtime dependencies: `data.table`, `dplyr`, `tidyr`, `ggplot2`, `scales`,
`magrittr`, `rlang`. Optional (feature-specific): `strawr` (.hic),
`rtracklayer` (bigWig), `gggenes`/`ggrepel` (gene track), `cowplot`/`ggpubr`
(stacking), `colorRamps` (matlab palette), `readxl` (xlsx config).

## Quick start (hic200 matrix)

```r
library(hicdraw)

hic <- createHiCmap()
hic <- addHiCdata(hic, "HiC_Single-MHM_200bp.txt.gz",
                  resolution = "200bp", dataName = "m1", format = "hic200")
hic <- addBoundary(hic, "HiC_Single-MHM_BS.txt", dataName = "b1")
hic <- setRegion(hic, "II", 800000, 900000)

p <- drawStack(hic, hic_name = "m1", boundary_name = "b1",
               col = "matlab", zmax = 0.99, unit = 1e3)
print(p)
```

## Configuration-driven batch drawing

```r
draw_from_config(
  config_file  = "config/Request_microdomain.xlsx",
  sample_sheet = "sample_sheet.tsv",
  matrix_dir   = "data/200bp_matrix",
  boundary_dir = "data/200bp_boundary",
  out_dir      = "out"
)
```

Each row of the configuration table describes one figure: `type`,
`data_format`, `input_id1`/`input_id2` (referencing `sample_id` in the sample
sheet), `resolution`, `chr`, `start`, `end`, `max_distance`, colour-scale
thresholds, `with_boundary`, `scale` and `output_name`. See
`?draw_from_config` for the full column list.

## Notes on the hic200 bin index

hic200-cpp concatenates chromosomes in the order given by `chrom_sizes`, so a
bin index is *genome-wide*. For _S. pombe_ (`pombe_chrom_sizes`) at 200 bp the
offsets are chrI = 0, chrII = 27896, chrIII = 50596, e.g.
`chrII:800,000` → bin `27896 + 800000/200 = 31896`. Matrices are stored as the
upper triangle only; `slice_hic200()` symmetrises on read.

## License

MIT © 2026 Hideki Noma
