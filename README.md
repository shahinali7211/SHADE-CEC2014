# SHADE-CEC2014
SHADE (Success-History Adaptive DE) for CEC-2014 (F1–F30). C++ with Boost RNG. Current-to-pbest/1, external archive, memory-based F/CR adaptation. 51 runs, 30D, 300K FEs. Auto CSV summary (mean ± std). Research-ready, clean code. Ideal for students &amp; EA researchers. MIT license

# SHADE-CEC2014

**Success-History based Adaptive Differential Evolution (SHADE)** applied to the **CEC-2014 real-parameter single-objective benchmark** (30 functions).

A clean, efficient, and well-documented C++ implementation with:
- `current-to-pbest/1` mutation
- External archive
- Historical memory for `F` and `CR`
- Lehmer mean update for `F`, weighted mean for `CR`
- Boost RNG (Mersenne Twister + normal/cauchy)
- Automatic CSV summary: `results_50_dim.csv`

---

## Algorithm Parameters
| Parameter | Value |
|---------|-------|
| Population size | 100 |
| Dimension (`dim`) | 30 |
| Max FEs | 300,000 (`10000 * dim`) |
| Runs per function | 51 |
| `p` (p-best ratio) | 0.1 |
| Memory size `H` | `dim` |
| Archive size | `2 × Pop` |
| Initial `F`, `CR` | 0.5 |

---

> **Note**: `input_data/` folder is **not included** due to GitHub size limits (>2 GB).  
> Download from: https://github.com/P-N-Suganthan/CEC2014  
> Place inside `CEC2010/` before running.

---

## Compile & Run

### Windows (MinGW / MSVC)
```bash
g++ -O3 -std=c++11 main.cpp -o SHADE_CEC2014 -lboost_random-mt
SHADE_CEC2014.exe
