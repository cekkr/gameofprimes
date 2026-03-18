# Game of Primes

Prime-factor chemistry simulation with two visual modes:

- `territory` (default): 2D chemistry-field view with resource map, extraction heat, movement vectors, and energy gradient.
- `3d`: original OpenGL molecule view.

## Features

- Prime-number based molecules with interaction and reaction rules.
- Index-based molecular mass: `mass = Σ(prime_index * exponent)`.
- Prime chemistry force field with:
  same-prime repulsion, even-index attraction, odd-index repulsion, mixed-parity orbital coupling, and recursive fractal polarization.
- Persistent inertial drive in a fluid-like medium (drag + jitter), so molecules keep moving even without immediate neighbors.
- Population dynamics with fusion/fission/decomposition.
- Territory resource extraction and regeneration.
- Dynamic growth zones that move, relocate, and change count over time.
- Per-agent extracted-energy tracking with real-time color gradient.
- Tunable growth profiles: `slow`, `normal`, `fast`.

## Requirements

- Python 3.9+
- Install dependencies:

```bash
pip install -r requirements.txt
```

## Run

Default run (`territory` view, `normal` growth):

```bash
python3 main.py
```

Select growth profile:

```bash
python3 main.py --growth slow
python3 main.py --growth normal
python3 main.py --growth fast
```

Use legacy 3D view:

```bash
python3 main.py --view 3d
```

Adjust FPS:

```bash
python3 main.py --fps 60
```

Run denser or wider worlds:

```bash
python3 main.py --size 120 --molecules 200
python3 main.py --size 180 --molecules 500
```

## Controls

### Territory View (default)

- `SPACE`: pause/resume
- `RIGHT`: single step (when paused)
- `A`: toggle auto-step
- `V`: toggle velocity vectors
- `I`: toggle info panel
- `1`: slow growth profile
- `2`: normal growth profile
- `3`: fast growth profile

### 3D View

- `SPACE`: pause/resume
- `RIGHT`: single step (when paused)
- `A`: toggle auto-step
- `V`: toggle velocity vectors
- `I`: toggle info overlay
- `R`: reset camera
- Mouse wheel: zoom
- Left drag: rotate
- Right drag: pan
