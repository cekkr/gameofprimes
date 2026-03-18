import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
import numpy as np
import random
import math
import argparse
import os
from collections import defaultdict, deque
from dataclasses import dataclass
from typing import Callable, List, Dict, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor


class UIFont:
    """Unified wrapper over pygame.freetype / pygame.font render APIs."""

    def __init__(self, backend: str, font_object, size: int):
        self.backend = backend
        self.font = font_object
        self.size = size

    def render(self, text: str, antialias: bool, color: Tuple[int, int, int]) -> pygame.Surface:
        if self.backend == 'pygame_freetype':
            rendered_surface, _ = self.font.render(text, fgcolor=color, size=self.size)
            return rendered_surface
        if self.backend == 'bitmap':
            return self.font.render(text, color)
        return self.font.render(text, antialias, color)

    def line_height(self) -> int:
        if self.backend == 'pygame_freetype':
            return max(1, int(self.font.get_sized_height()))
        if self.backend == 'bitmap':
            return self.font.line_height()
        if hasattr(self.font, 'get_linesize'):
            return int(self.font.get_linesize())
        return max(1, int(self.size))


class BitmapFont:
    """Minimal built-in bitmap font fallback for environments without pygame font modules."""

    GLYPHS = {
        ' ': ("000", "000", "000", "000", "000"),
        '?': ("111", "001", "010", "000", "010"),
        '0': ("111", "101", "101", "101", "111"),
        '1': ("010", "110", "010", "010", "111"),
        '2': ("111", "001", "111", "100", "111"),
        '3': ("111", "001", "111", "001", "111"),
        '4': ("101", "101", "111", "001", "001"),
        '5': ("111", "100", "111", "001", "111"),
        '6': ("111", "100", "111", "101", "111"),
        '7': ("111", "001", "010", "100", "100"),
        '8': ("111", "101", "111", "101", "111"),
        '9': ("111", "101", "111", "001", "111"),
        'A': ("010", "101", "111", "101", "101"),
        'B': ("110", "101", "110", "101", "110"),
        'C': ("011", "100", "100", "100", "011"),
        'D': ("110", "101", "101", "101", "110"),
        'E': ("111", "100", "110", "100", "111"),
        'F': ("111", "100", "110", "100", "100"),
        'G': ("011", "100", "101", "101", "011"),
        'H': ("101", "101", "111", "101", "101"),
        'I': ("111", "010", "010", "010", "111"),
        'J': ("111", "001", "001", "101", "111"),
        'K': ("101", "101", "110", "101", "101"),
        'L': ("100", "100", "100", "100", "111"),
        'M': ("101", "111", "101", "101", "101"),
        'N': ("101", "111", "111", "111", "101"),
        'O': ("111", "101", "101", "101", "111"),
        'P': ("111", "101", "111", "100", "100"),
        'Q': ("111", "101", "101", "111", "001"),
        'R': ("111", "101", "111", "110", "101"),
        'S': ("111", "100", "111", "001", "111"),
        'T': ("111", "010", "010", "010", "010"),
        'U': ("101", "101", "101", "101", "111"),
        'V': ("101", "101", "101", "101", "010"),
        'W': ("101", "101", "101", "111", "101"),
        'X': ("101", "101", "010", "101", "101"),
        'Y': ("101", "101", "010", "010", "010"),
        'Z': ("111", "001", "010", "100", "111"),
        ':': ("000", "010", "000", "010", "000"),
        '.': ("000", "000", "000", "000", "010"),
        ',': ("000", "000", "000", "010", "100"),
        '-': ("000", "000", "111", "000", "000"),
        '+': ("000", "010", "111", "010", "000"),
        '/': ("001", "001", "010", "100", "100"),
        '|': ("010", "010", "010", "010", "010"),
        '(': ("010", "100", "100", "100", "010"),
        ')': ("010", "001", "001", "001", "010"),
        '%': ("101", "001", "010", "100", "101"),
        '<': ("001", "010", "100", "010", "001"),
        '>': ("100", "010", "001", "010", "100"),
        '=': ("000", "111", "000", "111", "000"),
        '&': ("010", "101", "010", "101", "011"),
        '^': ("010", "101", "000", "000", "000"),
        '_': ("000", "000", "000", "000", "111"),
        '!': ("010", "010", "010", "000", "010")
    }

    def __init__(self, size: int):
        self.scale = max(1, size // 6)
        self.cache: Dict[Tuple[str, Tuple[int, int, int]], pygame.Surface] = {}

    def line_height(self) -> int:
        return 5 * self.scale

    def render(self, text: str, color: Tuple[int, int, int]) -> pygame.Surface:
        normalized = text.upper()
        key = (normalized, color)
        cached = self.cache.get(key)
        if cached is not None:
            return cached

        glyph_width = 3 * self.scale
        glyph_height = 5 * self.scale
        spacing = self.scale
        total_width = max(1, len(normalized) * (glyph_width + spacing) - spacing)
        surface = pygame.Surface((total_width, glyph_height), pygame.SRCALPHA)

        for index, char in enumerate(normalized):
            glyph = self.GLYPHS.get(char, self.GLYPHS['?'])
            x_origin = index * (glyph_width + spacing)
            for row_index, row in enumerate(glyph):
                for col_index, value in enumerate(row):
                    if value == '1':
                        pixel_rect = (
                            x_origin + col_index * self.scale,
                            row_index * self.scale,
                            self.scale,
                            self.scale
                        )
                        surface.fill(color, pixel_rect)

        if len(self.cache) > 1024:
            self.cache.clear()
        self.cache[key] = surface
        return surface


def load_default_font(size: int) -> UIFont:
    """Try freetype first, then pygame.font."""
    # Prefer pygame.freetype because it works even on some builds where pygame.font fails.
    try:
        import pygame.freetype as pygame_freetype
        if not pygame_freetype.get_init():
            pygame_freetype.init()
        return UIFont('pygame_freetype', pygame_freetype.Font(None, size), size)
    except Exception:
        pass

    try:
        import pygame.ftfont as pygame_ftfont
        if not pygame_ftfont.get_init():
            pygame_ftfont.init()
        return UIFont('pygame_font', pygame_ftfont.Font(None, size), size)
    except Exception:
        pass

    try:
        import pygame.font as pygame_font
        if not pygame_font.get_init():
            pygame_font.init()
        return UIFont('pygame_font', pygame_font.Font(None, size), size)
    except Exception:
        pass

    return UIFont('bitmap', BitmapFont(size), size)


@dataclass
class InteractionRule:
    """Defines how molecules interact based on their prime factors"""
    name: str
    priority: float
    condition: Callable[[Dict[int, int], Dict[int, int]], bool]
    force_function: Callable[[np.ndarray, float, float, float], np.ndarray]
    strength: float = 1.0
    description: str = ""


@dataclass
class ReactionRule:
    """Defines conditions and effects of molecular reactions"""
    name: str
    priority: float
    condition: Callable[[Dict[int, int], Dict[int, int]], bool]
    effect: Callable[['PrimeMolecule', 'PrimeMolecule'], List['PrimeMolecule']]
    probability: float = 0.5
    description: str = ""
    base_probability: Optional[float] = None


GROWTH_PROFILES: Dict[str, Dict[str, float]] = {
    'slow': {
        'time_scale': 0.007,
        'reaction_multiplier': 0.65,
        'damping': 0.965,
        'extraction_multiplier': 0.75
    },
    'normal': {
        'time_scale': 0.01,
        'reaction_multiplier': 1.0,
        'damping': 0.95,
        'extraction_multiplier': 1.0
    },
    'fast': {
        'time_scale': 0.018,
        'reaction_multiplier': 1.45,
        'damping': 0.93,
        'extraction_multiplier': 1.35
    }
}


class SimulationRules:
    """Configuration class for simulation rules and behaviors"""

    def __init__(self):
        self.interaction_rules: List[InteractionRule] = []
        self.reaction_rules: List[ReactionRule] = []
        self.constants = {
            'G': 0.1,  # Gravitational constant
            'k': 1.0,  # Coulomb-like constant
            'min_distance': 0.5,
            'max_force': 10.0,
            'damping': 0.95,
            'temperature_factor': 1.0,
            'quantum_strength': 0.2,
            'time_scale': 1.0,  # Added time scale constant, default is 1.0 (normal speed)
            'reaction_distance': 2.0,
            'threading_enabled': 1.0,
            'thread_workers': float(max(2, min(4, (os.cpu_count() or 4) - 1))),
            'threading_threshold': 240.0,
            'threading_pairwork_threshold': 180000.0,
            'region_partitions': 4.0,
            'interaction_region_range': 1.0,
            'interaction_distance': 0.0,
            'growth_zone_speed': 8.0,
            'growth_zone_jitter': 1.6,
            'growth_zone_radius': 9.0,
            'growth_zone_strength': 0.22,
            'growth_zone_relocation_chance': 0.004,
            'growth_zone_spawn_chance': 0.010,
            'growth_zone_despawn_chance': 0.007,
            'growth_zone_min_count': 3.0,
            'growth_zone_max_count': 10.0,
            'resource_mobility': 0.8,
            'proximity_repulsion_distance': 1.25,
            'proximity_repulsion_strength': 0.55
        }

    def add_interaction_rule(self, rule: InteractionRule):
        self.interaction_rules.append(rule)
        self.interaction_rules.sort(key=lambda x: x.priority, reverse=True)

    def add_reaction_rule(self, rule: ReactionRule):
        if rule.base_probability is None:
            rule.base_probability = rule.probability
        self.reaction_rules.append(rule)
        self.reaction_rules.sort(key=lambda x: x.priority, reverse=True)

    def set_constant(self, name: str, value: float):
        self.constants[name] = value

    def get_constant(self, name: str) -> float:
        return self.constants.get(name, 0.0)

    @staticmethod
    def default_rules() -> 'SimulationRules':
        """Creates a set of default rules for the simulation"""
        rules = SimulationRules()

        # Basic gravitational attraction
        def gravity_condition(factors1, factors2):
            return True  # Always applies

        def gravity_force(direction, distance, mass1, mass2):
            G = rules.get_constant('G')
            return G * mass1 * mass2 * direction / (distance ** 2)

        rules.add_interaction_rule(InteractionRule(
            name="Gravity",
            priority=1.0,
            condition=gravity_condition,
            force_function=gravity_force,
            strength=1.0,
            description="Basic gravitational attraction between all molecules"
        ))

        # Prime factor resonance
        def resonance_condition(factors1, factors2):
            return bool(set(factors1.keys()) & set(factors2.keys()))

        def resonance_force(direction, distance, charge1, charge2):
            k = rules.get_constant('k')
            quantum = rules.get_constant('quantum_strength')
            # Create orbital motion for resonant molecules
            orbital = np.cross(direction, np.array([0, 1, 0]))
            if np.linalg.norm(orbital) > 0:
                orbital = orbital / np.linalg.norm(orbital)
                return quantum * orbital / math.sqrt(distance)
            return np.zeros(3)

        rules.add_interaction_rule(InteractionRule(
            name="Quantum Resonance",
            priority=2.0,
            condition=resonance_condition,
            force_function=resonance_force,
            strength=0.5,
            description="Orbital interactions between molecules with shared prime factors"
        ))

        # Fusion reaction
        def fusion_condition(factors1, factors2):
            shared = set(factors1.keys()) & set(factors2.keys())
            return len(shared) >= 2

        def fusion_effect(mol1, mol2):
            new_number = mol1.number * mol2.number
            new_pos = (mol1.position * mol1.mass + mol2.position * mol2.mass) / (mol1.mass + mol2.mass)
            new_vel = (mol1.velocity * mol1.mass + mol2.velocity * mol2.mass) / (mol1.mass + mol2.mass)
            new_mol = PrimeMolecule(new_number, new_pos)
            new_mol.velocity = new_vel
            return [new_mol]

        rules.add_reaction_rule(ReactionRule(
            name="Fusion",
            priority=1.0,
            condition=fusion_condition,
            effect=fusion_effect,
            probability=0.3,
            description="Fusion of molecules with multiple shared prime factors"
        ))

        return rules


class PrimeMolecule:
    """Represents a molecule in the prime chemistry simulation"""

    def __init__(self, number: int, position: np.ndarray):
        self.number = number
        self.position = np.array(position, dtype=float)
        self.velocity = np.zeros(3)
        self.acceleration = np.zeros(3)
        self.prime_factors = self.factorize(number)
        self.mass = math.log2(number) * 2
        self.charge = self.calculate_charge()
        self.color = self.generate_color()
        self.wealth = max(0.05, self.mass * random.uniform(0.2, 0.9))

    def factorize(self, n):
        factors = defaultdict(int)
        d = 2
        while n > 1:
            while n % d == 0:
                factors[d] += 1
                n //= d
            d += 1
            if d * d > n:
                if n > 1:
                    factors[n] += 1
                break
        return dict(factors)

    def calculate_charge(self):
        charge = 0
        for prime, count in self.prime_factors.items():
            if prime == 2:
                charge += count * 2
            elif prime % 4 == 1:
                charge += count
            else:
                charge -= count
        return charge / (1 + math.log(self.number))

    def generate_color(self):
        if not self.prime_factors:
            return (0.5, 0.5, 0.5)

        color = np.zeros(3)
        max_prime = max(self.prime_factors.keys())

        for prime, count in self.prime_factors.items():
            # Use golden ratio to generate distinct colors
            hue = (prime * 0.618033988749895) % 1.0
            # Convert HSV to RGB (simplified)
            h = hue * 6.0
            c = count / (1 + math.log(max_prime))
            x = c * (1 - abs(h % 2 - 1))

            if 0 <= h < 1:
                rgb = (c, x, 0)
            elif 1 <= h < 2:
                rgb = (x, c, 0)
            elif 2 <= h < 3:
                rgb = (0, c, x)
            elif 3 <= h < 4:
                rgb = (0, x, c)
            elif 4 <= h < 5:
                rgb = (x, 0, c)
            else:
                rgb = (c, 0, x)

            color += np.array(rgb)

        # Normalize and make more vivid
        color = np.clip(color, 0, 1)
        color = color * 0.7 + 0.3

        return tuple(color)

class PrimeChemistry:
    """Main simulation class"""

    def __init__(self, rules: SimulationRules, size=5, molecule_count=10, max_number=1000):
        self.rules = rules
        self.size = size
        self.max_number = max_number
        self.molecules = []
        self.temperature = 1.0
        self.territory_resolution = 96
        self.territory_capacity = 1.0
        base_territory = np.random.normal(
            0.58,
            0.05,
            (self.territory_resolution, self.territory_resolution)
        )
        self.territory = np.clip(base_territory, 0.25, self.territory_capacity)
        grid_axis = np.arange(self.territory_resolution, dtype=float)
        self.territory_grid_x, self.territory_grid_z = np.meshgrid(grid_axis, grid_axis)
        min_zones = max(1, int(self.rules.get_constant('growth_zone_min_count')))
        max_zones = max(min_zones, int(self.rules.get_constant('growth_zone_max_count')))
        initial_zone_count = random.randint(min_zones, max_zones)
        self.growth_zone_positions = np.empty((0, 2), dtype=float)
        self.growth_zone_velocities = np.empty((0, 2), dtype=float)
        self.growth_zone_strengths = np.empty((0,), dtype=float)
        for _ in range(initial_zone_count):
            self._spawn_growth_zone()
        self.extraction_heat = np.zeros_like(self.territory)
        self.total_extracted = 0.0
        self.last_step_extraction = 0.0
        self.extraction_multiplier = 1.0
        self.profile_name = 'normal'
        self.step_count = 0
        self.population_history = deque(maxlen=360)
        self.extraction_history = deque(maxlen=360)
        self.thread_pool: Optional[ThreadPoolExecutor] = None
        self.thread_workers = 1
        self.threading_threshold = 240
        self.threading_pairwork_threshold = 180000
        self.region_partitions = 4
        self.last_pairwork = 0
        self.last_region_count = 0
        self.last_threaded_forces = False

        # Initialize random molecules
        for _ in range(molecule_count):
            pos = [random.uniform(-size / 2, size / 2) for _ in range(3)]
            number = random.randint(2, 100)
            self.molecules.append(PrimeMolecule(number, pos))
        self.population_history.append(len(self.molecules))
        self.extraction_history.append(0.0)
        self.refresh_parallel_configuration()

    @property
    def growth_zone_count(self) -> int:
        return int(self.growth_zone_positions.shape[0])

    def _random_growth_zone_position(self) -> np.ndarray:
        max_coord = self.territory_resolution - 1.0
        return np.random.uniform(0.0, max_coord, 2)

    @staticmethod
    def _random_growth_zone_velocity() -> np.ndarray:
        angle = np.random.uniform(0.0, 2.0 * math.pi)
        return np.array([math.cos(angle), math.sin(angle)], dtype=float)

    def _spawn_growth_zone(self):
        new_position = self._random_growth_zone_position().reshape(1, 2)
        new_velocity = self._random_growth_zone_velocity().reshape(1, 2)
        # Keep slight per-zone variance so zones are not visually identical.
        new_strength = np.array([np.random.uniform(0.75, 1.25)], dtype=float)
        self.growth_zone_positions = np.vstack((self.growth_zone_positions, new_position))
        self.growth_zone_velocities = np.vstack((self.growth_zone_velocities, new_velocity))
        self.growth_zone_strengths = np.concatenate((self.growth_zone_strengths, new_strength))

    def _despawn_growth_zone(self):
        if self.growth_zone_count <= 0:
            return
        remove_idx = random.randrange(self.growth_zone_count)
        self.growth_zone_positions = np.delete(self.growth_zone_positions, remove_idx, axis=0)
        self.growth_zone_velocities = np.delete(self.growth_zone_velocities, remove_idx, axis=0)
        self.growth_zone_strengths = np.delete(self.growth_zone_strengths, remove_idx, axis=0)

    def refresh_parallel_configuration(self):
        """Refresh threading settings from rules and create/update workers."""
        threading_enabled = self.rules.get_constant('threading_enabled') >= 0.5
        desired_workers = max(1, int(self.rules.get_constant('thread_workers')))
        desired_threshold = max(16, int(self.rules.get_constant('threading_threshold')))
        desired_pairwork_threshold = max(10000, int(self.rules.get_constant('threading_pairwork_threshold')))
        desired_partitions = max(1, int(self.rules.get_constant('region_partitions')))

        self.threading_threshold = desired_threshold
        self.threading_pairwork_threshold = desired_pairwork_threshold
        self.region_partitions = desired_partitions

        if not threading_enabled or desired_workers <= 1:
            if self.thread_pool is not None:
                self.thread_pool.shutdown(wait=False, cancel_futures=True)
                self.thread_pool = None
            self.thread_workers = 1
            return

        if self.thread_pool is None or self.thread_workers != desired_workers:
            if self.thread_pool is not None:
                self.thread_pool.shutdown(wait=False, cancel_futures=True)
            self.thread_pool = ThreadPoolExecutor(max_workers=desired_workers)
            self.thread_workers = desired_workers

    def close(self):
        if self.thread_pool is not None:
            self.thread_pool.shutdown(wait=False, cancel_futures=True)
            self.thread_pool = None

    def __del__(self):
        self.close()

    def set_growth_profile(self, profile_name: str, profile_settings: Dict[str, float]):
        self.profile_name = profile_name
        self.extraction_multiplier = profile_settings['extraction_multiplier']

    def _region_index_for_position(self, position: np.ndarray) -> Tuple[int, int]:
        normalized_x = np.clip((position[0] / self.size) + 0.5, 0.0, 0.999999)
        normalized_z = np.clip((position[2] / self.size) + 0.5, 0.0, 0.999999)
        region_x = min(self.region_partitions - 1, int(normalized_x * self.region_partitions))
        region_z = min(self.region_partitions - 1, int(normalized_z * self.region_partitions))
        return region_x, region_z

    def _build_regions(self) -> Dict[Tuple[int, int], List[int]]:
        regions: Dict[Tuple[int, int], List[int]] = defaultdict(list)
        for index, molecule in enumerate(self.molecules):
            regions[self._region_index_for_position(molecule.position)].append(index)
        return dict(regions)

    def _region_force_pairs(self, regions: Dict[Tuple[int, int], List[int]]) -> List[Tuple[List[int], List[int], bool]]:
        region_items = list(regions.items())
        interaction_region_range = max(0, int(self.rules.get_constant('interaction_region_range')))
        pairs: List[Tuple[List[int], List[int], bool]] = []

        for left_idx, (left_coord, left_region) in enumerate(region_items):
            for right_idx in range(left_idx, len(region_items)):
                right_coord, right_region = region_items[right_idx]
                if interaction_region_range > 0:
                    grid_distance = max(abs(left_coord[0] - right_coord[0]), abs(left_coord[1] - right_coord[1]))
                    if grid_distance > interaction_region_range:
                        continue
                pairs.append((left_region, right_region, left_idx == right_idx))

        return pairs

    def _compute_force_chunk(self, region_pair_chunk: List[Tuple[List[int], List[int], bool]]) -> Dict[int, np.ndarray]:
        local_forces: Dict[int, np.ndarray] = defaultdict(lambda: np.zeros(3))
        molecules = self.molecules
        interaction_distance = self.rules.get_constant('interaction_distance')
        distance_limit_sq = interaction_distance * interaction_distance if interaction_distance > 0.0 else 0.0
        repulsion_distance = self.rules.get_constant('proximity_repulsion_distance')
        repulsion_strength = self.rules.get_constant('proximity_repulsion_strength')
        repulsion_distance_sq = repulsion_distance * repulsion_distance if repulsion_distance > 0.0 else 0.0

        for left_indices, right_indices, same_region in region_pair_chunk:
            if same_region:
                for left_offset, i in enumerate(left_indices):
                    mol1 = molecules[i]
                    for j in right_indices[left_offset + 1:]:
                        delta = molecules[j].position - mol1.position
                        distance_sq = float(np.dot(delta, delta))
                        if distance_limit_sq > 0.0 and distance_sq > distance_limit_sq:
                            continue

                        force1, force2 = self.apply_rules(mol1, molecules[j])
                        if repulsion_distance_sq > 0.0 and 1e-8 < distance_sq < repulsion_distance_sq:
                            distance = math.sqrt(distance_sq)
                            direction = delta / distance
                            repulsion_scale = repulsion_strength * (1.0 - distance / repulsion_distance)
                            repulsion_vector = -direction * repulsion_scale
                            force1 += repulsion_vector
                            force2 -= repulsion_vector

                        local_forces[i] += force1
                        local_forces[j] += force2
            else:
                for i in left_indices:
                    mol1 = molecules[i]
                    for j in right_indices:
                        delta = molecules[j].position - mol1.position
                        distance_sq = float(np.dot(delta, delta))
                        if distance_limit_sq > 0.0 and distance_sq > distance_limit_sq:
                            continue

                        force1, force2 = self.apply_rules(mol1, molecules[j])
                        if repulsion_distance_sq > 0.0 and 1e-8 < distance_sq < repulsion_distance_sq:
                            distance = math.sqrt(distance_sq)
                            direction = delta / distance
                            repulsion_scale = repulsion_strength * (1.0 - distance / repulsion_distance)
                            repulsion_vector = -direction * repulsion_scale
                            force1 += repulsion_vector
                            force2 -= repulsion_vector

                        local_forces[i] += force1
                        local_forces[j] += force2

        return local_forces

    @staticmethod
    def _estimate_pairwork(region_pairs: List[Tuple[List[int], List[int], bool]]) -> int:
        pairwork = 0
        for left_indices, right_indices, same_region in region_pairs:
            if same_region:
                count = len(left_indices)
                pairwork += (count * (count - 1)) // 2
            else:
                pairwork += len(left_indices) * len(right_indices)
        return pairwork

    def _compute_forces_sequential(self, regions: Dict[Tuple[int, int], List[int]],
                                   region_pairs: Optional[List[Tuple[List[int], List[int], bool]]] = None) -> List[np.ndarray]:
        molecule_count = len(self.molecules)
        forces = [np.zeros(3) for _ in range(molecule_count)]
        if molecule_count <= 1:
            return forces

        if region_pairs is None:
            region_pairs = self._region_force_pairs(regions)

        local_forces = self._compute_force_chunk(region_pairs)
        for index, force in local_forces.items():
            forces[index] += force

        return forces

    def _compute_forces_threaded(self, regions: Dict[Tuple[int, int], List[int]],
                                 region_pairs: Optional[List[Tuple[List[int], List[int], bool]]] = None) -> List[np.ndarray]:
        molecule_count = len(self.molecules)
        forces = [np.zeros(3) for _ in range(molecule_count)]

        if self.thread_pool is None or molecule_count <= 1:
            return self._compute_forces_sequential(regions, region_pairs)

        if region_pairs is None:
            region_pairs = self._region_force_pairs(regions)
        if not region_pairs:
            return forces

        chunk_target = max(1, self.thread_workers * 2)
        chunk_size = max(1, math.ceil(len(region_pairs) / chunk_target))
        chunks = [region_pairs[start:start + chunk_size] for start in range(0, len(region_pairs), chunk_size)]

        futures = [self.thread_pool.submit(self._compute_force_chunk, chunk) for chunk in chunks]
        for future in futures:
            local_forces = future.result()
            for index, force in local_forces.items():
                forces[index] += force

        return forces

    def _collect_reaction_candidates(self, regions: Dict[Tuple[int, int], List[int]],
                                     reaction_distance: float) -> List[Tuple[float, int, int]]:
        if reaction_distance <= 0.0:
            return []

        candidates: List[Tuple[float, int, int]] = []
        neighbor_offsets = (-1, 0, 1)

        for region, left_indices in regions.items():
            region_x, region_z = region

            for dx in neighbor_offsets:
                for dz in neighbor_offsets:
                    neighbor_region = (region_x + dx, region_z + dz)
                    right_indices = regions.get(neighbor_region)
                    if right_indices is None:
                        continue

                    # Keep each region relation once, including inter-region neighbors.
                    if neighbor_region < region:
                        continue

                    if neighbor_region == region:
                        for left_offset, i in enumerate(left_indices):
                            for j in left_indices[left_offset + 1:]:
                                distance = np.linalg.norm(self.molecules[j].position - self.molecules[i].position)
                                if distance < reaction_distance:
                                    candidates.append((distance, i, j))
                    else:
                        for i in left_indices:
                            for j in right_indices:
                                if i >= j:
                                    continue
                                distance = np.linalg.norm(self.molecules[j].position - self.molecules[i].position)
                                if distance < reaction_distance:
                                    candidates.append((distance, i, j))

        candidates.sort(key=lambda value: value[0])
        return candidates

    def step(self):
        """Perform one step of the simulation"""
        self.refresh_parallel_configuration()
        new_molecules = []
        removed_molecules = set()

        # Apply temperature variation
        self.temperature = 1.0 + 0.1 * math.sin(pygame.time.get_ticks() / 1000)

        time_scale = max(self.rules.get_constant('time_scale'), 1e-4)
        reaction_distance = self.rules.get_constant('reaction_distance')
        step_extracted = 0.0

        # Thermal noise follows sqrt(dt) scaling for stable small timesteps.
        thermal_sigma = 0.01 * self.temperature * math.sqrt(time_scale)
        for molecule in self.molecules:
            molecule.velocity += np.random.normal(0, thermal_sigma, 3)

        regions = self._build_regions()
        self.last_region_count = len(regions)
        region_pairs = self._region_force_pairs(regions)
        pairwork = self._estimate_pairwork(region_pairs)
        self.last_pairwork = pairwork
        should_thread = (
            self.thread_pool is not None and
            len(self.molecules) >= self.threading_threshold and
            pairwork >= self.threading_pairwork_threshold
        )
        self.last_threaded_forces = should_thread

        if should_thread:
            forces = self._compute_forces_threaded(regions, region_pairs)
        else:
            forces = self._compute_forces_sequential(regions, region_pairs)

        # Reactions stay local to a reaction distance, including neighboring regions.
        for _, i, j in self._collect_reaction_candidates(regions, reaction_distance):
            if i in removed_molecules or j in removed_molecules:
                continue
            reaction_products = self.attempt_reactions(self.molecules[i], self.molecules[j], time_scale)
            if reaction_products:
                new_molecules.extend(reaction_products)
                removed_molecules.add(i)
                removed_molecules.add(j)

        # Update positions and velocities
        damping = self.rules.get_constant('damping')
        resource_mobility = self.rules.get_constant('resource_mobility')

        for i, molecule in enumerate(self.molecules):
            if i not in removed_molecules:
                # Integrate acceleration/velocity using the configured timestep.
                molecule.velocity += forces[i] * 0.1 * time_scale
                if resource_mobility > 0.0:
                    molecule.velocity += self.resource_gradient_vector(molecule) * resource_mobility * time_scale
                molecule.position += molecule.velocity * time_scale

                # Apply damping
                molecule.velocity *= damping

                # Boundary conditions
                for axis in range(3):
                    if abs(molecule.position[axis]) > self.size / 2:
                        molecule.position[axis] = np.sign(molecule.position[axis]) * self.size / 2
                        molecule.velocity[axis] *= -0.8

                step_extracted += self.extract_resources(molecule, time_scale)
                molecule.wealth *= max(0.0, 1.0 - 0.005 * time_scale)

        # Update molecule list
        self.molecules = [mol for i, mol in enumerate(self.molecules)
                         if i not in removed_molecules] + new_molecules
        self.last_step_extraction = step_extracted
        self.step_count += 1
        self.population_history.append(len(self.molecules))
        self.extraction_history.append(step_extracted)
        self.regenerate_territory(time_scale)

    def world_to_territory_index(self, position: np.ndarray) -> Tuple[int, int]:
        normalized_x = np.clip((position[0] / self.size) + 0.5, 0.0, 1.0)
        normalized_z = np.clip((position[2] / self.size) + 0.5, 0.0, 1.0)
        grid_x = int(normalized_x * (self.territory_resolution - 1))
        grid_z = int(normalized_z * (self.territory_resolution - 1))
        return grid_x, grid_z

    def _update_growth_zones(self, time_scale: float):
        zone_speed = self.rules.get_constant('growth_zone_speed')
        zone_jitter = self.rules.get_constant('growth_zone_jitter')
        relocation_chance = self.rules.get_constant('growth_zone_relocation_chance')
        spawn_chance = self.rules.get_constant('growth_zone_spawn_chance')
        despawn_chance = self.rules.get_constant('growth_zone_despawn_chance')
        min_zones = max(1, int(self.rules.get_constant('growth_zone_min_count')))
        max_zones = max(min_zones, int(self.rules.get_constant('growth_zone_max_count')))
        effective_dt = max(time_scale, 1e-4)

        if self.growth_zone_count < min_zones:
            for _ in range(min_zones - self.growth_zone_count):
                self._spawn_growth_zone()

        if self.growth_zone_count > max_zones:
            for _ in range(self.growth_zone_count - max_zones):
                self._despawn_growth_zone()

        # Scale chance by dt-like factor so behavior stays close across profiles.
        chance_scale = math.sqrt(effective_dt * 120.0)
        step_spawn_probability = np.clip(spawn_chance * chance_scale, 0.0, 1.0)
        step_despawn_probability = np.clip(despawn_chance * chance_scale, 0.0, 1.0)
        step_relocation_probability = np.clip(relocation_chance * chance_scale, 0.0, 1.0)

        if self.growth_zone_count < max_zones and random.random() < step_spawn_probability:
            self._spawn_growth_zone()
        if self.growth_zone_count > min_zones and random.random() < step_despawn_probability:
            self._despawn_growth_zone()

        if self.growth_zone_count <= 0:
            return

        random_jitter = np.random.normal(0.0, zone_jitter * math.sqrt(effective_dt), self.growth_zone_velocities.shape)
        self.growth_zone_velocities += random_jitter

        velocity_norms = np.linalg.norm(self.growth_zone_velocities, axis=1, keepdims=True)
        velocity_norms = np.where(velocity_norms < 1e-8, 1.0, velocity_norms)
        self.growth_zone_velocities /= velocity_norms

        self.growth_zone_positions += self.growth_zone_velocities * (zone_speed * effective_dt)
        self.growth_zone_strengths += np.random.normal(0.0, 0.03 * math.sqrt(effective_dt), self.growth_zone_strengths.shape)
        self.growth_zone_strengths = np.clip(self.growth_zone_strengths, 0.55, 1.45)

        max_coord = self.territory_resolution - 1.0
        for zone_idx in range(self.growth_zone_count):
            for axis in range(2):
                if self.growth_zone_positions[zone_idx, axis] < 0.0:
                    self.growth_zone_positions[zone_idx, axis] = 0.0
                    self.growth_zone_velocities[zone_idx, axis] *= -1.0
                elif self.growth_zone_positions[zone_idx, axis] > max_coord:
                    self.growth_zone_positions[zone_idx, axis] = max_coord
                    self.growth_zone_velocities[zone_idx, axis] *= -1.0

            if random.random() < step_relocation_probability:
                self.growth_zone_positions[zone_idx] = self._random_growth_zone_position()
                self.growth_zone_velocities[zone_idx] = self._random_growth_zone_velocity()
                self.growth_zone_strengths[zone_idx] = np.random.uniform(0.7, 1.35)

    def _apply_growth_zones(self, time_scale: float):
        if self.growth_zone_count <= 0:
            return
        zone_radius = max(1e-3, self.rules.get_constant('growth_zone_radius'))
        zone_strength = self.rules.get_constant('growth_zone_strength')
        radius_sq = zone_radius * zone_radius

        for zone_position, local_strength in zip(self.growth_zone_positions, self.growth_zone_strengths):
            delta_x = self.territory_grid_x - zone_position[0]
            delta_z = self.territory_grid_z - zone_position[1]
            influence = np.exp(-(delta_x * delta_x + delta_z * delta_z) / (2.0 * radius_sq))
            self.territory += influence * zone_strength * local_strength * time_scale
            self.extraction_heat += influence * 0.05 * time_scale

    def resource_gradient_vector(self, molecule: PrimeMolecule) -> np.ndarray:
        grid_x, grid_z = self.world_to_territory_index(molecule.position)
        left_x = max(0, grid_x - 1)
        right_x = min(self.territory_resolution - 1, grid_x + 1)
        up_z = max(0, grid_z - 1)
        down_z = min(self.territory_resolution - 1, grid_z + 1)

        gradient_x = self.territory[grid_z, right_x] - self.territory[grid_z, left_x]
        gradient_z = self.territory[down_z, grid_x] - self.territory[up_z, grid_x]
        return np.array([gradient_x, 0.0, gradient_z])

    def extract_resources(self, molecule: PrimeMolecule, time_scale: float) -> float:
        grid_x, grid_z = self.world_to_territory_index(molecule.position)
        available = self.territory[grid_z, grid_x]

        extraction_rate = (0.025 + 0.008 * math.log2(max(molecule.number, 2))) * self.extraction_multiplier
        extracted = min(available, extraction_rate * time_scale)

        if extracted <= 0:
            return 0.0

        self.territory[grid_z, grid_x] -= extracted
        self.extraction_heat[grid_z, grid_x] += extracted * 3.0
        molecule.wealth += extracted * 12.0
        self.total_extracted += extracted
        return extracted

    def regenerate_territory(self, time_scale: float):
        regeneration_rate = 0.14 * time_scale
        diffusion_rate = 0.06 * time_scale

        self.territory += (self.territory_capacity - self.territory) * regeneration_rate
        neighbors = (
            np.roll(self.territory, 1, axis=0) +
            np.roll(self.territory, -1, axis=0) +
            np.roll(self.territory, 1, axis=1) +
            np.roll(self.territory, -1, axis=1)
        ) * 0.25
        self.territory = self.territory * (1.0 - diffusion_rate) + neighbors * diffusion_rate
        self._update_growth_zones(time_scale)
        self._apply_growth_zones(time_scale)
        self.territory = np.clip(self.territory, 0.0, self.territory_capacity)
        self.extraction_heat *= max(0.0, 1.0 - 2.4 * time_scale)

    def wealth_stats(self) -> Tuple[float, float, float]:
        if not self.molecules:
            return 0.0, 0.0, 0.0

        wealth_values = np.array([molecule.wealth for molecule in self.molecules])
        return float(np.min(wealth_values)), float(np.mean(wealth_values)), float(np.max(wealth_values))

    def average_speed(self) -> float:
        if not self.molecules:
            return 0.0
        return float(np.mean([np.linalg.norm(molecule.velocity) for molecule in self.molecules]))

    def population_change(self, window_steps: int = 120) -> int:
        if len(self.population_history) < 2:
            return 0
        window = min(window_steps, len(self.population_history) - 1)
        return int(self.population_history[-1] - self.population_history[-1 - window])

    def growth_per_step(self, window_steps: int = 120) -> float:
        if len(self.population_history) < 2:
            return 0.0
        window = min(window_steps, len(self.population_history) - 1)
        if window <= 0:
            return 0.0
        return (self.population_history[-1] - self.population_history[-1 - window]) / window

    def extraction_per_step(self, window_steps: int = 120) -> float:
        if not self.extraction_history:
            return 0.0
        window = min(window_steps, len(self.extraction_history))
        return float(np.mean(list(self.extraction_history)[-window:]))

    def resource_depletion_ratio(self) -> float:
        depleted_mask = self.territory < (0.2 * self.territory_capacity)
        return float(np.mean(depleted_mask))

    def wealth_inequality(self) -> float:
        if not self.molecules:
            return 0.0
        wealth_values = np.array([molecule.wealth for molecule in self.molecules], dtype=float)
        mean_value = float(np.mean(wealth_values))
        if mean_value <= 1e-8:
            return 0.0
        return float(np.std(wealth_values) / mean_value)

    def apply_rules(self, mol1: PrimeMolecule, mol2: PrimeMolecule) -> Tuple[np.ndarray, np.ndarray]:
        """Apply all relevant interaction rules between two molecules"""
        direction = mol2.position - mol1.position
        distance = np.linalg.norm(direction)

        if distance < 0.0001:
            return np.zeros(3), np.zeros(3)

        direction_normalized = direction / distance
        total_force1 = np.zeros(3)
        total_force2 = np.zeros(3)

        for rule in self.rules.interaction_rules:
            if rule.condition(mol1.prime_factors, mol2.prime_factors):
                force = rule.force_function(
                    direction_normalized,
                    distance,
                    mol1.mass if 'mass' in rule.force_function.__code__.co_varnames else mol1.charge,
                    mol2.mass if 'mass' in rule.force_function.__code__.co_varnames else mol2.charge
                )
                force *= rule.strength
                total_force1 += force
                total_force2 -= force

        # Apply force limits
        max_force = self.rules.get_constant('max_force')
        total_force1 = np.clip(total_force1, -max_force, max_force)
        total_force2 = np.clip(total_force2, -max_force, max_force)

        return total_force1, total_force2

    def attempt_reactions(self, mol1: PrimeMolecule, mol2: PrimeMolecule, time_scale: float = 1.0) -> List[PrimeMolecule]:
        """Check and apply all relevant reaction rules"""
        for rule in self.rules.reaction_rules:
            if rule.condition(mol1.prime_factors, mol2.prime_factors):
                scaled_probability = 1.0 - (1.0 - rule.probability) ** max(time_scale, 0.0)
                if random.random() < scaled_probability:
                    products = rule.effect(mol1, mol2)
                    if not products:
                        return []

                    inherited_wealth = (mol1.wealth + mol2.wealth) * 0.85
                    new_entities = [product for product in products if product not in (mol1, mol2)]
                    if new_entities:
                        shared_wealth = inherited_wealth / len(new_entities)
                        for product in new_entities:
                            product.wealth += shared_wealth

                    return products

        return []

class Visualizer:
    """Handles the 3D visualization of the simulation"""

    def __init__(self, width=800, height=600):
        pygame.init()
        self.display = pygame.display.set_mode((width, height), DOUBLEBUF | OPENGL)
        pygame.display.set_caption("Prime Chemistry Simulation")

        # OpenGL initialization
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glEnable(GL_COLOR_MATERIAL)
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)

        # Lighting setup
        glLightfv(GL_LIGHT0, GL_POSITION, (5.0, 5.0, 5.0, 1.0))
        glLightfv(GL_LIGHT0, GL_AMBIENT, (0.2, 0.2, 0.2, 1.0))
        glLightfv(GL_LIGHT0, GL_DIFFUSE, (1.0, 1.0, 1.0, 1.0))

        # Camera setup
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45, (width / height), 0.1, 50.0)

        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glTranslatef(0.0, 0.0, -20)

        # View control
        self.rotation_x = 0
        self.rotation_y = 0
        self.pan_x = 0
        self.pan_y = 0
        self.camera_distance = 20

        # Simulation control
        self.paused = False
        self.auto_step = True
        self.show_vectors = False
        self.show_info = False

        # UI state
        self.selected_molecule = None
        self.font = load_default_font(36)
        if self.font is None:
            print("Warning: pygame font module unavailable. 3D text overlays are disabled.")

    def get_ray_from_mouse(self, mouse_pos):
        """Convert mouse position to a ray in world space"""
        # Get the viewport dimensions
        viewport = glGetIntegerv(GL_VIEWPORT)

        # Get the projection matrix
        projection = glGetDoublev(GL_PROJECTION_MATRIX)

        # Get the modelview matrix
        modelview = glGetDoublev(GL_MODELVIEW_MATRIX)

        # Get window y coordinate (OpenGL y coordinate is inverted)
        window_y = viewport[3] - mouse_pos[1]

        # Get the near and far points in normalized device coordinates
        near = gluUnProject(mouse_pos[0], window_y, 0.0, modelview, projection, viewport)
        far = gluUnProject(mouse_pos[0], window_y, 1.0, modelview, projection, viewport)

        # Calculate ray direction
        ray_dir = np.array(far) - np.array(near)
        ray_dir = ray_dir / np.linalg.norm(ray_dir)

        return np.array(near), ray_dir

    def distance_to_ray(self, point, ray_origin, ray_dir):
        """Calculate the distance from a point to a ray"""
        v = point - ray_origin
        t = np.dot(v, ray_dir)
        p = ray_origin + t * ray_dir
        return np.linalg.norm(point - p)

    def handle_molecule_selection(self):
        """Handle molecule selection when clicking"""
        mouse_pos = pygame.mouse.get_pos()
        ray_origin, ray_dir = self.get_ray_from_mouse(mouse_pos)

        # Transform ray origin based on current view
        glPushMatrix()
        glLoadIdentity()
        glTranslatef(self.pan_x, self.pan_y, -self.camera_distance)
        glRotatef(self.rotation_x, 1, 0, 0)
        glRotatef(self.rotation_y, 0, 1, 0)
        view_matrix = glGetDoublev(GL_MODELVIEW_MATRIX)
        glPopMatrix()

        # Apply view transformation to ray origin
        ray_origin = np.array([ray_origin[0], ray_origin[1], ray_origin[2], 1.0])
        ray_origin = np.dot(view_matrix, ray_origin)
        ray_origin = ray_origin[:3] / ray_origin[3]

        # Find closest molecule to ray
        closest_molecule = None
        min_distance = float('inf')
        for molecule in self.simulation.molecules:
            dist = self.distance_to_ray(molecule.position, ray_origin, ray_dir)
            # Consider molecule radius in selection
            radius = 0.2 + 0.1 * math.log(molecule.number, 2)
            if dist < min_distance and dist < radius + 0.5:  # Add some margin for easier selection
                min_distance = dist
                closest_molecule = molecule

        self.selected_molecule = closest_molecule

    def draw_molecule_info(self):
        """Draw information about the selected molecule"""
        if not self.selected_molecule:
            return

        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        glOrtho(0, 800, 600, 0, -1, 1)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()

        # Draw info box
        info = [
            f"Number: {self.selected_molecule.number}",
            "Prime factors: " + ", ".join(f"{p}^{c}" for p, c in self.selected_molecule.prime_factors.items()),
            f"Mass: {self.selected_molecule.mass:.2f}",
            f"Charge: {self.selected_molecule.charge:.2f}",
            f"Position: ({self.selected_molecule.position[0]:.2f}, "
            f"{self.selected_molecule.position[1]:.2f}, "
            f"{self.selected_molecule.position[2]:.2f})",
            f"Velocity: {np.linalg.norm(self.selected_molecule.velocity):.2f}"
        ]

        # Draw semi-transparent background
        x, y = 10, 10
        padding = 10
        line_height = 25
        box_height = len(info) * line_height + 2 * padding
        box_width = 300

        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glColor4f(0.0, 0.0, 0.0, 0.7)
        glBegin(GL_QUADS)
        glVertex2f(x, y)
        glVertex2f(x + box_width, y)
        glVertex2f(x + box_width, y + box_height)
        glVertex2f(x, y + box_height)
        glEnd()
        glDisable(GL_BLEND)

        # Draw text
        for i, text in enumerate(info):
            self.draw_text(text, (x + padding, y + padding + i * line_height))

        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)
        glPopMatrix()

    def draw(self, simulation: PrimeChemistry, runtime_info: Optional[Dict[str, float]] = None):
        """Updated draw method to include molecule info"""
        # Store simulation reference
        self.simulation = simulation

        # Clear buffers
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()

        # Set up camera
        glTranslatef(self.pan_x, self.pan_y, -self.camera_distance)
        glRotatef(self.rotation_x, 1, 0, 0)
        glRotatef(self.rotation_y, 0, 1, 0)

        # Draw coordinate axes
        self.draw_axes()

        # Draw molecules
        for molecule in simulation.molecules:
            # Highlight selected molecule
            if molecule == self.selected_molecule:
                glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT)
                glEnable(GL_BLEND)
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
                glColor4f(1.0, 1.0, 1.0, 0.3)
                self.draw_molecule(molecule)
                glPopAttrib()
            self.draw_molecule(molecule)

        # Draw boundaries
        self.draw_boundaries(simulation.size)

        # Draw UI overlays
        if self.show_info:
            self.draw_simulation_info(simulation, runtime_info)

        # Draw selected molecule info
        if self.selected_molecule:
            self.draw_molecule_info()

        pygame.display.flip()

    def draw_text(self, text, position, color=(255, 255, 255)):
        """Draw 2D text overlay"""
        if self.font is None:
            return
        text_surface = self.font.render(text, True, color)
        self.display.blit(text_surface, position)

    def handle_events(self):
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                return False

            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE:
                    self.paused = not self.paused
                elif event.key == pygame.K_a:
                    self.auto_step = not self.auto_step
                elif event.key == pygame.K_v:
                    self.show_vectors = not self.show_vectors
                elif event.key == pygame.K_i:
                    self.show_info = not self.show_info
                elif event.key == pygame.K_r:
                    self.reset_view()
                elif event.key == pygame.K_RIGHT and self.paused:
                    return 'step'

            elif event.type == pygame.MOUSEBUTTONDOWN:
                if event.button == 4:  # Mouse wheel up
                    self.camera_distance = max(5, self.camera_distance - 1)
                elif event.button == 5:  # Mouse wheel down
                    self.camera_distance = min(40, self.camera_distance + 1)
                elif event.button == 1:  # Left click
                    self.handle_molecule_selection()

            elif event.type == pygame.MOUSEMOTION:
                if event.buttons[0]:  # Left mouse button - rotate
                    self.rotation_x += event.rel[1]
                    self.rotation_y += event.rel[0]
                elif event.buttons[2]:  # Right mouse button - pan
                    self.pan_x += event.rel[0] * 0.1
                    self.pan_y -= event.rel[1] * 0.1

        return True

    def reset_view(self):
        """Reset camera to default position"""
        self.rotation_x = 0
        self.rotation_y = 0
        self.pan_x = 0
        self.pan_y = 0
        self.camera_distance = 20

    def draw_molecule(self, molecule: PrimeMolecule):
        """Draw a single molecule"""
        glPushMatrix()
        glTranslatef(*molecule.position)

        # Set color with material properties
        glColor3f(*molecule.color)

        # Draw sphere
        radius = 0.2 + 0.1 * math.log(molecule.number, 2)
        quad = gluNewQuadric()
        gluQuadricNormals(quad, GLU_SMOOTH)
        gluSphere(quad, radius, 32, 32)
        gluDeleteQuadric(quad)

        if self.show_vectors:
            # Draw velocity vector
            glBegin(GL_LINES)
            glColor3f(1.0, 1.0, 0.0)  # Yellow for velocity
            glVertex3f(0, 0, 0)
            glVertex3f(*molecule.velocity)
            glEnd()

        glPopMatrix()

    def draw_boundaries(self, size):
        """Draw simulation boundaries"""
        glColor3f(0.5, 0.5, 0.5)
        glBegin(GL_LINES)
        s = size / 2
        for x in [-s, s]:
            for y in [-s, s]:
                for z in [-s, s]:
                    if x == -s: glVertex3f(x, y, z); glVertex3f(s, y, z)
                    if y == -s: glVertex3f(x, y, z); glVertex3f(x, s, z)
                    if z == -s: glVertex3f(x, y, z); glVertex3f(x, y, s)
        glEnd()

    def draw_axes(self):
        """Draw coordinate axes"""
        glBegin(GL_LINES)
        # X axis (red)
        glColor3f(1.0, 0.0, 0.0)
        glVertex3f(0, 0, 0)
        glVertex3f(5, 0, 0)
        # Y axis (green)
        glColor3f(0.0, 1.0, 0.0)
        glVertex3f(0, 0, 0)
        glVertex3f(0, 5, 0)
        # Z axis (blue)
        glColor3f(0.0, 0.0, 1.0)
        glVertex3f(0, 0, 0)
        glVertex3f(0, 0, 5)
        glEnd()

    def draw_simulation_info(self, simulation: PrimeChemistry, runtime_info: Optional[Dict[str, float]] = None):
        """Draw simulation statistics and info"""
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        glOrtho(0, 800, 600, 0, -1, 1)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()

        fps_value = 0.0 if runtime_info is None else float(runtime_info.get('fps', 0.0))
        pop_delta = simulation.population_change(120)
        pop_delta_sign = "+" if pop_delta >= 0 else ""

        info_text = [
            f"Population: {len(simulation.molecules)} ({pop_delta_sign}{pop_delta} /120 steps)",
            f"Temperature: {simulation.temperature:.2f}",
            f"FPS: {fps_value:.1f}",
            f"{'PAUSED' if self.paused else 'RUNNING'}",
            f"{'AUTO' if self.auto_step else 'MANUAL'}"
        ]

        for i, text in enumerate(info_text):
            self.draw_text(text, (10, 30 * i + 10))

        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)
        glPopMatrix()


def create_custom_rules() -> SimulationRules:
    """Create a custom set of simulation rules with more dynamic reactions"""
    rules = SimulationRules()

    # ==== INTERACTION RULES ====

    # Base gravitational attraction (weak but always present)
    def gravity_condition(factors1, factors2):
        return True

    def gravity_force(direction, distance, mass1, mass2):
        G = rules.get_constant('G')
        return G * mass1 * mass2 * direction / (distance ** 2)

    rules.add_interaction_rule(InteractionRule(
        name="Gravity",
        priority=1.0,
        condition=gravity_condition,
        force_function=gravity_force,
        strength=0.1  # Reduced strength to prevent excessive clustering
    ))

    # Prime resonance orbital motion
    def resonance_condition(factors1, factors2):
        shared = set(factors1.keys()) & set(factors2.keys())
        return len(shared) >= 1

    def resonance_force(direction, distance, charge1, charge2):
        orbital = np.cross(direction, np.array([0, 1, 0]))
        if np.linalg.norm(orbital) > 0:
            orbital = orbital / np.linalg.norm(orbital)
            return orbital / (distance ** 0.5)  # Increased distance dependency
        return np.zeros(3)

    rules.add_interaction_rule(InteractionRule(
        name="Prime Resonance",
        priority=2.0,
        condition=resonance_condition,
        force_function=resonance_force,
        strength=0.3
    ))

    # ==== REACTION RULES ====

    # Fusion with emission
    def fusion_emission_condition(factors1, factors2):
        shared = set(factors1.keys()) & set(factors2.keys())
        return len(shared) >= 2 and (factors1[2] if 2 in factors1 else 0) + (factors2[2] if 2 in factors2 else 0) >= 1

    def fusion_emission_effect(mol1, mol2):
        # Calculate new main molecule
        new_number = mol1.number * mol2.number
        shared_primes = set(mol1.prime_factors.keys()) & set(mol2.prime_factors.keys())
        emitted_prime = min(shared_primes)  # Emit the smallest shared prime
        new_number = new_number // emitted_prime

        # Create main product
        new_pos = (mol1.position * mol1.mass + mol2.position * mol2.mass) / (mol1.mass + mol2.mass)
        main_product = PrimeMolecule(new_number, new_pos)

        # Create several emitted particles with high velocity
        particles = []
        for _ in range(3):  # Emit 3 particles
            emission_direction = np.random.normal(0, 1, 3)
            emission_direction = emission_direction / np.linalg.norm(emission_direction)
            emission_pos = new_pos + emission_direction * 0.5
            particle = PrimeMolecule(emitted_prime, emission_pos)
            particle.velocity = emission_direction * 3.0  # High velocity
            particles.append(particle)

        return [main_product] + particles

    rules.add_reaction_rule(ReactionRule(
        name="Fusion with Emission",
        priority=2.0,
        condition=fusion_emission_condition,
        effect=fusion_emission_effect,
        probability=0.4
    ))

    # Fission reaction for large molecules
    def fission_condition(factors1, factors2):
        # Calculate the total number for each molecule based on their prime factors
        num1 = math.prod(prime ** count for prime, count in factors1.items())
        num2 = math.prod(prime ** count for prime, count in factors2.items())
        # React if one molecule is significantly larger than the other
        return max(num1, num2) > 100 and min(num1, num2) < 50

    def fission_effect(mol1, mol2):
        # Split the larger molecule
        bigger = mol1 if mol1.number > mol2.number else mol2
        smaller = mol2 if mol1.number > mol2.number else mol1

        # Find prime factors
        factors = list(bigger.prime_factors.items())
        if len(factors) < 2:
            return []

        # Split into roughly equal parts
        product1 = 1
        product2 = 1
        for prime, count in factors:
            if product1 <= product2:
                product1 *= prime ** count
            else:
                product2 *= prime ** count

        # Create new molecules with opposite velocities
        vel = np.random.normal(0, 1, 3)
        vel = vel / np.linalg.norm(vel) * 2.0

        mol1 = PrimeMolecule(product1, bigger.position + np.random.normal(0, 0.1, 3))
        mol2 = PrimeMolecule(product2, bigger.position + np.random.normal(0, 0.1, 3))
        mol1.velocity = vel
        mol2.velocity = -vel

        # The smaller molecule gets deflected
        smaller.velocity = np.random.normal(0, 1, 3)

        return [mol1, mol2, smaller]

    rules.add_reaction_rule(ReactionRule(
        name="Fission",
        priority=1.5,
        condition=fission_condition,
        effect=fission_effect,
        probability=0.3
    ))

    # Catalytic decomposition
    def catalytic_condition(factors1, factors2):
        # One molecule must be prime (catalyst) and the other composite
        return (len(factors1) == 1 and len(factors2) > 1) or (len(factors2) == 1 and len(factors1) > 1)

    def catalytic_effect(mol1, mol2):
        # Identify catalyst and target
        catalyst = mol1 if len(mol1.prime_factors) == 1 else mol2
        target = mol2 if len(mol1.prime_factors) == 1 else mol1

        # Decompose target into prime factors
        new_molecules = []
        for prime, count in target.prime_factors.items():
            for _ in range(count):
                pos = target.position + np.random.normal(0, 0.5, 3)
                mol = PrimeMolecule(prime, pos)
                # Particles fly out in random directions
                vel = np.random.normal(0, 1, 3)
                mol.velocity = vel / np.linalg.norm(vel) * 2.0
                new_molecules.append(mol)

        # Catalyst remains unchanged but gets a small velocity boost
        catalyst.velocity += np.random.normal(0, 0.5, 3)
        new_molecules.append(catalyst)

        return new_molecules

    rules.add_reaction_rule(ReactionRule(
        name="Catalytic Decomposition",
        priority=1.0,
        condition=catalytic_condition,
        effect=catalytic_effect,
        probability=0.25
    ))

    # Domain decomposition defaults for threaded force computation.
    rules.set_constant('region_partitions', 5.0)
    rules.set_constant('interaction_region_range', 1.0)
    rules.set_constant('threading_threshold', 240.0)
    rules.set_constant('threading_pairwork_threshold', 180000.0)
    rules.set_constant('growth_zone_speed', 8.5)
    rules.set_constant('growth_zone_jitter', 1.7)
    rules.set_constant('growth_zone_radius', 8.0)
    rules.set_constant('growth_zone_strength', 0.24)
    rules.set_constant('growth_zone_relocation_chance', 0.0045)
    rules.set_constant('growth_zone_spawn_chance', 0.012)
    rules.set_constant('growth_zone_despawn_chance', 0.009)
    rules.set_constant('growth_zone_min_count', 3.0)
    rules.set_constant('growth_zone_max_count', 12.0)
    rules.set_constant('resource_mobility', 0.85)
    rules.set_constant('proximity_repulsion_distance', 1.3)
    rules.set_constant('proximity_repulsion_strength', 0.6)
    rules.set_constant('reaction_distance', 1.7)

    return rules


def apply_growth_profile(rules: SimulationRules, profile_name: str) -> Tuple[str, Dict[str, float]]:
    normalized_profile = profile_name.lower()
    if normalized_profile not in GROWTH_PROFILES:
        normalized_profile = 'normal'

    settings = GROWTH_PROFILES[normalized_profile]
    rules.set_constant('time_scale', settings['time_scale'])
    rules.set_constant('damping', settings['damping'])

    for rule in rules.reaction_rules:
        base_probability = rule.base_probability if rule.base_probability is not None else rule.probability
        rule.probability = min(0.95, max(0.001, base_probability * settings['reaction_multiplier']))

    return normalized_profile, settings


class TerritoryVisualizer:
    """2D territorial view: resources, extraction heat, and wealth-colored population."""

    def __init__(self, width=1280, height=800):
        pygame.init()
        self.display = pygame.display.set_mode((width, height))
        pygame.display.set_caption("Game of Primes: Territory Economy")
        self.width = width
        self.height = height

        self.paused = False
        self.auto_step = True
        self.show_vectors = True
        self.show_info = True
        self.profile_name = 'normal'

        self.font = load_default_font(18)
        self.title_font = load_default_font(26)

        # Bitmap fallback needs smaller nominal sizes to avoid dense overlap in side panels.
        if self.font.backend == 'bitmap':
            self.font = load_default_font(14)
        if self.title_font.backend == 'bitmap':
            self.title_font = load_default_font(20)

        territory_px = min(height - 40, width - 360)
        self.territory_rect = pygame.Rect(20, 20, territory_px, territory_px)
        self.sidebar_rect = pygame.Rect(self.territory_rect.right + 20, 20,
                                        width - self.territory_rect.right - 40, height - 40)

        self.body_line_step = self.font.line_height() + 4
        self.title_line_step = self.title_font.line_height() + 6

    def set_profile(self, profile_name: str):
        self.profile_name = profile_name

    def handle_events(self):
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                return False
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE:
                    self.paused = not self.paused
                elif event.key == pygame.K_a:
                    self.auto_step = not self.auto_step
                elif event.key == pygame.K_i:
                    self.show_info = not self.show_info
                elif event.key == pygame.K_v:
                    self.show_vectors = not self.show_vectors
                elif event.key == pygame.K_RIGHT and self.paused:
                    return 'step'
                elif event.key == pygame.K_1:
                    return ('profile', 'slow')
                elif event.key == pygame.K_2:
                    return ('profile', 'normal')
                elif event.key == pygame.K_3:
                    return ('profile', 'fast')

        return True

    @staticmethod
    def wealth_to_color(wealth: float, max_wealth: float) -> Tuple[int, int, int]:
        if max_wealth <= 0:
            return 40, 140, 220

        t = np.clip(wealth / max_wealth, 0.0, 1.0)
        low = np.array([40, 130, 220], dtype=float)
        mid = np.array([90, 220, 140], dtype=float)
        high = np.array([250, 205, 70], dtype=float)

        if t < 0.5:
            local_t = t / 0.5
            color = low + (mid - low) * local_t
        else:
            local_t = (t - 0.5) / 0.5
            color = mid + (high - mid) * local_t

        return tuple(color.astype(int))

    def world_to_screen(self, simulation: PrimeChemistry, position: np.ndarray) -> Tuple[int, int]:
        x_norm = np.clip((position[0] / simulation.size) + 0.5, 0.0, 1.0)
        y_norm = np.clip((position[2] / simulation.size) + 0.5, 0.0, 1.0)

        x_pos = int(self.territory_rect.left + x_norm * self.territory_rect.width)
        y_pos = int(self.territory_rect.top + y_norm * self.territory_rect.height)
        return x_pos, y_pos

    def build_territory_surface(self, simulation: PrimeChemistry) -> pygame.Surface:
        richness = np.clip(simulation.territory / simulation.territory_capacity, 0.0, 1.0)
        heat = np.clip(simulation.extraction_heat, 0.0, 1.0)

        red = np.clip(22 + richness * 70 + heat * 180, 0, 255)
        green = np.clip(30 + richness * 160 + heat * 40, 0, 255)
        blue = np.clip(20 + richness * 70, 0, 255)

        terrain_rgb = np.stack([red, green, blue], axis=-1).astype(np.uint8)
        terrain_surface = pygame.surfarray.make_surface(np.transpose(terrain_rgb, (1, 0, 2)))
        return pygame.transform.smoothscale(terrain_surface, self.territory_rect.size)

    def legend_rect(self) -> pygame.Rect:
        return pygame.Rect(self.sidebar_rect.left + 20, self.sidebar_rect.bottom - 80,
                           self.sidebar_rect.width - 40, 18)

    def draw_wealth_legend(self, max_wealth: float):
        legend_rect = self.legend_rect()
        for offset in range(legend_rect.width):
            ratio = offset / max(1, legend_rect.width - 1)
            color = self.wealth_to_color(ratio * max(max_wealth, 1.0), max(max_wealth, 1.0))
            pygame.draw.line(
                self.display,
                color,
                (legend_rect.left + offset, legend_rect.top),
                (legend_rect.left + offset, legend_rect.bottom)
            )
        pygame.draw.rect(self.display, (225, 225, 225), legend_rect, 1)
        if self.font is None:
            return

        self.display.blit(self.font.render("Low wealth", True, (220, 220, 220)),
                          (legend_rect.left, legend_rect.bottom + 6))
        high_text = self.font.render("High wealth", True, (220, 220, 220))
        self.display.blit(high_text, (legend_rect.right - high_text.get_width(), legend_rect.bottom + 6))

    def draw_sidebar(self, simulation: PrimeChemistry, runtime_info: Optional[Dict[str, float]] = None):
        pygame.draw.rect(self.display, (24, 24, 32), self.sidebar_rect)
        pygame.draw.rect(self.display, (72, 72, 88), self.sidebar_rect, 2)

        runtime_info = runtime_info or {}
        fps_value = float(runtime_info.get('fps', 0.0))
        status_text = 'PAUSED' if self.paused else 'RUNNING'
        stepping_text = 'AUTO' if self.auto_step else 'MANUAL'
        window_steps = min(120, max(1, simulation.step_count))
        pop_delta = simulation.population_change(window_steps)
        pop_delta_sign = "+" if pop_delta >= 0 else ""
        min_wealth, avg_wealth, max_wealth = simulation.wealth_stats()
        wealth_cv = simulation.wealth_inequality()
        growth_per_step = simulation.growth_per_step(window_steps)
        avg_extraction = simulation.extraction_per_step(window_steps)
        avg_speed = simulation.average_speed()
        depleted_ratio = simulation.resource_depletion_ratio()
        hotspot_heat = float(np.max(simulation.extraction_heat))
        richness = float(np.mean(simulation.territory))

        lines = [
            "Territory Economy Dashboard",
            f"Status: {status_text} | Step mode: {stepping_text} | FPS: {fps_value:.1f}",
            f"Profile: {self.profile_name.upper()} | Steps: {simulation.step_count}",
            f"Domain split: {simulation.region_partitions}x{simulation.region_partitions} | Active regions: {simulation.last_region_count}",
            f"Threaded forces: {'ON' if simulation.last_threaded_forces else 'OFF'} | Pair work: {simulation.last_pairwork}",
            f"Proximity repulsion: d={simulation.rules.get_constant('proximity_repulsion_distance'):.2f}",
            "",
            "Population & Motion",
            f"Population: {len(simulation.molecules)} ({pop_delta_sign}{pop_delta} over {window_steps} steps)",
            f"Growth velocity: {growth_per_step:+.3f} molecules/step",
            f"Average speed: {avg_speed:.3f} units/step",
            f"Temperature: {simulation.temperature:.2f}",
            "",
            "Economy & Extraction",
            f"Avg wealth: {avg_wealth:.2f} | Range: {min_wealth:.2f} -> {max_wealth:.2f}",
            f"Wealth inequality (CV): {wealth_cv:.3f}",
            f"Extraction current step: {simulation.last_step_extraction:.4f}",
            f"Extraction avg ({window_steps}): {avg_extraction:.4f}",
            f"Total extracted: {simulation.total_extracted:.2f}",
            "",
            "Territory Health",
            f"Richness mean: {richness:.3f}",
            f"Depleted area (<20%): {depleted_ratio * 100:.1f}%",
            f"Hotspot intensity: {hotspot_heat:.3f}",
            f"Moving growth zones: {simulation.growth_zone_count} "
            f"({int(simulation.rules.get_constant('growth_zone_min_count'))}-"
            f"{int(simulation.rules.get_constant('growth_zone_max_count'))})",
            "",
            "Controls:",
            "SPACE pause/resume",
            "RIGHT single-step when paused",
            "A toggle auto-step",
            "V toggle velocity vectors",
            "1 slow growth | 2 normal | 3 fast",
        ]

        if self.font is not None and self.title_font is not None:
            y_pos = self.sidebar_rect.top + 16
            max_text_bottom = self.legend_rect().top - 12
            previous_clip = self.display.get_clip()
            self.display.set_clip(self.sidebar_rect.inflate(-12, -12))
            for index, line in enumerate(lines):
                if index == 0:
                    label = self.title_font.render(line, True, (245, 245, 250))
                    step_height = self.title_line_step
                else:
                    label = self.font.render(line, True, (215, 215, 225))
                    step_height = self.body_line_step

                if y_pos + step_height > max_text_bottom:
                    overflow_label = self.font.render("...", True, (215, 215, 225))
                    self.display.blit(overflow_label, (self.sidebar_rect.left + 16, y_pos))
                    break

                self.display.blit(label, (self.sidebar_rect.left + 16, y_pos))
                y_pos += step_height
            self.display.set_clip(previous_clip)

        self.draw_wealth_legend(max_wealth)

    def draw(self, simulation: PrimeChemistry, runtime_info: Optional[Dict[str, float]] = None):
        self.display.fill((12, 14, 20))

        terrain_surface = self.build_territory_surface(simulation)
        self.display.blit(terrain_surface, self.territory_rect.topleft)
        pygame.draw.rect(self.display, (235, 235, 235), self.territory_rect, 2)

        zone_radius_cells = simulation.rules.get_constant('growth_zone_radius')
        zone_radius_px = max(4, int(zone_radius_cells * self.territory_rect.width / simulation.territory_resolution))
        for zone_position, zone_strength in zip(simulation.growth_zone_positions, simulation.growth_zone_strengths):
            zone_x = int(self.territory_rect.left + (zone_position[0] / (simulation.territory_resolution - 1))
                         * self.territory_rect.width)
            zone_y = int(self.territory_rect.top + (zone_position[1] / (simulation.territory_resolution - 1))
                         * self.territory_rect.height)
            tint = float(np.clip((zone_strength - 0.55) / 0.9, 0.0, 1.0))
            ring_color = (
                int(220 + 35 * tint),
                int(160 + 70 * tint),
                int(95 + 95 * tint)
            )
            core_radius = max(2, int(2 + 2.5 * tint))
            pygame.draw.circle(self.display, ring_color, (zone_x, zone_y), zone_radius_px, 1)
            pygame.draw.circle(self.display, (255, 235, 190), (zone_x, zone_y), core_radius)

        if simulation.molecules:
            wealth_cap = max(molecule.wealth for molecule in simulation.molecules)
        else:
            wealth_cap = 1.0

        for molecule in simulation.molecules:
            x_pos, y_pos = self.world_to_screen(simulation, molecule.position)
            radius = max(2, int(2 + 0.35 * math.log2(max(molecule.number, 2))))
            color = self.wealth_to_color(molecule.wealth, wealth_cap)

            if self.show_vectors:
                vector_scale = 22
                velocity_end = (
                    int(x_pos + molecule.velocity[0] * vector_scale),
                    int(y_pos + molecule.velocity[2] * vector_scale)
                )
                pygame.draw.line(self.display, (245, 245, 245), (x_pos, y_pos), velocity_end, 1)

            pygame.draw.circle(self.display, color, (x_pos, y_pos), radius)
            pygame.draw.circle(self.display, (16, 16, 20), (x_pos, y_pos), radius, 1)

        if self.show_info:
            self.draw_sidebar(simulation, runtime_info)

        pygame.display.flip()


def main(fps=60, growth_profile='normal', view_mode='territory'):
    """Main simulation loop"""
    rules = create_custom_rules()
    profile_name, profile_settings = apply_growth_profile(rules, growth_profile)

    simulation = PrimeChemistry(rules, size=1000, molecule_count=1000, max_number=10000)
    simulation.set_growth_profile(profile_name, profile_settings)

    if view_mode == '3d':
        visualizer = Visualizer(width=1024, height=768)
    else:
        visualizer = TerritoryVisualizer(width=1280, height=800)
        visualizer.set_profile(profile_name)

    clock = pygame.time.Clock()
    frame_time = 1.0 / fps
    accumulated_time = 0.0
    last_time = pygame.time.get_ticks() / 1000.0

    running = True
    while running:
        current_time = pygame.time.get_ticks() / 1000.0
        delta_time = current_time - last_time
        last_time = current_time
        accumulated_time += delta_time

        event_result = visualizer.handle_events()
        if event_result is False:
            running = False
            continue

        if isinstance(event_result, tuple) and event_result[0] == 'profile':
            profile_name, profile_settings = apply_growth_profile(rules, event_result[1])
            simulation.set_growth_profile(profile_name, profile_settings)
            if hasattr(visualizer, 'set_profile'):
                visualizer.set_profile(profile_name)

        if not visualizer.paused:
            while accumulated_time >= frame_time:
                if visualizer.auto_step:
                    simulation.step()
                accumulated_time -= frame_time
        elif event_result == 'step':
            simulation.step()

        runtime_info = {
            'fps': clock.get_fps()
        }
        visualizer.draw(simulation, runtime_info)
        clock.tick(fps)

    simulation.close()
    pygame.quit()


def parse_args():
    parser = argparse.ArgumentParser(description="Prime chemistry and territory economy simulation")
    parser.add_argument('--fps', type=int, default=30, help='Target frames per second')
    parser.add_argument('--growth', choices=sorted(GROWTH_PROFILES.keys()),
                        default='normal', help='Growth profile tuning')
    parser.add_argument('--view', choices=['territory', '3d'], default='territory',
                        help='Visualization mode (territory is default)')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(fps=args.fps, growth_profile=args.growth, view_mode=args.view)
