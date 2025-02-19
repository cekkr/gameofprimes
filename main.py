import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
import numpy as np
import random
import math
from collections import defaultdict
from dataclasses import dataclass
from typing import Callable, List, Dict, Optional, Tuple
import json


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
            'time_scale': 1.0  # Added time scale constant, default is 1.0 (normal speed)
        }

    def add_interaction_rule(self, rule: InteractionRule):
        self.interaction_rules.append(rule)
        self.interaction_rules.sort(key=lambda x: x.priority, reverse=True)

    def add_reaction_rule(self, rule: ReactionRule):
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
        self.accumulated_time = 0.0  # Added to track accumulated time for reactions

        # Initialize random molecules
        for _ in range(molecule_count):
            pos = [random.uniform(-size / 2, size / 2) for _ in range(3)]
            number = random.randint(2, 100)
            self.molecules.append(PrimeMolecule(number, pos))

    def step(self):
        """Perform one step of the simulation"""
        new_molecules = []
        removed_molecules = set()

        # Apply temperature variation
        self.temperature = 1.0 + 0.1 * math.sin(pygame.time.get_ticks() / 1000)

        # Calculate forces and check for reactions
        forces = defaultdict(lambda: np.zeros(3))
        time_scale = self.rules.get_constant('time_scale')

        # Update accumulated time
        self.accumulated_time += time_scale

        for i, mol1 in enumerate(self.molecules):
            if i in removed_molecules:
                continue

            # Random thermal motion - scaled by time_scale
            mol1.velocity += np.random.normal(0, 0.01 * self.temperature * time_scale, 3)

            for j, mol2 in enumerate(self.molecules[i + 1:], i + 1):
                if j in removed_molecules:
                    continue

                # Apply interaction rules - forces are scaled by time_scale
                force1, force2 = self.apply_rules(mol1, mol2)
                forces[i] += force1 * time_scale
                forces[j] += force2 * time_scale

                # Check for reactions - only attempt if enough time has accumulated
                if np.linalg.norm(mol2.position - mol1.position) < 2.0:
                    if self.accumulated_time >= 1.0:  # Only attempt reactions after accumulating enough time
                        reaction_products = self.attempt_reactions(mol1, mol2)
                        if reaction_products:
                            new_molecules.extend(reaction_products)
                            removed_molecules.add(i)
                            removed_molecules.add(j)
                            break

        # Reset accumulated time if it exceeds 1.0
        if self.accumulated_time >= 1.0:
            self.accumulated_time = 0.0

        # Update positions and velocities
        damping = self.rules.get_constant('damping')

        for i, molecule in enumerate(self.molecules):
            if i not in removed_molecules:
                # Apply forces - scaled movement by time_scale
                molecule.velocity += forces[i] * 0.1
                molecule.position += molecule.velocity * time_scale

                # Apply damping
                molecule.velocity *= damping

                # Boundary conditions
                for axis in range(3):
                    if abs(molecule.position[axis]) > self.size / 2:
                        molecule.position[axis] = np.sign(molecule.position[axis]) * self.size / 2
                        molecule.velocity[axis] *= -0.8

        # Update molecule list
        self.molecules = [mol for i, mol in enumerate(self.molecules)
                         if i not in removed_molecules] + new_molecules

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

    def attempt_reactions(self, mol1: PrimeMolecule, mol2: PrimeMolecule) -> List[PrimeMolecule]:
        """Check and apply all relevant reaction rules"""
        distance = np.linalg.norm(mol2.position - mol1.position)

        for rule in self.rules.reaction_rules:
            if rule.condition(mol1.prime_factors, mol2.prime_factors):
                if random.random() < rule.probability:
                    return rule.effect(mol1, mol2)

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
        self.font = pygame.font.Font(None, 36)

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

    def draw(self, simulation: PrimeChemistry):
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
            self.draw_simulation_info(simulation)

        # Draw selected molecule info
        if self.selected_molecule:
            self.draw_molecule_info()

        pygame.display.flip()

    def draw_text(self, text, position, color=(255, 255, 255)):
        """Draw 2D text overlay"""
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

    def draw_simulation_info(self, simulation: PrimeChemistry):
        """Draw simulation statistics and info"""
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        glOrtho(0, 800, 600, 0, -1, 1)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()

        info_text = [
            f"Molecules: {len(simulation.molecules)}",
            f"Temperature: {simulation.temperature:.2f}",
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
        num1 = np.prod([prime ** count for prime, count in factors1.items()])
        num2 = np.prod([prime ** count for prime, count in factors2.items()])
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

    return rules


def main(fps=60):
    """Main simulation loop"""
    # Create custom rules
    rules = create_custom_rules()
    rules.set_constant('time_scale', 0.01)

    # Initialize simulation
    simulation = PrimeChemistry(rules, size=10, molecule_count=50, max_number=1000)
    visualizer = Visualizer(width=1024, height=768)

    # Timing variables
    clock = pygame.time.Clock()
    frame_time = 1.0 / fps
    accumulated_time = 0
    last_time = pygame.time.get_ticks() / 1000.0

    running = True
    while running:
        current_time = pygame.time.get_ticks() / 1000.0
        delta_time = current_time - last_time
        last_time = current_time
        accumulated_time += delta_time

        # Handle events
        event_result = visualizer.handle_events()
        if event_result == False:
            running = False

        # Update simulation
        if not visualizer.paused:
            while accumulated_time >= frame_time:
                if visualizer.auto_step:
                    simulation.step()
                accumulated_time -= frame_time
        elif event_result == 'step':
            simulation.step()

        # Render
        visualizer.draw(simulation)
        clock.tick(fps)

    pygame.quit()


if __name__ == "__main__":
    main(fps=30)